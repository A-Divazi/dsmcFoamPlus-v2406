/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "boundBoxNew.H"
#include "PstreamReduceOps.H"
#include "tmp.H"
#include "plane.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::boundBoxNew Foam::boundBoxNew::greatBox
(
    point::uniform(-ROOTVGREAT),
    point::uniform(ROOTVGREAT)
);

const Foam::boundBoxNew Foam::boundBoxNew::invertedBox
(
    point::uniform(ROOTVGREAT),
    point::uniform(-ROOTVGREAT)
);

const Foam::faceList Foam::boundBoxNew::faces
({
    // Point and face order as per hex cellmodel
    face({0, 4, 7, 3}),  // 0: x-min, left
    face({1, 2, 6, 5}),  // 1: x-max, right
    face({0, 1, 5, 4}),  // 2: y-min, bottom
    face({3, 7, 6, 2}),  // 3: y-max, top
    face({0, 3, 2, 1}),  // 4: z-min, back
    face({4, 5, 6, 7})   // 5: z-max, front
});

const Foam::FixedList<Foam::vector, 6> Foam::boundBoxNew::faceNormals
({
    vector(-1,  0,  0), // 0: x-min, left
    vector( 1,  0,  0), // 1: x-max, right
    vector( 0, -1,  0), // 2: y-min, bottom
    vector( 0,  1,  0), // 3: y-max, top
    vector( 0,  0, -1), // 4: z-min, back
    vector( 0,  0,  1)  // 5: z-max, front
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::boundBoxNew::calculate(const UList<point>& points, const bool doReduce)
{
    if (points.empty())
    {
        min_ = point::zero;
        max_ = point::zero;

        if (doReduce && Pstream::parRun())
        {
            // Use values that get overwritten by reduce minOp, maxOp below
            min_ = point(VGREAT, VGREAT, VGREAT);
            max_ = point(-VGREAT, -VGREAT, -VGREAT);
        }
    }
    else
    {
        min_ = points[0];
        max_ = points[0];

        for (label i = 1; i < points.size(); i++)
        {
            min_ = ::Foam::min(min_, points[i]);
            max_ = ::Foam::max(max_, points[i]);
        }
    }

    // Reduce parallel information
    if (doReduce)
    {
        //reduce(min_, minOp<point>());
        //reduce(max_, maxOp<point>());
        reduce();
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundBoxNew::boundBoxNew(const UList<point>& points, bool doReduce)
:
    boundBoxNew()
{
    add(points);

    if (doReduce)
    {
        reduce();
    }
}


Foam::boundBoxNew::boundBoxNew(const tmp<pointField>& tpoints, bool doReduce)
:
    boundBoxNew()
{
    add(tpoints);

    if (doReduce)
    {
        reduce();
    }
}


Foam::boundBoxNew::boundBoxNew
(
    const UList<point>& points,
    const labelUList& indices,
    bool doReduce
)
:
    boundBoxNew()
{
    add(points, indices);

    if (doReduce)
    {
        reduce();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::boundBoxNew::points() const
{
    auto tpt = tmp<pointField>::New(8);
    auto& pt = tpt.ref();

    pt[0] = min_;                                   // min-x, min-y, min-z
    pt[1] = point(max_.x(), min_.y(), min_.z());    // max-x, min-y, min-z
    pt[2] = point(max_.x(), max_.y(), min_.z());    // max-x, max-y, min-z
    pt[3] = point(min_.x(), max_.y(), min_.z());    // min-x, max-y, min-z
    pt[4] = point(min_.x(), min_.y(), max_.z());    // min-x, min-y, max-z
    pt[5] = point(max_.x(), min_.y(), max_.z());    // max-x, min-y, max-z
    pt[6] = max_;                                   // max-x, max-y, max-z
    pt[7] = point(min_.x(), max_.y(), max_.z());    // min-x, max-y, max-z

    return tpt;
}

void Foam::boundBoxNew::inflate(const scalar s)
{
    vector ext = vector::one*s*mag();

    min_ -= ext;
    max_ += ext;
}

void Foam::boundBoxNew::expand(const scalar s)
{
    vector v = max_-min_;
    vector unitV = v/Foam::mag(v);

    min_ -= unitV*s*Foam::sqrt(3.0);
    max_ += unitV*s*Foam::sqrt(3.0);
}

void Foam::boundBoxNew::expandII(const scalar s)
{
//     Info << "(before) min_ " << min_ << endl;
//     Info << "(before) max_ " << max_ << endl;
    min_ -= vector(1,0,0)*s;
    min_ -= vector(0,1,0)*s;
    min_ -= vector(0,0,1)*s;    
    
    max_ += vector(1,0,0)*s;
    max_ += vector(0,1,0)*s;
    max_ += vector(0,0,1)*s;

//     Info << "(after) min_ " << min_ << endl;
//     Info << "(after) max_ " << max_ << endl;

}

void Foam::boundBoxNew::contractII(const scalar s)
{
//     Info << "(before) min_ " << min_ << endl;
//     Info << "(before) max_ " << max_ << endl;
    min_ += vector(1,0,0)*s;
    min_ += vector(0,1,0)*s;
    min_ += vector(0,0,1)*s;    
    
    max_ -= vector(1,0,0)*s;
    max_ -= vector(0,1,0)*s;
    max_ -= vector(0,0,1)*s;

//     Info << "(after) min_ " << min_ << endl;
//     Info << "(after) max_ " << max_ << endl;

}

Foam::tmp<Foam::pointField> Foam::boundBoxNew::faceCentres() const
{
    auto tpts = tmp<pointField>::New(6);
    auto& pts = tpts.ref();

    forAll(pts, facei)
    {
        pts[facei] = faceCentre(facei);
    }

    return tpts;
}


Foam::point Foam::boundBoxNew::faceCentre(const direction facei) const
{
    point pt = boundBoxNew::centre();

    if (facei > 5)
    {
        FatalErrorInFunction
            << "face should be [0..5]"
            << abort(FatalError);
    }

    switch (facei)
    {
        case 0: pt.x() = min().x(); break;  // 0: x-min, left
        case 1: pt.x() = max().x(); break;  // 1: x-max, right
        case 2: pt.y() = min().y(); break;  // 2: y-min, bottom
        case 3: pt.y() = max().y(); break;  // 3: y-max, top
        case 4: pt.z() = min().z(); break;  // 4: z-min, back
        case 5: pt.z() = max().z(); break;  // 5: z-max, front
    }

    return pt;
}

void Foam::boundBoxNew::reduce()
{
    Foam::reduce(min_, minOp<point>());
    Foam::reduce(max_, maxOp<point>());
}


bool Foam::boundBoxNew::intersect(const boundBoxNew& bb)
{
    min_ = ::Foam::max(min_, bb.min_);
    max_ = ::Foam::min(max_, bb.max_);

    return valid();
}


bool Foam::boundBoxNew::intersects(const plane& pln) const
{
    // Require a full 3D box
    if (nDim() != 3)
    {
        return false;
    }

    bool above = false;
    bool below = false;

    tmp<pointField> tpts(points());
    const auto& pts = tpts();

    for (const point& p : pts)
    {
        if (pln.sideOfPlane(p) == plane::FRONT)
        {
            above = true;
        }
        else
        {
            below = true;
        }
    }

    return (above && below);
}


bool Foam::boundBoxNew::contains(const UList<point>& points) const
{
    if (points.empty())
    {
        return true;
    }

    for (const point& p : points)
    {
        if (!contains(p))
        {
            return false;
        }
    }

    return true;
}


bool Foam::boundBoxNew::containsAny(const UList<point>& points) const
{
    if (points.empty())
    {
        return true;
    }

    for (const point& p : points)
    {
        if (contains(p))
        {
            return true;
        }
    }

    return false;
}


Foam::point Foam::boundBoxNew::nearest(const point& pt) const
{
    // Clip the point to the range of the bounding box
    const scalar surfPtx = Foam::max(Foam::min(pt.x(), max_.x()), min_.x());
    const scalar surfPty = Foam::max(Foam::min(pt.y(), max_.y()), min_.y());
    const scalar surfPtz = Foam::max(Foam::min(pt.z(), max_.z()), min_.z());

    return point(surfPtx, surfPty, surfPtz);
}

void Foam::boundBoxNew::resetBoundBoxNew
(
    const vector& startPoint,
    const vector& endPoint
)
{
    if(startPoint.x() < endPoint.x())
    {
        min_.x() = startPoint.x();
        max_.x() = endPoint.x();
    }
    else
    {
        min_.x() = endPoint.x();
        max_.x() = startPoint.x();
    }
    if(startPoint.y() < endPoint.y())
    {
        min_.y() = startPoint.y();
        max_.y() = endPoint.y();
    }
    else
    {
        min_.y() = endPoint.y();
        max_.y() = startPoint.y();
    }
    if(startPoint.z() < endPoint.z())
    {
        min_.z() = startPoint.z();
        max_.z() = endPoint.z();
    }
    else
    {
        min_.z() = endPoint.z();
        max_.z() = startPoint.z();
    }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const boundBoxNew& bb)
{
    if (os.format() == IOstream::ASCII)
    {
        os << bb.min_ << token::SPACE << bb.max_;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&bb.min_),
            sizeof(boundBoxNew)
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, boundBoxNew& bb)
{
    if (is.format() == IOstream::ASCII)
    {
        is >> bb.min_ >> bb.max_;
    }
    else
    {
        Detail::readContiguous<boundBoxNew>
        (
            is,
            reinterpret_cast<char*>(&bb.min_),
            sizeof(boundBoxNew)
        );
    }

    is.check(FUNCTION_NAME);
    return is;
}


// ************************************************************************* //

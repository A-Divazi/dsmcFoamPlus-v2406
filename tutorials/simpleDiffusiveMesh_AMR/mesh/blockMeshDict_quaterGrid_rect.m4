/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.0.1                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    `format'      ascii;

    root        "";
    case        "";
    instance    "";
    local       "";

    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// General m4 macros

changecom(//)changequote([,]) dnl>
define(calc, [esyscmd(perl -e 'use Math::Trig; print ($1)')]) dnl>
define(floor, [esyscmd(perl -e "printf ('%d', $1)")])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])

// 3D figures
define(slice, ($1 $4 $4 $1 $2 $5 $6 $3))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// User-defined parameters

scale 0.001;// 1:1

// General parameters
define(cosine,  calc(cos(pi/2.0)))			// 90. degrees at each side of the x-axis
define(sine,    calc(sin(pi/2.0)))
define(cosineH,  calc(cos(pi/4.0)))
define(sineH,    calc(sin(pi/4.0)))


//refine (dn):
define(dx, 1) // 1/10 mm cellsizes 

// Port dimensions
define(x_stag, 0)   // x-dimension stagnation
define(x_us, 12) //20     // x-dimension upstream
define(x_br, 1)    // x-dimension grid
define(x_ds, 12) //30    // x-dimension downstream

define(y_stag, 0)   // y-dimension stagnation
define(y_grd, 14)   // y-dimension grid
define(y_br, 1)     // y-dimension bracket
define(y_ts, 3)    // y-dimension topstream

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Derived parameters


// X-Direction:
define(x0, calc(x_stag-x_us))
define(x1, calc(x_stag + x_br + x_ds))


// Y-Direction:
define(y0, y_stag)
define(y1, calc(y_stag + y_grd + y_br + y_ts))


// Z-Direction: only for 90 degree quater
define(z0, -y0)
define(z1, -y1)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Parametric description

vertices
(
    // FRONT:
    (x0 y0 z0)  vlabel(FNT0)
    (x1 y0 z0)  vlabel(FNT1)
    (x1 y1 z0)  vlabel(FNT2)
    (x0 y1 z0)  vlabel(FNT3)


    // BACK:
    (x0 y0 z1)  vlabel(BCK0)
    (x1 y0 z1)  vlabel(BCK1)
    (x1 y1 z1)  vlabel(BCK2)
    (x0 y1 z1)  vlabel(BCK3)

);


blocks
(
    hex (FNT0 FNT3 FNT2 FNT1 BCK0 BCK3 BCK2 BCK1) dsmcZone1 (floor(calc((y1-(y0))/dx)) floor(calc((x1-(x0))/dx)) floor(calc((y1-(y0))/dx))) simpleGrading (1 1 1)
);


edges
(
);


patches
(    
    symmetry
        front
        (
            (FNT0 FNT1 FNT2 FNT3)
        )
    symmetry
        bottom
        (
            (BCK0 BCK1 FNT1 FNT0)
        )

    patch
        flow
        (
            (BCK0 FNT0 FNT3 BCK3)
            (BCK1 FNT1 FNT2 BCK2)
            (BCK3 FNT3 FNT2 BCK2)
            (BCK1 BCK0 BCK3 BCK2)
        )
);

boundary
(
);
mergePatchPairs
(
);

// ************************************************************************* //

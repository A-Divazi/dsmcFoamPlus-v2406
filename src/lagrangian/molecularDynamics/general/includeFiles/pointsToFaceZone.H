void pointsToFaceZone
(
    polyMesh& mesh,
    const labelList& pointList,
    const word& name
)
{
    DynamicList<label> faceList(0);

    forAll(pointList, p)
    {
        const labelList& faces = mesh.pointFaces()[pointList[p]];

        forAll(faces, f)
        {
            //if(findIndex(faceList, faces[f]) == -1)// is deprecated: Since 2017-10; use "UList find/found methods"
            if(faceList.find(faces[f]) == -1)
            {
                faceList.append(faces[f]);
            }
        }
    }



    createFaceZone
    (
        mesh,
        //faceList.shrink(),
		faceList,
        name
    );
}

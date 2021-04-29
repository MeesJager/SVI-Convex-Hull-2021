
/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

/*  A User Guide for this program can be found in the README file
    that has been added to the folder containing the source files   */

#include "imports.c"

// We execute this to end the program early.
int ExitProgram()
{
    printf("Ending program...\n");
    printf("-----------------------\n\n\n");
    return EXIT_FAILURE;
}

// This is the starting point of the program.
int main(int argc, char **argv)
{
    /*  PROGRAM START  */

    // Introduction
    printf("\n\n-----------------------\n");
    printf("Convex Hull 2021 -- SVI\n");
    printf("...\n");
    fflush(stdout);



    /*  READ INPUT  */

    // We read the necessary files from the given arguments
    // or ask the user.
    char *fileNames[8];
    int doComparison = 0;
    for (size_t i = 0; i < 8; i++)
        fileNames[i] = malloc(200 * sizeof(char));

    GetFileNames(argc, argv, fileNames, &doComparison);
    
    // Now we try to read all necessary data from files.
    // First we load the input size for the convex hull.
    printf("Loading convex hull data from file...\n");
    int trianglesCount;
    if (!ReadInputSize(fileNames[1], &trianglesCount))
        return ExitProgram();

    // Now we load the actual triangles of the convex hull.
    triangle *triangles = malloc(trianglesCount * sizeof(triangle));
    if (!ReadRawData(fileNames[0], fileNames[1], triangles, &trianglesCount))
        return ExitProgram();

    // If selected by the user, we load the data of the convex hull to
    // which we compare the original.
    int trianglesCountCompHull = 0;
    if (doComparison)
    {
        printf("Loading comparison convex hull data from file...\n");
        if (!ReadInputSize(fileNames[3], &trianglesCountCompHull))
            return ExitProgram();
    }

    triangle *trianglesCompHull = malloc(trianglesCountCompHull * sizeof(triangle));
    if (doComparison)
        if (!ReadRawData(fileNames[2], fileNames[3], trianglesCompHull, &trianglesCountCompHull))
            return ExitProgram();

    printf("Succesfully read all data from files...\n");
    printf("...\n");
    fflush(stdout);



    /*  BASIC PROPERTIES OF CONVEX HULL  */

    // We produce the dimensions of the bounding box of all triangles.
    // The boundingbox is snapped to the grid.
    // Ordering: (xMin, yMin, zMin), (xMax, yMax, zMax).
    point boundingBox[2] = {{0, 0, 0}, {0, 0, 0}};
    CalculateBoundingBox(triangles, trianglesCount, boundingBox);
    BoundingBoxToGrid(boundingBox);

    // Integrals that are calculated/approximated for use in multiple 
    // algorithms for calculating basic convex hull properties.
    // Order: 1, x, y, z, xˆ2, yˆ2, zˆ2, xy, yz, zx.
    long double integrals[11] = {0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0};
    // We calculate the normal values for each triangle as well, so will
    // be available for reuse. The normal at index i in this array
    // corresponds to the triangle at index i in the triangles array.
    point* normals = malloc(trianglesCount * sizeof(point));

    printf("Calculating basic properties of the convex hull...\n");
    fflush(stdout);
    // The calculations for the basic properties of the convex hull.
    CalculateIntegrals(triangles, trianglesCount, integrals, normals);

    // We create all variables containing the properties of the convex hull.
    long double volume = RetrieveVolume(integrals);
    long double surfaceArea = RetrieveSurfaceArea(integrals);
    point centerOfMass = RetrieveCenterOfMass(integrals);
    long double sphericity = CalculateSphericity(integrals);

    long double momentsOfInertia[3][3];
    long double axesLengths[] = {0, 0, 0};
    point axes[3];
    long double eigenvalues[3];
    CalculateMomentsOfInertia(integrals, momentsOfInertia);
    CalculatePrincipalAxes(momentsOfInertia, axes, eigenvalues);
    CalculateDimensionsAlongPrincipalAxes(triangles, trianglesCount, normals, axesLengths, axes, centerOfMass);



    /*  CONVEX HULL COMPARISON  */

    // We calculate the data for our comparison.
    long double similarity = 0;
    long double shortestDistance = 0;
    if (doComparison)
    {
        printf("Calculating comparison to a second convex hull...\n");
        
        long double integralsCompHull[11] = {0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0};
        point* normalsCompHull = malloc(trianglesCountCompHull * sizeof(point));
        CalculateIntegrals(trianglesCompHull, trianglesCountCompHull, integralsCompHull, normalsCompHull);
        long double volumeCompHull = RetrieveVolume(integralsCompHull);
        long double surfaceAreaCompHull = RetrieveSurfaceArea(integralsCompHull);
        point centerOfMassCompHull = RetrieveCenterOfMass(integralsCompHull);
  
        long double momentsOfInertiaCompHull[3][3];
        long double axesLengthsCompHull[] = {0, 0, 0};
        point axesCompHull[3];
        long double eigenValuesCompHull[3];
        CalculateMomentsOfInertia(integralsCompHull, momentsOfInertiaCompHull);
        CalculatePrincipalAxes(momentsOfInertiaCompHull, axesCompHull, eigenValuesCompHull);
        CalculateDimensionsAlongPrincipalAxes(trianglesCompHull, trianglesCountCompHull, normalsCompHull, axesLengthsCompHull, axesCompHull, centerOfMassCompHull);

        similarity = SimilarityTest(volume, volumeCompHull, surfaceArea, surfaceAreaCompHull, axesLengths, axesLengthsCompHull);
        shortestDistance = GJK(triangles, trianglesCompHull, trianglesCount, trianglesCountCompHull, centerOfMass, centerOfMassCompHull);
        
        free(normalsCompHull);
        free(trianglesCompHull);
    }



    /*  VOXELISATION  */

    // We translate the convex hull to the positive octant in space,
    // so that it can be used by our volume-voxelisation algorithm.
    point boundingBoxLocation; 
    boundingBoxLocation.x = boundingBox[0].x;
    boundingBoxLocation.y = boundingBox[0].y;
    boundingBoxLocation.z = boundingBox[0].z;
    point translation = scalar_point(&boundingBox[0], -1);
    centerOfMass = add_point(&centerOfMass, &translation);
    TranslateConvexHull(triangles, trianglesCount, boundingBox, translation);

    // We prepare the datastructures that will contain both voxelisations.
    int voxelCount = boundingBox[1].y * boundingBox[1].z * (boundingBox[1].x/8);
    unsigned char surfaceVoxels[voxelCount];
    unsigned char volumeVoxels[voxelCount];
    memset(surfaceVoxels, 0, voxelCount);
    memset(volumeVoxels, 0, voxelCount);

    // We produce both types of voxelisation of our object.
    printf("Building surface voxelisation of the convex hull...\n");
    fflush(stdout);
    SurfaceVoxelisation(triangles, trianglesCount, surfaceVoxels, boundingBox[1]);

    printf("Building volume voxelisation of the convex hull...\n\n");
    fflush(stdout);
    VolumeVoxelisation(triangles, trianglesCount, volumeVoxels, boundingBox[1]);


    /*  PROGRAM END  */

    // Write all output to file and the screen.
    WriteOutput(fileNames, surfaceVoxels, volumeVoxels, voxelCount, volume, surfaceArea, centerOfMass, sphericity, 
                momentsOfInertia, eigenvalues, axes, axesLengths, doComparison, similarity, shortestDistance);

    // We end the program.
    for (size_t i = argc - 1; i < 8; i++)
        free(fileNames[i]);
    
    free(normals);
    free(triangles);

    printf("Ending program...\n");
    printf("-----------------------\n\n\n");
    return EXIT_SUCCESS;
}



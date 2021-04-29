
/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

// We print all calculated output to files here, and the basic properties and comparison data
// is also printed on screen.
int WriteOutput(char **fileNames, unsigned char *surfaceVoxels, unsigned char *volumeVoxels, 
                int voxelCount, long double volume, long double surfaceArea, point centerOfMass, 
                long double sphericity, long double momentsOfInertia[3][3], long double eigenvalues[3], point *axes,
                long double *axesLength, int doComparison, long double similarity, long double shortestDistance)
{
    WritePropertiesOutput(fileNames[4], volume, surfaceArea, centerOfMass, sphericity, momentsOfInertia, eigenvalues, axes, axesLength);
    WriteVoxelisationOutput(fileNames[5], surfaceVoxels, voxelCount);
    WriteVoxelisationOutput(fileNames[6], volumeVoxels, voxelCount);
    if (doComparison)
        WriteComparisonOutput(fileNames[7], similarity, shortestDistance);

    return 1;
}

// All data corresponding to the basic properties of a convex hull are written to a text file.
int WritePropertiesOutput(char *fileName, long double volume, long double surfaceArea, 
                          point centerOfMass, long double sphericity, long double momentsOfInertia[3][3], 
                          long double eigenvalues[3], point *axes, long double *axesLength)
{
    FILE *file = fopen(fileName, "w");

    // Printing the results
    fprintf(file, "Properties of the convex hull:\n\n");
    printf("Properties of the convex hull:\n\n");
    fprintf(file, "Volume: %Lf\n---\n", volume);
    printf("Volume: %Lf\n---\n", volume);
    fprintf(file, "Surface Area: %Lf\n---\n", surfaceArea);
    printf("Surface Area: %Lf\n---\n", surfaceArea);
    fprintf(file, "Center Of Mass: (%Lf, %Lf, %Lf)\n---\n", centerOfMass.x, centerOfMass.y, centerOfMass.z);
    printf("Center Of Mass: (%Lf, %Lf, %Lf)\n---\n", centerOfMass.x, centerOfMass.y, centerOfMass.z);
    fprintf(file, "Sphericity: %Lf\n---\n", sphericity);
    printf("Sphericity: %Lf\n---\n", sphericity);
    fprintf(file, "Principal Moments of Inertia:\n| %Lf, %Lf, %Lf |\n---\n", eigenvalues[0], eigenvalues[1], eigenvalues[2]);
    printf("Principal Moments of Inertia:\n| %Lf, %Lf, %Lf |\n---\n", eigenvalues[0], eigenvalues[1], eigenvalues[2]);
    fprintf(file, "Moments-of-Inertia matrix:\n| %Lf, %Lf, %Lf |\n| %Lf, %Lf, %Lf |\n| %Lf, %Lf, %Lf |\n---\n", 
            momentsOfInertia[0][0], momentsOfInertia[0][1], momentsOfInertia[0][2],
            momentsOfInertia[1][0], momentsOfInertia[1][1], momentsOfInertia[1][2],
            momentsOfInertia[2][0], momentsOfInertia[2][1], momentsOfInertia[2][2]);
    printf("Moments-of-Inertia matrix:\n| %Lf, %Lf, %Lf |\n| %Lf, %Lf, %Lf |\n| %Lf, %Lf, %Lf |\n---\n", 
            momentsOfInertia[0][0], momentsOfInertia[0][1], momentsOfInertia[0][2],
            momentsOfInertia[1][0], momentsOfInertia[1][1], momentsOfInertia[1][2],
            momentsOfInertia[2][0], momentsOfInertia[2][1], momentsOfInertia[2][2]);
    fprintf(file, "Principal axes:\n(%Lf, %Lf, %Lf)\n(%Lf, %Lf, %Lf)\n(%Lf, %Lf, %Lf)\n---\n",
            axes[0].x, axes[0].y, axes[0].z, axes[1].x, axes[1].y, axes[1].z, axes[2].x, axes[2].y, axes[2].z);
    printf("Principal axes:\n(%Lf, %Lf, %Lf)\n(%Lf, %Lf, %Lf)\n(%Lf, %Lf, %Lf)\n---\n",
            axes[0].x, axes[0].y, axes[0].z, axes[1].x, axes[1].y, axes[1].z, axes[2].x, axes[2].y, axes[2].z);
    fprintf(file, "Dimensions along the principal axes:\n(%Lf, %Lf, %Lf)\n---\n\n", axesLength[0], axesLength[1], axesLength[2]);
    printf("Dimensions along the principal axes:\n(%Lf, %Lf, %Lf)\n---\n\n", axesLength[0], axesLength[1], axesLength[2]);


    fclose(file);
    return 1;
}

// We write the output of a voxelisation to a binary file.
int WriteVoxelisationOutput(char *fileName, unsigned char *voxels, int voxelCount)
{
    FILE *file = fopen(fileName, "wb");

    fwrite(voxels, sizeof(unsigned char), voxelCount, file);

    fclose(file);
    return 1;
}

// We write the output of convex hull comparisons to a text file.
int WriteComparisonOutput(char *fileName, long double similarity, long double shortestDistance)
{
    FILE *file = fopen(fileName, "w");

    fprintf(file, "Comparison with the second convex hull:\n\n");
    printf("Comparison with the second convex hull:\n\n");
    fprintf(file, "Similarity: %Lf\n---\n", similarity);
    printf("Similarity: %Lf\n---\n", similarity);
    fprintf(file, "ShortestDistance: %Lf\n---\n", shortestDistance);
    printf("ShortestDistance: %Lf\n---\n", shortestDistance);

    fclose(file);
    return 1;
}
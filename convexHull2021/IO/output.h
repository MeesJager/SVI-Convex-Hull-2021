
/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

int WriteOutput(char **fileNames, unsigned char *surfaceVoxels, unsigned char *volumeVoxels, 
                int voxelCount, long double volume, long double surfaceArea, point centerOfMass, 
                long double sphericity, long double momentsOfInertia[3][3], long double eigenvalues[3], point *axes,
                long double *axesLength, int doComparison, long double similarity, long double shortestDistance);

int WritePropertiesOutput(char *fileName, long double volume, long double surfaceArea, 
                          point centerOfMass, long double sphericity, long double momentsOfInertia[3][3], 
                          long double eigenvalues[3], point *axes, long double *axesLength);
                          
int WriteVoxelisationOutput(char *fileName, unsigned char *voxels, int voxelCount);

int WriteComparisonOutput(char *fileName, long double similarity, long double shortestDistance);
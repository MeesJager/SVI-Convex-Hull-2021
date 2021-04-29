
/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

int readBlock(void *buffer, size_t size, size_t count, FILE *file);

int GetFileNames(int argc, char **argv, char **fileNames, int *doComparison);

int ReadInputSize(char *indicesFileName, int *trianglesCount);

int ReadRawData(char *pointsFileName, char *indicesFileName, triangle *triangles, int *trianglesCount);
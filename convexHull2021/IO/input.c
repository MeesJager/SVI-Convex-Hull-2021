
/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

int readBlock(void *buffer, size_t size, size_t count, FILE *file)
{
    return (fread(buffer, size, count, file) == count) ? 1 : 0;
}

// We get all necessary file names from the user, either from the given
// arguments or directly from the user.
int GetFileNames(int argc, char **argv, char **fileNames, int *doComparison)
{
    // The point and index input files of the convex hull
    if (argc < 2)
    {
        printf("Points file: ");
        scanf("%150s", fileNames[0]);
        printf("\n");
    }
    else 
    {
        for (size_t i = 0; i < 150; i++)
            fileNames[0][i] = argv[1][i];
    }

    if (argc < 3)
    {
        printf("Indices file: ");
        scanf("%150s", fileNames[1]);
        printf("\n");
    }
    else 
    {
        for (size_t i = 0; i < 150; i++)
            fileNames[1][i] = argv[2][i];
    }

    // The output file that will contain the basic properties of the convex hull.
    if (argc < 4)
    {
        printf("Convex hull properties output file: "); 
        scanf("%150s", fileNames[4]);
        printf("\n");
    }
    else fileNames[4] = argv[3];

    // The output file that will contain the surface voxelisation of the convex hull.
    if (argc < 5)
    {
        printf("Surface voxelisation output file: "); 
        scanf("%150s", fileNames[5]);
        printf("\n");
    }
    else 
    {
        for (size_t i = 0; i < 150; i++)
            fileNames[5][i] = argv[4][i];
    }

    // The output file that will contain the volume voxelisation of the convex hull.
    if (argc < 6)
    {
        printf("Volume voxelisation output file: "); 
        scanf("%150s", fileNames[6]);
        printf("\n");
    }
    else 
    {
        for (size_t i = 0; i < 150; i++)
            fileNames[6][i] = argv[5][i];
    }

    // We ask the user if they want to do a comparison with a second convex hull,
    // and ask the necessary file names if so.
    char answer;
    while (answer != 'Y' && answer != 'N')
    {
        printf("Do you want to get a comparison with a second convex hull? [Y/N]: ");
        scanf("%s", &answer);
        printf("\n");

        if (answer == 'Y')
        {
            *doComparison = 1;
            // The input files for the comparison hull.
            if (argc < 7)
            {
                printf("Points file: ");
                scanf("%150s", fileNames[2]);
                printf("\n");
            }
            else 
            {
                for (size_t i = 0; i < 150; i++)
                    fileNames[2][i] = argv[6][i];
            }

            if (argc < 8)
            {
                printf("Indices file: ");
                scanf("%150s", fileNames[3]);
                printf("\n");
            }
            else 
            {
                for (size_t i = 0; i < 150; i++)
                    fileNames[3][i] = argv[7][i];
            }

            // The output file that will contain the comparison data.
            if (argc < 9)
            {
                printf("Comparison output file: ");
                scanf("%150s", fileNames[7]);
                printf("\n");
            }
            else 
            {
                for (size_t i = 0; i < 150; i++)
                    fileNames[7][i] = argv[8][i];
            }
        }
    }

    return 1;
}

// We read the input size, which corresponds to the number of triangles in the convex hull.
int ReadInputSize(char *indicesFileName, int *trianglesCount)
{
    // Open the file
    FILE *indexFile = fopen(indicesFileName, "rb");
    if(indexFile == NULL)
    {
        printf("The index file does not exist \"%s\".\n", indicesFileName);
        return 0;
    }
    // Read the input size.
    if (!readBlock(trianglesCount, sizeof(int), 1, indexFile))
    {
        printf("End of file reached.\n\n");
        fclose(indexFile);
        return 0;
    }

    printf("Succesfully read input size from file...\n");
    fclose(indexFile);
    return 1;
}

// We read from the given points file all points in the convex hull,
// and from the index file produced by the 2019 SVI project, which
// point are used in triangles, and we produce these triangles.
int ReadRawData(char *pointsFileName, char *indicesFileName, triangle *triangles, int *trianglesCount)
{
    // Open the points file.
    FILE *pointsFile = fopen(pointsFileName, "rb");
    if(pointsFile == NULL)
    {
        printf("The points file does not exist \"%s\".\n", pointsFileName);
        return 0;
    }

    // Read amount of points.
    int pointsCount = 0;
    if (!readBlock(&pointsCount, sizeof(int), 1, pointsFile))
    {
        printf("End of file reached.\n\n");
        fclose(pointsFile);
        return 0;
    }

    // Load all points into the program.
    point *points = malloc(pointsCount * sizeof(point));
    long double coordinates[3];
    for (int i = 0; i < pointsCount; i++)
    {
        if (!readBlock(&coordinates, sizeof(long double), 3, pointsFile))
        {
            free(points);
            printf("End of file reached.\n\n");
            fclose(pointsFile);
            return 0;
        }
        points[i].x = coordinates[0];
        points[i].y = coordinates[1];
        points[i].z = coordinates[2];
    }
    fclose(pointsFile);

    // Open the index file.
    FILE *indexFile = fopen(indicesFileName, "rb");
    if(indexFile == NULL)
    {
        printf("The index file does not exist \"%s\".\n", indicesFileName);
        return 0;
    }

    // Read amount of triangles.
    // This has already done before, but we need to step over
    // this integer.
    if (!readBlock(trianglesCount, sizeof(int), 1, indexFile))
    {
        printf("End of file reached.\n\n");
        fclose(indexFile);
        return 0;
    }

    // Load all indices into a single array, where each triple
    // of indices represents a triangle.
    int *indices = malloc(*trianglesCount * 3 * sizeof(int));
    int indexTriple[3];
    for (int i = 0; i < *trianglesCount; i++)
    {
        if (!readBlock(&indexTriple, sizeof(int), 3, indexFile))
        {
            free(indices);
            printf("End of file reached.\n\n");
            fclose(indexFile);
            return 0;
        }

        indices[3*i] = indexTriple[0];
        indices[3*i + 1] = indexTriple[1];
        indices[3*i + 2] = indexTriple[2];
    }
    printf("Succesfully read convex hull data from file...\n");
    fclose(indexFile);

    // Create all triangles.
    for (size_t i = 0; i < *trianglesCount; i++)
    {
        triangles[i].v0 = points[indices[3*i]];
        triangles[i].v1 = points[indices[3*i+1]];
        triangles[i].v2 = points[indices[3*i+2]];
    }

    return 1;
}
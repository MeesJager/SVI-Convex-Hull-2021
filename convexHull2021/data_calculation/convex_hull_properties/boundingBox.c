
/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

// We calculate a bounding box of a list of triangles that is axis-aligned.
int CalculateBoundingBox(triangle *triangles, int trianglesCount, point boundingBox[2])
{
    point min = {0, 0, 0};
    point max = {0, 0, 0};

    // Note that we do not parallelise this part to prevent errors
    // in finding the minimums and maximums.
    for (size_t i = 0; i < trianglesCount; i++)
    {
       min.x = fminl(min.x, triangles[i].v0.x);
       min.x = fminl(min.x, triangles[i].v1.x);
       min.x = fminl(min.x, triangles[i].v2.x);

       min.y = fminl(min.y, triangles[i].v0.y);
       min.y = fminl(min.y, triangles[i].v1.y);
       min.y = fminl(min.y, triangles[i].v2.y);

       min.z = fminl(min.z, triangles[i].v0.z);
       min.z = fminl(min.z, triangles[i].v1.z);
       min.z = fminl(min.z, triangles[i].v2.z);


       max.x = fmaxl(max.x, triangles[i].v0.x);
       max.x = fmaxl(max.x, triangles[i].v1.x);
       max.x = fmaxl(max.x, triangles[i].v2.x);

       max.y = fmaxl(max.y, triangles[i].v0.y);
       max.y = fmaxl(max.y, triangles[i].v1.y);
       max.y = fmaxl(max.y, triangles[i].v2.y);

       max.z = fmaxl(max.z, triangles[i].v0.z);
       max.z = fmaxl(max.z, triangles[i].v1.z);
       max.z = fmaxl(max.z, triangles[i].v2.z);
    }

    boundingBox[0].x = min.x;
    boundingBox[0].y = min.y;
    boundingBox[0].z = min.z;
    boundingBox[1].x = max.x;
    boundingBox[1].y = max.y;
    boundingBox[1].z = max.z;

    return 1;
}

// We use this to snap the bounding box of an object to the grid.
int BoundingBoxToGrid (point boundingBox[2])
{
    boundingBox[0].x = floorl(boundingBox[0].x);
    boundingBox[0].y = floorl(boundingBox[0].y);
    boundingBox[0].z = floorl(boundingBox[0].z);
    boundingBox[1].x = ceill(boundingBox[1].x);
    boundingBox[1].y = ceill(boundingBox[1].y);
    boundingBox[1].z = ceill(boundingBox[1].z);

    return 1;
}

// We use this to translate a complete convex hull by some given translation.
int TranslateConvexHull(triangle *triangles, int trianglesCount, point boundingBox[2], point translation)
{
    #pragma omp parallel for
    for (size_t i = 0; i < trianglesCount; i++)
    {
        triangles[i].v0 = add_point(&triangles[i].v0, &translation);
        triangles[i].v1 = add_point(&triangles[i].v1, &translation);
        triangles[i].v2 = add_point(&triangles[i].v2, &translation);
    }
    
    boundingBox[0] = add_point(&boundingBox[0], &translation);
    boundingBox[1] = add_point(&boundingBox[1], &translation);

    return 1;
}

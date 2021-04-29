
/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

// We perform a test for each triangle to check if it intersects with a certain voxel,
// and we turn on that voxel if that is the case.
int SurfaceVoxelisation (triangle *triangles, int number_of_triangles, unsigned char *voxels, point spaceBound)
{
    point voxelSize = {VOXELSIZE, VOXELSIZE, VOXELSIZE};
    
    #pragma omp parallel for
    for (size_t t = 0; t < number_of_triangles; t++)
        TriangleBoxIntersection(triangles[t], voxels, &voxelSize, spaceBound);

    return 1;
}

// Sets a bit from zero to one in the array of characters that represents the voxels.
int TurnOnVoxel (point voxel, unsigned char *voxels, point spaceBound)
{
    int index = (voxel.z * (spaceBound.x/8) * spaceBound.y) + (voxel.y * (spaceBound.x/8)) + (voxel.x/8);
    int subindex = fmod(voxel.x, 8);
    unsigned char mask = 1 << subindex | 0;
    #pragma omp critical
    {
        voxels[index] = (voxels[index] | mask);
    }
    return 1;
}


// We calculate the necessary data for a intersection test of a triangle and a voxel,
// and we return the result of the test.
int TriangleBoxIntersection(triangle triangle, unsigned char *voxels, point *boxSize, point spaceBound)
{
    int counter = 0;
    point v0 = triangle.v0;
    point v1 = triangle.v1;
    point v2 = triangle.v2;

    // Points representing the three edges of the triangle
    point e0 = subtract_point(&v1, &v0);
    point e1 = subtract_point(&v2, &v1);
    point e2 = subtract_point(&v0, &v2);

    // The normal vector of the triangle
    point n = cross_product(&e0, &e1);

    // Critical point c of this triangle
    point c;
    c.x = (n.x > 0) ? boxSize->x : 0;
    c.y = (n.y > 0) ? boxSize->y : 0;
    c.z = (n.z > 0) ? boxSize->z : 0;

    // Helper values d1 and d2
    point c_minus_v0 = subtract_point(&c, &v0);
    long double d1 = inner_product(&n, &c_minus_v0);
    point boxsize_minus_c = subtract_point(boxSize, &c);
    point boxsize_minus_c_minus_v0 = subtract_point(&boxsize_minus_c, &v0);
    long double d2 = inner_product(&n, &boxsize_minus_c_minus_v0);

    // Projection of the normal vector onto different planes
    point b;
    b.x = (n.x > 0) ? 1 : -1;
    b.y = (n.y > 0) ? 1 : -1;
    b.z = (n.z > 0) ? 1 : -1;

    point nxy[3] = {
        {-e0.y * b.z, e0.x * b.z, 0},
        {-e1.y * b.z, e1.x * b.z, 0},
        {-e2.y * b.z, e2.x * b.z, 0},
    };
    point nxz[3] = {
        {-e0.z * b.y, 0, e0.x * b.y},
        {-e1.z * b.y, 0, e1.x * b.y},
        {-e2.z * b.y, 0, e2.x * b.y},
    };
    point nyz[3] = {
        {0, -e0.z * b.x, e0.y * b.x},
        {0, -e1.z * b.x, e1.y * b.x},
        {0, -e2.z * b.x, e2.y * b.x},
    };

    /* 9 helper values */
    point dxy = {
        -inner_product(&nxy[0], &v0) + fmaxl(0, boxSize->x * nxy[0].x) + fmaxl(0, boxSize->y * nxy[0].y),
        -inner_product(&nxy[1], &v1) + fmaxl(0, boxSize->x * nxy[1].x) + fmaxl(0, boxSize->y * nxy[1].y),
        -inner_product(&nxy[2], &v2) + fmaxl(0, boxSize->x * nxy[2].x) + fmaxl(0, boxSize->y * nxy[2].y)
        };
    point dxz = {    
        -inner_product(&nxz[0], &v0) + fmaxl(0, boxSize->x * nxz[0].x) + fmaxl(0, boxSize->z * nxz[0].z),
        -inner_product(&nxz[1], &v1) + fmaxl(0, boxSize->x * nxz[1].x) + fmaxl(0, boxSize->z * nxz[1].z),
        -inner_product(&nxz[2], &v2) + fmaxl(0, boxSize->x * nxz[2].x) + fmaxl(0, boxSize->z * nxz[2].z)
        };
    point dyz = {    
        -inner_product(&nyz[0], &v0) + fmaxl(0, boxSize->y * nyz[0].y) + fmaxl(0, boxSize->z * nyz[0].z),
        -inner_product(&nyz[1], &v1) + fmaxl(0, boxSize->y * nyz[1].y) + fmaxl(0, boxSize->z * nyz[1].z),
        -inner_product(&nyz[2], &v2) + fmaxl(0, boxSize->y * nyz[2].y) + fmaxl(0, boxSize->z * nyz[2].z)
    };

    // We calculate the boundingbox of a triangle.
    point triangleBoundingBox[2] = {{0, 0, 0}, {0, 0, 0}};
    CalculateBoundingBox(&triangle, 1, triangleBoundingBox);
    BoundingBoxToGrid(triangleBoundingBox);

    // We only check the intersection with voxels that are in the bounding box of a triangle.
    point voxelSize = *boxSize;
    for (size_t x = triangleBoundingBox[0].x; x < triangleBoundingBox[1].x; x += voxelSize.x)
        for (size_t y = triangleBoundingBox[0].y; y < triangleBoundingBox[1].y; y += voxelSize.y)
            for (size_t z = triangleBoundingBox[0].z; z < triangleBoundingBox[1].z; z += voxelSize.z)
            {
                point boxlocation = {x, y, z};
                // We do the intersection test.

                if (!((inner_product(&n, &boxlocation) + d1) * (inner_product(&n, &boxlocation) + d2) <= 0))
                    continue;
                if (!((inner_product(&nxy[0], &boxlocation) + dxy.x >= 0) && (inner_product(&nxy[1], &boxlocation) + dxy.y >= 0) && (inner_product(&nxy[2], &boxlocation) + dxy.z >= 0)))
                    continue;
                if (!((inner_product(&nxz[0], &boxlocation) + dxz.x >= 0) && (inner_product(&nxz[1], &boxlocation) + dxz.y >= 0) && (inner_product(&nxz[2], &boxlocation) + dxz.z >= 0)))
                    continue;
                if (!((inner_product(&nyz[0], &boxlocation) + dyz.x >= 0) && (inner_product(&nyz[1], &boxlocation) + dyz.y >= 0) && (inner_product(&nyz[2], &boxlocation) + dyz.z >= 0)))
                    continue;
                counter++;

                TurnOnVoxel(boxlocation, voxels, spaceBound);
            }

    //printf("%d ", counter);
    //fflush(stdout);
    return 1;
}





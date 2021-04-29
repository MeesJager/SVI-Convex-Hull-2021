
/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

// Flip the value in the character array representing voxels, so
// turning a voxel on/off.
int FlipVoxel(point *voxel, unsigned char *voxels, point spaceBound)
{
    int index = voxel->z * spaceBound.y * spaceBound.x/8 + voxel->y * spaceBound.x/8 + voxel->x/8;
    int subindex = fmod(voxel->x, 8);
    unsigned char mask = 1 << subindex;
    #pragma omp critical
    {   
        voxels[index] = (voxels[index] ^ mask);
    }

    return 1;
}

// Check if a triangle is counterclockwise
int check_CCW(point2 *v0, point2 *v1, point2 *v2)
{
    point2 e0 = subtract_point2(v1, v0);
    point2 e1 = subtract_point2(v2, v0);
    double result = e0.x * e1.y - e0.y * e1.x;
    return (result > 0) ? 1 : 0;
}

// Lift 2D point to triangle, i.e. get its x coordinate
double get_x_coord(point *n, point *v0, point2 *p)
{
    return ((n->y * (v0->y - p->x) + n->z * (v0->z - p->y)) / n->x + v0->x);
}

// Check if an edge is a Top-Left edge.
// Since our triangle is counterclockwise, a left edge 'goes down'.
int TopLeftEdge(point2 *v0, point2 *v1)
{
    return ((v1->y < v0->y) || (v1->y == v0->y && v0->x > v1->x));
}

// Check if a point lies on edge 1,2 or 3 of triangle or inside
int check_point_triangle(point2 *v0, point2 *v1, point2 *v2, point2 *p)
{
    point2 PA = subtract_point2(p , v0);
    point2 PB = subtract_point2(p , v1);
    point2 PC = subtract_point2(p , v2);

    double t1 = PA.x * PB.y - PA.y * PB.x;
    if (fabsl(t1) < float_error && PA.x * PB.x <= 0 && PA.y * PB.y <= 0)
        return 1;

    double t2 = PB.x * PC.y - PB.y * PC.x;
    if (fabsl(t2) < float_error && PB.x * PC.x <= 0 && PB.y * PC.y <= 0)
        return 2;

    double t3 = PC.x * PA.y - PC.y * PA.x;
    if (fabsl(t3) < float_error && PC.x * PA.x <= 0 && PC.y * PA.y <= 0)
        return 3;

    return (t1 * t2 > 0 && t1 * t3 > 0) ? 0 : -1;
}

int VolumeVoxelisation(triangle *triangles, int number_of_triangles, unsigned char *voxels, point spaceBound)
{
    point voxelSize = {VOXELSIZE, VOXELSIZE, VOXELSIZE};

    #pragma omp parallel for
    for (size_t t = 0; t < number_of_triangles; t++) 
    {
        point v0 = triangles[t].v0;
        point v1 = triangles[t].v2;
        point v2 = triangles[t].v1;

        // Edge vectors
        point e0 = subtract_point(&v1, &v0);
        point e1 = subtract_point(&v2, &v1);
        point e2 = subtract_point(&v0, &v2);

        // normal vector
        point n = cross_product(&e0, &e1);

        // If normal vector has zero x-coord, then projecting to yz plane gives just a line segment
        if (fabsl(n.x) < float_error) {
            continue;
        }

        // Project points to yz plane
        point2 v0_yz = {v0.y, v0.z};
        point2 v1_yz = {v1.y, v1.z};
        point2 v2_yz = {v2.y, v2.z};

        if (!check_CCW(&v0_yz, &v1_yz, &v2_yz)) {
            point2 v3 = v1_yz;
            v1_yz = v2_yz;
            v2_yz = v3;
        }

        // Compute 2D bounding box of projection
        point2 a = max2(&v1_yz, &v2_yz);
        point2 bbox_max = max2(&v0_yz, &a);

        point2 b = min2(&v1_yz, &v2_yz);
        point2 bbox_min = min2(&v0_yz, &b);

        // compute the grid of (projections of) voxel centers contained in the bbox
        point2 bbox_max_grid = {
                floor(bbox_max.x / voxelSize.y - 0.5),
                floor(bbox_max.y / voxelSize.z - 0.5)
        };
        point2 bbox_min_grid = {
                ceil(bbox_min.x / voxelSize.y - 0.5),
                ceil(bbox_min.y / voxelSize.z - 0.5)
        };

        // We loop over every voxelcenter in the grid
        for (int y = bbox_min_grid.x; y <= bbox_max_grid.x; y++) {
            for (int z = bbox_min_grid.y; z <= bbox_max_grid.y; z++) {
                // compute coordinates of grid-point
                point2 voxel_center = {(y + 0.5) * voxelSize.y, (z + 0.5) * voxelSize.z};
                int checknum = check_point_triangle(&v0_yz, &v1_yz, &v2_yz, &voxel_center);

                // If the voxelcenter lies in the triangle or on a top or left edge,
                // find x-coord of intersection of the corresponding column of voxels with the triangle
                // and flip all voxels in the column between x=0 and the intersection point.
                if ((checknum == 1 && TopLeftEdge(&v0_yz, &v1_yz)) ||
                    (checknum == 2 && TopLeftEdge(&v1_yz, &v2_yz)) ||
                    (checknum == 3 && TopLeftEdge(&v2_yz, &v0_yz)) || (checknum == 0)) {

                    double intersection_x = get_x_coord(&n, &v0, &voxel_center);
                    int intersection_voxel_x = floor(intersection_x / voxelSize.x - 0.5);
                    for (int x = 0; x <= intersection_voxel_x; x++) {
                        point voxel = {x, y, z};
                        FlipVoxel(&voxel, voxels, spaceBound);
                    }
                }
            }
        }
    }

    return 1;
}

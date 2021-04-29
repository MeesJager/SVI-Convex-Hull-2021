
/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

int SurfaceVoxelisation(triangle *triangles, int number_of_triangles, unsigned char *voxels, point spaceBound);

int TurnOnVoxel(point voxel, unsigned char *voxels, point spaceBound);

int TriangleBoxIntersection(triangle triangle, unsigned char *voxels, point *boxsize, point spaceBound);
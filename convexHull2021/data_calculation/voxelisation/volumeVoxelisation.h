
/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

int FlipVoxel(point *voxel, unsigned char *voxels, point spaceBound);

int check_CCW(point2 *v0, point2 *v1, point2 *v2);

double get_x_coord(point *n, point *v0, point2 *p);

int TopLeftEdge(point2 *v0, point2 *v1);

int check_point_triangle(point2 *v0, point2 *v1, point2 *v2, point2 *p);

int VolumeVoxelisation(triangle *triangles, int number_of_triangles, unsigned char *voxels, point spaceBound);

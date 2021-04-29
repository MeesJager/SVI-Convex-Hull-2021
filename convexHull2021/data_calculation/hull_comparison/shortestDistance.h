
/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

int Sup(triangle *triangles, int trianglesCount, point vertex, point *sup);

int S3D(point tau[4], point dist[4], int *tauCounter, int *distCounter);

int S2D(point tau[4], point dist[4], int *tauCounter, int *distCounter);

int S1D(point tau[4], point dist[4], int *tauCounter, int *distCounter);

int Distance(point tau[4], point dist[4], int *tauCounter, int *distCounter);

long double GJK(triangle *triangles1, triangle *triangles2, int trianglesCount1, int trianglesCount2, point centerOfMass1, point centerOfMass2);

int CopyPointArray(point *copy, point *original, int *copyLength, int *originalLength);

int ErasePointArrayElement(point *points, int *counter, int index);
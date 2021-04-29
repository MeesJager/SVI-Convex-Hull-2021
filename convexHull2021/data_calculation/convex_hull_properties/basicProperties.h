
/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

int CalculateIntegrals(triangle *triangles, int trianglesCount, long double *integrals, point *normals);

int Subexpressions(long double w0, long double w1, long double w2, long double *fg);

long double RetrieveVolume(long double *integrals);

long double RetrieveSurfaceArea(long double *integrals);

point RetrieveCenterOfMass(long double *integrals);

long double CalculateSphericity(long double *integrals);

int CalculateMomentsOfInertia(long double* integrals, long double momentsOfInertia[3][3]);

int CalculatePrincipalAxes(long double momentsOfInertia[3][3], point *axes, long double* eigenvsalues);

int CalculateDimensionsAlongPrincipalAxes(triangle *triangles, int trianglesCount, point *normals, long double *axesLengths, point *axes, point centerOfMass);

/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

/* Based on the math.c file made for the convexHull2019 project for SVI */
/* Calculating the eigenvectors of a 3x3 matrix is adapted from the code
   given by https://gist.github.com/lh3/c280b2ac477e85c5c666*/

typedef struct point_tag {
    long double x, y, z;
} point;

typedef struct point2_tag {
    double x, y;
} point2;

typedef struct triangle_tag {
    point v0, v1, v2;
} triangle;

long double inner_product(point *v, point *w);
long double determinant_2d(long double x1, long double y1, long double x2, long double y2);
long double point_length(point *p);

point add_point(point *v, point *w);
point subtract_point(point *v, point *w);
point scalar_point(point *v, long double lambda);
point cross_product(point *v, point *w);

point2 subtract_point2(point2 *v, point2 *w);
point2 max2(point2 *v, point2 *w);
point2 min2(point2 *v, point2 *w);


double pythag(long double a, long double b);
long double Sign(long double a, long double b);
void tred2(long double a[3][3], long double d[3], long double e[3]);
void tqli(long double d[3], long double e[3], long double z[3][3]);
int EigenSymm(long double a[3][3], long double eval[3]);

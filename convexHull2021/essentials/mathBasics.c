
/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

/* Based on the math.c file made for the convexHull2019 project for SVI */
/* Calculating the eigenvectors of a 3x3 matrix is adapted from the code
   given by https://gist.github.com/lh3/c280b2ac477e85c5c666*/

long double inner_product(point *p, point *q) {return p->x * q->x + p->y * q->y + p->z * q->z;}
long double determinant_2d(long double x1, long double y1, long double x2, long double y2) {return x1 * y2 - x2 * y1;}
long double point_length(point *p) {long double r = inner_product(p, p); return sqrtl(r);}
int compare_points(point *p, point *q) {return ((p->x == q->x) && (p->y == q->y) && (p->z == q->z));}

point add_point(point *p, point *q) {point r = {p->x + q->x, p->y + q->y, p->z + q->z}; return r;}
point subtract_point(point *p, point *q) {point r = {p->x - q->x, p->y - q->y, p->z - q->z}; return r;}
point scalar_point(point *p, long double x) {point r = {x * p->x, x * p->y, x * p->z}; return r;}
point cross_product(point *p, point *q) {point r = {p->y * q->z - p->z * q->y, p->z * q->x - p->x * q->z, p->x * q->y - p->y * q->x}; return r;}

point2 subtract_point2(point2 *p, point2 *q) {point2 r = {p->x - q->x, p->y - q->y}; return r;}
point2 max2(point2 *v, point2 *w) 
{
    long double xMax = fmaxl(v->x, w->x);
    long double yMax = fmaxl(v->y, w->y);
    point2 maxPoint = {xMax, yMax};
    return maxPoint;
}

point2 min2(point2 *v, point2 *w) 
{
    long double xMin = fminl(v->x, w->x);
    long double yMin = fminl(v->y, w->y);
    point2 minPoint = {xMin, yMin};
    return minPoint;
}

// Methods to find the eigenvectors of a 3x3-matrix, adapted from
// https://gist.github.com/lh3/c280b2ac477e85c5c666.
double pythag(long double a, long double b)
{
	long double absa, absb;
	absa = fabsl(a);
	absb = fabsl(b);
	if (absa > absb) return absa * sqrt(1.0 + (absb / absa)*(absb / absa));
	else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + (absa / absb)*(absa / absb)));
}

long double Sign(long double a, long double b)
{
	return ((b) >= 0.0 ? fabsl(a) : -fabsl(a));
}

void tred2(long double a[3][3], long double d[3], long double e[3])
{
	int             l, k, j, i;
	long double     scale, hh, h, g, f;

	for (i = 2; i > 0; i--) {
		l = i - 1;
		h = scale = 0.0;
		if (l > 0) {
			for (k = 0; k < l + 1; k++)
				scale += fabsl(a[i][k]);
			if (scale == 0.0)
				e[i] = a[i][l];
			else {
				for (k = 0; k < l + 1; k++) {
					a[i][k] /= scale;
					h += a[i][k] * a[i][k];
				}
				f = a[i][l];
				g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i] = scale * g;
				h -= f * g;
				a[i][l] = f - g;
				f = 0.0;
				for (j = 0; j < l + 1; j++) {
					
					a[j][i] = a[i][j] / h;
					g = 0.0;
					for (k = 0; k < j + 1; k++)
						g += a[j][k] * a[i][k];
					for (k = j + 1; k < l + 1; k++)
						g += a[k][j] * a[i][k];
					e[j] = g / h;
					f += e[j] * a[i][j];
				}
				hh = f / (h + h);
				for (j = 0; j < l + 1; j++) {
					f = a[i][j];
					e[j] = g = e[j] - hh * f;
					for (k = 0; k < j + 1; k++)
						a[j][k] -= (f * e[k] + g * a[i][k]);
				}
			}
		} else
			e[i] = a[i][l];
		d[i] = h;
	}
	
	d[0] = 0.0;
	e[0] = 0.0;
	
	for (i = 0; i < 3; i++) {
		l = i;
		if (d[i] != 0.0) {
			for (j = 0; j < l; j++) {
				g = 0.0;
				for (k = 0; k < l; k++)
					g += a[i][k] * a[k][j];
				for (k = 0; k < l; k++)
					a[k][j] -= g * a[k][i];
			}
		}
		d[i] = a[i][i];
		a[i][i] = 1.0;
		for (j = 0; j < l; j++)
			a[j][i] = a[i][j] = 0.0;
	}
}

void tqli(long double d[3], long double e[3], long double z[3][3])
{
	int             m, l, iter, i, k;
	long double     s, r, p, g, f, dd, c, b;

	for (i = 1; i < 3; i++)
		e[i - 1] = e[i];
	e[2] = 0.0;
	for (l = 0; l < 3; l++) {
		iter = 0;
		do {
			for (m = l; m < 2; m++) {
				dd = fabsl(d[m]) + fabsl(d[m + 1]);
				if (fabsl(e[m]) + dd == dd)
					break;
			}
			if (m != l) {
				if (iter++ == 30) {
					fprintf(stderr, "[tqli] Too many iterations in tqli.\n");
					break;
				}
				g = (d[l + 1] - d[l]) / (2.0 * e[l]);
				r = pythag(g, 1.0);
				g = d[m] - d[l] + e[l] / (g + Sign(r, g));
				s = c = 1.0;
				p = 0.0;
				for (i = m - 1; i >= l; i--) {
					f = s * e[i];
					b = c * e[i];
					e[i + 1] = (r = pythag(f, g));
					if (r == 0.0) {
						d[i + 1] -= p;
						e[m] = 0.0;
						break;
					}
					s = f / r;
					c = g / r;
					g = d[i + 1] - p;
					r = (d[i] - g) * s + 2.0 * c * b;
					d[i + 1] = g + (p = s * r);
					g = c * r - b;
					
					for (k = 0; k < 3; k++) {
						f = z[k][i + 1];
						z[k][i + 1] = s * z[k][i] + c * f;
						z[k][i] = c * z[k][i] - s * f;
					}
				}
				if (r == 0.0 && i >= l)
					continue;
				d[l] -= p;
				e[l] = g;
				e[m] = 0.0;
			}
		} while (m != l);
	}
}

// This function is used to fill _a with its own eigenvectors as columns and
// fill eval with its eigenvalues.
int EigenSymm(long double a[3][3], long double eval[3])
{
	long double e[3];
	tred2(a, eval, e);
	tqli(eval, e, a);
	return 0;
}
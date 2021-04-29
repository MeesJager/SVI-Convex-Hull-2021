
/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

long double TrueSphericity(long double volume, long double surface);

long double KrumbeinSphericity(long double a, long double b, long double c);

long double ThirdSphericity(long double a, long double b, long double c);

long double TestCalculation(long double volume1, long double volume2, long double surface1, long double surface2, 
                              long double a1, long double a2, long double b1, long double b2, long double c1, long double c2);

long double SimilarityTest(long double volume1, long double volume2, long double surface1, long double surface2, 
                           long double *axesLength1, long double *axesLength2);

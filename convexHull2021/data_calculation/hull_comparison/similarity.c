
/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

// We can calculate the similarity between two convex hulls based on their volume,
// surface area and their dimensions along their principals axes.
long double SimilarityTest(long double volume1, long double volume2, long double surface1, long double surface2, 
                           long double *axesLength1, long double *axesLength2)
{
    // We first sort the values of the dimensions along the principal axes,
    // where a is the largest dimension and c is the smallest.
    long double a1, b1, c1, a2, b2, c2;
    long double temp1 = fmaxl(axesLength1[0], axesLength1[1]);
    long double temp2 = fminl(axesLength1[0], axesLength1[1]);
    if (axesLength1[2] > temp1)
    {
        a1 = axesLength1[2];
        b1 = temp1;
        c1 = temp2;
    }
    else if (axesLength1[2] > temp2)
    {
        a1 = temp1;
        b1 = axesLength1[2];
        c1 = temp2;
    }
    else
    {
        a1 = temp1;
        b1 = temp2;
        c1 = axesLength1[2];
    }

    long double temp3 = fmaxl(axesLength2[0], axesLength2[1]);
    long double temp4 = fminl(axesLength2[0], axesLength2[1]);
    if (axesLength2[2] > temp3)
    {
        a2 = axesLength2[2];
        b2 = temp3;
        c2 = temp4;
    }
    else if (axesLength2[2] > temp4)
    {
        a2 = temp3;
        b2 = axesLength2[2];
        c2 = temp4;
    }
    else
    {
        a2 = temp3;
        b2 = temp4;
        c2 = axesLength2[2];
    }

    // Then we perform the calculation.
    long double result = TestCalculation(volume1, volume2, surface1, surface2, a1, a2, b1, b2, c1, c2);
    return result; 
}

// We define three types of sphericities.
long double TrueSphericity(long double volume, long double surface)
{
    long double r = cbrt(6 * volume);
    return ((cbrt(M_PI) * cbrt(r * r))/ surface);
}

long double KrumbeinSphericity(long double a, long double b, long double c)
{
    return (cbrt((b*c)/(a*a)));
}

long double ThirdSphericity(long double a, long double b, long double c)
{
    return (cbrt((c*c)/(a * b)));
}

// The similarity between two convex hulls is determined by comparing
// them on their different types of sphericity.
long double TestCalculation(long double volume1, long double volume2, long double surface1, long double surface2, 
                              long double a1, long double a2, long double b1, long double b2, long double c1, long double c2)
{
    long double true1 = TrueSphericity(volume1, surface1);
    long double true2 = TrueSphericity(volume2, surface2);
    long double trueRatio = fminl(true1/true2 , true2/true1);

    long double krumbein1 = KrumbeinSphericity(a1,b1,c1);
    long double krumbein2 = KrumbeinSphericity(a2,b2,c2);
    long double krumbeinRatio = fminl(krumbein1/krumbein2, krumbein2/krumbein1);

    long double third1 = ThirdSphericity(a1,b1,c1);
    long double third2 = ThirdSphericity(a2,b2,c2);
    long double thirdRatio = fminl(third1/third2, third2/third1);

    long double min = fminl(trueRatio, krumbeinRatio);

    return (fminl(min, thirdRatio));
}



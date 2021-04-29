
/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

// We calculate 11 values that are returned in the integrals array,
// that are used to retrieve all basic properties of a convex hull.
int CalculateIntegrals(triangle *triangles, int trianglesCount, long double *integrals, point *normals)
{
    // These will be applied at the end to the values in integrals.
    long double multipliers[11] = {(1.0/6.0) , (1.0/24.0) , (1.0/24.0) , (1.0/24.0), (1.0/60.0) , (1.0/60.0) , (1.0/60.) , (1.0/120.0) , (1.0/120.0), (1.0/120.0), (1.0/2.0)};
       
    #pragma omp parallel for
    for (size_t i = 0; i < trianglesCount; i++)
    {
        // Get the vertices of the current triangle
        point v0 = triangles[i].v0;
        point v1 = triangles[i].v1;
        point v2 = triangles[i].v2;

        // Calculate points representing the edges of a triangle
        point e0 = subtract_point(&v1, &v0);
        point e1 = subtract_point(&v2, &v0);

        // Calculate the normal of a triangle and save it
        point d0 = cross_product(&e1, &e0);
        normals[i] = d0;

        // Calculate the area of a triangle
        long double magnitude = sqrtl(d0.x * d0.x + d0.y * d0.y + d0.z * d0.z);

        // Here we follow the exact calculations described in the paper.
        long double *fg = malloc(18 * sizeof(long double));
        Subexpressions(v0.x, v1.x, v2.x, fg);
        Subexpressions(v0.y, v1.y, v2.y, &fg[6]);
        Subexpressions(v0.z, v1.z, v2.z, &fg[12]);

        #pragma omp critical
        {
            integrals[0] += d0.x * fg[0];
            integrals[1] += d0.x * fg[1]; integrals[2] += d0.y * fg[7]; integrals[3] += d0.z * fg[13];
            integrals[4] += d0.x * fg[2]; integrals[5] += d0.y * fg[8]; integrals[6] += d0.z * fg[14];
            integrals[7] += d0.x * (v0.y * fg[3] + v1.y * fg[4] + v2.y * fg[5]);
            integrals[8] += d0.y * (v0.z * fg[9] + v1.z * fg[10] + v2.z * fg[11]);
            integrals[9] += d0.z * (v0.x * fg[15] + v1.x * fg[16] + v2.x * fg[17]);
            integrals[10] += magnitude;
        }

        free(fg);
    }

    // Apply the multipliers
    for (size_t i = 0; i < 11; i++)
        integrals[i] *= multipliers[i];

    return 1; 
}

// Specific calculations needed to calculate the integrals, as described in the paper.
int Subexpressions(long double w0, long double w1, long double w2, long double *fg)
{
    long double temp0 = w0 + w1;
    long double temp1 = w0 * w0;
    long double temp2 = temp1 + w1 * temp0;

    long double f0 = temp0 + w2;
    long double f1 = temp2 + w2 * f0;
    long double f2 = w0 * temp1 + w1 * temp2 + w2 * f1;

    long double g0 = f1 + w0 * (f0 + w0);
    long double g1 = f1 + w1 * (f0 + w1);
    long double g2 = f1 + w2 * (f0 + w2);

    fg[0] = f0; fg[1] = f1; fg[2] = f2; fg[3] = g0; fg[4] = g1; fg[5] = g2;

    return 1;
}

// We retrieve the volume of the convex hull from the calculated integrals.
// This is equal to the mass of the object, as we assume a uniform density of 1.
long double RetrieveVolume(long double *integrals)
{
    return integrals[0];
}

// We retrieve the surface area of the convex hull from the calculated integrals.
long double RetrieveSurfaceArea(long double *integrals)
{
    return integrals[10];
}

// We retrieve the center of mass of the convex hull from the calculated integrals.
point RetrieveCenterOfMass(long double *integrals)
{
    point centerOfMass;
    centerOfMass.x = integrals[1] / integrals[0];
    centerOfMass.y = integrals[2] / integrals[0];
    centerOfMass.z = integrals[3] / integrals[0];

    return centerOfMass;   
}

// We calculate the center of mass of the convex hull from the calculated integrals.
long double CalculateSphericity(long double *integrals)
{
    return (pow((double) M_PI, (double) 1/3) * pow(6 * integrals[0], (double) 2/3)) / integrals[10];   
}

// We fill a given 2-dimensional array with the values of a 3x3-matrix that
// represents the moments of inertia in row-major order.
int CalculateMomentsOfInertia(long double* integrals, long double momentsOfInertia[3][3])
{
    //hangt wel af van hoe we matrix willen hebben voor eigenvectoren vinden dus effe wachten met static weghalen
    point centerOfMass = RetrieveCenterOfMass(integrals); 
    
    momentsOfInertia[0][0] = integrals[5] + integrals[6] - integrals[0] * 
                             (centerOfMass.y * centerOfMass.y + centerOfMass.z * centerOfMass.z); //xx
    momentsOfInertia[0][1] = -(integrals[7] - integrals[0] * centerOfMass.x * centerOfMass.y); //xy
    momentsOfInertia[0][2] = -(integrals[9] - integrals[0] * centerOfMass.z * centerOfMass.x); //xz
    momentsOfInertia[1][0] = -(integrals[7] - integrals[0] * centerOfMass.x * centerOfMass.y); //yx
    momentsOfInertia[1][1] = integrals[4] + integrals[6] - integrals[0] * 
                             (centerOfMass.z * centerOfMass.z + centerOfMass.x * centerOfMass.x); //yy
    momentsOfInertia[1][2] = -(integrals[8] - integrals[0] * centerOfMass.y * centerOfMass.z); //yz
    momentsOfInertia[2][0] = -(integrals[9] - integrals[0] * centerOfMass.z * centerOfMass.x); //zx
    momentsOfInertia[2][1] = -(integrals[8] - integrals[0] * centerOfMass.y * centerOfMass.z); //zy
    momentsOfInertia[2][2] = integrals[4] + integrals[5] - integrals[0] * 
                             (centerOfMass.x * centerOfMass.x + centerOfMass.y * centerOfMass.y); //zz
    
    return 1;
}

// We fill the given array of points and long double to respectively contain 
// the eigenvectors and eigenvalues of the moments-of-inertia matrix.
int CalculatePrincipalAxes(long double momentsOfInertia[3][3], point axes[3], long double eigenvalues[3])
{
    // We make a copy of the momentsOfInertia-matrix that can be overwritten when calculating eigenvectors.
    long double copy[3][3];

    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 3; j++)
            copy[j][i] = momentsOfInertia[j][i];
    
    
    // We adapt the copy of the momentsOfInertia-matrix so that it contains eigenvectors
    // and we fill in the eigenvalues array. 
    EigenSymm(copy, eigenvalues);
    
    // We build the axes with the obtained eigenvectors in copy.
    for (size_t i = 0; i < 3; i++)
    {
        axes[i].x = copy[0][i];
        axes[i].y = copy[1][i];
        axes[i].z = copy[2][i];
    }

    return 1;
}

// We calculate the dimensions of the convex hull along each principle axis.
int CalculateDimensionsAlongPrincipalAxes(triangle *triangles, int trianglesCount, point *normals, long double *axesLengths, point *axes, point centerOfMass)
{
    // We make an array to save all values for each axis.
    long double *values0 = malloc(trianglesCount * sizeof(long double));
    long double *values1 = malloc(trianglesCount * sizeof(long double));
    long double *values2 = malloc(trianglesCount * sizeof(long double));

    // We place these three lists in a separate array so that we can call the corresponding array depending on what axis we are calculating t_i for
    long double *valuesOfAxes[3] = {values0, values1, values2};
    
    // We create counters for each list of values so that we can keep track of the end.
    int counters[3] = {0, 0, 0};
    
    #pragma omp parallel for
    for (int i = 0; i < trianglesCount; i++)
    {
        point v0 = triangles[i].v0;
    
        // We do calculations which are the same for every axis as described in the paper.
        point difference = subtract_point(&v0, &centerOfMass);
        
        long double innerProduct = inner_product(&difference, &normals[i]);  
            
        for(int j = 0; j < 3; j++)
        {
            point ac = subtract_point(&axes[j], &centerOfMass);
            long double D = inner_product(&ac, &normals[i]);
            
            if (D != 0)
            {                 
                long double t_i = innerProduct / D; 
                #pragma omp critical
                {
                    valuesOfAxes[j][counters[j]] = t_i;
                    counters[j] += 1;   
                }            
            }       
        }
    }
    
    // For each axis we collect the smallest positive value and largest negative value 
    // in all of our calculated values. They are combined to find the axis length.
    for(int j = 0; j < 3; j++)
    {
        long double positiveMinimum = 0;
        long double negativeMaximum = 0;
       
        for(int i = 0; i <= counters[j]; i++)
        {
            if(valuesOfAxes[j][i] > 0 && (valuesOfAxes[j][i] < positiveMinimum || positiveMinimum == 0))
                positiveMinimum = valuesOfAxes[j][i]; 

            else if(valuesOfAxes[j][i] < 0 && (valuesOfAxes[j][i] > negativeMaximum || negativeMaximum == 0))
                negativeMaximum = valuesOfAxes[j][i];
        }

        axesLengths[j] = positiveMinimum - negativeMaximum;
    }

    free(values0);
    free(values1);
    free(values2);
    return 1;
}

/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

int CopyPoint(point *copy, point *original)
{
	copy->x = original->x;
	copy->y = original->y;
	copy->z = original->z;
	return 1;
}

// We calculated the maximum innerproduct between a point and every vertex of a convex hull.
int Sup(triangle *triangles, int trianglesCount, point vertex, point *sup) 
{
	long double bigSup = inner_product(&vertex, &triangles[0].v0);
	point result = {0, 0, 0};

	#pragma omp parallel for
	for (int t = 0; t < trianglesCount; t++) 
	{
		long double supp0 = inner_product(&vertex, &triangles[t].v0);
		if (supp0 > bigSup) 
		{
			#pragma omp critical
			{
				bigSup = supp0;
				CopyPoint(&result, &triangles[t].v0);
			}
		}
		long double supp1 = inner_product(&vertex, &triangles[t].v1);
		if (supp1 > bigSup) 
		{
			#pragma omp critical
			{
				bigSup = supp1;
				CopyPoint(&result, &triangles[t].v1);
			}
		}
		long double supp2 = inner_product(&vertex, &triangles[t].v2);
		if (supp2 > bigSup) 
		{
			#pragma omp critical
			{
				bigSup = supp2;
				CopyPoint(&result, &triangles[t].v2);
			}
		}
	}

	CopyPoint(sup, &result);
	return 1;
}

// The three-dimensional subroutine. We refer to our paper for more info.
int S3D(point tau[4], point dist[4], int *tauCounter, int *distCounter)  
{
	printf("S3D REACHED\n");
	fflush(stdout);

	long double detM23 = tau[2].y * tau[3].z - tau[3].y * tau[2].z;
	long double detM13 = tau[1].y * tau[3].z - tau[3].y * tau[1].z;
	long double detM12 = tau[1].y * tau[2].z - tau[2].y * tau[1].z;
	long double detM03 = tau[0].y * tau[3].z - tau[3].y * tau[0].z;
	long double detM02 = tau[0].y * tau[2].z - tau[2].y * tau[0].z;
	long double detM01 = tau[0].y * tau[1].z - tau[1].y * tau[0].z;
	long double C[4] = {tau[1].x * detM23 - tau[2].x * detM13 + tau[3].x * detM12,
	-tau[0].x * detM23 + tau[2].x * detM03 - tau[3].x * detM02,
	tau[0].x * detM13 - tau[1].x * detM03 + tau[3].x * detM01,
	-tau[0].x * detM12 + tau[1].x * detM02 - tau[2].x * detM01};
	long double detM = 0;
	for (int i = 0; i<4; i++) 
		detM += C[i];

	int check = 0;
	for (int i = 0; i<4; i++) 
		if (detM * C[i] > 0)
			check++;

	if (check == 4) 
    {
		point v = {0, 0, 0};
		for (int i = 0; i < 4; i++) 
		{
			point addV = scalar_point(&tau[i], (C[i]/detM));
		    v = add_point(&v, &addV);
		}
		
		CopyPointArray(dist, tau, distCounter, tauCounter); 
        CopyPoint(&dist[*distCounter], &v); 
		*distCounter += 1;
	}
	else 
    {
		long double d;
		for (int i = 1; i < 4; i++) // i=2 according to paper
		{
			if (detM*C[i] < 0) // -C[i] according to paper
			{
				point tau0[4];
				int tau0Counter = 0;
				CopyPointArray(tau0, tau, &tau0Counter, tauCounter);
				
				ErasePointArrayElement(tau0, &tau0Counter, i);
                point dist0[4];
				int dist0Counter = 0;
				S2D(tau0, dist0, &tau0Counter, &dist0Counter);

				long double d0 = point_length(&dist0[dist0Counter - 1]);
				if (i==1) 
                {
					CopyPointArray(dist, dist0, distCounter, &dist0Counter);
					d = d0;
				}
				else 
                {
					if (d0 < d) 
                    {
						CopyPointArray(dist, dist0, distCounter, &dist0Counter);
						d = d0;
					}
				}
			}
		}
	}
	return 1;
}

// The two-dimensional subroutine. We refer to our paper for more info.
int S2D(point tau[4], point dist[4], int *tauCounter, int *distCounter)  
{
	printf("S2D REACHED\n");
	fflush(stdout);
	point subtau1 = subtract_point(&tau[1], &tau[2]);
	point subtau2 = subtract_point(&tau[0], &tau[2]);
	point n = {(tau[1].y*tau[0].z + tau[2].y*tau[1].z + tau[0].y*tau[2].z - tau[1].y*tau[2].z - tau[0].y*tau[1].z - tau[2].y*tau[0].z),
			   (tau[1].z*tau[0].x + tau[2].z*tau[1].x + tau[0].z*tau[2].x - tau[1].z*tau[2].x - tau[0].z*tau[1].x - tau[2].z*tau[0].x),// We took away a -.
			   (tau[1].x*tau[0].y + tau[2].x*tau[1].y + tau[0].x*tau[2].y - tau[1].x*tau[2].y - tau[0].x*tau[1].y - tau[2].x*tau[0].y)};
	point p = scalar_point(&n, inner_product(&tau[2], &n) / inner_product(&n, &n));

	long double mu;
	long double C[3];


	if (fabsl(n.x) < fabsl(n.y)) 
    {
		// We store the largest component of n in mu.
		if (fabsl(n.y) < fabsl(n.z)) 
        {
			mu = n.z;
			C[0] = -(p.x*tau[0].y + p.y*tau[1].x + tau[0].x*tau[1].y - p.x*tau[1].y - p.y*tau[0].x - tau[1].x*tau[0].y);
			C[1] =  (p.x*tau[0].y + p.y*tau[2].x + tau[0].x*tau[2].y - p.x*tau[2].y - p.y*tau[0].x - tau[2].x*tau[0].y);
			C[2] = -(p.x*tau[1].y + p.y*tau[2].x + tau[1].x*tau[2].y - p.x*tau[2].y - p.y*tau[1].x - tau[2].x*tau[1].y);
		}
		else 
        {
			mu = n.y;
			C[0] = -(p.x*tau[0].z + p.z*tau[1].x + tau[0].x*tau[1].z - p.x*tau[1].z - p.z*tau[0].x - tau[1].x*tau[0].z);
			C[1] =  (p.x*tau[0].z + p.z*tau[2].x + tau[0].x*tau[2].z - p.x*tau[2].z - p.z*tau[0].x - tau[2].x*tau[0].z);
			C[2] = -(p.x*tau[1].z + p.z*tau[2].x + tau[1].x*tau[2].z - p.x*tau[2].z - p.z*tau[1].x - tau[2].x*tau[1].z);
		}
	}
	else 
    {
		if (fabsl(n.x) < fabsl(n.z)) 
        {
			mu = n.z;
			C[0] = -(p.x*tau[0].y + p.y*tau[1].x + tau[0].x*tau[1].y - p.x*tau[1].y - p.y*tau[0].x - tau[1].x*tau[0].y);
			C[1] =  (p.x*tau[0].y + p.y*tau[2].x + tau[0].x*tau[2].y - p.x*tau[2].y - p.y*tau[0].x - tau[2].x*tau[0].y);
			C[2] = -(p.x*tau[1].y + p.y*tau[2].x + tau[1].x*tau[2].y - p.x*tau[2].y - p.y*tau[1].x - tau[2].x*tau[1].y);
		}
		else 
        {
			mu = n.x;
			C[0] = -(p.y*tau[0].z + p.z*tau[1].y + tau[0].y*tau[1].z - p.y*tau[1].z - p.z*tau[0].y - tau[1].y*tau[0].z);
			C[1] =  (p.y*tau[0].z + p.z*tau[2].y + tau[0].y*tau[2].z - p.y*tau[2].z - p.z*tau[0].y - tau[2].y*tau[0].z);
			C[2] = -(p.y*tau[1].z + p.z*tau[2].y + tau[1].y*tau[2].z - p.y*tau[2].z - p.z*tau[1].y - tau[2].y*tau[1].z);
		}
	}

	int check = 0;
	for (int i = 0; i < 3; i++)
		if (mu * C[i] > 0)
			check++;

	if (check == 0)
	{
		long double d = point_length(&tau[0]) + 1; // Must be higher to ensure dist0 is copied to dist at least once.
		for (int i = 1; i < 3; i++)  
        {
			if (mu * C[i] < 0) 
			{
				point tau0[4];
				int tau0Counter = 0;
				CopyPointArray(tau0, tau, &tau0Counter, tauCounter);
				
				ErasePointArrayElement(tau0, &tau0Counter, 2-i);

				point dist0[4];
				int dist0Counter = 0;
				S1D(tau0, dist0, &tau0Counter, &dist0Counter);

				long double d0 = point_length(&dist0[dist0Counter - 1]);
				if (i==1) 
				{
					CopyPointArray(dist, dist0, distCounter, &dist0Counter);
					d = d0;
				}
				else 
				{
					if (d0 < d) 
					{
						CopyPointArray(dist, dist0, distCounter, &dist0Counter);
						d = d0;
					}
				}
			}
		}
	}
	else if (check == 3) 
    {
		point v = {0, 0, 0};
		for (int i = 0; i<3; i++)
		{ 	
			point addV = scalar_point(&tau[i], C[i]/mu);
		    v = add_point(&v, &addV);
		}
		
		CopyPointArray(dist, tau, distCounter, tauCounter); 
        CopyPoint(&dist[3], &v); 
		*distCounter = 4;
	}
	else if (mu * C[2] < 0)
    {
		point tau0[4];
		int tau0Counter = 0;
		CopyPointArray(tau0, tau, &tau0Counter, tauCounter);
		
		ErasePointArrayElement(tau0, &tau0Counter, 2);

		point dist0[4];
		int dist0Counter = 0;
		S1D(tau0, dist0, &tau0Counter, &dist0Counter);

		CopyPointArray(dist, dist0, distCounter, &dist0Counter);
	}
	else if (mu * C[1] < 0)
    {
		point tau0[4];
		int tau0Counter = 0;
		CopyPointArray(tau0, tau, &tau0Counter, tauCounter);
		
		ErasePointArrayElement(tau0, &tau0Counter, 1);

		point dist0[4];
		int dist0Counter = 0;
		S1D(tau0, dist0, &tau0Counter, &dist0Counter);

		CopyPointArray(dist, dist0, distCounter, &dist0Counter);
	}
	else
	{
		point tau0[4];
		int tau0Counter = 0;
		CopyPointArray(tau0, tau, &tau0Counter, tauCounter);
		
		ErasePointArrayElement(tau0, &tau0Counter, 0);

		point dist0[4];
		int dist0Counter = 0;
		S1D(tau0, dist0, &tau0Counter, &dist0Counter);

		CopyPointArray(dist, dist0, distCounter, &dist0Counter);
	}
	return 1;
}

// The one-dimensional subroutine. We refer to our paper for more info.
int S1D(point tau[4], point dist[4], int *tauCounter, int *distCounter) 
{
	point n = subtract_point(&tau[0], &tau[1]);
	long double p; 

	long double mu;
	long double C[2];

	// We store the largest component of n in mu.
	if (fabsl(n.x) < fabsl(n.y)) 
    {	
		if (fabsl(n.y) < fabsl(n.z)) 
        {
			mu = -n.z;
			p = inner_product(&tau[0], &n) / inner_product(&n, &n) * mu + tau[0].z;
			C[0] = tau[1].z - p;
			C[1] = p - tau[0].z;
		}
		else 
        {
			mu = -n.y;
			p = inner_product(&tau[0], &n) / inner_product(&n, &n) * mu + tau[0].y;
			C[0] = tau[1].y - p;
			C[1] = p - tau[0].y;
		}
	}
	else 
    {
		if (fabsl(n.x) < fabsl(n.z)) 
        {
			mu = -n.z;
			p = inner_product(&tau[0], &n) / inner_product(&n, &n) * mu + tau[0].z;
			C[0] = tau[1].z - p;
			C[1] = p - tau[0].z;
		}
		else 
        {
			mu = -n.x;
			p = inner_product(&tau[0], &n) / inner_product(&n, &n) * mu + tau[0].x;
			C[0] = tau[1].x - p;
			C[1] = p - tau[0].x;
		}
	}

	int check = 0;
	for (int i = 0; i < 2; i++)
	{
		if (mu * C[i] > 0)
		{
			check++;
		}
	}

	if (check == 2) 
    {
		point v = {0, 0, 0};
		for (int i = 0; i < 2; i++)
		{
			point addV = scalar_point(&tau[i], C[i] / mu);
		    v = add_point(&v, &addV);
		}

		CopyPointArray(dist, tau, distCounter, tauCounter); 
        CopyPoint(&dist[2], &v); 
		*distCounter = 3;
	}
	else 
    {
		if (mu * C[0] > 0)
		{
			CopyPointArray(dist, tau, distCounter, tauCounter); 
       		CopyPoint(&dist[0], &tau[1]);
		}
		else
		{
			CopyPointArray(dist, tau, distCounter, tauCounter); 
        	CopyPoint(&dist[1], &tau[0]);
		}
	}
	return 1;
}

// We select the correct subroutine for calculating distance,
// based on the number of used elements in tau.
int Distance(point tau[4], point dist[4], int *tauCounter, int *distCounter) 
{
	if (*tauCounter == 4)
	{
		S3D(tau, dist, tauCounter, distCounter);
	}

	if (*tauCounter == 3)
	{
		S2D(tau, dist, tauCounter, distCounter);
	}
	
	if (*tauCounter == 2) 
	{
		S1D(tau, dist, tauCounter, distCounter);
	}

	if (*tauCounter == 1) 
    {
		CopyPointArray(dist, tau, distCounter, tauCounter); 
        CopyPoint(&dist[1], &tau[0]); 
		*distCounter = 2;
	}
	return 1;
}

// We calculate the shortest distance between two convex hulls based on
// the GJK algorithm described in our paper.
long double GJK(triangle *triangles1, triangle *triangles2, int trianglesCount1, int trianglesCount2, point centerOfMass1, point centerOfMass2) 
{
	point v = subtract_point(&centerOfMass1, &centerOfMass2);
	point w;

	// We keep an array of certain points and an integer representing 
	// the number of elements in the array that are used.
	point tau[4];
    int tauCounter = 0;

	int testCounter = 0;
	// We try to find the best distance until we encounter a problem.
	while (1) 
    {
        point minusV = scalar_point(&v, -1);
		point sup1;
		Sup(triangles1, trianglesCount1, minusV, &sup1);
		point sup2;
		Sup(triangles2, trianglesCount2, v, &sup2);

		w = subtract_point(&sup1, &sup2);

		
        for (size_t i = 0; i < tauCounter; i++)
        {
            long double lengthV = point_length(&v);

		    if (compare_points(&tau[i], &w) || fabsl(pow(lengthV, 2) - inner_product(&v, &w)) <= pow(10, -20) * pow(lengthV, 2))
			{
				long double test = pow(lengthV, 2) - inner_product(&v, &w);
				return lengthV;
			}
        }

		// We insert w in the first slot of tau
        for (size_t i = 0; i < tauCounter; i++)
            CopyPoint(&tau[tauCounter - i], &tau[tauCounter - i - 1]);
        CopyPoint(&tau[0], &w);
		tauCounter++;
        

		// We output all vertices of the simplex and v in one vector dist (v is always the last element)
		point dist[4];
		int distCounter = 0;

		//triangle *help = triangles1;
        Distance(tau, dist, &tauCounter, &distCounter);

		CopyPoint(&v, &dist[distCounter - 1]);
		tauCounter = 0;
		for (int i = 0; i < distCounter - 1; i++)
			CopyPoint(&tau[i], &dist[i]);
		tauCounter = distCounter - 1;

		float max = 0;
		for (int i = 0; i < tauCounter; i++) {
			float len = point_length(&tau[i]);
			if (len > max) {
				max = len;
			}
		}

		//check whether "Distance" doesn't have too many vertices and vertex v is not too close to the origin
		long double lengthV = point_length(&v);
		if (tauCounter == 4 || pow(lengthV, 2) <= pow(10, -10) * max)
		{
			return lengthV;
		}

		if (testCounter > 50)
			return -1;
		testCounter++;
	}
}

int CopyPointArray(point *copy, point *original, int *copyLength, int *originalLength)
{
	for (size_t i = 0; i < *originalLength; i++)
	{
		copy[i].x = original[i].x;
		copy[i].y = original[i].y;
		copy[i].z = original[i].z;
	}
	*copyLength = *originalLength;
	return 1;
}

int ErasePointArrayElement(point *points, int *counter, int index)
{
	for (size_t i = index; i < *counter - 1; i++)
		CopyPoint(&points[i], &points[i + 1]);

	*counter = *counter - 1;
	return 1;
}
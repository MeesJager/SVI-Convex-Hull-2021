
/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */
    
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define VOXELSIZE 1
#define float_error 0.000001

#include "essentials/mathBasics.h"
#include "data_calculation/convex_hull_properties/basicProperties.h"
#include "data_calculation/convex_hull_properties/boundingBox.h"
#include "data_calculation/voxelisation/surfaceVoxelisation.h"
#include "data_calculation/voxelisation/volumeVoxelisation.h"
#include "data_calculation/hull_comparison/similarity.h"
#include "data_calculation/hull_comparison/shortestDistance.h"
#include "IO/input.h"
#include "IO/output.h"

#include "essentials/mathBasics.c"
#include "data_calculation/convex_hull_properties/basicProperties.c"
#include "data_calculation/convex_hull_properties/boundingBox.c"
#include "data_calculation/voxelisation/surfaceVoxelisation.c"
#include "data_calculation/voxelisation/volumeVoxelisation.c"
#include "data_calculation/hull_comparison/similarity.c"
#include "data_calculation/hull_comparison/shortestDistance.c"
#include "IO/input.c"
#include "IO/output.c"


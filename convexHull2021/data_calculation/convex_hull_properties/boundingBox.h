
/*  This is a source file for the Convex Hull 2021 project for SVI 
    made as a part of the Mathematics for Industry course at Utrecht
    University by David de Best, Mees Jager, Luuk Lagendijk, Joost
    Mein, Joost Oliemans, Hannah Onverwagt, Andeos Rigas, Corijn 
    Rudrum, Mark Sterkenburg en Hannah van der Zande.               */

int CalculateBoundingBox(triangle *triangles, int trianglesCount, point boundingBox[2]);

int BoundingBoxToGrid (point boundingBox[2]);

int TranslateConvexHull(triangle *triangles, int trianglesCount, point boundingBox[2], point translation);

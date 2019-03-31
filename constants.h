#ifndef CONSTANTS_H_
#define CONSTANTS_H_


const double Pi = 3.14159;
const double t = 2.0;
const double t3 = 2.0;
const double tp = 2.0;
const double t3p = 2.0;
const double mp = 2 * 2.5;
const double Delta = 0.2;
const double DeltaOTE = 0.2;
const double hm = 0.0;

const bool STI = true;

const double xi = 8;
const double lmb = 8;

const int NPx = 24;
const int NPy = 24;
const int NPz = 1;
const int NP[3] = {NPx, NPy, NPz};
const int xcent = NPx / 2;
const int ycent = NPy / 2;

//Auxiliary definitions
const double sqrt2 = sqrt(2.0);

const double thresh = 1E-6;
const int neigs = 40;
const int nCore = 4;
const double shift = 0.0;
#endif

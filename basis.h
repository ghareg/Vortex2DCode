#ifndef BASIS_H_
#define BASIS_H_
#include <stdint.h>
#include <complex>
#include "constants.h"

typedef int8_t Occup;
typedef int Count;
typedef std::complex<double> MatType;

const MatType I(0.0, 1.0);

class State
{
public:
	State(Occup ph, Occup s, Occup orb, Occup xpos, Occup ypos, Occup zpos):
		mph(ph), ms(s), morb(orb), mpos{xpos, ypos, zpos} {}

	State(): mph(-1), ms(-1), morb(0), mpos{0, 0, 0} {} 

	Occup mph;
	Occup ms;
	Occup morb;
	Occup mpos[3];
};

struct Param
{
	double kz;
	double mu;
	double dmu;
	double hm;
	double DeltaS;
	double DeltaP;
};

void GenerateBasis(State*& basis, Count& bSize);

#endif

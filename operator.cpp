#include "operator.h"
#include "constants.h"
#include <cmath>

#define PBC 0
#define VRT_SIZE 2
#define FIELD_SIZE 2
#define OTE_SIZE 2

//Discretization with initial F combining first and second derivatives
MatType MatElement(Count i, Count j, const State* basis, const Param& pm)
{
	Occup ph = basis[i].mph;
	Occup s = basis[i].ms;
	Occup orb = basis[i].morb;
	Occup xpos = basis[i].mpos[0];
	Occup ypos = basis[i].mpos[1];
	Occup zpos = basis[i].mpos[2];
	
	Occup php = basis[j].mph;
	Occup sp = basis[j].ms;
	Occup orbp = basis[j].morb;
	Occup xposp = basis[j].mpos[0];
	Occup yposp = basis[j].mpos[1];
	Occup zposp = basis[j].mpos[2];
	MatType result = MatType(0.0, 0.0);

	if ((abs(xpos - xposp) > 1 && abs(xpos - xposp) != NPx - 1) || 
			(abs(ypos - yposp) > 1 && abs(ypos - yposp) != NPy - 1) || 
			(abs(zpos - zposp) > 1 && abs(zpos - zposp) != NPz - 1)) {
		return 0.0;
	}

	if (ph == php) {
		if (s == sp) {
			if (orb == orbp) {
#if FIELD_SIZE == 1
				double xdif = xpos - xcent - 0.5;
				double ydif = ypos - ycent - 0.5;
				double r = std::sqrt(xdif * xdif + ydif * ydif);
				result += (-ph * 1.0 * (-orb * (mp - t3 * std::cos(pm.kz)) - pm.mu) - s * pm.hm * std::exp(-r / lmb)) * Ident(basis[i], basis[j]) - ph * 1.0 *
					  orb * t * (Cosk(0, basis[i], basis[j]) + Cosk(1, basis[i], basis[j]));
#elif FIELD_SIZE == 2
				double xdif = xpos - xcent - 0.5;
				double ydif = ypos - ycent - 0.5;
				if (std::abs(xdif) < lmb && std::abs(ydif) < lmb) {
					result += (-ph * 1.0 * (-orb * (mp - t3 * std::cos(pm.kz)) - pm.mu) - s * pm.hm) * Ident(basis[i], basis[j]) - ph * 1.0 *
						orb * t * (Cosk(0, basis[i], basis[j]) + Cosk(1, basis[i], basis[j]));
				}
				else {
					result += -ph * 1.0 * (-orb * (mp - t3 * std::cos(pm.kz)) - pm.mu) * Ident(basis[i], basis[j]) - ph * 1.0 *
						orb * t * (Cosk(0, basis[i], basis[j]) + Cosk(1, basis[i], basis[j]));
				}
#else 
				result += (-ph * 1.0 * (-orb * (mp - t3 * std::cos(pm.kz)) - pm.mu) - s * pm.hm) * Ident(basis[i], basis[j]) - ph * 1.0 *
					  orb * t * (Cosk(0, basis[i], basis[j]) + Cosk(1, basis[i], basis[j]));
#endif
			}
			else {
				result += ph * s * 1.0 * t3p * std::sin(pm.kz) * Ident(basis[i], basis[j]);
			}
		}
		else {
			if (orb != orbp) {
				result += -ph * 1.0 * tp * Sink(0, basis[i], basis[j]);
				result += -ph * s * 1.0 * I * tp * Sink(1, basis[i], basis[j]);
			}	
		}
	}
	else {
		if (s == sp && orb == orbp) {
			result += Super(basis[i], basis[j], pm);
		}
		else if (s == sp && orb != orbp) {
			//result += -s * 1.0 * SuperOTEPPP(basis[i], basis[j], pm);
			result += -ph * 1.0 * I * SuperOTEFB(basis[i], basis[j], pm);
			//result += -ph * orb * s * 1.0 * SuperOTEFB(basis[i], basis[j], pm);
		}
	}

	return result;	
}

MatType Super(const State& lft, const State& rgt, const Param& pm)
{
	if (lft.mpos[0] == rgt.mpos[0] && lft.mpos[1] == rgt.mpos[1] &&
			lft.mpos[2] == rgt.mpos[2]) {
		//return pm.DeltaS;
		double xdif = lft.mpos[0] - xcent - 0.5;
		double ydif = lft.mpos[1] - ycent - 0.5;
		double ang = atan2(ydif, xdif);
#if VRT_SIZE == 1
		double r = std::sqrt(xdif * xdif + ydif * ydif);
		return pm.DeltaS * std::tanh(r / xi) * std::exp(-lft.mph * ang * I);
#elif VRT_SIZE == 2
		if (std::abs(xdif) > xi || std::abs(ydif) > xi) {
			return pm.DeltaS * std::exp(-lft.mph * ang * I);
		}
		else {
			return 0.0;
		}
#else 
		return pm.DeltaS * std::exp(-lft.mph * ang * I);
#endif
	}

	return 0.0;
}

MatType SuperOTEFB(const State& lft, const State& rgt, const Param& pm)
{
	if (lft.mpos[0] == rgt.mpos[0] && lft.mpos[1] == rgt.mpos[1] &&
			lft.mpos[2] == rgt.mpos[2]) {
		double xdif = lft.mpos[0] - xcent - 0.5;
		double ydif = lft.mpos[1] - ycent - 0.5;
		double ang = atan2(ydif, xdif);
	
#if OTE_SIZE == 1
		double r = std::sqrt(xdif * xdif + ydif * ydif);
		return pm.DeltaP * (1.0 - std::tanh(r / xi)) * std::exp(-lft.mph * ang * I);
#elif OTE_SIZE == 2
		if (std::abs(xdif) < xi && std::abs(ydif) < xi) {
			return pm.DeltaP * std::exp(-lft.mph * ang * I);
		}
		else {
			return 0.0;
		}
#else
		return pm.DeltaP * std::exp(-lft.mph * ang * I);
#endif
	}

	return 0.0;
}

MatType SuperOTEPPP(const State& lft, const State& rgt, const Param& pm)
{
	double ang1 = atan2(lft.mpos[1] - ycent - 0.5, lft.mpos[0] - xcent - 0.5);
	double ang2 = atan2(rgt.mpos[1] - ycent - 0.5, rgt.mpos[0] - xcent - 0.5);

#if OTE_SIZE == 1
	double r = 0.5 * (std::sqrt((lft.mpos[0] - xcent - 0.5) * (lft.mpos[0] - xcent - 0.5) + (lft.mpos[1] - ycent - 0.5) * (lft.mpos[1] - ycent - 0.5)) +
		std::sqrt((rgt.mpos[0] - xcent - 0.5) * (rgt.mpos[0] - xcent - 0.5) + (rgt.mpos[1] - ycent - 0.5) * (rgt.mpos[1] - ycent - 0.5)));

	return 0.5 * pm.DeltaP * (Sink(0, lft, rgt) -lft.mph * 1.0 * I * Sink(1, lft, rgt)) * (1.0 - std::tanh(r / xi)) *
		(std::exp(-lft.mph * ang1 * I) + std::exp(-lft.mph * ang2 * I));
#elif OTE_SIZE == 2
	if (std::abs(lft.mpos[0] - xcent - 0.5) < xi && std::abs(lft.mpos[1] - ycent - 0.5) < xi && 
			std::abs(rgt.mpos[0] - xcent - 0.5) < xi && std::abs(rgt.mpos[1] - ycent - 0.5) < xi) {
		return 0.5 * pm.DeltaP * (Sink(0, lft, rgt) -lft.mph * 1.0 * I * Sink(1, lft, rgt)) *
			(std::exp(-lft.mph * ang1 * I) + std::exp(-lft.mph * ang2 * I));
	}
	return 0.0;
#else 
	return 0.5 * pm.DeltaP * (Sink(0, lft, rgt) -lft.mph * 1.0 * I * Sink(1, lft, rgt)) *
		(std::exp(-lft.mph * ang1 * I) + std::exp(-lft.mph * ang2 * I));
#endif
}

MatType Cosk(int cord, const State& lft, const State& rgt)
{
	MatType matel = 0.0;
	for (int cd = 0; cd < 3; ++cd) {
		if (cd == cord) {
#if PBC
			if (lft.mpos[cord] == (rgt.mpos[cord] + 1) % NP[cord]) {
#else
			if (lft.mpos[cord] == rgt.mpos[cord] + 1) {
#endif
				matel = 0.5;
			}
#if PBC
			else if (rgt.mpos[cord] == (lft.mpos[cord] + 1) % NP[cord]) {
#else
			else if (rgt.mpos[cord] == lft.mpos[cord] + 1) {
#endif
				matel = 0.5;
			}
			else {
				return 0.0;
			}
		}
		else {
			if (lft.mpos[cd] != rgt.mpos[cd]) {
				return 0.0;
			}
		}
	}

	return matel;
}

MatType Sink(int cord, const State& lft, const State& rgt)
{
	MatType matel = 0.0;
	for (int cd = 0; cd < 3; ++cd) {
		if (cd == cord) {
#if PBC
			if (lft.mpos[cord] == (rgt.mpos[cord] + 1) % NP[cord]) {
#else
			if (lft.mpos[cord] == rgt.mpos[cord] + 1) {
#endif
				matel = -0.5 * I;
			}
#if PBC
			else if (rgt.mpos[cord] == (lft.mpos[cord] + 1) % NP[cord]) {
#else
			else if (rgt.mpos[cord] == lft.mpos[cord] + 1) {
#endif
				matel = 0.5 * I;
			}
			else {
				return 0.0;
			}
		}
		else {
			if (lft.mpos[cd] != rgt.mpos[cd]) {
				return 0.0;
			}
		}
	}

	return matel;
}

MatType Ident(const State& lft, const State& rgt)
{
	for (int cd = 0; cd < 3; ++cd) {
		if (lft.mpos[cd] != rgt.mpos[cd]) {
			return 0.0;
		}
	}

	return 1.0;
}

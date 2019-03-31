#include "basis.h"
#include "matrix.h"
#include <iostream>

void GenerateBasis(State*& basis, Count& bSize)
{
	bSize = 0;
	Count maxbSize = NPx * NPy * NPz * 2 * 2 * 2;
	basis = new State[maxbSize];
	for (Occup ph = -1; ph <= 1; ph += 2) {
		for (Occup s = -1; s <= 1; s += 2) {
			for (Occup orb = -1; orb <= 1; orb += 2) {
				for (Occup xpos = 0; xpos < NPx; ++xpos) {
					for (Occup ypos = 0; ypos < NPy; ++ypos) {
						for (Occup zpos = 0; zpos < NPz; ++zpos) {
							if (bSize == maxbSize) {
								maxbSize = 3 * maxbSize / 2;
								resize(basis, bSize, maxbSize);
							}
							basis[bSize++] = State(ph, s, orb, xpos, ypos, zpos);
						}
					}
				}
			}
		}
	}
}

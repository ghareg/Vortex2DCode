#include <cmath>
#include "matrix.h"
#include "constants.h"
#include "operator.h"
#include "diag.h"
#include <thread>
#include <stdio.h>
#include <algorithm>
#include <time.h>
#include <limits>

struct eigen {
	MatType value;
	const MatType *vector;
	
	bool operator<(eigen const &other) const {
		return value.real() < other.value.real();
	}
};

void GenerateMuMat(double*& muMat, Count& muSize);
void ProbDensity(int st, const Param& pm, const MatType* state, Count bSize, const State* basis);
void sortEValues(eigen* eval, const MatType* evalues, const MatType* evecs, Count bSize);
void FourierTransform(int st, int size, const Param& pm, const eigen* eval, Count bSize, const State* basis, MatType* Fstate);
void free(SparseMat& baseHam, MatType* evalues, MatType* evecs);
void free(SparseMat& baseHam);



int main(void)
{
	State* basis = nullptr;
	Count bSize;
	GenerateBasis(basis, bSize);
	
	double* muMat = nullptr;
	Count muSize;
	GenerateMuMat(muMat, muSize);

	SparseMat baseHam;
	bool zeroBased = false;
	MatType* evalues = new MatType[neigs];
	MatType* evecs = new MatType[neigs * bSize];
	eigen* eval = new eigen[neigs];
	GenMatProd op;
	bool matrixCons = false;
	int Ni = 0;

	Param pm;
	double kz = 0.0;
//	double mu = 1.603;

	FILE* evalFile;
	char str[100];
	sprintf (str, "Endepmu2DNx%dNy%dNz%d%sTIkz%sh%.2fxi%dlmb%dOTE%.2fFBD2ESE%.2fRB.dat", NPx, NPy, NPz, STI ? "S" : "W", kz > 0.1 ? "Pi" : "0.0", hm, (int) xi, (int) lmb, DeltaOTE, Delta);
//	sprintf (str, "Endepkz2DNx%dNy%dNz%d%sTImu%.2fh%.2fxi%dlmb%dOTE%.2fFBD2ESE%.2fRB.dat", NPx, NPy, NPz, STI ? "S" : "W", mu, hm, (int) xi, (int) lmb, DeltaOTE, Delta);
	
	//muSize = Ni + 1; 
	evalFile = fopen(str, "w");
	for (int muCount = Ni; muCount < muSize; ++muCount) {
//	for (double kz = -Pi; kz <= Pi; kz += 2 * Pi / 600) {
		fprintf(evalFile, "%.4f", muMat[muCount]);
		pm.kz = kz;
		pm.mu = muMat[muCount];
		pm.dmu = 0.0;
		pm.hm = hm;
		pm.DeltaS = Delta;
		pm.DeltaP = DeltaOTE;
		if (!matrixCons) {
			calcFullMat(baseHam, basis, bSize, pm, zeroBased);
			
			op.init(baseHam);
			matrixCons = true;
		}
		else {
			pm.dmu = muMat[muCount] - muMat[muCount - 1];
			changeFullUMat(baseHam, pm.dmu, zeroBased);
			fflush(evalFile);
			//op.restart(pm.dmu);
			op.restart(baseHam);
		}
		calcEValues(baseHam, op, evalues, evecs);
		sortEValues(eval, evalues, evecs, bSize);
/*
   		int st = 19;
		int size = 4;
		MatType* Fstate = new MatType[size * bSize];
		FourierTransform(st, size, pm, eval, bSize, basis, Fstate);
		for (int i = 0; i < size; ++i) {
			ProbDensity(st + i, pm, &Fstate[i * bSize], bSize, basis);
		}
		delete[] Fstate;
	

		ProbDensity(19, pm, eval[19].vector, bSize, basis);
		ProbDensity(20, pm, eval[20].vector, bSize, basis);
		ProbDensity(21, pm, eval[21].vector, bSize, basis);
		ProbDensity(22, pm, eval[22].vector, bSize, basis);
	
*/	
	//	free(baseHam);
		
		for (int i = 0; i < neigs; ++i) {
			fprintf(evalFile, "\t%.6f", eval[i].value.real());
		}
		fprintf(evalFile, "\n");	
	
		fflush(evalFile);
	}
	fclose(evalFile);

	free(baseHam, evalues, evecs);
	delete[] basis;
	delete[] muMat;
	delete[] eval;
}

void sortEValues(eigen* eval, const MatType* evalues, const MatType* evecs, Count bSize)
{
	for (int i = 0; i < neigs; ++i) {
		eval[i].value = evalues[i];
		eval[i].vector = &evecs[i * bSize];
	}
	std::sort(&eval[0], &eval[0] + neigs);
}

void GenerateMuMat(double*& muMat, Count& muSize)
{
	muSize = 0;
	int nmax = 600;
	muMat = new double[nmax + 1];
	double mu0 = 0.003;
	double murange = 4.0;
	for (int j = 0; j <= nmax; ++j) {
		muMat[muSize] = mu0 + murange * j / nmax;
		muSize++;
	}
}

void ProbDensity(int st, const Param& pm, const MatType* state, Count bSize, const State* basis)
{
	double* dens = new double[NPx * NPy * NPz];
	double sum = 0.0;
	double edge = 0.0;
	for (int i = 0; i < NPx * NPy * NPz; ++i) {
		dens[i] = 0.0;
	}
	for (Count i = 0; i < bSize; ++i) {
		dens[basis[i].mpos[0] * NPy * NPz + basis[i].mpos[1] * NPz + basis[i].mpos[2]] +=
			std::norm(state[i]);
		if ((basis[i].morb == -1 && ((basis[i].mph == -1 && basis[i].ms == -1) || (basis[i].mph == 1 && basis[i].ms == 1))) || 
			(basis[i].morb == 1 && ((basis[i].mph == -1 && basis[i].ms == 1) || (basis[i].mph == 1 && basis[i].ms == -1)))) {
			sum += std::norm(state[i]);
		}
		if (basis[i].mpos[0] > xcent + 0.5 + 10 || basis[i].mpos[1] > ycent + 0.5 + 10) {
			edge += std::norm(state[i]);
		}
	}
	
	char str[100];
	sprintf (str, "ProbDens2DNx%dNy%dNz%dkz%sh%.2fmu%.2fxi%dlmb%dOTE%.2fESE%.2fFBst%d.dat", NPx, NPy, NPz, pm.kz > 0.1 ? "Pi" : "0.0", pm.hm, pm.mu, (int) xi, (int) lmb, pm.DeltaP, pm.DeltaS, st);
	FILE* evalFile = fopen(str, "w");
	for (int ix = 0; ix < NPx; ++ix) {
		for (int iy = 0; iy < NPy; ++iy) {
			for (int iz = 0; iz < NPz; ++iz) {
				fprintf(evalFile, "%d\t%d\t%d\t%.4f\n", ix, iy, iz, dens[ix * NPy * NPz + iy * NPz + iz]);
			}
		}
		fprintf(evalFile, "\n");
	}

	fclose(evalFile);
	printf("%d\t%.2f\t%.2f\n", st, sum, edge);
	delete[] dens;
}

void FourierTransform(int st, int size, const Param& pm, const eigen* eval, Count bSize, const State* basis, MatType* Fstate)
{
	const int NP = NPx * NPy;
	DFTI_DESCRIPTOR_HANDLE fft_handle;
	MKL_LONG status;
	MKL_LONG length[2] = {NPx, NPy};
	status = DftiCreateDescriptor(&fft_handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, length);
	DftiSetValue(fft_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	DftiSetValue(fft_handle, DFTI_FORWARD_SCALE, 1.0 / sqrt(NP));
	DftiSetValue(fft_handle, DFTI_BACKWARD_SCALE, 1.0 / sqrt(NP));
	MKL_LONG is[3] = {0, NPy, 1};
	DftiSetValue(fft_handle, DFTI_INPUT_STRIDES, is);
	DftiSetValue(fft_handle, DFTI_OUTPUT_STRIDES, is);
	DftiSetValue(fft_handle, DFTI_NUMBER_OF_TRANSFORMS, 8);
	DftiSetValue(fft_handle, DFTI_INPUT_DISTANCE, NP);
	DftiSetValue(fft_handle, DFTI_OUTPUT_DISTANCE, NP);
	status = DftiCommitDescriptor(fft_handle);
	//	printf(DftiErrorMessage(status));
	for (int i = 0; i < size; ++i) {
		status = DftiComputeBackward(fft_handle, (void *) eval[st + i].vector, (void *) &Fstate[i * bSize]);
		//		printf(DftiErrorMessage(status));
	}
	status = DftiFreeDescriptor(&fft_handle);
}


void free(SparseMat& baseHam, MatType* evalues, MatType* evecs)
{
	delete[] baseHam.mat;
	delete[] baseHam.ia;
	delete[] baseHam.ja;

	baseHam.mat = nullptr;
	baseHam.ia = nullptr;
	baseHam.ja = nullptr;

	delete[] evalues;
	delete[] evecs;
}

void free(SparseMat& baseHam)
{
	delete[] baseHam.mat;
	delete[] baseHam.ia;
	delete[] baseHam.ja;
	baseHam.mat = nullptr;
	baseHam.ia = nullptr;
	baseHam.ja = nullptr;
}



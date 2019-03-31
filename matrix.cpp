#include "matrix.h"
#include "constants.h"
#include "operator.h"
#include <thread>
#include <algorithm>
#include <mkl_types.h>
#undef MKL_Complex16
#define MKL_Complex16 MatType
#include <mkl.h>
#include <cstdlib>
#include <fstream>
#include <stdio.h>

void calcMat(SparseMat& baseHam, const State* basis, Count bSize, const Param& pm)
{
	Count start =  baseHam.start;
	Count end = baseHam.end;
	Count valueSize  = 3 * bSize;
	baseHam.mat = new MatType[valueSize];
	baseHam.ia = new Count[end - start + 1];
	baseHam.ja = new Count[valueSize];

	MatType matel = 0.0;
	Count currHNum = 0;
	for (Count i = start; i < end; ++i) {
		for (Count j = i; j < bSize; ++j) {
			matel = MatElement(i, j, basis, pm);
			if (j == i) {
				baseHam.ia[i - start] = currHNum;
			}
			if (j == i || std::abs(matel) > thresh) {
				if (currHNum == valueSize) {
					valueSize *= 2;
					resize(baseHam.mat, currHNum, valueSize);
					resize(baseHam.ja, currHNum, valueSize);
				}
				baseHam.mat[currHNum] = matel;
				baseHam.ja[currHNum] = j;
				++currHNum;
			}
		}
	}
	baseHam.ia[end - start] = currHNum;
}

void initStartEnd(SparseMat* baseHamMat, Count bSize)
{
	double total = 0.5 * (1.0 * bSize * bSize + bSize);	
	Count end0 = (Count) (bSize + 0.5 - 0.5 * sqrt((2.0 * bSize + 1.0) * (2.0 * bSize + 1.0) - 
			   8.0 * total / nCore));
	baseHamMat[0].start = 0;
	baseHamMat[0].end = end0;
	for (int i = 1; i < nCore - 1; ++i) {
		Count oldstart = baseHamMat[i - 1].start;
		Count oldend = baseHamMat[i - 1].end;
		baseHamMat[i].start = oldend;
		Count endi = (Count) (bSize + 0.5 - 0.5 * sqrt((2.0 * bSize + 1.0) * (2.0 * bSize + 1.0) -
				16.0 * bSize * oldend - 8.0 * oldend + 8.0 * oldend * oldend + 8.0 * bSize * 
				oldstart - 4.0 * oldstart * oldstart + 4.0 * oldstart));
		//Count endi = (Count) (bSize + 0.5 - 0.5 * sqrt((2.0 * bSize - 2.0 * oldend + 1) * (2.0 * bSize - 2.0 * oldend + 1)
		//			- 8.0 * total / nCore));
		baseHamMat[i].end = endi;
	}
	if (nCore > 1) {
		baseHamMat[nCore - 1].start = baseHamMat[nCore - 2].end;
	}
	baseHamMat[nCore - 1].end = bSize;
/*	for (int i = 0; i < nCore - 1; ++i) {
		baseHamMat[i].start = bSize * i / nCore;
		baseHamMat[i].end = bSize * (i + 1) / nCore;
	}
	baseHamMat[nCore - 1].start = bSize * (nCore - 1) / nCore;
	baseHamMat[nCore - 1].end = bSize;*/
}

void joinMat(SparseMat& MatFull, const SparseMat* MatArray, Count bSize, bool zeroBased)
{
	Count valuesSize = 0;
	for (int i = 0; i < nCore; ++i) {
		valuesSize += MatArray[i].ia[MatArray[i].end - MatArray[i].start];
	}
	MatFull.mat = new MatType[valuesSize];
	MatFull.ja = new Count[valuesSize];
	MatFull.ia = new Count[bSize + 1];
	Count jadd = zeroBased ? 0 : 1;
	Count indCount = 0;
	Count iaIndCount = 0;
	Count currTot = zeroBased ? 0 : 1;
	for (int i = 0; i < nCore; ++i) {
		for (Count j = 0; j < MatArray[i].ia[MatArray[i].end - MatArray[i].start]; ++j) {
			MatFull.mat[indCount] = MatArray[i].mat[j];
			MatFull.ja[indCount] = MatArray[i].ja[j] + jadd;
			indCount++;
		}
		for (Count j = 0; j < MatArray[i].end - MatArray[i].start; ++j) {
			MatFull.ia[iaIndCount] = currTot + MatArray[i].ia[j];
			++iaIndCount;
		}
		currTot += MatArray[i].ia[MatArray[i].end - MatArray[i].start];
		
		delete[] MatArray[i].mat;
		delete[] MatArray[i].ja;
		delete[] MatArray[i].ia;
	}
	
	MatFull.ia[bSize] = currTot; 
	MatFull.start = 0;
	MatFull.end = bSize;
}

void calcFullMat(SparseMat& baseHam, const State* basis, Count bSize, const Param& pm, bool zeroBased)
{
	SparseMat baseHamMat[nCore];
	initStartEnd(baseHamMat, bSize);
	std::thread t[nCore - 1];
	for (int i = 0; i < nCore - 1; ++i) {
		t[i] = std::thread(calcMat, std::ref(baseHamMat[i]), basis, bSize, std::cref(pm));
	}
	
	calcMat(baseHamMat[nCore - 1], basis, bSize, pm);

	for (int i = 0; i < nCore - 1; ++i) {
		t[i].join();
	}

	joinMat(baseHam, baseHamMat, bSize, zeroBased);
}

void writeMatrixToFile(const SparseMat& baseHam, bool zeroBased)
{
	std::ofstream matfile;
	std::ofstream iafile;
	std::ofstream jafile;
	matfile.open("mat.dat");
	iafile.open("ia.dat");
	jafile.open("ja.dat");
	Count add = zeroBased ? 0 : 1;
	for (Count i = 0; i < baseHam.ia[baseHam.end - baseHam.start] - add; ++i) {
		matfile << baseHam.mat[i].real() << "\t" << baseHam.mat[i].imag() << std::endl;
		jafile << baseHam.ja[i] - add << std::endl;
	}

	for (Count i = 0; i <= baseHam.end - baseHam.start; ++i) {
		iafile << baseHam.ia[i] - add << std::endl;
	}
	matfile.close();
	iafile.close();
	jafile.close();
}

void changeFullUMat(SparseMat& baseHam, double dmu, bool zeroBased)
{
	Count Size = baseHam.end - baseHam.start;
	Count add = zeroBased ? 0 : 1;
	#pragma omp parallel for
	for (Count i = 0; i < Size / 2; ++i) {
		baseHam.mat[baseHam.ia[i] - add] += MatType(-dmu, 0.0);
	}
	#pragma omp parallel for
	for (Count i = Size / 2; i < Size; ++i) {
		baseHam.mat[baseHam.ia[i] - add] += MatType(dmu, 0.0);
	}
}

Count DiagValueIndex(Count rowNum, const SparseMat& baseHam, bool zeroBased)
{
	Count add = zeroBased ? 0 : 1;
	Count rowStart = baseHam.ia[rowNum] - add;
	Count rowEnd = baseHam.ia[rowNum + 1] - add;
	Count* posStart = &baseHam.ja[rowStart];
	Count* posEnd = &baseHam.ja[rowEnd];
	Count* pos = std::lower_bound(posStart, posEnd, rowNum + add);
	if (*pos == rowNum + add && pos != posEnd) {
		return rowStart + (pos - posStart);
	}

	return -1;
}

void changeFullMat(SparseMat& baseHam, double dmu, bool zeroBased)
{
	Count Size = baseHam.end - baseHam.start;
	#pragma omp parallel for
	for (Count i = 0; i < Size / 2; ++i) {
		baseHam.mat[DiagValueIndex(i, baseHam, zeroBased)] += MatType(-dmu, 0.0);
	}
	#pragma omp parallel for
	for (Count i = Size / 2; i < Size; ++i) {
		baseHam.mat[DiagValueIndex(i, baseHam, zeroBased)] += MatType(dmu, 0.0);
	}
}

void HermitianTransform(SparseMat& MatFull, bool zeroBased, bool& isHermitian)
{
	isHermitian = true;
	const int add = zeroBased ? 0 : 1;
	const Count nimbsize = MatFull.end - MatFull.start;
	const Count Size = MatFull.ia[nimbsize] - add;
	MatType* HCmat = nullptr;
	Count* HCja = nullptr;
	Count* HCia = new Count[nimbsize + 1];
	const char trans = 'C';
	MKL_INT request = 1;
	MatType beta = MatType(1.0, 0.0);
	const MKL_INT nzmax = Size;
	MKL_INT sort = 0;
	MKL_INT n = nimbsize;
	MKL_INT info = 0;
	mkl_zcsradd (&trans, &request, &sort, &n , &n , MatFull.mat, MatFull.ja, MatFull.ia, &beta, MatFull.mat, MatFull.ja, MatFull.ia, HCmat, HCja, HCia, &nzmax, &info);
	const Count NSize = HCia[nimbsize] - add;
	HCmat = new MatType[NSize];
	HCja = new Count[NSize];
	request = 2;
	mkl_zcsradd (&trans, &request, &sort, &n , &n , MatFull.mat, MatFull.ja, MatFull.ia, &beta, MatFull.mat, MatFull.ja, MatFull.ia, HCmat, HCja, HCia, &nzmax, &info);
	if (NSize == Size) {
		for (Count i = 0; i < Size; ++i) {
			HCmat[i] *= 0.5;
			if (std::abs(HCmat[i] - MatFull.mat[i]) > thresh) {
				isHermitian = false;
			}
		}
	}
	else {
		isHermitian = false;
		for (Count i = 0; i < NSize; ++i) {
			HCmat[i] *= 0.5;
		}
	}

	delete[] MatFull.mat;
	delete[] MatFull.ia;
	delete[] MatFull.ja;
	MatFull.mat = HCmat;
	MatFull.ja = HCja;
	MatFull.ia = HCia;
}

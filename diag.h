#ifndef DIAG_H_
#define DIAG_H_
#include "matrix.h"
#include <mkl_types.h>
#undef MKL_Complex16
#define MKL_Complex16 MatType
#include <mkl.h>
//#include <Eigen/Sparse>
//#define NDEBUG

class GenMatProd
{
public:
	GenMatProd();
	~GenMatProd();

	void init(SparseMat baseHam);
	void restart(SparseMat baseHam);

	Count rows() {return qsize_;}
	Count cols() {return qsize_;}
	void set_shift(double sigma) {sigma_ = sigma;}
	void perform_op(const MatType* xIn, MatType* yOut);
private:
	SparseMat baseHam_;
	Count qsize_;
	double sigma_;
	_MKL_DSS_HANDLE_t pt_;
	MKL_INT mtype_;
	MKL_INT* iparm_;
	MKL_INT* perm_;
	MKL_INT phase_;
};
/*
class GenMatProdEig
{
public:
	typedef Eigen::SparseMatrix<MatType, Eigen::ColMajor, long int> SpMatType; 

	void init(SparseMat baseHam, bool zeroBased);
	void restart(double dmu);

	Count rows() {return SpMat_.rows();}
	Count cols() {return SpMat_.cols();}
	void set_shift(double sigma) {sigma_ = sigma;}
	void perform_op(const MatType* xIn, MatType* yOut);

private:
	SpMatType SpMat_;
	Eigen::SparseLU<SpMatType, Eigen::COLAMDOrdering<Eigen::Index> > solver_;
	double sigma_;
};

*/
void calcEValues(const SparseMat& baseHam, GenMatProd& op, MatType* evalues, MatType* evecs);
#endif

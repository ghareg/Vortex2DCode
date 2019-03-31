#ifndef OPERATOR_H_
#define OPERATOR_H_
#include "matrix.h"

MatType MatElement(Count i, Count j, const State* basis, const Param& pm);
MatType Super(const State& lft, const State& rgt, const Param& pm);
MatType SuperOTEFB(const State& lft, const State& rgt, const Param& pm);
MatType SuperOTEPPP(const State& lft, const State& rgt, const Param& pm);
MatType Cosk(int cord, const State& lft, const State& rgt);
MatType Sink(int cord, const State& lft, const State& rgt);
MatType Ident(const State& lft, const State& rgt);


#endif

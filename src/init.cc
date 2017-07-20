#include "qtlmt.h"
#include <R_ext/Rdynload.h> //R_CMethodDef
//#include <Rdefines.h>

extern void runifc(double*, int&, long*);
extern void svdc(double*, int&, int&, double*, double*);
extern void cholc(double*, int&, double*);
extern void cholsolve(double*, int&, double*, double*);
extern void lusolve(double*, int&, double*, double*);
extern void sureEstc(double*, int&, int&, double*, int&, int*, int*, double*, double*,
   double&, int&, int&, double&);
extern void sureStepc(double*, int&, int&, double*, int&, int*, int*, int*, int*,
   double&, int&, int*, double*, int&, int&, int&, double&);
extern void mimEstc(double*, double*, double*, double*, int&, int&, int&, int&,
   double*, double*, double&, double&, int&, int&, double&);
extern void mtcmimEstc(double*, int&, int&, double*, int&, double*, int&, int*,
   int*, double*, int&, int*, int*, double*, double*, double*, double&, int&, int&, double&);
extern void fPc(int*, int&, int&, int*, int&, int&, double*, int*, int*, double*,
   int*, int&, double*, int&);

static const R_CMethodDef cMethods[] = {
    {"runifc",       (DL_FUNC) &runifc,          3},
    {"svdc",         (DL_FUNC) &svdc,            5},
    {"cholc",        (DL_FUNC) &cholc,           3},
    {"cholsolve",    (DL_FUNC) &cholsolve,       4},
    {"lusolve",      (DL_FUNC) &lusolve,         4},
    {"sureEstc",     (DL_FUNC) &sureEstc,       13},
    {"sureStepc",    (DL_FUNC) &sureStepc,      17},
    {"mimEstc",      (DL_FUNC) &mimEstc,        15},
    {"mtcmimEstc",   (DL_FUNC) &mtcmimEstc,     20},
    {"fPc",          (DL_FUNC) &fPc,            14},
    {NULL, NULL, 0}
};

void R_init_qtlmt(DllInfo *dll)
{
    R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}


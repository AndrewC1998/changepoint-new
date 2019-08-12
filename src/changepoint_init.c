#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void binseg(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void CptReg_Normal_AMOC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Free_CptReg_Normal_AMOC(void *);
extern void FreeNPPELT(void *);
extern void FreePELT(void *);
extern void NPPELT(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void PELT(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"binseg",                  (DL_FUNC) &binseg,                  17},
    {"CptReg_Normal_AMOC",      (DL_FUNC) &CptReg_Normal_AMOC,      17},
    {"Free_CptReg_Normal_AMOC", (DL_FUNC) &Free_CptReg_Normal_AMOC,  1},
    {"FreeNPPELT",              (DL_FUNC) &FreeNPPELT,               1},
    {"FreePELT",                (DL_FUNC) &FreePELT,                 1},
    {"NPPELT",                  (DL_FUNC) &NPPELT,                  12},
    {"PELT",                    (DL_FUNC) &PELT,                    18},
    {NULL, NULL, 0}
};

void R_init_changepoint(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

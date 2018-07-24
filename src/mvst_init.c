#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "Utils.h"

/* .C calls */
extern void rvST(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rvSTX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rvT(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rvTX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rzSN(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rzSNX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rzST(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rzSTX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"rvST",  (DL_FUNC) &rvST,  11},
    {"rvSTX", (DL_FUNC) &rvSTX, 13},
    {"rvT",   (DL_FUNC) &rvT,    9},
    {"rvTX",  (DL_FUNC) &rvTX,  11},
    {"rzSN",  (DL_FUNC) &rzSN,  10},
    {"rzSNX", (DL_FUNC) &rzSNX, 12},
    {"rzST",  (DL_FUNC) &rzST,  11},
    {"rzSTX", (DL_FUNC) &rzSTX, 13},
    {NULL, NULL, 0}
};

void R_init_mvst(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}


#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "precision.h"
#include "power_wrappers.h"

extern void pghyperR(void *, void *, void *, void *, void *, void *);
extern void ughyperR(void *, void *, void *, void *, void *, void *);

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CMethodDef CEntries[] = {

    CALLDEF(fecrt_pee_wrap, 13),
    CALLDEF(fecrt_pvals, 23),
    CALLDEF(fecrt_power_comparison, 23),

    CALLDEF(pghyperR, 6),
    CALLDEF(ughyperR, 6),

    CALLDEF(pbnb_lower_wrap, 5),
    CALLDEF(pbnb_upper_wrap, 5),

    CALLDEF(waavp_ci_wrap, 9),
    CALLDEF(dobson_ci_wrap, 7),
    CALLDEF(conjbeta_ci_wrap, 12),
    CALLDEF(asymptotic_ci_wrap, 10),

    CALLDEF(precision_count, 12),
    CALLDEF(precision_reduction, 13),

    {NULL, NULL, 0}
};

void R_init_bayescount(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}

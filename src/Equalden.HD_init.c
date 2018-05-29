#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void pdf_sum_unif(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"pdf_sum_unif", (DL_FUNC) &pdf_sum_unif, 4},
    {NULL, NULL, 0}
};

void R_init_Equalden_HD(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

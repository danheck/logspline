#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void logcensor(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void nlogcensor(void *, void *, void *, void *, void *, void *);
extern void nlogcensorx(void *);
extern void pqlsd(void *, void *, void *, void *, void *, void *, void *, void *);
extern void rpqlsd(void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"logcensor",   (DL_FUNC) &logcensor,   13},
    {"nlogcensor",  (DL_FUNC) &nlogcensor,   6},
    {"nlogcensorx", (DL_FUNC) &nlogcensorx,  1},
    {"pqlsd",       (DL_FUNC) &pqlsd,        8},
    {"rpqlsd",      (DL_FUNC) &rpqlsd,       7},
    {NULL, NULL, 0}
};

void R_init_logspline(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

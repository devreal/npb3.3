#ifndef HAVE_PRINT_RESULTS_H
#define HAVE_PRINT_RESULTS_H

#ifdef __cplusplus
extern "C" {
#endif

extern
void c_print_results( const char   *name,
                      const char   CLASS,
                      int    n1,
                      int    n2,
                      int    n3,
                      int    niter,
                      double t,
                      double mops,
                      const char   *optype,
                      int    passed_verification,
                      const char   *npbversion,
                      const char   *compiletime,
                      const char   *cc,
                      const char   *clink,
                      const char   *c_lib,
                      const char   *c_inc,
                      const char   *cflags,
                      const char   *clinkflags );

static
void print_results( const char   *name,
                    const char   CLASS,
                    int    n1,
                    int    n2,
                    int    n3,
                    int    niter,
                    double t,
                    double mops,
                    const char   *optype,
                    int    passed_verification,
                    const char   *npbversion,
                    const char   *compiletime,
                    const char   *cc,
                    const char   *clink,
                    const char   *c_lib,
                    const char   *c_inc,
                    const char   *cflags,
                    const char   *clinkflags )
{
  c_print_results(name, CLASS, n1, n2, n3, niter, t, mops, optype,
                  passed_verification, npbversion, compiletime, cc,
                  clink, c_lib, c_inc, cflags, clinkflags);
}

#ifdef __cplusplus
}
#endif

#endif // HAVE_PRINT_RESULTS_H

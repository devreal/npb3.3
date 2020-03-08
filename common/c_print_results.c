/*****************************************************************/
/******     C  _  P  R  I  N  T  _  R  E  S  U  L  T  S     ******/
/*****************************************************************/
#include <stdlib.h>
#include <stdio.h>

#include "print_results.h"

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
                      const char   *clinkflags )
{
    printf( "\n\n %s Benchmark Completed\n", name );

    printf( " Class           =                        %c\n", CLASS );

    if( n2 == 0 && n3 == 0 )
        printf( " Size            =             %12d\n", n1 );   /* as in IS */
    else
        printf( " Size            =              %3dx %3dx %3d\n", n1,n2,n3 );

    printf( " Iterations      =             %12d\n", niter );

    printf( " Time in seconds =             %12.2f\n", t );

    printf( " Mop/s total     =             %12.2f\n", mops );

    printf( " Operation type  = %24s\n", optype);

    if( passed_verification )
        printf( " Verification    =               SUCCESSFUL\n" );
    else
        printf( " Verification    =             UNSUCCESSFUL\n" );

    printf( " Version         =             %12s\n", npbversion );

    printf( " Compile date    =             %12s\n", compiletime );

    printf( "\n Compile options:\n" );

    printf( "    CC           = %s\n", cc );

    printf( "    CLINK        = %s\n", clink );

    printf( "    C_LIB        = %s\n", c_lib );

    printf( "    C_INC        = %s\n", c_inc );

    printf( "    CFLAGS       = %s\n", cflags );

    printf( "    CLINKFLAGS   = %s\n", clinkflags );
#ifdef SMP
    evalue = getenv("MP_SET_NUMTHREADS");
    printf( "   MULTICPUS = %s\n", evalue );
#endif

    printf( "\n\n" );
    printf( " Please send all errors/feedbacks to:\n\n" );
    printf( " NPB Development Team\n" );
    printf( " npb@nas.nasa.gov\n\n" );
/*    printf( " Please send the results of this run to:\n\n" );
    printf( " NPB Development Team\n" );
    printf( " Internet: npb@nas.nasa.gov\n \n" );
    printf( " If email is not available, send this to:\n\n" );
    printf( " MS T27A-1\n" );
    printf( " NASA Ames Research Center\n" );
    printf( " Moffett Field, CA  94035-1000\n\n" );
    printf( " Fax: 650-604-3957\n\n" ); */
}

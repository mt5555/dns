/*
//   FORTRAN is used if the original name does not have an underscore
//   FORTRAN2 is used if the original name does have an underscore
*/
#ifdef F_NO_UNDERSCORE
#define FORTRAN(A) A
#define FORTRAN2(A) A
#else
#ifdef G77_UNDERSCORE
#define FORTRAN(A) A##_
#define FORTRAN2(A) A##__
#else
#define FORTRAN(A) A##_
#define FORTRAN2(A) A##_
#endif
#endif







#include        <stdio.h>
#include        <signal.h>
#include        <errno.h>
#include        <stdlib.h>
#include        <time.h>

int     sig_received=0;


/************************************************************************/
/*                                                                      */
/************************************************************************/
extern  void    ExceptionHandler( int nSignal )
{


        fprintf( stderr, "\nCaught signal %d\n", nSignal );


        /*      Reset the Signal so we can catch it again               */


        if( signal( nSignal, ExceptionHandler ) == SIG_ERR )
        {
                fprintf( stderr, "**WARNING** Can not reset signal %d\n", nSignal );
        }


        /*      Indicate that the signal has been caught        */
        /*if (nSignal==SIGUSR1) sig_received=1;*/
        if (nSignal==SIGURG) sig_received=1;


}

void  FORTRAN2(set_sighandler)(void) {

        if( signal( SIGURG, ExceptionHandler ) == SIG_ERR ) {
                fprintf( stderr, "** WARNING** Can't catch signal SIGURG" );
        }
        if( signal( SIGUSR1, ExceptionHandler ) == SIG_ERR ) {
                fprintf( stderr, "** WARNING** Can't catch signal SIGUSR1" );
        }

}

void FORTRAN2(caught_sig)(int *i) {

*i=sig_received;

}





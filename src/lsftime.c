/*
Copyright 2007.  Los Alamos National Security, LLC. This material was
produced under U.S. Government contract DE-AC52-06NA25396 for Los
Alamos National Laboratory (LANL), which is operated by Los Alamos
National Security, LLC for the U.S. Department of Energy. The
U.S. Government has rights to use, reproduce, and distribute this
software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY
FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so
as not to confuse it with the version available from LANL.

Additionally, this program is free software; you can redistribute it
and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version. Accordingly, this
program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.
*/



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





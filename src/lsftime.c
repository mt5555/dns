#if (defined AIX || defined HPUX || defined OSF1)
#define FORTRAN(A) A
#else
#define FORTRAN(A) A##_
#endif





#if 0
To use this in fortran:

program tester

integer iret,lsf_time_remaining,itimeleft
iret=lsf_time_remaining(itimeleft)
if (iret==0) then
   print *,'timeleft = ',itimeleft
else
   print *,'error getting lsf time'
endif






Here is the subroutine we talked about.  Its name is lsf_time_remaining
and is imbedded in the program I call lsfTTRL.  The usage for
lsf_time_remaining is:
ret=lsf_time_remaining(&timeleft);
where: ret = 0 is successful and ret != 0 is failure.
If ret=0, then "timeleft" will be an integer with the time left to
RUNLIMIT in seconds.


Jerry



#include        <stdio.h>
#include        <signal.h>
#include        <errno.h>
#include        <stdlib.h>
#include        <time.h>
#include        <lsf/lsbatch.h>

int     SigURG=0;

/************************************************************************/
/*      Prototypes:                                                     */
/************************************************************************/


extern  void    CleanUp( void );
extern  void    ExceptionHandler( int nSignal );
extern  int     lsf_time_remaining(int *TTRL);




/************************************************************************/
/*      Function:       main                                            */
/************************************************************************/


extern  void    main(int argc,char **argv)
{
        auto    int     i;
        auto    int     timeLeft,ret;


        /*      Make arragements to capture the following excetpions    */
        /*                                                              */
        /*      SIGINT          Signal 2                                */
        /*      SIGTERM         Signal 15                               */
        /*      SIGURG          Signal 21                               */
        /*                                                              */


        if( signal( SIGINT, ExceptionHandler ) == SIG_ERR ) {
                fprintf( stderr, "Can't catch signal SIGINT" );
                exit( 1 );
        }


        if( signal( SIGTERM, ExceptionHandler ) == SIG_ERR ) {
                fprintf( stderr, "Can't catch signal SIGTERM" );
                exit( 1 );
        }


        if( signal( SIGURG, ExceptionHandler ) == SIG_ERR ) {
                fprintf( stderr, "Can't catch signal SIGURG" );
                exit( 1 );
        }


/*      start with SigInt off  */
        SigInt = 0;


        /*  process the data  */
        for(;;) {
                ret = lsf_time_remaining(&timeLeft);
                if( SigInt != 0 ) { break; }
                if ( ret == 0 ) {
                        printf("Time remaining to RUNLIMIT is %d sec \n",timeLeft);
                } else {
                        printf("lsf_time_remaining failed\n");
                }
                sleep(10);
        }


        CleanUp();
        exit(0);
}


/************************************************************************/
/*                                                                      */
/************************************************************************/


extern  void    CleanUp( void )
{
        fprintf( stderr, "\nWe are now terminating ...\n" );
}
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
        sig_received=1;


}

void  FORTRAN(set_sighandler)(void) {

        if( signal( SIGURG, ExceptionHandler ) == SIG_ERR ) {
                fprintf( stderr, "** WARNING** Can't catch signal SIGURG" );
        }
        if( signal( SIGUSR1, ExceptionHandler ) == SIG_ERR ) {
                fprintf( stderr, "** WARNING** Can't catch signal SIGUSR1" );
        }

}

void FORTRAN(caught_sig)(int *i) {

*i=sig_received;

}




#ifdef USE_LSFTIME
#include        <lsf/lsbatch.h>

/************************************************************************/
/*            lsf_time_remaining                                        */
/*                                                                      */
/* usage: ret=lsf_time_remaiing(&timeLeft)                              */
/*    On return, timeLeft is the time left                              */
/*               in seconds until RUNLIMIT.                             */
/*    ret = 0 signals successful execution                              */
/*        > 0 signals unsuccessful execution                            */
/*                                                                      */
/* Written by Jerry Melendez, CIC-7                                     */
/************************************************************************/

extern  int     FORTRAN(lsf_time_remaining)(int *TTRL) {
    /* typedef enum    { FALSE = 0, TRUE = 1 } Bool_t;*/
typedef struct {
                int  initialized;
                time_t  start_time;
                time_t  end_time;
        } lsf_data_t;
static  lsf_data_t lsf_data = {0, NULL, NULL};
auto    char    *lsf_jobid_str;
auto    int     lsf_jobid, more, numJobs;
auto    struct  jobInfoEnt      *job;
auto    time_t  Ctime;

        if (! lsf_data.initialized ) {
                lsf_jobid_str = getenv("LSB_JOBID");
                if ( lsf_jobid_str == NULL ) { return 1; };
                lsf_jobid = (int)strtol(lsf_jobid_str,NULL,10);


                if (lsb_init(NULL) < 0) {
                        return 1;
                }


                numJobs = lsb_openjobinfo(lsf_jobid, NULL, NULL, NULL, NULL, 0);
                if (numJobs < 0) {
                        lsb_closejobinfo();
                        return 1;
                }
                job = lsb_readjobinfo(&more);
                lsf_data.start_time = job->startTime;
                lsf_data.end_time = lsf_data.start_time +
                                (time_t) job->submit.rLimits[LSF_RLIMIT_RUN];
                lsf_data.initialized=1;
                lsb_closejobinfo();
        }


        time(&Ctime);
        *TTRL = (int) (lsf_data.end_time - Ctime);
        return 0;
}


#else

extern  int     FORTRAN(lsf_time_remaining)(int *TTRL){
return 1;

}

#endif

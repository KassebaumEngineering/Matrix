//  TimeValue.C
//
//  C++ Implementation of struct timeval operators
//
//  ----------
//  See TimeValue.H for the Interface description and Usage directions.
//  ----------
//
// Notes: 
//
//  This file and it's contents are my unique creation, and you
//  can't have it.  I retain all copy rights.  I may from time to time 
//  grant license to others to use it, but I retain ownership of
//  MY software.  
//
// History:
/*    $Log: TimeValue.C,v $
/*    Revision 1.4  1993/11/24 03:45:02  jak
/*    Added a division operator for time calculation for benchmarking
/*    purposes.   -jak
/*
 * Revision 1.3  1993/11/20  21:53:21  jak
 * Fixed a bug in the Linked_List_Template to allow it to be correctly
 * included and used in a library situation.  -jak
 *
 * Revision 1.2  1993/11/20  06:10:02  jak
 * Bug fixes and optimization turned on.   -jak
 *
 * Revision 1.1  1993/11/20  02:19:45  jak
 * Added Time and resource usage programs.  Also, the class is now
 * built into a library (libMatrix.a).  The Linked_List now has
 * reference counts and is correctly copied and deleted by the new
 * inc and dec ,methods for the reference count.  -jak
 **/
//  ************************************************************
static char rcsid_TimeValue_C[] = "$Id: TimeValue.C,v 1.4 1993/11/24 03:45:02 jak Exp $";

#pragma implementation

#include "TimeValue.H"
#include <sys/time.h>
#include <math.h>

//
//  Constructors and Destructors
//
TimeValue::TimeValue()
{
    struct timeval now;
    gettimeofday( &now, (struct timezone *)0);

    seconds  = now.tv_sec;
    micro_seconds = now.tv_usec;
};

TimeValue::TimeValue( const TimeValue& atime )
{
    seconds  = atime.seconds;
    micro_seconds = atime.micro_seconds;
};

TimeValue::TimeValue( struct timeval & atime )
{
    seconds  = atime.tv_sec;
    micro_seconds = atime.tv_usec;
};

TimeValue::TimeValue( float time ) // seconds
{
    seconds = (long) floor( fabs( (double) time ));
    if( time < 0.0 ) seconds = 0 - seconds;
    micro_seconds = (long)((time - (float)seconds)*1000000.0);
};

TimeValue::TimeValue( double time) // seconds
{
    seconds = (long) floor( fabs( time ));
    if( time < 0.0 ) seconds = 0 - seconds;
    micro_seconds = (long)( (time - (double)seconds) * 1000000.0 );
};

TimeValue::TimeValue( long secs, long mseconds )
{
    seconds  = secs;
    micro_seconds = mseconds;
};

//
//  Destructor
//
TimeValue::~TimeValue()
{
};

//
//  Member Operations
//
TimeValue& TimeValue::operator=(const TimeValue& atimeval)
{
    seconds  = atimeval.seconds;
    micro_seconds = atimeval.micro_seconds;

    return *this;
};

// -------------------------------------------------
// Utility operators
//

TimeValue operator+( const TimeValue& t1, const TimeValue& t2 )
{
    double result, op2;
    //TimeValue result;

    //result.seconds  = t1.secs() + t2.secs();
    //result.micro_seconds  = t1.msecs() + t2.msecs();
    //while( result.micro_seconds >= 1000000 ){
    //    result.micro_seconds -= 1000000;
    //    result.seconds += 1;
    //};

    result = (double) t1.secs() + ((double) t1.msecs() / 1000000.0 );
    op2    = (double) t2.secs() + ((double) t2.msecs() / 1000000.0 );

    result += op2;

    return TimeValue( result );
};

TimeValue  operator*(const TimeValue& t, double scale )
{
    double result;

    result = (double) t.seconds + ((double) t.micro_seconds / 1000000.0 );
    return TimeValue( result*scale );
};

TimeValue operator-( const TimeValue& t1, const TimeValue& t2 )
{
    double result, op2;

    result = (double) t1.secs() + ((double) t1.msecs()) / 1000000.0 ;
    op2    = (double) t2.secs() + ((double) t2.msecs()) / 1000000.0 ;

    result -= op2;

    return TimeValue( result );
};

int operator>=( const TimeValue& t1, const TimeValue& t2 )
{
    long seconds, micro_seconds;

    seconds =  t1.seconds - t2.seconds;
    micro_seconds =  t1.micro_seconds - t2.micro_seconds;

    if ( seconds < 0 ){
        return 0;
    } else if ( seconds == 0){
        if( micro_seconds < 0)
            return 0;
        else
            return 1;
    } else /* seconds > 0 */ {
        return 1;
    }
};

int operator<=( const TimeValue& t1, const TimeValue& t2 )
{
    long seconds, micro_seconds;

    seconds =  t1.seconds - t2.seconds;
    micro_seconds =  t1.micro_seconds - t2.micro_seconds;

    if ( seconds > 0 ){
        return 0;
    } else if ( seconds == 0){
        if( micro_seconds > 0)
            return 0;
        else
            return 1;
    } else /* seconds < 0 */ {
        return 1;
    }
};

TimeValue& time_abs( TimeValue& atime )
{
    if( !(atime >= TimeValue( 0,0 )) ){
        return TimeValue( -atime.seconds, -atime.micro_seconds);
    } else
        return atime;
};

ostream & operator << (ostream &cbuf, const TimeValue& atime)
{
    double time;
    time = ((double) atime.secs()) + ((double) atime.msecs())/1000000.0 ;
    cbuf << time << " seconds ";
    return cbuf;
};


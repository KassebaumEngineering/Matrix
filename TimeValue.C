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
    Revision 1.6  1993/11/30 05:36:37  jak
    A Small fix to the operator /()  -jak
    
 * Revision 1.5  1993/11/30  05:24:46  jak
 * Changed the implementation of the TimeValue Class.  -jak
 *
 * Revision 1.4  1993/11/24  03:45:02  jak
 * Added a division operator for time calculation for benchmarking
 * purposes.   -jak
 *
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
static char rcsid_TimeValue_C[] = "$Id: TimeValue.C,v 1.6 1993/11/30 05:36:37 jak Exp $";

// #pragma implementation // Removed obsolete pragma

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

    the_time = (double) now.tv_sec + ((double) now.tv_usec / 1000000.0 );
};

TimeValue::TimeValue( const TimeValue& atime )
{
    the_time = atime.the_time;
};

TimeValue::TimeValue( struct timeval & atime )
{
    the_time = (double) atime.tv_sec + ((double) atime.tv_usec / 1000000.0 );
};

TimeValue::TimeValue( float time ) // seconds
{
    the_time = time;
};

TimeValue::TimeValue( double time) // seconds
{
    the_time = time;
};

TimeValue::TimeValue( long secs, long mseconds )
{
    the_time = (double) secs + ((double) mseconds / 1000000.0 );
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
    the_time = atimeval.the_time;

    return *this;
};

// -------------------------------------------------
// Utility operators
//

TimeValue operator+( const TimeValue& t1, const TimeValue& t2 )
{
    double result;
    result = t1.the_time + t2.the_time;

    return TimeValue( result );
};

TimeValue  operator*(const TimeValue& t, double scale )
{
    double result;

    result = t.the_time * scale;
    return TimeValue( result );
};

double  TimeValue::operator/( const TimeValue& t2 )
{
    double result;

    result = the_time / t2.the_time;
    return result;
};

TimeValue  operator/(const TimeValue& t1, double scale )
{
    double result;

    result = t1.the_time / scale;
    return TimeValue( result );
};

TimeValue operator-( const TimeValue& t1, const TimeValue& t2 )
{
    double result;

    result = t1.the_time - t2.the_time;
    return TimeValue( result );
};

int operator>=( const TimeValue& t1, const TimeValue& t2 )
{
    return ( t1.the_time >= t2.the_time );
};

int operator<=( const TimeValue& t1, const TimeValue& t2 )
{
    return ( t1.the_time <= t2.the_time );
};

TimeValue time_abs( TimeValue& atime ) // Return by value
{
    return TimeValue( fabs(atime.the_time) );
};

std::ostream & operator << (std::ostream &cbuf, const TimeValue& atime) // Qualify ostream
{
    cbuf << atime.the_time << " seconds ";
    return cbuf;
};


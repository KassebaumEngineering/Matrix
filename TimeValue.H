//  TimeValue.H
//  C++ Interface Definition for struct timeval operators
//
//  *************************  Description  **************************
//
//  Author(s) and Credits:
//  John Kassebaum.
//
//  This file and it's contents are my unique creation, and you
//  can't have it.  I retain all copy rights.  I may from time to time 
//  grant license to others to use it, but I retain ownership of
//  MY software.  
//
//  History:
/*    $Log: TimeValue.H,v $
/*    Revision 1.3  1993/11/20 21:53:22  jak
/*    Fixed a bug in the Linked_List_Template to allow it to be correctly
/*    included and used in a library situation.  -jak
/*
 * Revision 1.2  1993/11/20  06:10:03  jak
 * Bug fixes and optimization turned on.   -jak
 *
 * Revision 1.1  1993/11/20  02:19:47  jak
 * Added Time and resource usage programs.  Also, the class is now
 * built into a library (libMatrix.a).  The Linked_List now has
 * reference counts and is correctly copied and deleted by the new
 * inc and dec ,methods for the reference count.  -jak
 **/
//  ******************************************************************
#ifndef _TimeValue_H
#define _TimeValue_H

static char rcsid_TimeValue_H[] = "$Id: TimeValue.H,v 1.3 1993/11/20 21:53:22 jak Exp $";

#pragma interface 

#include <iostream.h>
#include <sys/types.h>
#include <sys/time.h>
//FYI-   struct timeval { long tv_sec;  long tv_usec; };

class TimeValue {
public:
    long    seconds;       // seconds since 1970
    long    micro_seconds; // in micro-seconds

//  Constructors and Destructors
    TimeValue( void );
    TimeValue( const TimeValue& );
    TimeValue( struct timeval & );
    TimeValue( float ); // seconds
    TimeValue( double ); // seconds
    TimeValue( long /*seconds*/, long = 0 /*micro-seconds*/);
    ~TimeValue( void );

    inline long msecs( void ) const {
        return micro_seconds;
    };
    inline long secs( void ) const {
        return seconds;
    };
    TimeValue&  operator=( const TimeValue& );

    friend TimeValue   operator+ ( const TimeValue&, const TimeValue&);
    friend TimeValue   operator- ( const TimeValue&, const TimeValue&);
    friend int         operator>=( const TimeValue&, const TimeValue&);
    friend int         operator<=( const TimeValue&, const TimeValue&);

    friend TimeValue&   time_abs ( TimeValue& );

    friend ostream & operator << (ostream &, const TimeValue&);

};

#endif
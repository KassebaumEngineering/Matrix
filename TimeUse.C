// 
// Implementation of TimeUse Usage Class
//
//  Author: John Kassebaum
//  
// 
//  This file and it's contents are my unique creation, and you
//  can't have it.  I retain all copy rights.  I may from time to time 
//  grant license to others to use it, but I retain ownership of
//  MY software. 
//
//  Revision:
/*  $Id: TimeUse.C,v 1.3 1993/11/20 21:53:18 jak Exp $
*/
//  History:
/*  $Log: TimeUse.C,v $
/*  Revision 1.3  1993/11/20 21:53:18  jak
/*  Fixed a bug in the Linked_List_Template to allow it to be correctly
/*  included and used in a library situation.  -jak
/*
 * Revision 1.2  1993/11/20  06:09:59  jak
 * Bug fixes and optimization turned on.   -jak
 *
 * Revision 1.1  1993/11/20  02:19:42  jak
 * Added Time and resource usage programs.  Also, the class is now
 * built into a library (libMatrix.a).  The Linked_List now has
 * reference counts and is correctly copied and deleted by the new
 * inc and dec ,methods for the reference count.  -jak
 **/
// =====================================

static char rcsid_TimeUse_C[] =  "$Id: TimeUse.C,v 1.3 1993/11/20 21:53:18 jak Exp $";

#pragma implementation

#include "TimeUse.H"
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>

extern void Abort( char * );

TimeUse:: TimeUse( int ident )
{
    start( ident );
};

TimeUse:: ~TimeUse( )
{
};
	
TimeValue TimeUse:: start( int ident )
{
    struct rusage rusg;
	TimeValue *rtn;
	
	if ( getrusage( RUSAGE_SELF, &rusg ) == -1) {
		Abort("TimeUse:: start( int ): Could not get TimeUse Info!");
    }
	
 // new TimeValue is eventually deleted by milestones list
	rtn = new TimeValue( rusg.ru_utime ); 

	milestones.add( rtn , ident );
	
	return *rtn;
};

TimeValue TimeUse:: restart( int ident )
{
    stop( ident );
	return start( ident );
};

TimeValue TimeUse:: look_at( int ident )
{
    struct rusage rusg;
    TimeValue now, *then;
	
	if ( getrusage( RUSAGE_SELF, &rusg ) == -1){
		Abort("TimeUse:: start( int ): Could not get TimeUse Info!");
    }
	
	now = TimeValue( rusg.ru_utime ); 
	then = milestones.find( ident );


    now = now - (*then);

    return now;
};

TimeValue TimeUse:: stop( int ident )
{
    TimeValue rtn;
	
	rtn = look_at(ident);
    milestones.del( ident );
	
	return rtn;
};

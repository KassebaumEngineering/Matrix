//
//  Interface to TimeUse Usage Class
//
//  Author: John Kassebaum
//  
//  This file and it's contents are my unique creation, and you
//  can't have it.  I retain all copy rights.  I may from time to time 
//  grant license to others to use it, but I retain ownership of
//  MY software.  
//
//  Revision:
/*  $Id: TimeUse.H,v 1.1 1993/11/20 02:19:44 jak Exp $
 */
//  History:
/*  $Log: TimeUse.H,v $
/*  Revision 1.1  1993/11/20 02:19:44  jak
/*  Added Time and resource usage programs.  Also, the class is now
/*  built into a library (libMatrix.a).  The Linked_List now has
/*  reference counts and is correctly copied and deleted by the new
/*  inc and dec ,methods for the reference count.  -jak
/**/
// =====================================

#ifndef _TimeUse_H
#define _TimeUse_H
static char rcsid_TimeUse_H[] =  "$Id: TimeUse.H,v 1.1 1993/11/20 02:19:44 jak Exp $";

#include "TimeValue.H"
#include "Linked_List_Template.H"

class TimeUse {
private:
    Linked_List<TimeValue>	milestones;

public:
    TimeUse( int = 0 );
    ~TimeUse( void );
	
	TimeValue start( int ident );
	TimeValue restart( int ident );
	TimeValue look_at( int ident );
	TimeValue stop( int ident );
	
};


#endif
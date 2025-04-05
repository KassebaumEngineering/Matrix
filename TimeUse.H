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
/*  $Id: TimeUse.H,v 1.4 1993/11/27 00:20:28 jak Exp $
 */
//  History:
/*  $Log: TimeUse.H,v $
 *  Revision 1.4  1993/11/27 00:20:28  jak
 *  Matrix Class has been ported for use with the AT&T cfront compiler version 3
 *  (with templates).   -jak
 *
 * Revision 1.3  1993/11/20  21:53:20  jak
 * Fixed a bug in the Linked_List_Template to allow it to be correctly
 * included and used in a library situation.  -jak
 *
 * Revision 1.2  1993/11/20  06:10:01  jak
 * Bug fixes and optimization turned on.   -jak
 *
 * Revision 1.1  1993/11/20  02:19:44  jak
 * Added Time and resource usage programs.  Also, the class is now
 * built into a library (libMatrix.a).  The Linked_List now has
 * reference counts and is correctly copied and deleted by the new
 * inc and dec ,methods for the reference count.  -jak
 **/
// =====================================

#ifndef _TimeUse_H
#define _TimeUse_H
static char rcsid_TimeUse_H[] =  "$Id: TimeUse.H,v 1.4 1993/11/27 00:20:28 jak Exp $";


#include "TimeValue.H"
#include "Linked_List_Template.H"

class TimeUse {
private:
    Linked_List<TimeValue>	milestones;

public:
    TimeUse( int = 0 );
    ~TimeUse( void );
	
	TimeValue start( int = 0 );
	TimeValue restart( int  = 0  );
	TimeValue look_at( int  = 0  );
	TimeValue stop( int  = 0  );
	
};


#endif


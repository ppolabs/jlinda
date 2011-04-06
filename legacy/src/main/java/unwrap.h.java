public class GlobalMembersUnwrap
{





	// ______ Use stanford software through unix calls _____
//C++ TO JAVA CONVERTER TODO TASK: The implementation of the following method could not be found:
	//void unwraptreeframon(input_gen generalinput, input_unwrap unwrapinput, productinfo interferogram);

	// ______ Use snaphu software through unix calls _____
//C++ TO JAVA CONVERTER TODO TASK: The implementation of the following method could not be found:
	//void snaphu_unwrap(input_gen generalinput, input_unwrap unwrapinput, productinfo interferogram, slcimage master, slcimage slave, RefObject<orbit> masterorbit, RefObject<orbit> slaveorbit, input_ell ellips);

	//#endif // UNWRAP_H
}
//
// * Copyright (c) 1999-2009 Delft University of Technology, The Netherlands
// *
// * This file is part of Doris, the Delft o-o radar interferometric software.
// *
// * Doris program is free software; you can redistribute it and/or modify
// * it under the terms of the GNU General Public License as published by
// * the Free Software Foundation; either version 2 of the License, or
// * (at your option) any later version.
// *
// * Doris is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// * GNU General Public License for more details.
// *
// * You should have received a copy of the GNU General Public License
// * along with this program; if not, write to the Free Software
// * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// *
// * Publications that contain results produced by the Doris software should
// * contain an acknowledgment. (For example: The interferometric processing
// * was performed using the freely available Doris software package developed
// * by the Delft Institute of Earth Observation and Space Systems (DEOS), Delft
// * University of Technology, or include a reference to: Bert Kampes and
// * Stefania Usai. \"Doris: The Delft Object-oriented Radar Interferometric
// * software.\" In: proceedings 2nd ITC ORS symposium, August 1999. (cdrom)).
// *
// 
//***************************************************************
// * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/unwrap.hh,v $
// * $Revision: 3.8 $
// * $Date: 2005/08/24 10:03:18 $
// * $Author: kampes $
// *
// * Declaration of routines for unwrapping 
// ***************************************************************



//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
//#if ! UNWRAP_H
//#define UNWRAP_H


// Jia defined this for compilation under windows
// Bert Kampes, 24-Aug-2005
//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
//#if _MSC_VER > 1000
//#endif // _MSC_VER > 1000




//----------------------------------------------------------------------------------------
//	Copyright © 2006 - 2008 Tangible Software Solutions Inc.
//
//	This class is used to simulate the ability to pass arguments by reference in Java.
//----------------------------------------------------------------------------------------
final class RefObject<T>
{
	T argvalue;
	RefObject(T refarg)
	{
		argvalue = refarg;
	}
}
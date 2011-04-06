public class GlobalMembersExceptions
{
	//
	// @file   excpetions.cc exception handling for Doris InSAR processor
	// @brief  exception handling for Doris InSAR processor
	//
	//
	// * Copyright (c) 1999-2005 Bert Kampes
	// * Copyright (c) 1999-2005 Delft University of Technology, The Netherlands
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




	// ====== Globals to throw everywhere, e.g., throw(some_error) ======
	public static SOME_ERROR some_error; // can be thrown from all programs
	public static INPUT_ERROR input_error; // can be thrown from all programs
	public static FILE_ERROR file_error; // can be thrown from all programs
	public static MEMORY_ERROR memory_error; // can be thrown from all programs
	public static UNHANDLED_CASE_ERROR unhandled_case_error; // can be thrown from all programs
	public static ARGUMENT_ERROR argument_error; // can be thrown from all programs
	public static KEYWORD_ERROR keyword_error; // can be thrown from all programs
	public static USAGE_ERROR usage_error; // can be thrown from all programs




	//********************************************************************
	// * @brief exception handler for floating point exception
	// ********************************************************************
	// function: CatchSignals()
	// * ------------------------
	// * Traps common signals that by default cause the program to abort.  
	// * Sets (pointer to function) Handler as the signal handler for all.
	// * Note that SIGKILL usually cannot be caught.  No return value.
	//C++ TO JAVA CONVERTER NOTE: Beginning of line comments are not maintained by C++ to Java Converter
	//ORIGINAL LINE: */   // following is from snaphu:
	 // following is from snaphu:
//C++ TO JAVA CONVERTER TODO TASK: There are no simple equivalents to function pointers in Java:
	void CatchSignals(void (*SigHandler)(int))
	{
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if ! WIN32
	  signal(SIGHUP,SigHandler); // Hang-up signal
	  signal(SIGQUIT,SigHandler);
	  signal(SIGPIPE,SigHandler);
	  signal(SIGALRM,SigHandler);
	  signal(SIGBUS,SigHandler);
	//#endif
	  signal(SIGINT,SigHandler);
	  signal(SIGILL,SigHandler);
	  signal(SIGABRT,SigHandler);
	  signal(SIGFPE,SigHandler); // floating point exception
	  signal(SIGSEGV,SigHandler); // segmentation fault: introduces when compiled with -O in gcc4?
	  signal(SIGTERM,SigHandler);
	}



	// following is based on snaphu code: but needs some work.
	// function: SetDump()
	// * -------------------
	// * Set the global variable dumpresults_global to TRUE if SIGINT or SIGHUP
	// * signals recieved.  Also sets requestedstop_global if SIGINT signal 
	// * received.  This function should only be called via signal() when 
	// * a signal is caught.
	// 
	public static void handle_signal(int signum)
	  {
	  switch (signum)
		{
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if ! WIN32
		case SIGHUP:
		  System.out.print("Caught SIGHUP: Hang-up signal.");
		  System.out.print("\n");
		  cerr << "Caught SIGHUP: Hang-up signal." << "\n";
		  break;
		case SIGQUIT:
		  System.out.print("Caught SIGQUIT: Quit signal.");
		  System.out.print("\n");
		  cerr << "Caught SIGQUIT: Quit signal." << "\n";
		  break;
		case SIGPIPE:
		  System.out.print("Caught SIGPIPE: ? signal.");
		  System.out.print("\n");
		  cerr << "Caught SIGPIPE: ? signal." << "\n";
		  break;
		case SIGALRM:
		  System.out.print("Caught SIGALRM: Alarm signal.");
		  System.out.print("\n");
		  cerr << "Caught SIGALRM: Alarm signal." << "\n";
		  break;
		case SIGBUS:
		  System.out.print("Caught SIGBUS: Bus error (accessing memory incorrectly)?");
		  System.out.print("\n");
		  cerr << "Caught SIGBUS: Bus error (accessing memory incorrectly)?" << "\n";
		  break;
	//#endif
		case SIGINT:
		  System.out.print("Caught SIGINT: User interupt signal.");
		  System.out.print("\n");
		  cerr << "Caught SIGINT: User interupt signal." << "\n";
		  exit(1);
		  break;
		case SIGFPE:
		  System.out.print("Caught SIGFPE: floating point exception, zero division, etc.");
		  System.out.print("\n");
		  cerr << "Caught SIGFPE: floating point exception, zero division, etc." << "\n";
		  break;
		case SIGILL:
		  System.out.print("Caught SIGILL: ? signal.");
		  System.out.print("\n");
		  cerr << "Caught SIGILL: ? signal." << "\n";
		  break;
		case SIGABRT:
		  System.out.print("Caught SIGABRT: Abort signal.");
		  System.out.print("\n");
		  cerr << "Caught SIGABRT: Abort signal." << "\n";
		  break;
		case SIGSEGV:
		  System.out.print("Caught SIGSEGV: Segmentation fault.");
		  System.out.print("\n");
		  cerr << "Caught SIGSEGV: Segmentation fault." << "\n";
		  exit(1);
		  break;
		case SIGTERM:
		  System.out.print("Caught SIGTERM: ? signal.");
		  System.out.print("\n");
		  cerr << "Caught SIGTERM: ? signal." << "\n";
		  break;
		default:
		  System.out.print("Caught an unknown signal.  Signum = ");
		  System.out.print(signum);
		  System.out.print("\n");
		  cerr << "Caught an unknown signal.  Signum = " << signum << "\n";
		}
	  }
}





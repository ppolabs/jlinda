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
// * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/productinfo.cc,v $        *
// * $Revision: 3.13 $                                            *
// * $Date: 2005/08/01 15:51:19 $                                 *
// * $Author: kampes $                                            *
// *                                                              *
// * implementation of product info class.                        *
// * - data filling/updating.                                     * 
// * - reading in matrices.                                       *
// * - etc.                                                       *
// #%// BK 01-Sep-2000
// ***************************************************************





//***************************************************************
// *    productinfo::fillproductinfo                              *
// * Fills data of productinfo object, read from resultfile after *
// * steps that produced a product, identified by 'iden'.         *
// * 13 lines are checked for info, first occurence counts.       *
// * exit if 8 are found, or 15 lines are checked.                *
// *  and read information upto string ":_NORMAL" (i.e., only read  *
// *  appropriate section                                         *
// *                                                              *
// *  "Data_output_file:"                                         *
// *  "Data_output_format:"                                       *
// *  "Multilookfactor_azimuth_direction:"                        *
// *  "Multilookfactor_range_direction:"                          *
// *  "First_line"  (w.r.t. original):                            *
// *  "Last_line"  (w.r.t. original):                             *
// *  "First_pixel"  (w.r.t. original):                           *
// *  "Last_pixel"  (w.r.t. original):                            *
// *                                                              *
// * input:                                                       *
// *  - resultfilename                                            *
// *  - identifier                                                *
// * output:                                                      *
// *  - (updated) productinfo                                     *
// *                                                              *
// *    Bert Kampes, 22-Dec-1998                                  *
// *  Dont know why found1 is used more then ones?                *
// *  BK 13-april-2000                                            *
// #%// BK 01-Sep-2000 class implementation                       *
// * added check for filesize                                     *
// #%// BK 08-Mar-2001                                            *
// * added stop after ":_NORMAL" found                            *
// #%// Bert Kampes, 01-Aug-2005
// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void productinfo::fillproductinfo(String resultfile, String iden)


//***************************************************************
// * productinfo::readphase                                       *
// *  read data from file in a real4 matrix                       *
// *  phase is extracted from hgt, not unwrapped set to NaN       *
// *  real4 is read 'as is'.                                      *
// *  complex real4: angle is taken.                              *
// *  hgt: phase is read, set to NaN if amplitude==0              *
// * cutout window in own system starting at line 1,              *
// * (no offset win)                                              *
// * This function is useful for reading products: complex interf.*
// * unwrapped interf. and to get phase for s2h etc.              *
// #%// BK 22-Sep-2000                                            *
// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<real4> productinfo::readphase(window cutout) const
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
 // window to be read

//***************************************************************
// * productinfo::readdata                                        *
// *  read data from file in a complex real4 matrix               *
// *  complex real4                                               *
// * cutout window in own system starting at line 1,              *
// * (no offset win)                                              *
// * This function is useful for reading products: complex interf.*
// * unwrapped interf. and to get phase for s2h etc.              *
// #%// BK 22-Sep-2000 (readphase)                                *
// #%// Batu 30-JUL-2007                                          *
// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<complr4> productinfo::readdata(window cutout) const
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
 // window to be read

//***************************************************************
// * productinfo::readdatar4                                      *
// *  read data from file in a complex real4 matrix               *
// *  complex real4                                               *
// * cutout window in own system starting at line 1,              *
// * (no offset win)                                              *
// * This function is useful for reading products: complex interf.*
// * unwrapped interf. and to get phase for s2h etc.              *
// #%// BK 22-Sep-2000 (readphase)                                *
// #%// Batu 30-JUL-2007                                          *
// #%// MA   19-NOV-2008 [TODO] template                          *
// ***************************************************************

//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<real4> productinfo::readdatar4(window cutout) const
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
 // window to be read
//----------------------------------------------------------------------------------------
//	Copyright © 2006 - 2008 Tangible Software Solutions Inc.
//
//	This class provides the ability to simulate various classic C string functions
//	which don't have exact equivalents in the Java framework.
//----------------------------------------------------------------------------------------
final class SimulateStringFunctions
{
	//------------------------------------------------------------------------------------
	//	This method simulates the classic C string function 'isxdigit' (and 'iswxdigit').
	//------------------------------------------------------------------------------------
	static boolean isXDigit(char character)
	{
		if (Character.isDigit(character))
			return true;
		else if ("ABCDEFabcdef".indexOf(character) > -1)
			return true;
		else
			return false;
	}

	//------------------------------------------------------------------------------------
	//	This method simulates the classic C string function 'strchr' (and 'wcschr').
	//------------------------------------------------------------------------------------
	static String strChr(String stringtosearch, char chartofind)
	{
		int index = stringtosearch.indexOf(chartofind);
		if (index > -1)
			return stringtosearch.substring(index);
		else
			return null;
	}

	//------------------------------------------------------------------------------------
	//	This method simulates the classic C string function 'strrchr' (and 'wcsrchr').
	//------------------------------------------------------------------------------------
	static String strRChr(String stringtosearch, char chartofind)
	{
		int index = stringtosearch.lastIndexOf(chartofind);
		if (index > -1)
			return stringtosearch.substring(index);
		else
			return null;
	}

	//------------------------------------------------------------------------------------
	//	This method simulates the classic C string function 'strstr' (and 'wcsstr').
	//------------------------------------------------------------------------------------
	static String strStr(String stringtosearch, String stringtofind)
	{
		int index = stringtosearch.indexOf(stringtofind);
		if (index > -1)
			return stringtosearch.substring(index);
		else
			return null;
	}

	//------------------------------------------------------------------------------------
	//	This method simulates the classic C string function 'strtok' (and 'wcstok').
	//------------------------------------------------------------------------------------
	private static String activestring;
	private static int activeposition;
	static String strTok(String stringtotokenize, String delimiters)
	{
		if (stringtotokenize != null)
		{
			activestring = stringtotokenize;
			activeposition = -1;
		}

		//the stringtotokenize was never set:
		if (activestring == null)
			return null;

		//all tokens have already been extracted:
		if (activeposition == activestring.length())
			return null;

		//bypass delimiters:
		activeposition++;
		while (activeposition < activestring.length() && delimiters.indexOf(activestring[activeposition]) > -1)
		{
			activeposition++;
		}

		//only delimiters were left, so return null:
		if (activeposition == activestring.length())
			return null;

		//get starting position of string to return:
		int startingposition = activeposition;

		//read until next delimiter:
		do
		{
			activeposition++;
		} while (activeposition < activestring.length() && delimiters.indexOf(activestring[activeposition]) == -1);

		return activestring.substring(startingposition, activeposition - startingposition);
	}
}

public class productinfo
{
	public void fillproductinfo(String resultfile, String iden)
	  {
	  TRACE_FUNCTION("fillproductinfo (BK 01-Sep-2000)")
	  String word =" ";
	  String dummyline = new String(new char[ONE27]);
	
	  // ______Open file______
	  ifstream resfile = new ifstream(resultfile, ios.in);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(resfile,resultfile,__FILE__,__LINE__);
	
	  boolean foundiden = false;
	  while (resfile != null)
		{
		if (strcmp(word,iden)) // Lookfor identifier
		  {
		  resfile.getline(dummyline,ONE27,'\n'); // next line
		  resfile >> word; // read word
		  }
		else
		  {
		  foundiden = true;
		  break;
		  }
		}
	
	// ______Check if section has been found______
	  if (!foundiden)
		{
		ERROR << "(fillproductinfo). identifier: \"" << iden << "\" not found in file: " << resultfile;
		PRINT_ERROR(ERROR.get_str())
		throw(file_error);
		}
	
	// ======Extract info (file name/format and window size)======
	  int32 numlineschecked = 0; // number of lines checked after iden
	  boolean foundfilename = false;
	  boolean foundfileformat = false;
	  boolean foundfirstline = false;
	  boolean foundlastline = false;
	  boolean foundfirstpixel = false;
	  boolean foundlastpixel = false;
	  boolean foundmlL = false;
	  boolean foundmlP = false;
	
	  while (resfile != null) // i.e. rest of file
		{
		resfile.getline(dummyline,ONE27,'\n'); // next line
		TRACE << "dummyline: " << dummyline;
		TRACE.print();
		// ___ Check if we are at end of section already ____
		String pch;
		pch = SimulateStringFunctions.strStr(dummyline,":_NORMAL"); // section terminator
		if (!pch.equals(null))
		  {
		  DEBUG.print("Section terminator found (string \":_NORMAL\").");
		  break; // section ends
		  }
		// ___ Check if all parameters are found ______
		if (foundfilename && foundfileformat && foundfirstline && foundlastline && foundfirstpixel && foundlastpixel && foundmlL && foundmlP)
			break;
		//if (numlineschecked == 13) break;
		// report that 13 was too little for some output sections
		// so increased it to 15, BK 07-05-2002
		if (numlineschecked == 15)
			break;
		numlineschecked++;
	
		// ___ Check for strings with parameters ___
		resfile >> word; // read word
		if (!strcmp(word,"Data_output_file:"))
		  {
		  if (!foundfilename)
			{
			foundfilename = true;
			resfile >> file;
			DEBUG << "String: \"Data_output_file:\", \t\t\tvalue: " << file;
			DEBUG.print();
			}
		  else
			WARNING.print("String: \"Data_output_file:\" found again (ignored).");
		  }
		else if (!strcmp(word,"Data_output_format:"))
		  {
		  if (!foundfileformat)
			{
			foundfileformat = true;
			resfile >> word;
			if (!strcmp(word,"complex_short"))
			  formatflag = FORMATCI2;
			else if (!strcmp(word,"complex_real4"))
			  formatflag = FORMATCR4;
			else if (!strcmp(word,"real4"))
			  formatflag = FORMATR4;
			else if (!strcmp(word,"hgt"))
			  formatflag = FORMATHGT;
			else
			  {
			  PRINT_ERROR("wrong format specifier (impossible?)")
			  throw(unhandled_case_error);
			  }
			DEBUG << "String: \"Data_output_format:\", \t\tvalue: " << formatflag;
			DEBUG.print();
			}
		  else
			WARNING.print("String: \"Data_output_format:\" found again (ignored).");
		  }
	
		else if (!strcmp(word,"Multilookfactor_azimuth_direction:"))
		  {
		  if (!foundmlL)
			{
			foundmlL = true;
			resfile >> multilookL;
			DEBUG << "String: \"Multilookfactor_azimuth_direction:\", value: " << multilookL;
			DEBUG.print();
			}
		  else
			WARNING.print("String: \"Multilookfactor_azimuth_direction:\" found again (ignored).");
		  }
		else if (!strcmp(word,"Multilookfactor_range_direction:"))
		  {
		  if (!foundmlP)
			{
			foundmlP = true;
			resfile >> multilookP;
			DEBUG << "String: \"Multilookfactor_range_direction:\", \tvalue: " << multilookP;
			DEBUG.print();
			}
		  else
			WARNING.print("String: \"Multilookfactor_range_direction:\" found again (ignored).");
		  }
	
		else if (!strcmp(word,"First_line")) // (w.r.t. original):
		  {
		  if (!foundfirstline)
			{
			foundfirstline = true;
			resfile >> word >> word >> win.linelo;
			DEBUG << "String: \"First_line:\", \t\t\tvalue: " << win.linelo;
			DEBUG.print();
			}
		  else
			WARNING.print("String: \"First_line:\" found again (ignored).");
		  }
	
		else if (!strcmp(word,"Last_line")) // (w.r.t. original):
		  {
		  if (!foundlastline)
			{
			foundlastline = true;
			resfile >> word >> word >> win.linehi;
			DEBUG << "String: \"Last_line:\", \t\t\tvalue: " << win.linehi;
			DEBUG.print();
			}
		  else
			WARNING.print("String: \"Last_line:\" found again (ignored).");
		  }
	
		else if (!strcmp(word,"First_pixel")) // (w.r.t. original):
		  {
		  if (!foundfirstpixel)
			{
			foundfirstpixel = true;
			resfile >> word >> word >> win.pixlo;
			DEBUG << "String: \"First_pixel:\", \t\t\tvalue: " << win.pixlo;
			DEBUG.print();
			}
		  else
			WARNING.print("String: \"First_pixel:\" found again (ignored).");
		  }
	
		else if (!strcmp(word,"Last_pixel")) // (w.r.t. original):
		  {
		  if (!foundlastpixel)
			{
			foundlastpixel = true;
			resfile >> word >> word >> win.pixhi;
			DEBUG << "String: \"Last_pixel:\", \t\t\tvalue: " << win.pixhi;
			DEBUG.print();
			}
		  else
			WARNING.print("String: \"Last_pixel:\" found again (ignored).");
		  }
		} // while resultfile
	
	
	  // ______ Check filesize with format/dimensions (BK 08-Mar-2001) ______
	  // BK 08-Mar-2001
	  //ifstream tmpfile(file, ios::in | ios::nocreate);
	  ifstream tmpfile = new ifstream(file, ios.in);
	  if (tmpfile != null)
		{
		tmpfile.seekg(0,ios.end); // internal filesize, normal one exists if not exists
		// uint filesizetrue  = tmpfile.tellg();
		final streamoff filesizetrue = tmpfile.tellg(); // [MA] file > 4GB support, this fix eliminates wrong warning
		int16 bytesperelem = 4;
		if (formatflag==FORMATCI2)
			bytesperelem =4;
		if (formatflag==FORMATCR4)
			bytesperelem =8;
		if (formatflag==FORMATR4)
			bytesperelem =4;
		if (formatflag==FORMATI2)
			bytesperelem =2;
		if (formatflag==FORMATI2_BIGENDIAN)
			bytesperelem =2;
		if (formatflag==FORMATR8)
			bytesperelem =8;
		if (formatflag==FORMATHGT)
			bytesperelem =8;
		// int32 filesizecomp = int32(win.lines()/multilookL) *
		uint64 filesizecomp = (uint64)(win.lines()/multilookL) * (win.pixels()/multilookP) * bytesperelem;
		DEBUG << "Checking format/dimensions file=" << file;
		DEBUG.print();
		if (filesizecomp != filesizetrue)
		  {
		  WARNING << "File: \'" << file << "\' has wrong fileformat or dimensions" << ": bytesperpix=" << bytesperelem << ", #l=" << win.lines()/multilookL << ", #p=" << win.pixels()/multilookP << "; size_computed=" << filesizecomp << "B" << " v." << " size_ondisk=" << filesizetrue << "B";
		  WARNING.print();
		  }
		else
		  {
		  DEBUG.print("Fileformat and dimensions are checked and ok.");
		  }
		tmpfile.close();
		} // stream ok
	  else
		{
		WARNING << "File: " << file << " does not seem to exist (may not be a problem).";
		WARNING.print();
		}
	
	
	
	  // ______Tidy up______
	  DEBUG.print("");
	  resfile.close();
	
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if __DEBUG
	  DEBUG.print("finished fillproductinfo");
	  DEBUG.print("content of struct:");
	  showdata();
	//#endif
	  } // END fillproductinfo
	public matrix<real4> readphase(window cutout)
	  {
	  TRACE_FUNCTION("productinfo::readphase (BK 22-Sep-2000)")
	  // ______ Check if cutout < #lines/pixels ______
	  final int32 linesondisk = win.lines()/multilookL;
	  final int32 pixelsondisk = win.pixels()/multilookP;
	
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
		if (cutout.linelo < 1)
		  {
		  PRINT_ERROR("readphase: cutout window larger than what's on disk!")
		  throw(input_error);
		  }
		if (cutout.linehi>linesondisk)
		  {
		  PRINT_ERROR("readphase: cutout window larger than what's on disk!")
		  throw(input_error);
		  }
		if (cutout.pixlo < 1)
		  {
		  PRINT_ERROR("readphase: cutout window larger than what's on disk!")
		  throw(input_error);
		  }
		if (cutout.pixhi>pixelsondisk)
		  {
		  PRINT_ERROR("readphase: cutout window larger than what's on disk!")
		  throw(input_error);
		  }
	
	//  #endif
	  // ______ Open file ______
	//#ifdef __NO_IOS_BINARY__
	//  ifstream ifile(file, ios::in);
	//#else
	//  ifstream ifile(file, ios::in | ios::binary);
	//#endif
	  ifstream ifile;
	  openfstream(ifile,file);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(ifile,file,__FILE__,__LINE__);
	
	  matrix<real4> Result = new matrix(cutout.lines(),cutout.pixels());
	
	  // ====== Actual read data ======
	  switch (formatflag)
		{
		case FORMATR4:
		  {
		  // why not use readfile function?
		  DEBUG.print("reading real4 from file");
		  matrix<real4> LINE = new matrix(1,cutout.pixels()); // phase
		  for (int32 line =cutout.linelo; line<=cutout.linehi; ++line)
			{
			uint64 start = (uint64)(cutout.pixlo-1 + pixelsondisk*(line-1)) * sizeof(real4);
			ifile.seekg(start,ios.beg);
			ifile >> LINE;
			Result.setrow(line-cutout.linelo, LINE); // not the fastest way...
															// but cannot write directly
															// to private data ... (bk)
			}
		  break;
		  }
	
		case FORMATCR4:
		  {
		  DEBUG.print("reading complex real4 from file (get phase)");
		  final window dummyoffset = new window(1,99999,1,99999); // only 1's used -> no offset
		  matrix<complr4> TMP = new matrix(cutout.lines(),cutout.pixels());
		  readfile(TMP,file,linesondisk,cutout,dummyoffset);
		  Result = angle(TMP);
		  break;
		  }
	
		// ______ hgt: first amplitude band, then phase, then next line ______
		case FORMATHGT:
		  {
		  DEBUG.print("reading hgt (band interleaved) from file (get phase)");
		  // ______ Read in one line, set to NaN, fill result ______
		  matrix<real4> LINE = new matrix(1,win.pixels()*2); // disk phase/ampl band interleaved
		  for (int32 line =cutout.linelo; line<=cutout.linehi; ++line)
			{
			final int32 start = (line-1) * pixelsondisk * 2 * sizeof(real4);
			ifile.seekg(start,ios.beg); // file pointer to 1st pix of line
			ifile >> LINE;
			for (int32 pix =cutout.pixlo; pix<=cutout.pixhi; ++pix)
			  {
			  // ______ Check amplitude defined as 0 if not unwrapped ok. ______
			  if (LINE(0,pix-1) == 0.) // ampl. band
				Result(line-cutout.linelo,pix-cutout.pixlo) = NaN;
			  else // phase band
				Result(line-cutout.linelo,pix-cutout.pixlo) = LINE(0,pix-1+pixelsondisk);
			  }
			}
		  break;
		  }
	
		default:
		  PRINT_ERROR("readphase::not correct format on file.")
		  throw(unhandled_case_error);
		}
	
	  ifile.close();
	  return Result;
	  } // END readphase
	public matrix<complr4> readdata(window cutout)
	  {
	  TRACE_FUNCTION("productinfo::readdata (Batu 30-JUL-2007)")
	  // ______ Check if cutout < #lines/pixels ______
	  final int32 linesondisk = win.lines()/multilookL;
	  final int32 pixelsondisk = win.pixels()/multilookP;
	
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
		if (cutout.linelo < 1)
		  {
		  PRINT_ERROR("readdata: cutout window larger than what's on disk!")
		  throw(input_error);
		  }
		if (cutout.linehi>linesondisk)
		  {
		  PRINT_ERROR("readdata: cutout window larger than what's on disk!")
		  throw(input_error);
		  }
		if (cutout.pixlo < 1)
		  {
		  PRINT_ERROR("readdata: cutout window larger than what's on disk!")
		  throw(input_error);
		  }
		if (cutout.pixhi>pixelsondisk)
		  {
		  PRINT_ERROR("readdata: cutout window larger than what's on disk!")
		  throw(input_error);
		  }
	
	//  #endif
	  // ______ Open file ______
	//#ifdef __NO_IOS_BINARY__
	//  ifstream ifile(file, ios::in);
	//#else
	//  ifstream ifile(file, ios::in | ios::binary);
	//#endif
	  ifstream ifile;
	  openfstream(ifile,file);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(ifile,file,__FILE__,__LINE__);
	
	  matrix<complr4> Result = new matrix(cutout.lines(),cutout.pixels());
	
	  // ====== Actual read data ======
		  DEBUG.print("reading complex real4 from file (readdata)");
		  final window dummyoffset = new window(1,99999,1,99999); // only 1's used -> no offset
		  readfile(Result,file,linesondisk,cutout,dummyoffset);
	
	  ifile.close();
	  return Result;
	  } // END readdata
	public matrix<real4> readdatar4(window cutout)
	  {
	  TRACE_FUNCTION("productinfo::readdatar4 (MA & BO 19-NOV-2008)")
	  // --- Log debug info ---
	  DEBUG << "Reading file:     " << file;
	  DEBUG.print();
	  DEBUG << "Formatflag:       " << formatflag;
	  DEBUG.print();
	  DEBUG << "Currentwindow:    "; // appends to DEBUG and prints
	  win.disp();
													  // [MA] tricked since productinfo lacks currentwindow()
	  DEBUG << "Window from file: "; // appends to DEBUG and prints
	  cutout.disp();
	  //
	
	  // ______ Check if cutout < #lines/pixels ______
	  final int32 linesondisk = win.lines()/multilookL;
	  final int32 pixelsondisk = win.pixels()/multilookP;
	
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
		if (cutout.linelo < 1)
		  {
		  PRINT_ERROR("readdata: cutout window larger than what's on disk!")
		  throw(input_error);
		  }
		if (cutout.linehi>linesondisk)
		  {
		  PRINT_ERROR("readdata: cutout window larger than what's on disk!")
		  throw(input_error);
		  }
		if (cutout.pixlo < 1)
		  {
		  PRINT_ERROR("readdata: cutout window larger than what's on disk!")
		  throw(input_error);
		  }
		if (cutout.pixhi>pixelsondisk)
		  {
		  PRINT_ERROR("readdata: cutout window larger than what's on disk!")
		  throw(input_error);
		  }
	
	//  #endif
	  ifstream ifile;
	  openfstream(ifile,file);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(ifile,file,__FILE__,__LINE__);
	  ifile.close();
	
	  //matrix<real4> Result(cutout.lines(),cutout.pixels());
	
	  // ====== Actual read data ======
	  switch (formatflag)
		{
		case FORMATR4: // [MA] bigendian support
		  {
		  matrix<real4> Result = new matrix(cutout.lines(),cutout.pixels());
		  DEBUG.print("reading real4 from file (readdata)");
		  //const window dummyoffset(1,99999,1,99999);      // only 1's used -> no offset
		  //readfile(Result,file,linesondisk,cutout,dummyoffset);
		  readfile(Result,file,linesondisk,cutout,win);
		  return Result;
		  break;
		  }
	//   case FORMATR8:
	//     {
	//     matrix<real8> Result(cutout.lines(),cutout.pixels());
	//     DEBUG.print("reading real8 from file (readdata)");
	//     readfile(Result,file,linesondisk,cutout,win);
	//     return Result;
	//     break;
	//     }
	//   case FORMATCR4:
	//     {
	//     matrix<complr4> Result(cutout.lines(),cutout.pixels());
	//     DEBUG.print("reading complex real4 from file (readdata)");
	//     readfile(Result,file,linesondisk,cutout,win);
	//     return Result;
	//     break;
	//     }
		default:
		  PRINT_ERROR("readdata::not correct format on file.")
		  throw(file_error);
	   }
	  //return Result;
	  } // END readdatar4
}
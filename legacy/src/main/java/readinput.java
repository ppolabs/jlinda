public class GlobalMembersReadinput
{
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
	// * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/readinput.cc,v $
	// * $Revision: 3.36 $
	// * $Date: 2009/01/09 11:09:20 $
	// * $Author: TUDelft $
	// *
	// * implementation of readinput.
	// ***************************************************************





	// ______ displevel used in ioroutines.h, changed here ______
	public static String[] WARNS = new String[6]; // remember 6 last warnings in WARNS
	public static int32 beeplevel =1; // global variable for beeping.
									// 0:  nobeep; BEEP OFF
									// -1: beep on error exit ; BEEP ERROR
									// 1:  beep on error, warnings; BEEP WARNING
									// 2:  beep on error, warnings, progress;
									//     BEEP PROGRESS, BEEP [ON]
	public static int32 displevel =30000; // controls level of screen output

//***************************************************************
// *    checkgeneral                                              *
// *                                                              *
// * Checks general cards.                                        *
// * Repair process card                                          *
// *                                                              *
// *    Bert Kampes, 06-Sep-1999                                  *
// ***************************************************************
									// -100 only errors
									// 0:     warnings and errors
									// 10000: progress, warn and err
									// 20000: info, pro, warn and err
									// 30000: debug, info, pro, warn, err


	// ______ Checks input prototypes (see below in this file) ______
	public static void checkgeneral(RefObject<input_gen> generalinput, int16 onlyprocess)
	  {
	  TRACE_FUNCTION("checkgeneral (BK 06-Sep-1999)")
//C++ TO JAVA CONVERTER NOTE: 'register' variable declarations are not supported in Java:
//ORIGINAL LINE: register int32 i;
	  int32 i;
	  int32 cs_processcard = 0;

	// ______ Repair process control if ONLYPROCESS card was present ______
	  for (i =0; i<NUMPROCESSES; i++)
		cs_processcard += generalinput.argvalue.process[i];
	  if (onlyprocess != -1) // initialized to -1;
		{
		generalinput.argvalue.interactive=false;
		INFO.print("ONLYPROCESS card present, batch processing.");
		if (cs_processcard == 1)
		  WARNING.print("PROCESS card ignored due to presence ONLYPROCESS card.");
		if (cs_processcard > 1)
		  WARNING.print("PROCESS cards ignored due to presence ONLYPROCESS card.");
		for (i =0; i<NUMPROCESSES; i++) // do not process anything...
		  generalinput.argvalue.process[i]=0;
		generalinput.argvalue.process[onlyprocess]=1; // exept this one.
		}

	  else // check process card presence
		{
		if (cs_processcard == 0) // no cards
		  {
		  PRINT_ERROR("code 303: No (ONLY)PROCESS card present, exiting.")
		  throw(keyword_error);
		  }
		}

	  INFO.print("\n\t*** General input cards ***");
	  INFO << "MEMORY: \tAvailable to Doris [MB]: \t" << generalinput.argvalue.memory/1e6;
	  INFO.print();

	  INFO << "M_RESFILE: \tResultfile for master: \t\t" << generalinput.argvalue.m_resfile;
	  INFO.print();
	  INFO << "S_RESFILE: \tResultfile for slave: \t\t" << generalinput.argvalue.s_resfile;
	  INFO.print();
	  INFO << "I_RESFILE: \tResultfile for products: \t" << generalinput.argvalue.i_resfile;
	  INFO.print();
	  INFO << "LOGFILE: \tOut file for logging: \t\t" << generalinput.argvalue.logfile;
	  INFO.print();
	  INFO << "ORB_INTERP: \tmethod selector value:\t\t" << generalinput.argvalue.orb_interp;
	  INFO.print();
	  INFO << "DUMPBASELINE: evaluation grid for baseline: \t" << generalinput.argvalue.dumpbaselineL << " lines x " << generalinput.argvalue.dumpbaselineP << " pixels: ";
	  INFO.print();
	  INFO << "HEIGHT: \taverage terrain height:\t\t" << generalinput.argvalue.terrain_height;
	  INFO.print();
	  INFO << "TIEPOINT: \tlat/lon/hei: " << generalinput.argvalue.tiepoint.x << " " << generalinput.argvalue.tiepoint.y << " " << generalinput.argvalue.tiepoint.z;
	  INFO.print();

	  if (!strcmp(generalinput.argvalue.m_resfile,generalinput.argvalue.s_resfile))
		{
		PRINT_ERROR("same name master and slave resultfile not allowed.")
		throw(keyword_error);
		}
	  if (!strcmp(generalinput.argvalue.m_resfile,generalinput.argvalue.i_resfile))
		{
		PRINT_ERROR("same name master and interferogram resultfile not allowed.");
		throw(keyword_error);
		}
	  if (!strcmp(generalinput.argvalue.m_resfile,generalinput.argvalue.logfile))
		{
		PRINT_ERROR("same name master resultfile and logfile not allowed.");
		throw(keyword_error);
		}
	  if (!strcmp(generalinput.argvalue.i_resfile,generalinput.argvalue.logfile))
		{
		PRINT_ERROR("same name interferogram resultfile and logfile not allowed.");
		throw(keyword_error);
		}
	  } // END checkgeneral

//***************************************************************
// *    checkreadfiles                                            *
// *                                                              *
// * Checks cards for step readfiles.                             *
// *                                                              *
// *    Bert Kampes, 06-Sep-1999                                  *
// ***************************************************************
	public static void checkreadfiles(input_readfiles readfilesinput, int16 id)
	  {
	  TRACE_FUNCTION("checkreadfiles (BK 06-Sep-1999)")
	  switch (id)
		{
		case MASTERID:
		  INFO.print("\n\t*** Input for step M_READFILES (master) ***");
		  INFO << "M_IN_METHOD: \tmethod selected for master: \t\t" << readfilesinput.sensor_id;
		  INFO.print();
		  INFO << "M_IN_VOL:    \tVolumefile of master: \t\t" << readfilesinput.volfile;
		  INFO.print();
		  INFO << "M_IN_LEA:    \tLeaderfile of master: \t\t" << readfilesinput.leaderfile;
		  INFO.print();
		  INFO << "M_IN_NULL:   \tNullfile of master:  \t\t" << readfilesinput.nullfile;
		  INFO.print();
		  INFO << "M_IN_DAT:    \tDatfile of master:    \t\t" << readfilesinput.datfile;
		  INFO.print();
		  if (readfilesinput.sensor_id == SLC_ERS)
			{
			if (!specified(readfilesinput.volfile))
			  {
			  PRINT_ERROR("M_IN_VOL not defined");
			  throw(keyword_error);
			  }
			if (!specified(readfilesinput.leaderfile))
			  {
			  PRINT_ERROR("M_IN_LEA not defined");
			  throw(keyword_error);
			  }
			if (!specified(readfilesinput.nullfile))
			  WARNING.print("M_IN_NULL not defined");
			}
		  // ___ use datfile for asar input file ___
		  if (!specified(readfilesinput.datfile))
			{
			PRINT_ERROR("M_IN_DAT not defined");
			throw(keyword_error);
			}
		  if (readfilesinput.sensor_id == SLC_TSX) // [MA] TSX
			{
			if (!specified(readfilesinput.leaderfile)) // if changed, update also processor for leaderfile
			  {
			  PRINT_ERROR("M_IN_LEA not defined");
			  throw(keyword_error);
			  }
			}
		  break;

		case SLAVEID:
		  INFO.print("\n\t*** Input for step S_READFILES (slave) ***");
		  INFO << "S_IN_METHOD: \tmethod selected for slave: \t\t" << readfilesinput.sensor_id;
		  INFO.print();
		  INFO << "S_IN_VOL:    \tvolumefile of slave:  \t\t" << readfilesinput.volfile;
		  INFO.print();
		  INFO << "S_IN_LEA:    \tleaderfile of slave:  \t\t" << readfilesinput.leaderfile;
		  INFO.print();
		  INFO << "S_IN_NULL:   \tnullfile of slave:   \t\t" << readfilesinput.nullfile;
		  INFO.print();
		  INFO << "S_IN_DAT:    \tdatfile of slave:     \t\t" << readfilesinput.datfile;
		  INFO.print();
		  if (readfilesinput.sensor_id == SLC_ERS)
			{
			if (!specified(readfilesinput.volfile))
			  {
			  PRINT_ERROR("S_IN_VOL not defined");
			  throw(keyword_error);
			  }
			if (!specified(readfilesinput.leaderfile))
			  {
			  PRINT_ERROR("S_IN_LEA not defined");
			  throw(keyword_error);
			  }
			if (!specified(readfilesinput.nullfile))
			  WARNING.print("S_IN_NULL not defined");
			}
		  // ___ use datfile for asar input file ___
		  if (!specified(readfilesinput.datfile))
			{
			PRINT_ERROR("S_IN_DAT not defined");
			throw(keyword_error);
			}
		  if (readfilesinput.sensor_id == SLC_TSX) // [MA] TSX
			{
			if (!specified(readfilesinput.leaderfile)) // if changed, update also processor.cc for the [leader]file
			  {
			  PRINT_ERROR("M_IN_LEA not defined");
			  throw(keyword_error);
			  }
			}
		  break;

		default:
		  PRINT_ERROR("panic: impossible");
		  throw(keyword_error);
		} // switch
	  if (readfilesinput.sensor_id == SLC_ERS)
		{
		if (!strcmp(readfilesinput.volfile,readfilesinput.leaderfile))
		  {
		  PRINT_ERROR("same file name volume and leader file not allowed.");
		  throw(keyword_error);
		  }
		if (!strcmp(readfilesinput.volfile,readfilesinput.nullfile))
		  {
		  PRINT_ERROR("same file name volume and null file not allowed.");
		  throw(keyword_error);
		  }
		if (!strcmp(readfilesinput.volfile,readfilesinput.datfile))
		  {
		  PRINT_ERROR("same file name volume and data file not allowed.");
		  throw(keyword_error);
		  }
		if (!strcmp(readfilesinput.leaderfile,readfilesinput.nullfile))
		  {
		  PRINT_ERROR("same file name leader and null file not allowed.");
		  throw(keyword_error);
		  }
		if (!strcmp(readfilesinput.leaderfile,readfilesinput.datfile))
		  {
		  PRINT_ERROR("same file name leader and data file not allowed.");
		  throw(keyword_error);
		  }
		if (!strcmp(readfilesinput.nullfile,readfilesinput.datfile))
		  {
		  PRINT_ERROR("same file name null and data file not allowed.");
		  throw(keyword_error);
		  }
		}
	  } // END checkreadfiles

//***************************************************************
// *    checkcrop                                                 *
// *                                                              *
// * Checks cards for step crop.                          *
// *                                                              *
// *    Bert Kampes, 06-Sep-1999                                  *
// ***************************************************************
	public static void checkcrop(input_crop cropinput, int16 id)
	  {
	  TRACE_FUNCTION("checkcrop (BK 06-Sep-1999)")
	  switch (id)
		{
		case MASTERID:
		  INFO.print("\n\t*** Input for step M_CROP (master) ***");
		  INFO << "M_IDCROP: \tidentifier master write slc data to raster format: \t" << cropinput.idcrop;
		  INFO.print();
		  INFO << "M_CROP_IN: \tslc data inputfile for master:   \t" << cropinput.filein1;
		  INFO.print();
		  INFO << "M_CROP_OUT: \traw data outputfile for master: \t" << cropinput.fileout1;
		  INFO.print();
		  INFO << "M_DBOW: \tProcessing master line " << cropinput.dbow.linelo << " to " << cropinput.dbow.linehi << ". pixel " << cropinput.dbow.pixlo << " to " << cropinput.dbow.pixhi << ".";
		  INFO.print();
		  if (cropinput.dbow_geo.pixhi != 0)
			{
			INFO.print("M_DBOW_GEO: overrides M_DBOW! processing: ");
			INFO << "center latitude " << cropinput.dbow_geo.linelo/1e6-360.0 << "; center longitude " << cropinput.dbow_geo.linehi/1e6-360.0 << "; height, width: " << cropinput.dbow_geo.pixlo << ", " << cropinput.dbow_geo.pixhi;
			INFO.print();
			}
		  break;

		case SLAVEID:
		  INFO.print("\n\t*** Input for step S_CROP (slave) ***");
		  INFO << "S_IDCROP: \tidentifier slave write slc data to raster format: \t" << cropinput.idcrop;
		  INFO.print();
		  INFO << "S_CROP_IN: \tslc data inputfile for slave:   \t" << cropinput.filein1;
		  INFO.print();
		  INFO << "S_CROP_OUT: \traw data outputfile for slave: \t" << cropinput.fileout1;
		  INFO.print();
		  INFO << "S_DBOW: \tProcessing slave line " << cropinput.dbow.linelo << " to " << cropinput.dbow.linehi << ". pixel " << cropinput.dbow.pixlo << " to " << cropinput.dbow.pixhi << ".";
		  INFO.print();
		  if (cropinput.dbow_geo.pixhi != 0)
			{
			INFO.print("S_DBOW_GEO: overrides S_DBOW! processing: ");
			INFO << "center latitude " << cropinput.dbow_geo.linelo/1e6-360.0 << "; center longitude " << cropinput.dbow_geo.linehi/1e6-360.0 << "; height, width: " << cropinput.dbow_geo.pixlo << ", " << cropinput.dbow_geo.pixhi;
			INFO.print();
			}
		  break;

		default:
		  PRINT_ERROR("panic: impossible");
		  throw(keyword_error);
		}
	  } // END checkcrop

//____RaffaeleNutricato START MODIFICATION SECTION 11
//***************************************************************
// *    checkoversample                                           *
// *                                                              *
// * Checks cards for step oversample.                            *
// *                                                              *
// *    Raffaele Nutricato, 12-Jan-2004                           *
// ***************************************************************
	//____RaffaeleNutricato START MODIFICATION SECTION 1
	public static void checkoversample(input_oversample oversampleinput, int16 id)
	  {
	  TRACE_FUNCTION("checkoversample (Raffaele Nutricato 12-Jan-2004)")
	  switch (id)
		{
		case MASTERID:
		  INFO.print("\n\t*** Input for step M_OVS (master) ***");
		  INFO << "M_OVS_OUT: \t\tData output file for ovs master: " << oversampleinput.fileoutovs; // Output file for the oversampled master
		  INFO.print();
		  INFO << "M_OVS_FACT_RNG: \tOversampling ratio in the master range direction: " << oversampleinput.OsrRange; // Oversampling ratio in the range direction.
		  INFO.print();
		  INFO << "M_OVS_FACT_AZI: \tOversampling ratio in the master azimuth direction: " << oversampleinput.OsrAzimuth; // Oversampling ratio in the azimuth direction.
		  INFO.print();
		  INFO << "M_OVS_KERNELSIZE: \tKernel length for the master oversampling: " << oversampleinput.FilterSize; // Length of the interpolation kernel in range.
		  INFO.print();
		  INFO << "M_OVS_OUT_FORMAT: \tOutput data format for the oversampled master: ";
		  if (oversampleinput.oformatflag==FORMATCR4)
			INFO << "complex_real4";
		  if (oversampleinput.oformatflag==FORMATCI2)
			INFO << "complex_short";
		  INFO.print();
		  break;

		case SLAVEID:
		  INFO.print("\n\t*** Input for step S_OVS (slave) ***");
		  INFO << "S_OVS_OUT: \t\tData output file for ovs slave: " << oversampleinput.fileoutovs; // Output file for the oversampled slave
		  INFO.print();
		  INFO << "S_OVS_FACT_RNG: \tOversampling ratio in the slave range direction: " << oversampleinput.OsrRange; // Oversampling ratio in the range direction.
		  INFO.print();
		  INFO << "S_OVS_FACT_AZI: \tOversampling ratio in the slave azimuth direction: " << oversampleinput.OsrAzimuth; // Oversampling ratio in the azimuth direction.
		  INFO.print();
		  INFO << "S_OVS_KERNELSIZE: \tKernel length for the slave oversampling: " << oversampleinput.FilterSize; // Length of the interpolation kernel in range.
		  INFO.print();
		  INFO << "S_OVS_OUT_FORMAT: \tOutput data format for the oversampled slave: ";
		  if (oversampleinput.oformatflag==FORMATCR4)
			INFO << "complex_real4";
		  if (oversampleinput.oformatflag==FORMATCI2)
			INFO << "complex_short";
		  INFO.print();
		  break;

		default:
		  PRINT_ERROR("panic: impossible");
		  throw(keyword_error);
		}
	  } // END checkoversample

//____RaffaeleNutricato END MODIFICATION SECTION 11


//***************************************************************
// *    checkporbits                                              *
// *                                                              *
// * Checks cards for step porbits.                               *
// *                                                              *
// *    Bert Kampes, 06-Sep-1999                                  *
// ***************************************************************
	//____RaffaeleNutricato END MODIFICATION SECTION 1
	public static void checkporbits(input_pr_orbits porbitsinput, int16 id)
	  {
	  TRACE_FUNCTION("checkporbits (BK 06-Sep-1999)")
	  switch (id)
		{
		case MASTERID:
		  INFO.print("\n\t*** Input for step M_PORBITS (master) ***");
		  INFO << "M_ORBDIR: \tPrecise orbits master in: " << porbitsinput.m_orbdir;
		  INFO.print();
		  if (!specified(porbitsinput.m_orbdir))
			{
			PRINT_ERROR("M_ORBDIR: no directory specified.");
			throw(keyword_error);
			}
		  break;

		case SLAVEID:
		  INFO.print("\n\t*** Input for step S_PORBITS (slave) ***");
		  INFO << "S_ORBDIR: \tPrecise orbits slave in: " << porbitsinput.s_orbdir;
		  INFO.print();
		  if (!specified(porbitsinput.s_orbdir))
			{
			PRINT_ERROR("S_ORBDIR: no directory specified.");
			throw(keyword_error);
			}
		  break;

		default:
		  PRINT_ERROR("panic: impossible");
		  throw(keyword_error);
		}
	  INFO << "ORB_INTERVAL: \ttime between ephemerides: \t" << porbitsinput.timeinterval;
	  INFO.print();
	  INFO << "ORB_EXTRATIME: \ttime before first, after last line: \t" << porbitsinput.timebefore;
	  INFO.print();
	  if (!porbitsinput.dumpmasterorbit<0)
		INFO.print("dumping masterorbit to file.");
	  if (!porbitsinput.dumpslaveorbit<0)
		INFO.print("dumping slaveorbit to file.");
	  } // END checkporbits

//***************************************************************
// *    checksimamp                                               *
// *                                                              *
// * Checks cards for step simamp.                                *
// *                                                              *
// *    Mahmut Arikan, 30-Oct-2008                                *
// ***************************************************************
	public static void checksimamp(input_simamp simampinput)
	  {
	  TRACE_FUNCTION("checksimamp (MA 30-oct-2008)")
	  INFO.print("\n\t*** Input for step SIMAMP ***");
	  INFO << "SAM_IN_DEM:    \t" << simampinput.firefdem;
	  INFO.print();
	  if (specified(simampinput.fodem))
		{
		INFO << "SAM_OUT_DEM:   \t" << simampinput.fodem << "; output requested of input DEM.";
		INFO.print();
		if (!strcmp(simampinput.fosimamp,simampinput.fodem))
		  {
		  PRINT_ERROR("OUT_DEM, OUT_FILE: Same output file name.");
		  throw(keyword_error);
		  }
		}
	  if (specified(simampinput.fosimamp))
		{
		INFO << "SAM_OUT_FILE:   \t" << simampinput.fosimamp << "; output requested of simulated amplitude.";
		INFO.print();
		if (!strcmp(simampinput.fosimamp,simampinput.fodem))
		  {
		  PRINT_ERROR("OUT_DEM, OUT_FILE: Same output file name.");
		  throw(keyword_error);
		  }
		}
	  INFO << "SAM_IN_SIZE:   \t" << simampinput.demrows << " " << simampinput.demcols << "; number of rows (latitude), columns (lon) in DEM.";
	  INFO.print();
	  INFO << "SAM_IN_UL:     \t" << rad2deg(simampinput.demlatleftupper) << " " << rad2deg(simampinput.demlonleftupper) << "; coordinates of upper left corner (first row/col).";
	  INFO.print();
	  INFO << "SAM_IN_DELTA:  \t" << rad2deg(simampinput.demdeltalat) << " " << rad2deg(simampinput.demdeltalon);
	  INFO.print();
	  INFO << "SAM_IN_NODATA:  \t" << simampinput.demnodata << "; this number in DEM will be set to 0 reference phase.";
	  INFO.print();

	  INFO << "SAM_IN_FORMAT: \tinput format DEM: \t";
	  switch (simampinput.iformatflag)
		{
		case FORMATR4:
		  INFO << "real4.";
		  break;
		case FORMATR8:
		  INFO << "real8.";
		  break;
		case FORMATI2:
		  INFO << "signed short.";
		  break;
		case FORMATI2_BIGENDIAN:
		  INFO << "signed short big endian.";
		  break;
		default:
		  PRINT_ERROR("totally impossible, checked input.");
		  throw(keyword_error);
		}
	  INFO.print();



	// ______ Check some things ______
	  if (!existed(simampinput.firefdem))
		{
		ERROR << "SAM_IN_DEM:   \t" << simampinput.firefdem << " can not be opened.";
		PRINT_ERROR(ERROR.get_str())
		throw(keyword_error);
		}
	  if (rad2deg(simampinput.demdeltalat)<.0002) // [MA] 1 arc secs = 1/3600=0.00027 : allow finer resolutions
		{
		WARNING << "SAM_IN_DELTA: \t" << rad2deg(simampinput.demdeltalat) << " [deg] seems very small (input in decimal degrees).";
		WARNING.print();
		}
	  if (simampinput.demrows<1 || simampinput.demrows>100000)
		{
		WARNING << "SAM_DEM_SIZE: numrows: \t" << simampinput.demrows << " seems to be wrong.";
		WARNING.print();
		}
	  if (simampinput.demcols<1 || simampinput.demcols>100000)
		{
		WARNING << "SAM_DEM_SIZE: numcols: \t" << simampinput.demcols << " seems to be wrong.";
		WARNING.print();
		}
	  if (rad2deg(simampinput.demdeltalon)<.0002)
		{
		WARNING << "SAM_IN_DELTA: \t" << simampinput.demdeltalon << " seems very small (it should be in decimal degrees).";
		WARNING.print();
		}
	  if (rad2deg(simampinput.demlatleftupper) < -90. || rad2deg(simampinput.demlatleftupper) > 90.)
		{
		ERROR << "SAM_IN_LU:    \t" << rad2deg(simampinput.demlatleftupper) << " out of range (-90:90).";
		PRINT_ERROR(ERROR.get_str())
		throw(keyword_error);
		}
	  if (rad2deg(simampinput.demlonleftupper) < -180. || rad2deg(simampinput.demlonleftupper) > 180.)
		{
		WARNING << "SAM_IN_LU:    \t" << rad2deg(simampinput.demlonleftupper) << " out of range (-180:180).";
		WARNING.print();
		}
	  } // END simamp [MA]

//***************************************************************
// *    checkmtiming                                              *
// *                                                              *
// * Checks cards for step simamp corr.                           *
// *                                                              *
// *    Mahmut Arikan, 11-Nov-2008                                *
// ***************************************************************
	public static void checkmtiming(input_mtiming mtiminginput)
	  {
	  TRACE_FUNCTION("checkmtiming (MA, BO 11-Nov-2008)")
	  INFO.print("\n\t*** Input for step MASTER TIMING ***");
	  switch (mtiminginput.method)
		{
		case cc_magfft:
		  INFO.print("MTE_METHOD: \tMAGFFT is used for MASTER TIMING ERROR estimation.");
		  break;
		case cc_magspace:
		  INFO.print("MTE_METHOD: \tMAGSPACE is used for MASTER TIMING ERROR estimation.");
		  break;
		default:
		  PRINT_ERROR("panic.");
		  throw(keyword_error);
		}
	  if (specified(mtiminginput.ifpositions))
		{
		INFO << "MTE_IN_POS: \tFile: " << mtiminginput.ifpositions << " containing " << mtiminginput.Nwin << " positions is used.";
		INFO.print();
		}
	  else // no input file
		{
		INFO << "MTE_NWIN: \tDistributing " << mtiminginput.Nwin << " windows for master timing error estimation (correlation).";
		INFO.print();
		}
	  } // END checkmtiming

//***************************************************************
// *    checkslant2h                                              *
// *                                                              *
// * Checks cards for step slant2h.                               *
// *                                                              *
// *    Bert Kampes, 29-Sep-1999                                  *
// ***************************************************************
	public static void checkslant2h(input_slant2h slant2hinput)
	  {
	  TRACE_FUNCTION("checkslant2h (BK 29-Sep-1999)")
	  INFO.print("\n\t*** Input for step SLANT2H ***");
	  switch (slant2hinput.method)
		{
		case s2h_schwabisch:
		  INFO.print("Method schwabisch is used for slant2height conversion.");

		  INFO << "S2H_NPOINTS: \tNumber of points used in computation:   " << slant2hinput.Npoints;
		  INFO.print();
		  INFO << "S2H_DEGREE1D: \tDegree of 1d polynomial at each point: " << slant2hinput.degree1d;
		  INFO.print();
		  INFO << "S2H_DEGREE2D: \tDegree of 2d polynomial to compute     " << "coefficients of 1D polynomial: " << slant2hinput.degree2d;
		  INFO.print();
		  INFO << "S2H_NHEIGHTS: \t#heights evaluation ref. phase for 1d polynomial: " << slant2hinput.Nheights;
		  INFO.print();
		  break;

		case s2h_ambiguity:
		  INFO.print("Method exact is used for slant2height conversion.");
		  INFO << "S2H_OUT_LAM: \tData output file for lambda: " << slant2hinput.folam;
		  INFO.print();
		  INFO << "S2H_OUT_PHI: \tData output file for phi:    " << slant2hinput.fophi;
		  INFO.print();
		  if (!strcmp(slant2hinput.folam,slant2hinput.fophi))
			{
			PRINT_ERROR("Same filename S2H_OUT_LAM and S2H_OUT_PHI not allowed.");
			throw(keyword_error);
			}
		  if (!strcmp(slant2hinput.fohei,slant2hinput.fophi))
			{
			PRINT_ERROR("Same filename S2H_OUT_HEI and S2H_OUT_PHI not allowed.");
			throw(keyword_error);
			}
		  if (!strcmp(slant2hinput.folam,slant2hinput.fohei))
			{
			PRINT_ERROR("Same filename S2H_OUT_LAM and S2H_OUT_HEI not allowed.");
			throw(keyword_error);
			}
		  break;

		case s2h_rodriguez:
		  DEBUG.print("do some checks here as well.");
		  break;

		default:
		  PRINT_ERROR("impossible, unknown method s2h");
		  throw(keyword_error);
		}

		INFO << "S2H_OUT_HEI: Data output file for height:   " << slant2hinput.fohei;
		INFO.print();
	  } // END checkslant2h

//***************************************************************
// *    checkunwrap                                               *
// *                                                              *
// * Checks cards for step unwrap.                                *
// *                                                              *
// *    Bert Kampes, 29-Sep-1999                                  *
// ***************************************************************
	public static void checkunwrap(input_unwrap unwrapinput)
	  {
	  TRACE_FUNCTION("checkunwrap (BK 29-Sep-1999)")
	  INFO.print("\n\t*** Input for step UNWRAP ***");
	  INFO << "UW_OUT_FILE: \tOutput data file for unwrapped interferogram: " << unwrapinput.fouint;
	  INFO.print();
	  INFO << "UW_OUT_FORMAT: \tOutput data format for unwrapped interferogram: " << unwrapinput.oformatflag;
	  INFO.print();
	  switch (unwrapinput.method)
		{
		case uw_method1:
		  INFO.print("Method 1 is used for unwrapping (treef_ramon unix calls).");
		  INFO << "UW_SEEDS: ";
		  if (Character.isWhitespace(unwrapinput.seedfile[0])) // no file specified
			INFO << "Delta seed in line/pixel direction: " << unwrapinput.deltaLseed << " " << unwrapinput.deltaPseed;
		  else
			INFO << " Input file with seeds for unwrapping: " << unwrapinput.seedfile;
		  INFO.print();
		  INFO << "UW_OUT_REGIONS: \tOutput data file with region numbers: " << unwrapinput.foregions;
		  INFO.print();
		  break;

		case uw_method2:
		  INFO.print("Method 2: SNAPHU is used for unwrapping.");
		  INFO.print("Please make sure snaphu is installed.  check results carefuly.");
		  INFO << "UW_SNAPHU_LOG: \tOutput log file of snaphu: " << unwrapinput.snaphu_log;
		  INFO.print();
		  INFO << "UW_SNAPHU_INIT: \tinitialization method of snaphu: " << unwrapinput.snaphu_init;
		  INFO.print();
		  INFO << "UW_SNAPHU_COH: \tinput coherence file for snaphu: " << unwrapinput.snaphu_coh;
		  INFO.print();
		  INFO << "UW_SNAPHU_MODE: \tsnaphu mode: " << unwrapinput.snaphu_mode;
		  INFO.print();
		  break;

		case uw_method3:
		  INFO.print("Method 3 is used for unwrapping (?).");
		  PRINT_ERROR("only for delft, standalone application.");
		  throw(keyword_error);
		  break;

		default:
		  PRINT_ERROR("probably forgot to update this, should not be possible.");
		  throw(keyword_error);
		}
	  } // END checkunwrap

//***************************************************************
// *    checkgeocode                                              *
// *                                                              *
// * Checks cards for step geocode.                               *
// *                                                              *
// *    Bert Kampes, 29-Sep-1999                                  *
// ***************************************************************
	public static void checkgeocode(input_geocode geocodeinput)
	  {
	  TRACE_FUNCTION("checkgeocode (BK 29-Sep-1999)")
	  INFO.print("\n\t*** Input for step GEOCODE ***");
	  INFO << "GEO_OUT_PHI: \tOutput file for latitude: \t" << geocodeinput.fophi;
	  INFO.print();
	  INFO << "GEO_OUT_PHI: \tOutput file for longitude: \t" << geocodeinput.folam;
	  INFO.print();
	  if (!strcmp(geocodeinput.fophi,geocodeinput.folam))
		{
		ERROR << "Same file name GEO_OUT_PHI and GEO_OUT_LAM not allowed.";
		PRINT_ERROR(ERROR.get_str())
		throw(keyword_error);
		}
	  } // END checkgeocode

//***************************************************************
// *    checkcoarsecorr                                           *
// *                                                              *
// * Checks cards for step coarsecorr.                            *
// *                                                              *
// *    Bert Kampes, 29-Sep-1999                                  *
// ***************************************************************
	public static void checkcoarsecorr(input_coarsecorr coarsecorrinput)
	  {
	  TRACE_FUNCTION("checkcoarsecorr (BK 29-Sep-1999)")
	  INFO.print("\n\t*** Input for step COARSE_CORR ***");
	  switch (coarsecorrinput.method)
		{
		case cc_magfft:
		  INFO.print("CC_METHOD: \tMAGFFT is used for coarse correlation.");
		  break;
		case cc_magspace:
		  INFO.print("CC_METHOD: \tMAGSPACE is used for coarse correlation.");
		  break;
		default:
		  PRINT_ERROR("panic.");
		  throw(keyword_error);
		}
	  if (specified(coarsecorrinput.ifpositions))
		{
		INFO << "CC_IN_POS: \tFile: " << coarsecorrinput.ifpositions << " containing " << coarsecorrinput.Nwin << " positions is used.";
		INFO.print();
		}
	  else // no input file
		{
		INFO << "CC_NWIN: \tDistributing " << coarsecorrinput.Nwin << " windows for coarse coregistration based on correlation.";
		INFO.print();
		}
	  } // END checkcoarsecorr

//***************************************************************
// *    checkfine                                                 *
// *                                                              *
// * Checks cards for step fine.                                  *
// *                                                              *
// *    Bert Kampes, 29-Sep-1999                                  *
// ***************************************************************
	public static void checkfine(input_fine fineinput)
	  {
	  TRACE_FUNCTION("checkfine (BK 29-Sep-1999)")
	  INFO.print("\n\t*** Input for step FINE ***");
	  switch (fineinput.method)
		{
	//    case fc_cmplxfft:
	//      INFO.print("FC_METHOD: \tCMPLXFFT is used for fine correlation.");
	//      break;
	//    case fc_cmplxspace:
	//      INFO.print("FC_METHOD: \tCMPLXSPACE is used for fine correlation.");
	//      break;
		case fc_magfft:
		  INFO.print("FC_METHOD: \tMAGFFT is used for fine correlation.");
		  break;
		case fc_magspace:
		  INFO.print("FC_METHOD: \tMAGSPACE is used for fine correlation.");
		  break;
		case fc_oversample:
		  INFO.print("FC_METHOD: \tOVERSAMPLE is used for fine correlation.");
		  break;
		default:
		  PRINT_ERROR("panic.");
		  throw(keyword_error);
		}

	  if (specified(fineinput.ifpositions))
		{
		INFO << "FC_IN_POS: File: " << fineinput.ifpositions << " containing " << fineinput.Nwin << " positions is used.";
		INFO.print();
		}
	  else // no input file
		{
		INFO << "FC_NWIN: \tDistributing " << fineinput.Nwin << " windows for fine coregistration.";
		INFO.print();
		}
	  if (fineinput.plotoffsets)
		{
		INFO << "FC_PLOT " << fineinput.plotthreshold << ": Plotting estimated offsets with correlation > " << fineinput.plotthreshold;
		INFO.print();
		}
	  else
		{
		WARNING.print("It is highly recommended to use FC_PLOT to plot the computed FINE offset vectors.");
		}
	  if (fineinput.plotmagbg)
		INFO.print("FC_PLOT: magnitude will be in background of plot.");

	  INFO << "FC_WINSIZE for correlation: (l,p) = (" << fineinput.MasksizeL << ", " << fineinput.MasksizeP << ").";
	  INFO.print();
	  INFO.print("Max. estimable offset is half the window size");
	  INFO.print("If the offset varies a lot over the image, use larger winsize");
	  INFO << "FC_WINSIZE for oversampling (2*Acc): (l,p) = (" << 2 *fineinput.AccL << ", " << 2 *fineinput.AccP << ").";
	  INFO.print();
	  INFO << "FC_OSFACTOR: " << fineinput.osfactor << ". Oversampling factor for fine coregistration.";
	  INFO.print();
	  if (!ispower2(fineinput.osfactor))
		{
		ERROR << " no power of 2.";
		PRINT_ERROR(ERROR.get_str())
		throw(keyword_error);
		}
	  } // END checkfine

// ____end added by MA ____


// ____start added by FvL ____

//***************************************************************
// *    checkreltiming                                            *
// *                                                              *
// * Checks cards for step timing                                 *
// *                                                              *
// *    Freek van Leijen, 22-Aug-2007                             *
// ***************************************************************
	public static void checkreltiming(input_reltiming timinginput)
	  {
	  TRACE_FUNCTION("checkreltiming (FvL 22-aug-2007)")
	  INFO.print("\n\t*** Input for step REL_TIMING ***");
	  INFO << "RTE_THRESHOLD: \tThreshold correlation for model: \t" << timinginput.threshold;
	  INFO.print();
	  INFO << "RTE_MAXITER: \tNumber of points to remove (max): \t" << timinginput.maxiter;
	  INFO.print();
	  INFO << "RTE_K_ALPHA: \tCritical value for outlier test: \t" << timinginput.k_alpha;
	  INFO.print();
	  } // END checkreltiming

//***************************************************************
// *    checkdemassist                                            *
// *                                                              *
// * Checks cards for step demassist.                             *
// *                                                              *
// *    Freek van Leijen, 22-Aug-2007                             *
// ***************************************************************
	public static void checkdemassist(input_demassist demassistinput)
	  {
	  TRACE_FUNCTION("checkdemassist (FvL 22-aug-2007)")
	  INFO.print("\n\t*** Input for step DEMASSIST ***");
	  INFO << "DAC_IN_DEM:    \t" << demassistinput.firefdem;
	  INFO.print();
	  if (specified(demassistinput.forefdemhei))
		{
		INFO << "DAC_OUT_DEM_LP: \t" << demassistinput.forefdemhei << "; output requested of DEM [m] in radar coordinates.";
		INFO.print();
		}
	  if (specified(demassistinput.fodem))
		{
		INFO << "DAC_OUT_DEM:   \t" << demassistinput.fodem << "; output requested of input DEM.";
		INFO.print();
		if (!strcmp(demassistinput.fodemi,demassistinput.fodem))
		  {
		  PRINT_ERROR("OUT_DEM, OUT_DEMI: Same output file name.");
		  throw(keyword_error);
		  }
		}
	  if (specified(demassistinput.fodemi))
		{
		INFO << "DAC_OUT_DEMI:   \t" << demassistinput.fodemi << "; output requested of interpolated DEM.";
		INFO.print();
		if (!strcmp(demassistinput.fodemi,demassistinput.fodem))
		  {
		  PRINT_ERROR("OUT_DEM, OUT_DEMI: Same output file name.");
		  throw(keyword_error);
		  }
		}
	  INFO << "DAC_IN_SIZE:   \t" << demassistinput.demrows << " " << demassistinput.demcols << "; number of rows (latitude), columns (lon) in DEM.";
	  INFO.print();
	  INFO << "DAC_IN_UL:     \t" << rad2deg(demassistinput.demlatleftupper) << " " << rad2deg(demassistinput.demlonleftupper) << "; coordinates of upper left corner (first row/col).";
	  INFO.print();
	  INFO << "DAC_IN_DELTA:  \t" << rad2deg(demassistinput.demdeltalat) << " " << rad2deg(demassistinput.demdeltalon);
	  INFO.print();
	  INFO << "DAC_IN_NODATA:  \t" << demassistinput.demnodata << "; this number in DEM will be set to 0 reference phase.";
	  INFO.print();

	  INFO << "DAC_IN_FORMAT: \tinput format DEM: \t";
	  switch (demassistinput.iformatflag)
		{
		case FORMATR4:
		  INFO << "real4.";
		  break;
		case FORMATR8:
		  INFO << "real8.";
		  break;
		case FORMATI2:
		  INFO << "signed short.";
		  break;
		case FORMATI2_BIGENDIAN:
		  INFO << "signed short big endian.";
		  break;
		default:
		  PRINT_ERROR("totally impossible, checked input.");
		  throw(keyword_error);
		}
	  INFO.print();



	// ______ Check some things ______
	  if (!existed(demassistinput.firefdem))
		{
		ERROR << "DAC_IN_DEM:   \t" << demassistinput.firefdem << " can not be opened.";
		PRINT_ERROR(ERROR.get_str())
		throw(keyword_error);
		}
	  if (rad2deg(demassistinput.demdeltalat)<.0002) //[MA] 1 arc secs = 1/3600=0.00027 : allow finer resolutions
		{
		WARNING << "DAC_IN_DELTA: \t" << rad2deg(demassistinput.demdeltalat) << " [deg] seems very small (input in decimal degrees).";
		WARNING.print();
		}
	  if (demassistinput.demrows<1 || demassistinput.demrows>100000)
		{
		WARNING << "DAC_DEM_SIZE: numrows: \t" << demassistinput.demrows << " seems to be wrong.";
		WARNING.print();
		}
	  if (demassistinput.demcols<1 || demassistinput.demcols>100000)
		{
		WARNING << "DAC_DEM_SIZE: numcols: \t" << demassistinput.demcols << " seems to be wrong.";
		WARNING.print();
		}
	  if (rad2deg(demassistinput.demdeltalon)<.0002)
		{
		WARNING << "DAC_IN_DELTA: \t" << demassistinput.demdeltalon << " seems very small (it should be in decimal degrees).";
		WARNING.print();
		}
	  if (rad2deg(demassistinput.demlatleftupper) < -90. || rad2deg(demassistinput.demlatleftupper) > 90.)
		{
		ERROR << "DAC_IN_LU:    \t" << rad2deg(demassistinput.demlatleftupper) << " out of range (-90:90).";
		PRINT_ERROR(ERROR.get_str())
		throw(keyword_error);
		}
	  if (rad2deg(demassistinput.demlonleftupper) < -180. || rad2deg(demassistinput.demlonleftupper) > 180.)
		{
		WARNING << "DAC_IN_LU:    \t" << rad2deg(demassistinput.demlonleftupper) << " out of range (-180:180).";
		WARNING.print();
		}
	  } // END demassist

// ____end added by FvL ____


//***************************************************************
// *    checkcoregpm                                              *
// *                                                              *
// * Checks cards for step coregpm.                               *
// *                                                              *
// *    Bert Kampes, 29-Sep-1999                                  *
// ***************************************************************
	public static void checkcoregpm(input_coregpm coregpminput)
	  {
	  TRACE_FUNCTION("checkcoregpm (BK 29-Sep-1999)")
	  INFO.print("\n\t*** Input for step COREGPM ***");
	  INFO << "CPM_THRESHOLD: \tThreshold correlation for model: \t" << coregpminput.threshold;
	  INFO.print();
	  INFO << "CPM_DEGREE: \tPolynomial for coregistration: \t" << coregpminput.degree;
	  INFO.print();
	  INFO << "CPM_WEIGHT: \tData weighting option: \t" << coregpminput.weightflag << " (0 none, 1 linear, 2 quadratic, 3 bamler paper) is used.";
	  INFO.print();
	  INFO << "CPM_MAXITER: \tNumber of points to remove (max): \t" << coregpminput.maxiter;
	  INFO.print();
	  INFO << "CPM_K_ALPHA: \tCritical value for outlier test: \t" << coregpminput.k_alpha;
	  INFO.print();
	  if (coregpminput.dumpmodel)
		INFO.print("CPM_DUMP: dumping model to files (see INFO).");
	  if (coregpminput.plot)
		INFO.print("CPM_PLOT: error vectors etc. will be plotted.");
	  else
		WARNING.print("It is higly recommended to use CPM_PLOT to review the estimated model.");
	  if (coregpminput.plotmagbg)
		INFO.print("CPM_PLOT: magnitude will be in background.");
	  } // END checkcoregpm

//***************************************************************
// *    checkcomprefpha                                           *
// *                                                              *
// * Checks cards for step comprefpha.                            *
// *                                                              *
// *    Bert Kampes, 29-Sep-1999                                  *
// ***************************************************************
	public static void checkcomprefpha(input_comprefpha comprefphainput)
	  {
	  TRACE_FUNCTION("checkcomprefpha (BK 29-Sep-1999)")
	  INFO.print("\n\t*** Input for step COMPREFPHA ***");
	  switch (comprefphainput.method)
		{
		case fe_porbits:
		  INFO.print("FE_METHOD: method for flatearth correction: PORBITS");
		  // ______ Check points from file or random distributed ______
		  if (specified(comprefphainput.ifpositions))
			{
			INFO << "FE_IN_POS: file: " << comprefphainput.ifpositions << " containing " << comprefphainput.Npoints << " positions used for ref. phase estimation.";
			INFO.print();
			}
		  else // no input file
			{
			INFO << "FE_NPOINTS: Using " << comprefphainput.Npoints << " (random like) distributed points for estimation of refpha polynomial.";
			INFO.print();
			}
		  if (comprefphainput.Npoints > 5000)
			WARNING.print("FE_NPOINTS: Too many points requested?");
		  INFO << "FE_DEGREE: Using " << comprefphainput.degree << " order polynomial for flat Earth correction.";
		  INFO.print();
		  if (comprefphainput.degree > 10)
			WARNING.print("FE_DEGREE: degree > 10?");
		  break;
		case fe_method2:
		  INFO.print("FE_METHOD: method for flatearth correction: method2 :NOT IMPLEMENTED");
		  break;
		default:
		  PRINT_ERROR("impossible, method name is checked while reading cards.");
		  throw(keyword_error);
		}
	  } // END checkcomprefpha

//***************************************************************
// *    checksubtrrefpha                                          *
// *                                                              *
// * Checks cards for step subtrrefpha.                           *
// *                                                              *
// *    Bert Kampes, 09-Feb-2000                                  *
// ***************************************************************
	public static void checksubtrrefpha(input_subtrrefpha subtrrefphainput)
	  {
	  TRACE_FUNCTION("checksubtrrefpha (BK 09-Feb-2000)")
	  INFO.print("\n\t*** Input for step SUBTRREFPHA ***");
	  // ______ Print info ______
	  switch (subtrrefphainput.method)
		{
		case srp_polynomial:
		  INFO.print("SRP_METHOD: \tpolynomial: \tPolynomial from COMPREFPHA used.");
		  break;
		case srp_exact:
		  INFO.print("SRP_METHOD: \texact:      \treference phase computed foreach pixel.");
		  break;
		default:
		  PRINT_ERROR("impossible, checked above.");
		  throw(keyword_error);
		}
	  if (subtrrefphainput.dumponlyrefpha==false)
		{
		INFO << "SRP_OUT_CINT: \toutputfile complex interferogram: \t" << subtrrefphainput.focint;
		INFO.print();
		}
	  INFO << "SRP_MULTILOOK: \tFactors (line pixel): \t" << subtrrefphainput.multilookL << " " << subtrrefphainput.multilookP;
	  INFO.print();
	  if (subtrrefphainput.dumponlyrefpha==true)
		{
		WARNING.print("SRP_DUMPREFPHA: Only dumping reference phase, no subtraction.");
		INFO << "SRP_DUMPREFPHA: only dump refphase: " << subtrrefphainput.dumponlyrefpha;
		INFO.print();
		INFO << "SRP_OUT_REFPHA: Output file reference phase: " << subtrrefphainput.forefpha;
		INFO.print();
		}
		// __________ added by FvL
	  if (specified(subtrrefphainput.foh2ph))
		{
		INFO << "SRP_OUT_H2PH:   \t" << subtrrefphainput.foh2ph << "; output requested of height-to-phase constant per resolution cell.";
		INFO.print();
		}
	  // ____________ end added by FvL

	  // ______ Check ______ // [MA] increased from 100 to 1000
	  if (subtrrefphainput.multilookL > 1000 || subtrrefphainput.multilookP > 1000 || subtrrefphainput.multilookL < 1 || subtrrefphainput.multilookP < 1)
		WARNING.print("SRP_MULTILOOK: multilookfactor seems very high.");
	  } // END checksubtrrefpha

//***************************************************************
// *    checkresample                                             *
// *                                                              *
// * Checks cards for step resample.                              *
// *                                                              *
// *    Bert Kampes, 29-Sep-1999                                  *
// ***************************************************************
	public static void checkresample(input_resample resampleinput)
	  {
	  TRACE_FUNCTION("checkresample (BK 29-Sep-1999)")
	  INFO.print("\n\t*** Input for step RESAMPLE ***");
	  INFO << "RS_OUT_FILE: \toutput filename: \t\t" << resampleinput.fileout;
	  INFO.print();
	  INFO << "RS_OUT_FORMAT: output format: \t\t";
	  switch (resampleinput.oformatflag)
		{
		case FORMATCR4:
		  INFO << "complex_real4.";
		  break;
		case FORMATCI2:
		  INFO << "complex_short.";
		  break;
		default:
		  PRINT_ERROR("totally impossible, checked input.");
		  throw(keyword_error);
		}
	  INFO.print();
	  if (resampleinput.dbow.linehi != 0 || resampleinput.dbow.pixhi != 0)
		{
		INFO << "RS_DBOW: \tOutput window: \t\t\t" << resampleinput.dbow.linelo << " " << resampleinput.dbow.linehi << " " << resampleinput.dbow.pixlo << " " << resampleinput.dbow.pixhi;
		INFO.print();
		}
	  INFO << "RS_SHIFTAZI: \tshift azimuth spectrum: \t" << resampleinput.shiftazi;
	  INFO.print();
	  } // END checkresample

//***************************************************************
// *    checkinterfero                                            *
// *                                                              *
// * Checks cards for step interfero.                             *
// *                                                              *
// *    Bert Kampes, 29-Sep-1999                                  *
// ***************************************************************
	public static void checkinterfero(input_interfero interferoinput)
	  {
	  TRACE_FUNCTION("checkinterfero (BK 29-Sep-1999)")
	  INFO.print("\n\t*** Input for step INTERFERO ***");
	  boolean filespecified = false;
	  boolean samefilename = false;
	  if (specified(interferoinput.foint))
		{
		filespecified = true;
		INFO << "INT_OUT_INT: \tOutputfile interferogram: \t" << interferoinput.foint;
		INFO.print();
		if (!(strcmp(interferoinput.foint,interferoinput.focint)))
		  samefilename = true;
		}
	  if (specified(interferoinput.focint))
		{
		filespecified = true;
		INFO << "INT_OUT_CINT: Outfile complex interferogram: \t" << interferoinput.focint;
		INFO.print();
		if (!(strcmp(interferoinput.focint,interferoinput.foint)))
		  samefilename = true;
		}
	//  if (strcmp(interferoinput.foflatearth," "))         // something is specified
	//    {
	//    filespecified = true;
	//    INFO << "INT_OUT_FE: data outputfile reference phase: \t"
	//         <<  interferoinput.foflatearth;
	//    INFO.print();
	//    if (!(strcmp(interferoinput.foflatearth,interferoinput.foint))  ||
	//        !(strcmp(interferoinput.foflatearth,interferoinput.focint)))
	//      samefilename = true;
	//    }

	  INFO << "INT_MULTILOOK: \tFactor (line pixel): \t" << interferoinput.multilookL << " " << interferoinput.multilookP;
	  INFO.print();
	  // [MA] increased from 100 to 1000
	  if (interferoinput.multilookL > 1000 || interferoinput.multilookP > 1000)
		{
		PRINT_ERROR("code ???: INT_MULTILOOK: > 1000.");
		throw(keyword_error);
		}
	  if (interferoinput.multilookL < 1 || interferoinput.multilookP < 1)
		{
		PRINT_ERROR("code ???: INT_MULTILOOK: < 1.");
		throw(keyword_error);
		}

	  if (!filespecified)
		{
		PRINT_ERROR("code ???: INT_OUT_*: at least one output file must be specified.");
		throw(keyword_error);
		}
	  if (samefilename)
		{
		PRINT_ERROR("code ???: INT_OUT_*: same name output files.");
		throw(keyword_error);
		}
	  } // END checkinterfero

//***************************************************************
// *    checkcoherence                                            *
// *                                                              *
// * Checks cards for step coherence.                             *
// *                                                              *
// *    Bert Kampes, 29-Sep-1999                                  *
// ***************************************************************
	public static void checkcoherence(input_coherence coherenceinput)
	  {
	  TRACE_FUNCTION("checkcoherence (BK 29-Sep-1999)")
	  INFO.print("\n\t*** Input for step COHERENCE ***");
	  boolean filespecified = false;
	  boolean samefilename = false;
	  if (specified(coherenceinput.foccoh))
		{
		filespecified = true;
		INFO << "COH_OUT_CCOH: \tOutfile complex coherence: \t" << coherenceinput.foccoh;
		INFO.print();
		if (!(strcmp(coherenceinput.foccoh,coherenceinput.focoh)))
		  samefilename = true;
		}

	  if (specified(coherenceinput.focoh))
		{
		filespecified = true;
		INFO << "COH_OUT_COH: \tOutputfile coherence image: " << coherenceinput.focoh;
		INFO.print();
		if (!(strcmp(coherenceinput.focoh,coherenceinput.foccoh)))
		  samefilename = true;
		}

	  INFO << "COH_MULTILOOK: \tFactor (line pixel): \t" << coherenceinput.multilookL << " " << coherenceinput.multilookP;
	  INFO.print();
	  INFO << "COH_WINSIZE: \t window size coh. estimation (l/p): \t" << coherenceinput.cohsizeL << " " << coherenceinput.cohsizeP;
	  INFO.print();
	   // [MA] increased from 100 to 1000
	  if (coherenceinput.multilookL > 1000 || coherenceinput.multilookP > 1000)
		{
		PRINT_ERROR("code ???: COH_MULTILOOK: > 1000.");
		throw(keyword_error);
		}
	  if (coherenceinput.multilookL < 1 || coherenceinput.multilookP < 1)
		{
		PRINT_ERROR("code ???: COH_MULTILOOK: < 1.");
		throw(keyword_error);
		}
	  if (coherenceinput.cohsizeL > 500 || coherenceinput.cohsizeP > 500)
		{
		PRINT_ERROR("code ???: COH_WINSIZE: > 500.");
		throw(keyword_error);
		}
	  if (coherenceinput.cohsizeL < 1 || coherenceinput.cohsizeP < 1)
		{
		PRINT_ERROR("code ???: COH_WINSIZE: < 1.");
		throw(keyword_error);
		}

	  if (!filespecified)
		{
		PRINT_ERROR("code ???: COH_OUT_*: at least one output file must be specified.");
		throw(keyword_error);
		}
	  if (samefilename)
		{
		PRINT_ERROR("code ???: COH_OUT_*: same name output files.");
		throw(keyword_error);
		}
	  } // END checkcoherence

//***************************************************************
// *    checkcomprefdem                                           *
// *                                                              *
// * Checks cards for step comprefdem.                            *
// *                                                              *
// *    Bert Kampes, 14-Feb-2000                                  *
// ***************************************************************
	public static void checkcomprefdem(input_comprefdem comprefdeminput)
	  {
	  TRACE_FUNCTION("checkcomprefdem (BK 14-Feb-2000)")
	  INFO.print("\n\t*** Input for step COMPREFDEM ***");
	//  switch (comprefdeminput.method)
	//    {
	//    case crd_nearest:
	//      INFO.print("NEAREST_NEIGHBOR, use DENSE=2 or so.");
	//      break;
	//    case crd_trilinear:
	//      INFO.print("TRI_LINEAR; use DENSE=0.2 for speed.");
	//      break;
	//    default:
	//      PRINT_ERROR("totally impossible, checked input.");
	//      throw(keyword_error);
	//    }
	  INFO << "CRD_IN_DEM:    \t" << comprefdeminput.firefdem;
	  INFO.print();
	  INFO << "CRD_OUT_FILE:  \t" << comprefdeminput.forefdem;
	  INFO.print();
	  if (specified(comprefdeminput.forefdemhei))
		{
		INFO << "CRD_OUT_DEM_LP: \t" << comprefdeminput.forefdemhei << "; output requested of DEM [m] in radar coordinates.";
		INFO.print();
		if (!strcmp(comprefdeminput.forefdem,comprefdeminput.forefdemhei))
		  {
		  PRINT_ERROR("CRD_OUT_FILE, CRD_OUT_DEM_LP: Same output file name.");
		  throw(keyword_error);
		  }
		}
	  if (specified(comprefdeminput.fodem))
		{
		INFO << "CRD_OUT_DEM:   \t" << comprefdeminput.fodem << "; output requested of input DEM.";
		INFO.print();
		if (!strcmp(comprefdeminput.fodemi,comprefdeminput.fodem))
		  {
		  PRINT_ERROR("OUT_DEM, OUT_DEMI: Same output file name.");
		  throw(keyword_error);
		  }
		}
	  if (specified(comprefdeminput.fodemi))
		{
		INFO << "CRD_OUT_DEMI:   \t" << comprefdeminput.fodemi << "; output requested of interpolated DEM.";
		INFO.print();
		if (!strcmp(comprefdeminput.fodemi,comprefdeminput.fodem))
		  {
		  PRINT_ERROR("OUT_DEM, OUT_DEMI: Same output file name.");
		  throw(keyword_error);
		  }
		}
		// __________ added by FvL
	  if (specified(comprefdeminput.foh2ph))
		{
		INFO << "CRD_OUT_H2PH:   \t" << comprefdeminput.foh2ph << "; output requested of height-to-phase constant per resolution cell.";
		INFO.print();
		}
	  // ____________ end added by FvL
	  INFO << "CRD_IN_SIZE:   \t" << comprefdeminput.demrows << " " << comprefdeminput.demcols << "; number of rows (latitude), columns (lon) in DEM.";
	  INFO.print();
	  INFO << "CRD_IN_UL:     \t" << rad2deg(comprefdeminput.demlatleftupper) << " " << rad2deg(comprefdeminput.demlonleftupper) << "; coordinates of upper left corner (first row/col).";
	  INFO.print();
	  INFO << "CRD_IN_DELTA:  \t" << rad2deg(comprefdeminput.demdeltalat) << " " << rad2deg(comprefdeminput.demdeltalon);
	  INFO.print();
	  INFO << "CRD_IN_NODATA:  \t" << comprefdeminput.demnodata << "; this number in DEM will be set to 0 reference phase.";
	  INFO.print();
	//  INFO << "CRD_DENSE:      \t" << comprefdeminput.extradense
	//       << "; this is the factor for oversampling DEM more than minimum.";
	//  INFO.print();
	  if (comprefdeminput.includerefpha)
		INFO << "CRD_INCLUDE_FE: \tref. dem is computed including flat earth.";
	  else
		INFO << "CRD_INCLUDE_FE: \tref. dem is computed w.r.t. ellipsoid (topo only).";
	  INFO.print();

	  INFO << "CRD_IN_FORMAT: \tinput format DEM: \t";
	  switch (comprefdeminput.iformatflag)
		{
		case FORMATR4:
		  INFO << "real4.";
		  break;
		case FORMATR8:
		  INFO << "real8.";
		  break;
		case FORMATI2:
		  INFO << "signed short.";
		  break;
		case FORMATI2_BIGENDIAN:
		  INFO << "signed short big endian.";
		  break;
		default:
		  PRINT_ERROR("totally impossible, checked input.");
		  throw(keyword_error);
		}
	  INFO.print();



	// ______ Check some things ______
	  if (!existed(comprefdeminput.firefdem))
		{
		ERROR << "CRD_IN_DEM:   \t" << comprefdeminput.firefdem << " can not be opened.";
		PRINT_ERROR(ERROR.get_str())
		throw(keyword_error);
		}
	  if (rad2deg(comprefdeminput.demdeltalat)<.0002) //[MA] 1 arc secs = 1/3600=0.00027 : allow finer resolutions
		{
		WARNING << "CRD_IN_DELTA: \t" << rad2deg(comprefdeminput.demdeltalat) << " [deg] seems very small (input in decimal degrees).";
		WARNING.print();
		}
	  if (comprefdeminput.demrows<1 || comprefdeminput.demrows>100000)
		{
		WARNING << "CRD_DEM_SIZE: numrows: \t" << comprefdeminput.demrows << " seems to be wrong.";
		WARNING.print();
		}
	  if (comprefdeminput.demcols<1 || comprefdeminput.demcols>100000)
		{
		WARNING << "CRD_DEM_SIZE: numcols: \t" << comprefdeminput.demcols << " seems to be wrong.";
		WARNING.print();
		}
	//  if (comprefdeminput.extradense>5.0)
	//    {
	//    WARNING << "CRD_DENSE:    \t" << comprefdeminput.extradense
	//         << " seems to be quite large.";
	//    WARNING.print();
	//    }
	//  if (comprefdeminput.extradense<0.2)
	//    {
	//    WARNING << "CRD_DENSE:    \t" << comprefdeminput.extradense
	//         << " seems too small.";
	//    WARNING.print();
	//    }
	  if (rad2deg(comprefdeminput.demdeltalon)<.0002)
		{
		WARNING << "CRD_IN_DELTA: \t" << comprefdeminput.demdeltalon << " seems very small (it should be in decimal degrees).";
		WARNING.print();
		}
	  if (rad2deg(comprefdeminput.demlatleftupper) < -90. || rad2deg(comprefdeminput.demlatleftupper) > 90.)
		{
		ERROR << "CRD_IN_LU:    \t" << rad2deg(comprefdeminput.demlatleftupper) << " out of range (-90:90).";
		PRINT_ERROR(ERROR.get_str())
		throw(keyword_error);
		}
	  if (rad2deg(comprefdeminput.demlonleftupper) < -180. || rad2deg(comprefdeminput.demlonleftupper) > 180.)
		{
		WARNING << "CRD_IN_LU:    \t" << rad2deg(comprefdeminput.demlonleftupper) << " out of range (-180:180).";
		WARNING.print();
		}
	  if (!strcmp(comprefdeminput.fodem,comprefdeminput.forefdem))
		{
		PRINT_ERROR("OUT_DEM, OUT_FILE: Same output file name.");
		throw(keyword_error);
		}
	  if (!strcmp(comprefdeminput.firefdem,comprefdeminput.forefdem))
		{
		PRINT_ERROR("OUT_FILE, IN_DEM: Same file name.");
		throw(keyword_error);
		}
	  } // END comprefdem

//***************************************************************
// *    checksubtrrefdem                                          *
// *                                                              *
// * Checks cards for step subtrrefdem.                           *
// *                                                              *
// *    Bert Kampes, 14-Feb-2000                                  *
// ***************************************************************
	public static void checksubtrrefdem(input_subtrrefdem subtrrefdeminput)
	  {
	  TRACE_FUNCTION("checksubtrrefdem (BK 14-Feb-2000)")
	  INFO.print("\n\t*** Input for step SUBTRREFDEM ***");
	  INFO << "SRD_OUT_CINT:   \t" << subtrrefdeminput.focint;
	  INFO.print();
	  INFO << "SRD_OFFSET:     \t" << subtrrefdeminput.offsetL << " " << subtrrefdeminput.offsetP;
	  INFO.print();
	  if (Math.abs(subtrrefdeminput.offsetL)>5)
		WARNING.print("Apply offset in azimuth larger than 5 pixels?");
	  if (Math.abs(subtrrefdeminput.offsetP)>5)
		WARNING.print("Apply offset in range larger than 5 pixels?");
	  } // END checksubtrrefdem

//***************************************************************
// *    checkfiltrange                                            *
// * Checks cards for step filtrange.                             *
// *    Bert Kampes, 14-Feb-2000                                  *
// ***************************************************************
	public static void checkfiltrange(input_filtrange filtrangeinput)
	  {
	  TRACE_FUNCTION("checkfiltrange (BK 14-Feb-2000)")
	  INFO.print("\n\t*** Input for step FILTRANGE ***");
	  // ______ Give info ______
	  switch (filtrangeinput.method)
		{
		case rf_adaptive:
		  INFO.print("RF_METHOD:        ADAPTIVE \t(estimate fringe freq.)");
		  INFO << "RF_NLMEAN:     " << filtrangeinput.nlmean;
		  INFO.print();
		  INFO << "RF_THRESHOLD:  " << filtrangeinput.SNRthreshold;
		  INFO.print();
		  INFO << "RF_OVERSAMPLE: " << filtrangeinput.oversample;
		  INFO.print();
		  INFO << "RF_WEIGHTCORR: " << filtrangeinput.doweightcorrel;
		  INFO.print();
		  INFO << "RF_OVERLAP:    " << filtrangeinput.overlap;
		  INFO.print();
		  if (filtrangeinput.nlmean > 51)
			WARNING.print("RF_NLMEAN:     mean over more than 51 lines (?)");
		  if (filtrangeinput.SNRthreshold<1.99)
			WARNING.print("RF_THRESHOLD:  < 2");
		  if (filtrangeinput.SNRthreshold>10.01)
			WARNING.print("RF_THRESHOLD:  > 10 ?");
		  if (filtrangeinput.oversample<=1)
			WARNING.print("RF_OVERSAMPLE: no oversampling.");
		  if (filtrangeinput.oversample>8)
			WARNING.print("RF_OVERSAMPLE: >8 ?");
		  if (!ispower2(filtrangeinput.oversample))
			WARNING.print("RF_OVERSAMPLE: not power of two.");
		  if (filtrangeinput.doweightcorrel==true)
			WARNING.print("RF_WEIGHTCORR: weighting, not sure it has effect.");
		  if (filtrangeinput.fftlength > 1024)
			WARNING.print("RF_FFTLENGTH:  adaptive filterlength > 1024 ?");
		  if (filtrangeinput.SNRthreshold<0.)
			{
			PRINT_ERROR("RF_THRESHOLD:  < 0.");
			throw(keyword_error);
			}
		  break;
		case rf_porbits:
		  INFO.print("RF_METHOD:        PORBITS  \t(based on orbits.)");
		  INFO << "RF_SLOPE:      " << rad2deg(filtrangeinput.terrainslope) << "\t[deg] terrainslope.";
		  INFO.print();
		  if (filtrangeinput.fftlength < 256)
			WARNING.print("RF_FFTLENGTH:  porbits filterlength < 256 (?)");
		  break;
		default:
		  PRINT_ERROR("totally impossible, checked input.");
		  throw(keyword_error);
		}
	  // ______ Both methods cards ______
	  INFO << "RF_FFTLENGTH:  " << filtrangeinput.fftlength;
	  INFO.print();
	  INFO << "RF_HAMMING:    " << filtrangeinput.hammingalpha;
	  INFO.print();
	  INFO << "RF_OUT_MASTER: " << filtrangeinput.fomaster;
	  INFO.print();
	  INFO << "RF_OUT_SLAVE:  " << filtrangeinput.foslave;
	  INFO.print();
	  INFO << "RF_OUT_FORMAT: output format: ";
	  switch (filtrangeinput.oformatflag)
		{
		case FORMATCR4:
		  INFO << "complex_real4.";
		  break;
		case FORMATCI2:
		  INFO << "complex_short.";
		  break;
		default:
		  PRINT_ERROR("totally impossible, checked input.");
		  throw(keyword_error);
		}
	  INFO.print();

	  // ______ Check input ______
	  if (filtrangeinput.hammingalpha>0.999)
		WARNING.print("RF_HAMMING: no hamming filtering.");
	  if (existed(filtrangeinput.fomaster))
		WARNING.print("RF_OUT_MASTER: file exists.");
	  if (existed(filtrangeinput.foslave))
		WARNING.print("RF_OUT_SLAVE: file exists.");
	  if (!ispower2(filtrangeinput.fftlength))
		{
		PRINT_ERROR("RF_FFTLENGTH: not power of 2.");
		throw(keyword_error);
		}
	  if (filtrangeinput.overlap >= 0.5 *filtrangeinput.fftlength)
		{
		PRINT_ERROR("RF_OVERLAP >= 0.5*RF_FFTLENGTH");
		throw(keyword_error);
		}
	  if (filtrangeinput.hammingalpha>1. || filtrangeinput.hammingalpha<0.)
		{
		PRINT_ERROR("RF_HAMMING: not e[0,1].");
		throw(keyword_error);
		}
	  } // END checkfiltrange

//***************************************************************
// *    checkdinsar                                               *
// *                                                              *
// * Checks cards for step dinsar.                                *
// *                                                              *
// #%// BK 25-Sep-2000
// ***************************************************************
	public static void checkdinsar(input_dinsar dinsarinput)
	  {
	  TRACE_FUNCTION("checkdinsar (BK 25-Sep-2000)")
	  INFO.print("\n\t*** Input for step DINSAR ***");

	  if (!specified(dinsarinput.topomasterresfile))
		{
		INFO.print("Using 3 pass differential (for 4 pass, see DI_IN_TOPOMASTER card).");
		}
	  else
		{
		INFO << "DI_IN_TOPOMASTER: \t" << dinsarinput.topomasterresfile << " (4 pass)";
		INFO.print();
		}
	  INFO << "DI_IN_TOPOSLAVE: \t" << dinsarinput.toposlaveresfile;
	  INFO.print();
	  INFO << "DI_IN_TOPOINT:   \t" << dinsarinput.topointresfile;
	  INFO.print();
	  INFO << "DI_OUT_FILE:     \t" << dinsarinput.fodinsar;
	  INFO.print();
	  if (!specified(dinsarinput.foscaleduint))
		INFO << "DI_OUT_SCALED: \tNo (debug) output requested scaled topography interf.";
	  else
		INFO << "DI_OUT_SCALED: \t" << dinsarinput.foscaleduint << "; debug output requested.";
	  INFO.print();
	  if (!specified(dinsarinput.toposlaveresfile))
		{
		PRINT_ERROR("DI_IN_TOPOSLAVE: result file topo slave not specified.");
		throw(keyword_error);
		}
	  if (!specified(dinsarinput.topointresfile))
		{
		PRINT_ERROR("DI_IN_TOPOINT: result file topo interferogram not specified.");
		throw(keyword_error);
		}
	  if (!strcmp(dinsarinput.toposlaveresfile,dinsarinput.topointresfile))
		{
		PRINT_ERROR("IN_TOPOSLAVE, IN_TOPOINT: Same input file name.");
		throw(keyword_error);
		}
	  } // END checkdinsar

//***************************************************************
// *    checkfiltphase                                            *
// *                                                              *
// * Checks cards for step filtphase.                             *
// *                                                              *
// #%// BK 25-Sep-2000
// ***************************************************************
	public static void checkfiltphase(input_filtphase filtphaseinput)
	  {
	  TRACE_FUNCTION("checkfiltphase (BK 25-Sep-2000)")
	  INFO.print("\n\t*** Input for step FILTPHASE ***");
	  if (specified(filtphaseinput.fifiltphase))
		{
		INFO << "PF_IN_FILE: \t" << filtphaseinput.fifiltphase << " " << filtphaseinput.finumlines << " (this cr4 file will be filtered)";
		INFO.print();
		if (!existed(filtphaseinput.fifiltphase))
		  WARNING.print("Impossible? PF input file does not exist?");
		}
	  INFO << "PF_OUT_FILE: \t" << filtphaseinput.fofiltphase << " (output filename).";
	  INFO.print();

	  // ______ Method goldstein filter ______
	  if (filtphaseinput.method==fp_goldstein)
		{
		INFO.print("FILTPHASE: Method goldstein.");
		INFO << "PF_ALPHA: \t" << filtphaseinput.alpha << " (weigthing parameter for spectrum).";
		INFO.print();
		INFO << "PF_BLOCKSIZE: " << filtphaseinput.blocksize << " (size of block to perform filtering on).";
		INFO.print();
		INFO << "PF_OVERLAP: \t" << filtphaseinput.overlap << " (half overlap between consequetive blocks).";
		INFO.print();

		// ______ Use 1d kernel to smooth amplitude, e.g. 12321 ______
		// ______ Which is normalized by me ______
		INFO << "PF_KERNEL: \t";
		for (int32 ii =0; ii<filtphaseinput.kernel.pixels(); ++ii)
		  INFO << " " << filtphaseinput.kernel(0,ii);
		INFO << " (smooth |spectrum| with this).";
		INFO.print();
		if (filtphaseinput.kernel.pixels()==1)
		  INFO.print("No smoothing of amplitude spectrum!");

		// ______ Check errors _____
		if (filtphaseinput.alpha<0. || filtphaseinput.alpha>1.)
		  WARNING.print("PF_ALPHA not 0<a<1");
		if (filtphaseinput.blocksize>64 || filtphaseinput.blocksize<16)
		  WARNING.print("PF_BLOCKSIZE very small or large?");
		if (filtphaseinput.kernel.pixels()>11)
		  WARNING.print("smoothing kernel > 11:  very large?");
		if (filtphaseinput.overlap<0)
		  {
		  PRINT_ERROR("PF_OVERLAP < 0");
		  throw(keyword_error);
		  }
		if (2 *filtphaseinput.overlap>filtphaseinput.blocksize)
		  {
		  PRINT_ERROR("2*PF_OVERLAP > PF_BLOCKSIZE");
		  throw(keyword_error);
		  }
		if (filtphaseinput.kernel.pixels()>filtphaseinput.blocksize)
		  {
		  PRINT_ERROR("smoothing kernel > PF_BLOCKSIZE");
		  throw(keyword_error);
		  }
		if (!ispower2(filtphaseinput.blocksize))
		  {
		  PRINT_ERROR("PF_BLOCKSIZE not a power of 2");
		  throw(keyword_error);
		  }
		}

	  // ______ Method spatial convolution ______
	  else if (filtphaseinput.method==fp_spatialconv)
		{
		INFO.print("FILTPHASE: Method spatial convolution.");
		if (!specified(filtphaseinput.fikernel2d))
		  {
		  INFO.print("Using 1d kernel for spatial convolution (no PF_IN_KERNEL2D).");
		  INFO << "PF_KERNEL: used: \t";
		  for (int32 ii =0; ii<filtphaseinput.kernel.pixels(); ++ii)
			INFO << " " << filtphaseinput.kernel(0,ii);
		  INFO.print();
		  }
		else
		  {
		  INFO.print("Using 2d kernel for spatial convolution.");
		  INFO << "PF_IN_KERNEL2D: \t" << filtphaseinput.fikernel2d << " (ascii input file with 2d kernel).";
		  INFO.print();
		  INFO.print("PF_IN_KERNEL2D: \t(input file has 1 line header: numrows numcols scale");
		  if (filtphaseinput.kernel.size()!=0)
			WARNING.print("PF_KERNEL card ignored due to card PF_IN_KERNEL2D.");
		  if (!existed(filtphaseinput.fikernel2d))
			WARNING.print("PF_IN_KERNEL2D infile cannot be found.");
		  }
		}
	  else if (filtphaseinput.method==fp_spectral)
		{
		INFO.print("FILTPHASE: Method spectral filter with 2D kernel.");
		INFO << "PF_BLOCKSIZE: " << filtphaseinput.blocksize << " (size of block to perform filtering on).";
		INFO.print();
		INFO << "PF_OVERLAP: \t" << filtphaseinput.overlap << " (half overlap between consequetive blocks).";
		INFO.print();
		if (filtphaseinput.kernel.size()!=0)
		  WARNING.print("PF_KERNEL card ignored for method spectral.");
		if (filtphaseinput.overlap<0)
		  {
		  PRINT_ERROR("PF_OVERLAP < 0");
		  throw(keyword_error);
		  }
		if (2 *filtphaseinput.overlap>filtphaseinput.blocksize)
		  {
		  PRINT_ERROR("2*PF_OVERLAP > PF_BLOCKSIZE");
		  throw(keyword_error);
		  }
		if (!ispower2(filtphaseinput.blocksize))
		  {
		  PRINT_ERROR("PF_BLOCKSIZE not a power of 2");
		  throw(keyword_error);
		  }
		if (!specified(filtphaseinput.fikernel2d))
		  {
		  PRINT_ERROR("method spectral needs card PF_IN_KERNEL2D");
		  throw(keyword_error);
		  }
		}
	  else
		{
		PRINT_ERROR("Method phasefiltering != {goldstein,spatialconv,spectral}.");
		throw(keyword_error);
		}
	  } // END checkfiltphase

//***************************************************************
// *    checkfiltazi                                              *
// * Checks cards for step filtazi.                               *
// *    Bert Kampes, 02-Nov-2000                                  *
// ***************************************************************
	public static void checkfiltazi(input_filtazi filtaziinput, int16 id)
	  {
	  TRACE_FUNCTION("checkfiltazi (BK 02-Nov-2000)")
	  INFO.print("\n\t*** Input for step FILTAZI ***");
	  INFO << "AF_BLOCKSIZE:   \t" << filtaziinput.fftlength;
	  INFO.print();
	  INFO << "AF_OVERLAP:     \t" << filtaziinput.overlap;
	  INFO.print();
	  INFO << "AF_HAMMING:     \t" << filtaziinput.hammingalpha;
	  INFO.print();
	  if (filtaziinput.oformatflag==FORMATCR4)
		INFO.print("AF_OUT_FORMAT:  \tcomplex_real4");
	  else if (filtaziinput.oformatflag==FORMATCI2)
		INFO.print("AF_OUT_FORMAT:  \tcomplex_short");
	  else
		{
		PRINT_ERROR("formatflag not ok for output.");
		throw(keyword_error);
		}

	  if (id!=SLAVEID)
		{
		INFO << "AF_OUT_MASTER:  \t" << filtaziinput.fomaster;
		INFO.print();
		}
	  if (id!=MASTERID)
		{
		INFO << "AF_OUT_SLAVE:   \t" << filtaziinput.foslave;
		INFO.print();
		}

	  if (filtaziinput.fftlength<256)
		WARNING.print("AF_BLOCKSIZE < 256 (too small?)");
	  if (2 *filtaziinput.overlap>0.5 *filtaziinput.fftlength)
		WARNING.print("2*AF_OVERLAP > .5*AF_BLOCKSIZE (very large?)");
	  if (filtaziinput.hammingalpha<0 || filtaziinput.hammingalpha>1)
		{
		PRINT_ERROR("AF_HAMMING not e[0,1]");
		throw(keyword_error);
		}
	  if (filtaziinput.overlap<0)
		{
		PRINT_ERROR("AF_BLOCKSIZE < 0");
		throw(keyword_error);
		}
	  if (2 *filtaziinput.overlap>filtaziinput.fftlength)
		{
		PRINT_ERROR("AF_BLOCKSIZE < 2*AF_BLOCKSIZE");
		throw(keyword_error);
		}
	  if (!ispower2(filtaziinput.fftlength))
		{
		PRINT_ERROR("AF_BLOCKSIZE must be power of 2.");
		throw(keyword_error);
		}
	  } // END checkfiltazi
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>

	//***************************************************************
	// *    writearg                                                  *
	// * echo arg to screen if debug defined                          *
	// #%// BK 13-Jul-2000                                            *
	// ***************************************************************
	public static <Type> void writearg(Type argument)
	  {
	  TRACE_FUNCTION("writearg");
	  DEBUG << "Read argument: " << argument;
	  DEBUG.print();
	  } // END writearg



	//***************************************************************
	// *    readinput                                                 *
	// *                                                              *
	// * Read and interpret                                           *
	// * "inputoptionsfile" (file in variable: logfile).              *
	// * Write to logfile (via "scratchreadinput")                    *
	// * mandatory input is checked on presence by checksums (cs).    *
	// * If there are more methods, process card switches default     *
	// *  methodselector (if present) is used later to correct        *
	// *  see for example unwrap methodeselector                      *
	// * linecounter only used in case of errors                      *
	// *                                                              *
	// * input:                                                       *
	// *  - struct: general input options                             *
	// *  - struct: master readfiles input options                    *
	// *  - struct: slave readfiles input options                     *
	// *  - struct: master crop input options                         *
	// *  - struct: slave crop input options                          *
	// *  - ...                                                       *
	// * output:                                                      *
	// *  - (updated input structs)                                   *
	// *  - (file copy of input)                                      *
	// *                                                              *
	// *    Bert Kampes,      11-Dec-1998                             *
	// *    Mahmut Arikan,    09-Jan-2009 (code to fix linecnt bug    *
	// *    G.J. van Zwieten, 09-Jan-2009  and line extend         )  *
	// ***************************************************************
	//____RaffaeleNutricato START MODIFICATION SECTION 2
	//____RaffaeleNutricato END MODIFICATION SECTION 2
	//____RaffaeleNutricato START MODIFICATION SECTION 3
	//____RaffaeleNutricato END MODIFICATION SECTION 3
	public static void readinput(RefObject<input_gen> generalinput, RefObject<input_ell> ellipsinput, RefObject<input_pr_orbits> porbitsinput, RefObject<input_readfiles> m_readfilesinput, RefObject<input_crop> m_cropinput, RefObject<input_oversample> m_oversample, RefObject<input_simamp> simampinput, RefObject<input_mtiming> mtiminginput, RefObject<input_readfiles> s_readfilesinput, RefObject<input_crop> s_cropinput, RefObject<input_oversample> s_oversample, RefObject<input_filtazi> filtaziinput, RefObject<input_coarsecorr> coarsecorrinput, RefObject<input_fine> fineinput, RefObject<input_reltiming> reltiminginput, RefObject<input_demassist> demassistinput, RefObject<input_coregpm> coregpminput, RefObject<input_resample> resampleinput, RefObject<input_filtrange> filtrangeinput, RefObject<input_interfero> interferoinput, RefObject<input_coherence> coherenceinput, RefObject<input_comprefpha> comprefphainput, RefObject<input_subtrrefpha> subtrrefphainput, RefObject<input_comprefdem> comprefdeminput, RefObject<input_subtrrefdem> subtrrefdeminput, RefObject<input_filtphase> filtphaseinput, RefObject<input_dinsar> dinsarinput, RefObject<input_unwrap> unwrapinput, RefObject<input_slant2h> slant2hinput, RefObject<input_geocode> geocodeinput)
	  {
	  //TRACE_FUNCTION("readinput (BK 11-Dec-1998)");
	  TRACE_FUNCTION("readinput rev.4 (TUDelft 09-Jan-2009)");

	  // ______ Set ids ______
	  m_readfilesinput.argvalue.fileid = MASTERID;
	  m_cropinput.argvalue.fileid = MASTERID;
	  s_readfilesinput.argvalue.fileid = SLAVEID;
	  s_cropinput.argvalue.fileid = SLAVEID;

	  // ______ Misuse generalinput.logfile to store name: open file here! ______
	  //ifstream optionsfile(generalinput.logfile, ios::in | ios::nocreate);
	  ifstream optionsfile = new ifstream(generalinput.argvalue.logfile, ios.in);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(optionsfile, generalinput.argvalue.logfile, __FILE__, __LINE__);
	  String inputoptionsfile = generalinput.argvalue.logfile; // [MA] keep input filename
														// later variable is updated
														// as logfile name.

	  int16 onlyprocess = -1; // flag for ONLYPROCESS card
	  int16 linecnt = 0; // counter
	  final int16 BASE10 = 10; // [MA] base 10 defined for strtol
	  //char                  keyword[EIGHTY];
	  //char                  filename[4*ONE27];            // string for filenames // [MA] changed EIGHTY --> 2*ONE27, due to comments line gets longer
	  String eachline = new String(new char[4 *ONE27]); // assuming maximum char lenght of the line is 4*ONE27. It should be sufficient.

	  // ______ Check (multiple) occurence of cards ______
	  boolean priorscreen = false; // no screen card present
	  boolean priormemory = false; // check if present for info
	  boolean priorbatch = false; // check if present for info
	  boolean prioroverwrite = false; // check if present for info
	  boolean priorlistinput = false; // check if present for info
	  boolean priorrs_fileout = false;


	// ====== Initialization, defaults ======
	  input_ell WGS84;
	  boolean listinput = true; // default copy input to log
	  boolean ellipsoid = false; // set default if no card present
	  ellipsinput.argvalue = WGS84; // default

//C++ TO JAVA CONVERTER NOTE: 'register' variable declarations are not supported in Java:
//ORIGINAL LINE: register int32 i;
	  int32 i;
	  for (i =0; i<NUMPROCESSES; i++)
		generalinput.argvalue.process[i] = 0; // default no processing step
	  generalinput.argvalue.interactive = true; // default interactive mode
	  generalinput.argvalue.overwrit = false; // default no overwriting
	  generalinput.argvalue.memory = 500000000; // default memory (500MB)
	  generalinput.argvalue.logfile = "log.out"; // default logfile
	  generalinput.argvalue.m_resfile = "master_result.out"; // default resultfile
	  generalinput.argvalue.s_resfile = "slave_result.out"; // default resultfile
	  generalinput.argvalue.i_resfile = "interferogram.out"; // default interf_out
	  generalinput.argvalue.orb_interp = ORB_DEFAULT; // default polyfit
	  generalinput.argvalue.dumpbaselineL = 0; // default no dump
	  generalinput.argvalue.dumpbaselineP = 0; // default no dump
	  generalinput.argvalue.preview = 0; // default no preview
	  generalinput.argvalue.terrain_height = 0.0; // above ellipsoid
	  generalinput.argvalue.tiepoint.x = 0.0; // default none
	  generalinput.argvalue.tiepoint.y = 0.0; // default none
	  generalinput.argvalue.tiepoint.z = 0.0; // default none

	  //m_readfilesinput.method     = readfiles_ers;        // default
	  m_readfilesinput.argvalue.sensor_id = SLC_ERS; // default
	  setunspecified(m_readfilesinput.argvalue.volfile); // check later, then set default
	  setunspecified(m_readfilesinput.argvalue.leaderfile); // check later, then set default
	  setunspecified(m_readfilesinput.argvalue.nullfile); // check later, then set default
	  setunspecified(m_readfilesinput.argvalue.datfile); // check later, then set default
	  m_readfilesinput.argvalue.rg_timing_error = 0.0; // default
	  m_readfilesinput.argvalue.az_timing_error = 0.0; // default
	  //s_readfilesinput.sensor_id          = readfiles_ers;        // default
	  s_readfilesinput.argvalue.sensor_id = SLC_ERS; // default
	  setunspecified(s_readfilesinput.argvalue.volfile); // check later, then set default
	  setunspecified(s_readfilesinput.argvalue.leaderfile); // check later, then set default
	  setunspecified(s_readfilesinput.argvalue.nullfile); // check later, then set default
	  setunspecified(s_readfilesinput.argvalue.datfile); // check later, then set default
	  s_readfilesinput.argvalue.rg_timing_error = 0.0; // default
	  s_readfilesinput.argvalue.az_timing_error = 0.0; // default

	  setunspecified(porbitsinput.argvalue.m_orbdir); // check later, then set default
	  setunspecified(porbitsinput.argvalue.s_orbdir); // check later, then set default
	  porbitsinput.argvalue.timeinterval = 1; // default time interval
	  porbitsinput.argvalue.timebefore = 4 *porbitsinput.argvalue.timeinterval; // default 4 extra datapoints
	  porbitsinput.argvalue.dumpmasterorbit = -1.; // default no dump
	  porbitsinput.argvalue.dumpslaveorbit = -1.; // default no dump

	  m_cropinput.argvalue.idcrop = "master step01"; // default identifier
	  s_cropinput.argvalue.idcrop = "slave step01"; // default identifier
	  m_cropinput.argvalue.fileout1 = "master.raw"; // default output filename
	  s_cropinput.argvalue.fileout1 = "slave.raw"; // default output filename
	  m_cropinput.argvalue.dbow.linelo = 0; // min. line coord. initialization
	  m_cropinput.argvalue.dbow.linehi = 0; // max. line coord. initialization
	  m_cropinput.argvalue.dbow.pixlo = 0; // min. pixel coord. initialization
	  m_cropinput.argvalue.dbow.pixhi = 0; // max. pixel coord. initialization
	  s_cropinput.argvalue.dbow.linelo = 0; // min. line coord. initialization
	  s_cropinput.argvalue.dbow.linehi = 0; // max. line coord. initialization
	  s_cropinput.argvalue.dbow.pixlo = 0; // min. pixel coord. initialization
	  s_cropinput.argvalue.dbow.pixhi = 0; // max. pixel coord. initialization
	  m_cropinput.argvalue.dbow_geo.linelo = 0; // min. line coord. initialization
	  m_cropinput.argvalue.dbow_geo.linehi = 0; // max. line coord. initialization
	  m_cropinput.argvalue.dbow_geo.pixlo = 0; // min. pixel coord. initialization
	  m_cropinput.argvalue.dbow_geo.pixhi = 0; // max. pixel coord. initialization
	  s_cropinput.argvalue.dbow_geo.linelo = 0; // min. line coord. initialization
	  s_cropinput.argvalue.dbow_geo.linehi = 0; // max. line coord. initialization
	  s_cropinput.argvalue.dbow_geo.pixlo = 0; // min. pixel coord. initialization
	  s_cropinput.argvalue.dbow_geo.pixhi = 0; // max. pixel coord. initialization

	//____RaffaeleNutricato START MODIFICATION SECTION 4
	  m_oversample.argvalue.OsrRange = 1; // Default for oversampling ratio in range.
	  m_oversample.argvalue.OsrAzimuth = 1; // Default for oversampling ratio in azimuth.
	  m_oversample.argvalue.FilterSize = 16; // Default for length of the interpolation kernel in range.
	  m_oversample.argvalue.oformatflag = FORMATCI2; // Default for output format.
	  m_oversample.argvalue.fileoutovs = "master_ovs.raw"; // Default output filename.
	  s_oversample.argvalue.OsrRange = 1; // Default for oversampling ratio in range.
	  s_oversample.argvalue.OsrAzimuth = 1; // Default for oversampling ratio in azimuth.
	  s_oversample.argvalue.FilterSize = 16; // Default for length of the interpolation kernel in range.
	  s_oversample.argvalue.oformatflag = FORMATCI2; // Default for output format.
	  s_oversample.argvalue.fileoutovs = "slave_ovs.raw"; // Default output filename.
	//____RaffaeleNutricato END MODIFICATION SECTION 4

	  //____ added by MA ____

	  setunspecified(simampinput.argvalue.firefdem); // check later, mandatory
	  //setunspecified(simampinput.fodemi);                                // check later, then set default
	  //setunspecified(simampinput.forefdemhei);                           // check later, then set default
	  simampinput.argvalue.fodem = "demcrop_sam.raw"; // default name
	  setunspecified(simampinput.argvalue.fosimamp); // check later, then set default
	  simampinput.argvalue.iformatflag = FORMATI2; // default gtopo30
	  simampinput.argvalue.demrows = 6000; // default gtopo30
	  simampinput.argvalue.demcols = 4800; // default gtopo30
	  simampinput.argvalue.demnodata = -9999; // default gtopo30
	  simampinput.argvalue.demdeltalat = deg2rad(0.00833333333333333333); // default gtopo30
	  simampinput.argvalue.demdeltalon = deg2rad(0.00833333333333333333); // default gtopo30
	  simampinput.argvalue.demlatleftupper = deg2rad(89.995833333333333333); // w020n90.DEM
	  simampinput.argvalue.demlonleftupper = deg2rad(-19.995833333333333333); // w020n90.DEM
	  simampinput.argvalue.fosimamp = "master.sam"; // default name


	  final int32 def_mte_nwin = 16; // default #windows
	  setunspecified(mtiminginput.argvalue.ifpositions); // check later, then set default
	  mtiminginput.argvalue.MasksizeL = 256; // default correlation size
	  //mtiminginput.MasksizeP     = mtiminginput.MasksizeL;                // default correlation size
	  mtiminginput.argvalue.MasksizeP = 128; // default correlation size
	  mtiminginput.argvalue.AccL = 32; // default searching limit
	  mtiminginput.argvalue.AccP = mtiminginput.argvalue.AccL; // default searching limit
	  mtiminginput.argvalue.initoffsetL = 0; // default initial offset
	  mtiminginput.argvalue.initoffsetP = 0; // default initial offset

	  // ____ end added by MA ____

	  filtaziinput.argvalue.oformatflag = FORMATCR4; // default
	  filtaziinput.argvalue.fftlength = 1024; // default
	  filtaziinput.argvalue.overlap = -1; // default to fftlength/8
	  filtaziinput.argvalue.hammingalpha = 0.75; // default (slc)
	  filtaziinput.argvalue.fomaster = "master.afilter"; // default
	  filtaziinput.argvalue.foslave = "slave.afilter"; // default

	  final int32 def_cc_nwin = 11; // default #windows
	  setunspecified(coarsecorrinput.argvalue.ifpositions); // check later, then set default
	  coarsecorrinput.argvalue.MasksizeL = 64; // default correlation size
	  coarsecorrinput.argvalue.MasksizeP = coarsecorrinput.argvalue.MasksizeL; // default correlation size
	  coarsecorrinput.argvalue.AccL = 8; // default searching limit
	  coarsecorrinput.argvalue.AccP = coarsecorrinput.argvalue.AccL; // default searching limit
	  coarsecorrinput.argvalue.initoffsetL = 0; // default initial offset
	  coarsecorrinput.argvalue.initoffsetP = 0; // default initial offset

	  setunspecified(fineinput.argvalue.ifpositions); // check later, then set default
	  final int32 def_fc_nwin = 601; // default #windows
	  fineinput.argvalue.MasksizeL = 64; // default correlation size
	  fineinput.argvalue.MasksizeP = fineinput.argvalue.MasksizeL; // default correlation size
	  fineinput.argvalue.AccL = 8; // default searching limit
	  fineinput.argvalue.AccP = fineinput.argvalue.AccL; // default searching limit
	  fineinput.argvalue.initoffsetL = 0; // default initial offset
	  fineinput.argvalue.initoffsetP = 0; // default initial offset
	  fineinput.argvalue.osfactor = 32; // default oversampling factor
	  fineinput.argvalue.plotoffsets = false; // default no plotting
	  fineinput.argvalue.plotmagbg = false; // default no plotting
	  fineinput.argvalue.plotthreshold = 0.3; // default no plotting

	  //____ added by FvL ____
	  reltiminginput.argvalue.threshold = 0.4; // default threshold data
	  reltiminginput.argvalue.maxiter = 10000; // default max. 10000 outliers removed
	  reltiminginput.argvalue.k_alpha = 1.97; // critical value for outliers

	  setunspecified(demassistinput.argvalue.firefdem); // check later, mandatory
	  setunspecified(demassistinput.argvalue.fodemi); // check later, then set default
	  setunspecified(demassistinput.argvalue.forefdemhei); // check later, then set default
	  demassistinput.argvalue.fodem = "demcrop.raw"; // default name
	  demassistinput.argvalue.iformatflag = FORMATI2; // default gtopo30
	  demassistinput.argvalue.demrows = 6000; // default gtopo30
	  demassistinput.argvalue.demcols = 4800; // default gtopo30
	  demassistinput.argvalue.demnodata = -9999; // default gtopo30
	  demassistinput.argvalue.demdeltalat = deg2rad(0.00833333333333333333); // default gtopo30
	  demassistinput.argvalue.demdeltalon = deg2rad(0.00833333333333333333); // default gtopo30
	  demassistinput.argvalue.demlatleftupper = deg2rad(89.995833333333333333); // w020n90.DEM
	  demassistinput.argvalue.demlonleftupper = deg2rad(-19.995833333333333333); // w020n90.DEM
	  // ____ end added by FvL ____

	  coregpminput.argvalue.threshold = 0.2; // default threshold data
	  coregpminput.argvalue.degree = 1; // default degree polynomial
	  coregpminput.argvalue.weightflag = 0; // default no weighting of data
	  coregpminput.argvalue.maxiter = 10000; // default max. 10000 outliers removed, changed by [FvL]
	  coregpminput.argvalue.k_alpha = 1.97; // critical value for outliers
	  coregpminput.argvalue.dumpmodel = false; // default no files
	  coregpminput.argvalue.plot = false; // default no plots
	  coregpminput.argvalue.plotmagbg = false; // default no plots

	  filtrangeinput.argvalue.method = rf_adaptive; // default
	  filtrangeinput.argvalue.terrainslope = 0.0; // default porbits
	  filtrangeinput.argvalue.fftlength = -999; // set default later
	  filtrangeinput.argvalue.overlap = 0; // default no overlap
	  filtrangeinput.argvalue.nlmean = 15; // default
	  filtrangeinput.argvalue.hammingalpha = 0.75; // default
	  filtrangeinput.argvalue.SNRthreshold = 5.0; // default
	  filtrangeinput.argvalue.oversample = 2; // default
	  filtrangeinput.argvalue.doweightcorrel = false; // default
	  filtrangeinput.argvalue.fomaster = "master.rfilter"; // default
	  filtrangeinput.argvalue.foslave = "slave.rfilter"; // default
	  filtrangeinput.argvalue.oformatflag = FORMATCR4; // default

	  setunspecified(comprefphainput.argvalue.ifpositions); // check later, then set default
	  final int32 def_fe_Npoints = 501; // default
	  final int32 def_fe_degree = 5; // default

	  resampleinput.argvalue.fileout = "s_resampled.raw"; // default
	  resampleinput.argvalue.oformatflag = FORMATCR4; // default
	  resampleinput.argvalue.dbow.linelo = 0; // min. line coord. initialization
	  resampleinput.argvalue.dbow.linehi = 0; // max. line coord. initialization
	  resampleinput.argvalue.dbow.pixlo = 0; // min. pixel coord. initialization
	  resampleinput.argvalue.dbow.pixhi = 0; // max. pixel coord. initialization
	  resampleinput.argvalue.dbow_geo.linelo = 0; // min. line coord. initialization
	  resampleinput.argvalue.dbow_geo.linehi = 0; // max. line coord. initialization
	  resampleinput.argvalue.dbow_geo.pixlo = 0; // min. pixel coord. initialization
	  resampleinput.argvalue.dbow_geo.pixhi = 0; // max. pixel coord. initialization
	  resampleinput.argvalue.shiftazi = true; // default apply shift

	  interferoinput.argvalue.method = int_oldmethod; // default method
	  setunspecified(interferoinput.argvalue.focint); // check later, then set default
	  setunspecified(interferoinput.argvalue.foint); // use later if specified
	  //setunspecified(interferoinput.foflatearth);         // check later, then set default
	  interferoinput.argvalue.multilookL = 5; // default multilookfactor
	  interferoinput.argvalue.multilookP = 1; // default multilookfactor

	  coherenceinput.argvalue.method = coh_oldmethod; // default method
	  setunspecified(coherenceinput.argvalue.focoh); // check later, then set default
	  setunspecified(coherenceinput.argvalue.foccoh); // check later, then set default
	  coherenceinput.argvalue.multilookL = 10; // default multilookfactor
	  coherenceinput.argvalue.multilookP = 2; // default multilookfactor
	  coherenceinput.argvalue.cohsizeL = coherenceinput.argvalue.multilookL; // default windowsize.
	  coherenceinput.argvalue.cohsizeP = coherenceinput.argvalue.multilookP; // default windowsize.

	  subtrrefphainput.argvalue.method = srp_polynomial; // default method
	  subtrrefphainput.argvalue.forefpha = "refphase.raw"; // default name
	  subtrrefphainput.argvalue.focint = "cint.minrefpha.raw"; // default
	  subtrrefphainput.argvalue.multilookL = 1; // default multilookfactor
	  subtrrefphainput.argvalue.multilookP = 1; // default multilookfactor
	  subtrrefphainput.argvalue.dumponlyrefpha = false; // default not
	  //_____ added by FvL
	  setunspecified(subtrrefphainput.argvalue.foh2ph); // check later, then set default
	  // ____ end added by FvL

	  filtphaseinput.argvalue.method = fp_goldstein; // default method
	  filtphaseinput.argvalue.alpha = 0.2; // default
	  // ______ 32 blocks default for goldstein, kernel as large as possible, ______
	  // ______ 2Dkernel then specify ______
	  filtphaseinput.argvalue.blocksize = 32; // default
	  filtphaseinput.argvalue.overlap = 3; // default
	  // set default later with alpha in name
	  setunspecified(filtphaseinput.argvalue.fofiltphase); // check later, then set default
	  setunspecified(filtphaseinput.argvalue.fifiltphase); // if specified, use it
	  setunspecified(filtphaseinput.argvalue.fikernel2d); // if specified, use it
	  filtphaseinput.argvalue.finumlines = 0; // numlines.

	  dinsarinput.argvalue.fodinsar = "differentialinterf.raw"; // default
	  setunspecified(dinsarinput.argvalue.foscaleduint); // default no output
	  // ______ set mastertopo file to same as master resfile if not specified ______
	  setunspecified(dinsarinput.argvalue.topomasterresfile); // if specified, use 4 pass
	  setunspecified(dinsarinput.argvalue.toposlaveresfile); // check later, mandatory
	  setunspecified(dinsarinput.argvalue.topointresfile); // check later, mandatory

	  setunspecified(comprefdeminput.argvalue.firefdem); // check later, mandatory
	  setunspecified(comprefdeminput.argvalue.fodemi); // check later, then set default
	  //_____ added by FvL
	  setunspecified(comprefdeminput.argvalue.foh2ph); // check later, then set default
	  // ____ end added by FvL
	  setunspecified(comprefdeminput.argvalue.forefdemhei); // check later, then set default
	  comprefdeminput.argvalue.forefdem = "refdem.raw"; // default name
	  comprefdeminput.argvalue.fodem = "demcrop.raw"; // default name [FvL]
	  comprefdeminput.argvalue.iformatflag = FORMATI2; // default gtopo30
	//  comprefdeminput.method        = crd_trilinear;      // default method
	//  comprefdeminput.extradense    = 0.5;                        // default (now interpolated in l,p)
	  comprefdeminput.argvalue.demrows = 6000; // default gtopo30
	  comprefdeminput.argvalue.demcols = 4800; // default gtopo30
	  comprefdeminput.argvalue.demnodata = -9999; // default gtopo30
	  comprefdeminput.argvalue.demdeltalat = deg2rad(0.00833333333333333333); // default gtopo30
	  comprefdeminput.argvalue.demdeltalon = deg2rad(0.00833333333333333333); // default gtopo30
	  comprefdeminput.argvalue.demlatleftupper = deg2rad(89.995833333333333333); // w020n90.DEM
	  comprefdeminput.argvalue.demlonleftupper = deg2rad(-19.995833333333333333); // w020n90.DEM

	  subtrrefdeminput.argvalue.focint = "cint.minrefdem.raw"; // default name
	  subtrrefdeminput.argvalue.offsetL = 0; // default no offset
	  subtrrefdeminput.argvalue.offsetP = 0; // default no offset

	  unwrapinput.argvalue.fouint = "unwrapped_interferogram.raw"; // default name
	  unwrapinput.argvalue.foregions = "regions_unwrapped.raw"; // default name
	  setunspecified(unwrapinput.argvalue.seedfile); // check later, then set default
	  unwrapinput.argvalue.deltaLseed = 100; // default 100 pixels;
	  unwrapinput.argvalue.deltaPseed = unwrapinput.argvalue.deltaLseed; // default 100 pixels;
	  setunspecified(unwrapinput.argvalue.snaphu_log); // " "
	  setunspecified(unwrapinput.argvalue.snaphu_coh); // " "
	  unwrapinput.argvalue.snaphu_mode = "DEFO"; // default to DEFO from TOPO
	  unwrapinput.argvalue.snaphu_init = "MST"; // default method
	  unwrapinput.argvalue.snaphu_verbose = "TRUE"; // default verbose
	  unwrapinput.argvalue.oformatflag = FORMATR4; // default to
															// REAL4 from FORMATHGT

	  // block for Tile.CONTROL: only for SNAPHU; defaults for single CPU
	  unwrapinput.argvalue.ntilerow = 1; // number of tiles in range
	  unwrapinput.argvalue.ntilecol = 1; // number of tiles in azimuth
	  unwrapinput.argvalue.rowovrlp = 0; // overlap between tiles in rng
	  unwrapinput.argvalue.colovrlp = 0; // overlap between tiles in az
	  unwrapinput.argvalue.nproc = 1; // no.cpus or nodes on load
										 // balancing cluster
	  unwrapinput.argvalue.tilecostthresh = 500; // cost threshold boundaries of reliable regions


	  slant2hinput.argvalue.Npoints = 200; // default 100 pixels;
	  slant2hinput.argvalue.degree1d = 2; // default 2
	  slant2hinput.argvalue.degree2d = 5; // default 5
	  slant2hinput.argvalue.Nheights = slant2hinput.argvalue.degree1d+1; // minimum
	  slant2hinput.argvalue.fohei = "hei.raw"; // default
	  slant2hinput.argvalue.fophi = "phi.raw"; // default
	  slant2hinput.argvalue.folam = "lam.raw"; // default

	  geocodeinput.argvalue.fophi = "geo_phi.raw"; // default name
	  geocodeinput.argvalue.folam = "geo_lambda.raw"; // default name


	// ====== (default) Methods ======
	  final int16 def_mte_method = cc_magspace;
	  final int16 def_cc_method = cc_magfft;
	  final int16 def_fc_method = fc_magfft;
	  final int16 def_fe_method = fe_porbits;
	  final int16 def_rs_method = rs_cc4p;
	  final int16 def_uw_method = uw_method2; // snaphu, system call, see www
	  final int16 def_s2h_method = s2h_ambiguity;

	// ______ To check later if default should be used ______
	  mtiminginput.argvalue.method = def_mte_method - 999; // default method (repair later)
	  mtiminginput.argvalue.Nwin = def_mte_nwin + 999; // default #windows
	  coarsecorrinput.argvalue.method = def_cc_method - 999; // default method (repair later)
	  coarsecorrinput.argvalue.Nwin = def_cc_nwin + 999; // default #windows
	  fineinput.argvalue.method = def_fc_method - 999; // default method (repair later)
	  fineinput.argvalue.Nwin = def_fc_nwin + 999; // default #windows
	  unwrapinput.argvalue.method = def_uw_method - 999; // default method (repair later)
	  comprefphainput.argvalue.method = def_fe_method - 999; // default method (repair later)
	  comprefphainput.argvalue.degree = def_fe_degree - 999; // default degree polynomial
	  comprefphainput.argvalue.Npoints = def_fe_Npoints - 999; // default degree polynomial
	  resampleinput.argvalue.method = def_rs_method - 999; // default method (repair later)
	  slant2hinput.argvalue.method = def_s2h_method - 999; // default method (repair later)





	// ====== Process "inputoptionsfile" ======
	  boolean continuereading = true;
	  //while (! optionsfile.eof())                               // read file until end of file.
	  //while (continuereading || ! optionsfile.eof())              // read file until STOP card or if fails exit at end of the file
	  while (continuereading) // read file until STOP card
		{
		linecnt++;
		optionsfile.getline(eachline,4 *ONE27,'\n'); // get line. [MA] if the line is longer than 4*ONE27 then loops forever.
		//cerr << "line: " << linecnt << " : " << eachline << endl;

		if (optionsfile.eof()) // && strcmp(keyword,"STOP") ) // [MA] STOP is missing
		  {
		  ERROR << "STOP:   \t is missing" << " in file '" << inputoptionsfile << "'";
		  PRINT_ERROR(ERROR.get_str())
		  WARNING.print("Please make sure the inputfile line: STOP has ended with [eol] character");
		  throw(keyword_error);
		  continuereading = false; // break while loop
		  }

		final int maxwords = 8; // [GJ, MA]
		String[] word = new String[maxwords];

//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
		char *c = eachline;
		for (int i = 0; i < maxwords; i++)
		{
		  while (*c == ' ' || *c == '\t') // get rid of leadin white space
		  {
			c++; // next address
		  }
		  word = StringHelper.changeCharacter(word, i, c); // pass first char address of a word

		  while (*c != ' ' && *c != '\t' && *c != '\0')
		  {
		   c++;
		  }

		  if (*c != '\0') // at last char
		  {
			*c = '\0';
			c++;
		  }
		}

		String keyword = word.charAt(0); // start off with the card.
									// word[1] --> first argument
									// word[2] --> second argument and so on.
		Character.toUpperCase(keyword);

		DEBUG << linecnt << ": Read keyword: " << keyword;
		DEBUG.print();

	// *******************************************************************
	// *** GENERAL
	// *******************************************************************
		if (!strcmp(keyword,"COMMENT") || !strcmp(keyword,"C"))
		  {
		  ; // comment: no action
		  }
		else if (!strncmp(keyword,"//",2) || !strncmp(keyword,"#",1)) // with // or '#' as delimiter
		  { // but no blank after keyword
		  ; // comment: no action
		  }
		//else if (!strncmp(keyword,'\0',1))          // empty line? crashes
		else if (!keyword.length()) // empty line
		  {
		  ; // comment: no action
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"BEEP")) // level of beep output (THE CARD)
		  { // /progress/[warning]/error/ON/OFF
		  //keyword =  word[1] ;    // pass keyword                 // argument
		  //writearg((char*)keyword);
		  keyword = word.charAt(1); // pass next word (the argument)
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"ERROR"))
			{
			TRACE.bellrings(0);
			DEBUG.bellrings(0);
			INFO.bellrings(0);
			PROGRESS.bellrings(0);
			WARNING.bellrings(0);
			ERROR.bellrings(1);
			beeplevel = -1;
			INFO.print("BEEP: \tbeeping enabled at level: \tERROR");
			}
		  else if (!strcmp(keyword,"PROGRESS"))
			{
			TRACE.bellrings(0);
			DEBUG.bellrings(0);
			INFO.bellrings(0);
			PROGRESS.bellrings(1);
			WARNING.bellrings(2);
			ERROR.bellrings(3);
			beeplevel = 2;
			INFO.print("BEEP: \tbeeping enabled at level: \tPROGRESS");
			}
		  else if (!strcmp(keyword,"WARNING"))
			{
			TRACE.bellrings(0);
			DEBUG.bellrings(0);
			INFO.bellrings(0);
			PROGRESS.bellrings(0);
			WARNING.bellrings(1);
			ERROR.bellrings(2);
			beeplevel = 1;
			INFO.print("BEEP: \tbeeping enabled at level: \tWARNING");
			}
		  else if (!strcmp(keyword,"OFF"))
			{
			TRACE.bellrings(0);
			DEBUG.bellrings(0);
			INFO.bellrings(0);
			PROGRESS.bellrings(0);
			WARNING.bellrings(0);
			ERROR.bellrings(0);
			beeplevel = 0;
			INFO.print("BEEP: \tbeeping disabled");
			}
			 //!strcmp(keyword,\'\\0\'))
		  else if (!strcmp(keyword,"ON") || !strncmp(keyword,"//",2) || !strncmp(keyword,"#",1) || !(keyword.charAt(0) == '\0')) // no keyword
	//             !strcmp(keyword,""))             // no keyword
			{
			TRACE.bellrings(0);
			DEBUG.bellrings(0);
			INFO.bellrings(0);
			PROGRESS.bellrings(1);
			WARNING.bellrings(2);
			ERROR.bellrings(3);
			beeplevel = 2;
			INFO.print("BEEP: \tbeeping enabled for all levels: \tON");
			}
		  else
			{
			beeplevel = 1;
			TRACE.bellrings(0);
			DEBUG.bellrings(0);
			INFO.bellrings(0);
			PROGRESS.bellrings(0);
			WARNING.bellrings(1);
			ERROR.bellrings(2);
			WARNING << "BEEP:   line " << linecnt << ": Argument " << keyword << " not recognized." << " [error/warning/progress/on/off] I used WARNING.";
			WARNING.print();
			}
		  } // BEEP key

	// **********************************************************************
		else if (!strcmp(keyword,"SCREEN")) // level of screen output
		  { // debug/info/progress/warning or error
		  switch (priorscreen)
			{
			case true:
			  WARNING << "SCREEN: line " << linecnt << ": stdout: " << " ignored due to prior occurence.";
			  WARNING.print();
			  break;

			default:
			  priorscreen = true;
			  //keyword =  word[1] ;        // pass keyword                         // argument
			  keyword = word.charAt(1); // argument
			  writearg(keyword);
			  Character.toUpperCase(keyword);
			  if (!strcmp(keyword,"INFO"))
				{
				TRACE.doprint(0);
				DEBUG.doprint(0);
				INFO.doprint(1);
				PROGRESS.doprint(1);
				WARNING.doprint(1);
				ERROR.doprint(1);
				displevel = 20000 + displevel%10000; // for cnt #warnings
				INFO.print("SCREEN: \tverboseness: \t\t\tINFO");
				}
			  else if (!strcmp(keyword,"PROGRESS"))
				{
				TRACE.doprint(0);
				DEBUG.doprint(0);
				INFO.doprint(1);
				PROGRESS.doprint(1);
				WARNING.doprint(1);
				ERROR.doprint(1);
				displevel = 10000 + displevel%10000; // for cnt #warnings
				INFO.print("SCREEN: \tverboseness: \t\t\tPROGRESS");
				}
			  else if (!strcmp(keyword,"DEBUG"))
				{
				TRACE.doprint(0);
				DEBUG.doprint(1);
				INFO.doprint(1);
				PROGRESS.doprint(1);
				WARNING.doprint(1);
				ERROR.doprint(1);
				displevel = 30000 + displevel%10000; // for cnt #warnings
				INFO.print("SCREEN: \tverboseness: \t\t\tDEBUG");
				}
			  else if (!strcmp(keyword,"TRACE"))
				{
				TRACE.doprint(1);
				DEBUG.doprint(1);
				INFO.doprint(1);
				PROGRESS.doprint(1);
				WARNING.doprint(1);
				ERROR.doprint(1);
				}
			  else if (!strcmp(keyword,"WARNING"))
				{
				TRACE.doprint(0);
				DEBUG.doprint(0);
				INFO.doprint(0);
				PROGRESS.doprint(0);
				WARNING.doprint(1);
				ERROR.doprint(1);
				displevel = 0 + displevel%10000; // for cnt #warnings;
				INFO.print("SCREEN: \tverboseness: \t\t\tWARNING");
				}
			  else if (!strcmp(keyword,"ERROR"))
				{
				TRACE.doprint(0);
				DEBUG.doprint(0);
				INFO.doprint(0);
				PROGRESS.doprint(0);
				WARNING.doprint(0);
				ERROR.doprint(1);
				displevel = -100 + displevel%10000; // for cnt #warnings
				INFO.print("SCREEN: \tverboseness: \t\t\tERROR");
				}
			  else
				{
				TRACE.doprint(0);
				DEBUG.doprint(1);
				INFO.doprint(1);
				PROGRESS.doprint(1);
				WARNING.doprint(1);
				ERROR.doprint(1);
				WARNING << "SCREEN: line " << linecnt << ": Argument " << keyword << " not recognized." << " [error/warning/progress/info/debug] I used DEBUG.";
				WARNING.print();
				displevel = 30000 + displevel%10000; // for cnt #warnings
				}
			} // switch
		  } // SCREEN key

	// **********************************************************************
		else if (!strcmp(keyword,"MEMORY")) // available mem in MB
		  {
		  switch (priormemory)
			{
			case true:
			  WARNING << "MEMORY: line " << linecnt << ": ignored due to prior occurence.";
			  WARNING.print();
			  break;

			default:
			  priormemory = true;
			  //generalinput.memory =  word[1] ;    // pass keyword
			  keyword = word.charAt(1);
			  String pLast = null;
			  generalinput.argvalue.memory = uint(strtod(keyword, pLast)); // [MA] try strtoul( keyword, &pLast, 10 ) for uint
			  if (pLast.equals(keyword)) // fail
			  {
				ERROR << "memory argument: " << keyword << " is not valid.";
				PRINT_ERROR(ERROR.get_str())
				throw(keyword_error);
			  }
			  writearg(generalinput.argvalue.memory);
			  if (generalinput.argvalue.memory > 10000)
				WARNING.print("MEMORY: > 10000 MB seems unlikely.");
			  generalinput.argvalue.memory *= 1000000; // in B
			} // switch
		  } // MEMORY card

	// **********************************************************************
		else if (!strcmp(keyword,"BATCH")) // overrides interactive mode
		  {
		  switch (priorbatch)
			{
			case true:
			  WARNING << "BATCH: line: " << linecnt << ": " << "ignored due to prior occurence.";
			  WARNING.print();
			  break;
			default:
			  priorbatch = true; // flag for occurence
			 // keyword =  word[1] ;        // pass keyword                 // argument
			  keyword = word.charAt(1); // pass next word (the argument)
			  writearg(keyword);
		Character.toUpperCase(keyword);
			  if (!strcmp(keyword,"OFF"))
				generalinput.argvalue.interactive = true;
			  else if (!strcmp(keyword,"ON") || !strncmp(keyword,"//",2) || !strncmp(keyword,"#",1) || !(keyword.charAt(0) == '\0')) // no keyword
				generalinput.argvalue.interactive = false;
			  else
				{
				generalinput.argvalue.interactive = true;
				WARNING << "BATCH: line: " << linecnt << ": " << "argument: " << keyword << " not recognized, interactive processing.";
				WARNING.print();
				}
			} // switch
		  } // BATCH key

	// **********************************************************************
		else if (!strcmp(keyword,"OVERWRITE"))
		  {
		  switch (prioroverwrite)
			{
			case true:
			  WARNING << "OVERWRITE: line: " << linecnt << ": " << "ignored due to prior occurence.";
			  WARNING.print();
			  break;
			default:
			  prioroverwrite = true; // flag for occurence
			  keyword = word.charAt(1); // pass next word (the argument)
			  writearg(keyword);
			  Character.toUpperCase(keyword);
			  if (!strcmp(keyword,"OFF"))
				generalinput.argvalue.overwrit = false;
			  else if (!strcmp(keyword,"ON") || !strncmp(keyword,"//",2) || !strncmp(keyword,"#",1) || !(keyword.charAt(0) == '\0')) // no keyword
				generalinput.argvalue.overwrit = true;
			  else
				{
				generalinput.argvalue.overwrit=false; // don't overwrite files
				WARNING << "OVERWRITE: line " << linecnt << ": argument: " << keyword << " not recognized, existing files are not overwritten.";
				WARNING.print();
				}
			} // switch
		  } // OVERWRITE key

	// **********************************************************************
		else if (!strcmp(keyword,"LISTINPUT"))
		  {
		  switch (priorlistinput)
			{
			case true:
			  WARNING << "LISTINPUT: line: " << linecnt << ": " << "ignored due to prior occurence.";
			  WARNING.print();
			  break;
			default:
			  priorlistinput = true; // flag for occurence
			  keyword = word.charAt(1); // pass next word (the argument)
			  writearg(keyword);
			  Character.toUpperCase(keyword);
			  if (!strcmp(keyword,"OFF"))
				listinput = false;
			  else if (!strcmp(keyword,"ON") || !strncmp(keyword,"//",2) || !strncmp(keyword,"#",1) || !(keyword.charAt(0) == '\0')) // no keyword
				listinput = true;
			  else
				{
				listinput = true; // default list input
				WARNING << "LISTINPUT: line " << linecnt << ": argument: " << keyword << " not recognized, input will be appended to logfile.";
				WARNING.print();
				}
			} // switch
		  } // LISTINPUT key

	// **********************************************************************
		else if (!strcmp(keyword,"ONLYPROCESS")) // process only one step
		  {
		  // Read argument and set onlyprocess to value for filling
		  //  the flag input array 'process[NUMPROCESSES]' after reading reset input
		  //  to avoid interference with PROCESS cards (ONLYPROCESS overrides)
		  //
		  if (onlyprocess == -1) // check multiple occurences
			{
			  keyword = word.charAt(1); // pass next word (the argument)
			writearg(keyword);
		Character.toUpperCase(keyword);
			if (!strcmp(keyword,"M_READFILES"))
			  onlyprocess =pr_m_readfiles;
			else if (!strcmp(keyword,"M_CROP"))
			  onlyprocess =pr_m_crop;
	//____RaffaeleNutricato START MODIFICATION SECTION 5
			else if (!strcmp(keyword,"M_OVS"))
			  onlyprocess =pr_m_oversample;
	//____RaffaeleNutricato END MODIFICATION SECTION 5
			else if (!strcmp(keyword,"M_PORBITS"))
			  onlyprocess =pr_m_porbits;
			else if (!strcmp(keyword,"M_SIMAMP")) // [MA]
			  onlyprocess =pr_m_simamp;
			else if (!strcmp(keyword,"M_TIMING"))
			  onlyprocess =pr_m_mtiming;
			else if (!strcmp(keyword,"M_FILTAZI"))
			  onlyprocess =pr_m_filtazi;
			else if (!strcmp(keyword,"FILTRANGE"))
			  {
			  onlyprocess =pr_m_filtrange;
			  onlyprocess =pr_s_filtrange;
			  }
			else if (!strcmp(keyword,"M_EXTRA"))
			  onlyprocess =pr_m_EXTRA;

			else if (!strcmp(keyword,"S_READFILES"))
			  onlyprocess =pr_s_readfiles;
			else if (!strcmp(keyword,"S_CROP"))
			  onlyprocess =pr_s_crop;
	//____RaffaeleNutricato START MODIFICATION SECTION 6
			else if (!strcmp(keyword,"S_OVS"))
			  onlyprocess =pr_s_oversample;
	//____RaffaeleNutricato END MODIFICATION SECTION 6
			else if (!strcmp(keyword,"S_PORBITS"))
			  onlyprocess =pr_s_porbits;
			else if (!strcmp(keyword,"S_FILTAZI"))
			  onlyprocess =pr_s_filtazi;
			else if (!strcmp(keyword,"RESAMPLE"))
			  onlyprocess =pr_s_resample;
			else if (!strcmp(keyword,"S_EXTRA"))
			  onlyprocess =pr_s_EXTRA;

			else if (!strcmp(keyword,"COARSEORB"))
			  onlyprocess =pr_i_coarse;
			else if (!strcmp(keyword,"COARSECORR"))
			  onlyprocess =pr_i_coarse2;
			else if (!strcmp(keyword,"FINE"))
			  onlyprocess =pr_i_fine; // see methodselector
			else if (!strcmp(keyword,"RELTIMING")) // [FvL]
			  onlyprocess =pr_i_timing;
			else if (!strcmp(keyword,"DEMASSIST")) // [FvL]
			  onlyprocess =pr_i_demassist;
			else if (!strcmp(keyword,"COREGPM"))
			  onlyprocess =pr_i_coregpm;
			else if (!strcmp(keyword,"INTERFERO"))
			  onlyprocess =pr_i_interfero;
			else if (!strcmp(keyword,"COHERENCE"))
			  onlyprocess =pr_i_coherence;
			else if (!strcmp(keyword,"FILTPHASE"))
			  onlyprocess =pr_i_filtphase;
			else if (!strcmp(keyword,"COMPREFPHA"))
			  onlyprocess =pr_i_comprefpha;
			else if (!strcmp(keyword,"SUBTRREFPHA"))
			  onlyprocess =pr_i_subtrrefpha;
			else if (!strcmp(keyword,"COMPREFDEM"))
			  onlyprocess =pr_i_comprefdem;
			else if (!strcmp(keyword,"SUBTRREFDEM"))
			  onlyprocess =pr_i_subtrrefdem;
			else if (!strcmp(keyword,"UNWRAP"))
			  onlyprocess =pr_i_unwrap; // see methodselector
			else if (!strcmp(keyword,"SLANT2H"))
			  onlyprocess =pr_i_slant2h;
			else if (!strcmp(keyword,"GEOCODE"))
			  onlyprocess =pr_i_geocoding;
			else if (!strcmp(keyword,"DINSAR"))
			  onlyprocess =pr_i_dinsar;
			else if (!strcmp(keyword,"I_EXTRA2"))
			  onlyprocess =pr_i_EXTRA2;
			else
			  {
			  ERROR << "ONLYPROCESS: line " << linecnt << ": Argument " << keyword << " not recognized.";
			  PRINT_ERROR(ERROR.get_str())
			  throw(keyword_error);
			  }
			INFO << "ONLYPROCESS: \tonly processing step: \t\t" << keyword;
			INFO.print();
			}
		  else
			{
			WARNING << "ONLYPROCESS: more than one occurence of card, ignored line: " << linecnt << ".";
			WARNING.print();

			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"PROCESS")) // which routine to run
		  {
		  if (onlyprocess+1) // initialized to -1;
			{
			WARNING << "PROCESS card on line " << linecnt << " ignored due to presence of ONLYPROCESS card.";
			WARNING.print();
			}
		  else
			{
			  keyword = word.charAt(1); // pass next word (the argument)
			writearg(keyword);
		Character.toUpperCase(keyword);
			if (!strcmp(keyword,"M_READFILES"))
			  generalinput.argvalue.process[pr_m_readfiles] = 1;
			else if (!strcmp(keyword,"M_CROP"))
			  generalinput.argvalue.process[pr_m_crop] = 1;
	//____RaffaeleNutricato START MODIFICATION SECTION 7
			else if (!strcmp(keyword,"M_OVS"))
			  generalinput.argvalue.process[pr_m_oversample] = 1;
	//____RaffaeleNutricato END MODIFICATION SECTION 7
			else if (!strcmp(keyword,"M_PORBITS"))
			  generalinput.argvalue.process[pr_m_porbits] = 1;
			else if (!strcmp(keyword,"M_SIMAMP")) // [MA]
			  generalinput.argvalue.process[pr_m_simamp] = 1;
			else if (!strcmp(keyword,"M_TIMING"))
			  generalinput.argvalue.process[pr_m_mtiming] = 1;
			else if (!strcmp(keyword,"M_FILTAZI"))
			  generalinput.argvalue.process[pr_m_filtazi] = 1;
			else if (!strcmp(keyword,"FILTRANGE"))
			  {
			  generalinput.argvalue.process[pr_m_filtrange] = 1; // use for s_ as well
			  generalinput.argvalue.process[pr_s_filtrange] = 1;
			  }
			else if (!strcmp(keyword,"M_EXTRA"))
			  generalinput.argvalue.process[pr_m_EXTRA] = 1;

			else if (!strcmp(keyword,"S_READFILES"))
			  generalinput.argvalue.process[pr_s_readfiles] = 1;
			else if (!strcmp(keyword,"S_CROP"))
			  generalinput.argvalue.process[pr_s_crop] = 1;
	//____RaffaeleNutricato START MODIFICATION SECTION 8
			else if (!strcmp(keyword,"S_OVS"))
			  generalinput.argvalue.process[pr_s_oversample] = 1;
	//____RaffaeleNutricato END MODIFICATION SECTION 8
			else if (!strcmp(keyword,"S_PORBITS"))
			  generalinput.argvalue.process[pr_s_porbits] = 1;
			else if (!strcmp(keyword,"S_FILTAZI"))
			  generalinput.argvalue.process[pr_s_filtazi] = 1;
			else if (!strcmp(keyword,"RESAMPLE"))
			  generalinput.argvalue.process[pr_s_resample] = 1;
			else if (!strcmp(keyword,"S_EXTRA"))
			  generalinput.argvalue.process[pr_s_EXTRA] = 1;

			else if (!strcmp(keyword,"COARSEORB"))
			  generalinput.argvalue.process[pr_i_coarse] = 1;
			else if (!strcmp(keyword,"COARSECORR"))
			  generalinput.argvalue.process[pr_i_coarse2] = 1;
			else if (!strcmp(keyword,"FINE"))
			  generalinput.argvalue.process[pr_i_fine] = 1;
			else if (!strcmp(keyword,"RELTIMING")) // [FvL]
			  generalinput.argvalue.process[pr_i_timing] = 1;
			else if (!strcmp(keyword,"DEMASSIST")) // [FvL]
			  generalinput.argvalue.process[pr_i_demassist] = 1;
			else if (!strcmp(keyword,"COREGPM"))
			  generalinput.argvalue.process[pr_i_coregpm] = 1;
			else if (!strcmp(keyword,"COMPREFPHA"))
			  generalinput.argvalue.process[pr_i_comprefpha] = 1;
			else if (!strcmp(keyword,"SUBTRREFPHA"))
			  generalinput.argvalue.process[pr_i_subtrrefpha] = 1;
			else if (!strcmp(keyword,"COMPREFDEM"))
			  generalinput.argvalue.process[pr_i_comprefdem] = 1;
			else if (!strcmp(keyword,"SUBTRREFDEM"))
			  generalinput.argvalue.process[pr_i_subtrrefdem] = 1;
			else if (!strcmp(keyword,"INTERFERO"))
			  generalinput.argvalue.process[pr_i_interfero] = 1;
			else if (!strcmp(keyword,"COHERENCE"))
			  generalinput.argvalue.process[pr_i_coherence] = 1;
			else if (!strcmp(keyword,"FILTPHASE"))
			  generalinput.argvalue.process[pr_i_filtphase] = 1;
			else if (!strcmp(keyword,"UNWRAP"))
			  generalinput.argvalue.process[pr_i_unwrap] = 1;
			else if (!strcmp(keyword,"SLANT2H"))
			  generalinput.argvalue.process[pr_i_slant2h] = 1;
			else if (!strcmp(keyword,"GEOCODE"))
			  generalinput.argvalue.process[pr_i_geocoding] = 1;
			else if (!strcmp(keyword,"DINSAR"))
			  generalinput.argvalue.process[pr_i_dinsar] = 1;

			else if (!strcmp(keyword,"I_EXTRA2"))
			  generalinput.argvalue.process[pr_i_EXTRA2] = 1;
			else
			  {
			  ERROR << "PROCESS: line " << linecnt << ": Argument " << keyword << " not recognized.";
			  PRINT_ERROR(ERROR.get_str())
			  throw(keyword_error);
			  }
			INFO << "PROCESS: \tI will process step: \t\t" << keyword;
			INFO.print();
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"ELLIPSOID")) // ref. system
		  { // inputoptionsfile
		  ellipsoid = true; // use below
		  keyword = word.charAt(1);
		  String keyword2 = word.charAt(2);
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  writearg(keyword2);
		  if (!strcmp(keyword,"WGS84"))
			{
			ellipsinput.argvalue = WGS84; // default
			}
		  else if (!strcmp(keyword,"GRS80"))
			{
			WARNING.print("ELLIPS: not ok, sat. ephemerides should be in this system.");
				  input_ell GRS80 = new input_ell(6378137.0,6356752.3);
				  GRS80.set_name("GRS80");
				  ellipsinput.argvalue = GRS80; // copy
			}
		  else if (!strcmp(keyword,"BESSEL"))
			{
			WARNING.print("ELLIPS: not ok, sat. ephemerides should be in this system.");
				  input_ell BESSEL = new input_ell(6377397.155,6356078.963);
				  BESSEL.set_name("BESSEL");
				  ellipsinput.argvalue = BESSEL; // copy
			}
		  else if (Character.isDigit(keyword2.charAt(0))) // likely to be a,b
			{
			WARNING.print("ELLIPS: not ok, sat. ephemerides should be in this system.");
			input_ell ELL_USER_DEFINED = new input_ell(Double.parseDouble(keyword),Double.parseDouble(keyword2)); // a,b
				  ELL_USER_DEFINED.set_name("user_defined");
			ellipsinput.argvalue = ELL_USER_DEFINED; // copy
			if (ellipsinput.argvalue.a<ellipsinput.argvalue.b || ellipsinput.argvalue.b<EPS)
			  {
			  ERROR << "ELLIPSOID keyword (real8A real8B): B==0 or A<B: " << ellipsinput.argvalue.a << "<" << ellipsinput.argvalue.b << " at line " << linecnt << ".";
					PRINT_ERROR(ERROR.get_str())
					throw(keyword_error);
			  }
			}
		  else
			{
			PRINT_ERROR("unknown argument for ellipsoid card.")
			throw(keyword_error);
			}
		  INFO << "ELLIPSOID: \tsemimajor=" << ellipsinput.argvalue.a << ", semiminor=" << ellipsinput.argvalue.b << "; line " << linecnt << ".";
		  INFO.print();
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"M_RESFILE")) // resultfile filename
		  {
		  //generalinput.m_resfile =  word[1] ;     // pass keyword
		  keyword = word.charAt(1);
		  // strcpy(generalinput.m_resfile, keyword);
		  generalinput.argvalue.m_resfile = keyword;
		  writearg(generalinput.argvalue.m_resfile);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S_RESFILE")) // resultfile filename
		  {
		  keyword = word.charAt(1);
		  generalinput.argvalue.s_resfile = keyword;
		  writearg(generalinput.argvalue.s_resfile);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"LOGFILE")) // logfile filename
		  {
		  keyword = word.charAt(1);
		  generalinput.argvalue.logfile = keyword;
		  writearg(generalinput.argvalue.logfile);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"I_RESFILE")) // interferogram.out
		  {
		  keyword = word.charAt(1);
		  generalinput.argvalue.i_resfile = keyword;
		  writearg(generalinput.argvalue.i_resfile);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"ORB_INTERP")) // orbit
		  {
		  //keyword =  word[1] ;    // pass keyword                       // argument
		  keyword = word.charAt(1);
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"POLYFIT"))
			{
			INFO.print("ORB_INTERP:  polynomial fit for interpolation");
			generalinput.argvalue.orb_interp = ORB_DEFAULT; // depends on number of points
			// ___ Check second optional argument with degree ___
			// keyword =  word[1] ;  // pass keyword
			String keyword2 = word.charAt(2);
			int32 degree = Integer.parseInt(keyword2);
			if (degree > 0) // atoi returns 0 if not convertible
			  {
			  generalinput.argvalue.orb_interp = degree;
			  INFO << "ORB_INTERP:  second argument read: degree = " << degree;
			  INFO.print();
			  }
			}
		  else if (!strcmp(keyword,"SPLINE"))
			{
			INFO.print("ORB_INTERP:  natural cubic splines used fit interpolation");
			generalinput.argvalue.orb_interp = ORB_SPLINE; // natural cubic splines
			}
		  else
			{
			WARNING.print("argument ORB_INTERP not reconized, using polyfit");
			generalinput.argvalue.orb_interp = ORB_DEFAULT; // depends on number of points
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"DUMPBASELINE")) // eval on grid
		  {
	//      generalinput.dumpbaselineL >> generalinput.dumpbaselineP =  word[1] ;   // pass keyword
			  keyword = word.charAt(1);
			  String keyword2 = word.charAt(2); // in local scope
			  generalinput.argvalue.dumpbaselineL = Integer.parseInt(keyword);
			  generalinput.argvalue.dumpbaselineP = Integer.parseInt(keyword2);
			  writearg(generalinput.argvalue.dumpbaselineL);
			  writearg(generalinput.argvalue.dumpbaselineP);
		  if (generalinput.argvalue.dumpbaselineL==0 || generalinput.argvalue.dumpbaselineP==0)
			{
			ERROR << "DUMPBASELINE: " << generalinput.argvalue.dumpbaselineL << " " << generalinput.argvalue.dumpbaselineP << " line: " << linecnt << ": ==0.\n";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"PREVIEW")) // system call to cpxfiddle to get SUNraster
		  {
		  //keyword =  word[1] ;    // pass keyword                 // argument
		  keyword = word.charAt(1); // argument
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"OFF"))
			{
			generalinput.argvalue.preview = 0;
			INFO.print("PREVIEW: \tOFF: generation of SUNraster files disabled.");
			}
		  else if (!strcmp(keyword,"XV"))
			{
			generalinput.argvalue.preview = 2;
			INFO.print("PREVIEW: \tON: generation of SUNraster files enabled + XV sytem call.");
			}
		  else if (!strcmp(keyword,"ON") || !strncmp(keyword,"//",2) || !strncmp(keyword,"#",1) || !(keyword.charAt(0) == '\0')) // no keyword
			{
			generalinput.argvalue.preview = 1;
			INFO.print("PREVIEW: \tON: generation of SUNraster files enabled.");
			}
		  else
			{
			generalinput.argvalue.preview = 0;
			WARNING << "PREVIEW: line: " << linecnt << ": " << "argument: " << keyword << " not recognized, no preview generated.";
			WARNING.print();
			}
		  } // PREVIEW key

	// **********************************************************************
		else if (!strcmp(keyword,"HEIGHT")) // mean height or for CROP
		  {
		  //generalinput.terrain_height =  word[1] ;        // pass keyword
			  keyword = word.charAt(1);
			  String pLast = null;
			  generalinput.argvalue.terrain_height = strtod(keyword, pLast);
			  if (pLast.equals(keyword)) // fail
			  {
				ERROR << "Height argument: " << keyword << " is not valid.";
				PRINT_ERROR(ERROR.get_str())
				throw(keyword_error);
			  }
		  writearg(generalinput.argvalue.terrain_height);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"TIEPOINT")) // lat/lon/hei dec.degrees
		  {

			keyword = word.charAt(1);
			String pLast;
			String pLast2;
			String pLast3 = null;
			generalinput.argvalue.tiepoint.x = strtod(keyword, pLast);
			generalinput.argvalue.tiepoint.y = strtod(word.charAt(2), pLast2); // 2nd arg
						generalinput.argvalue.tiepoint.z = strtod(word.charAt(3), pLast3); // 3rd arg
			if (pLast.equals(keyword) || pLast2.equals(word.charAt(2)) || pLast3.equals(word.charAt(3))) // fail
			{
			  ERROR << "Tiepoints: " << keyword << " " << word.charAt(2) << " " << word.charAt(3) << " are not valid.";
			  PRINT_ERROR(ERROR.get_str())
			  throw(keyword_error);
			}
		   writearg(generalinput.argvalue.tiepoint.x);
		   writearg(generalinput.argvalue.tiepoint.y);
		   writearg(generalinput.argvalue.tiepoint.z);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"STOP")) // STOP interpreting input
		  {
		  INFO.print("STOP:   \tEncountered.");
		  DEBUG << "STOP card encountered at line " << linecnt << "\n";
		  DEBUG.print();
		  continuereading = false; // break while loop
		  }

	// *******************************************************************
	// *** ?_READFILES
	// *******************************************************************
		else if (!strcmp(keyword,"M_IN_METHOD")) // ERS or ASAR ENVISAT
		  {
		  //filename =  word[1] ;   // pass keyword
		  keyword = word.charAt(1);
		  writearg(keyword);
		  Character.toUpperCase(keyword);

		  if (!strcmp(keyword,"ERS"))
			m_readfilesinput.argvalue.sensor_id=SLC_ERS; // default ers
		  else if (!strcmp(keyword,"ERS-1")) // ers
			m_readfilesinput.argvalue.sensor_id=SLC_ERS;
		  else if (!strcmp(keyword,"ERS1")) // ers
			m_readfilesinput.argvalue.sensor_id=SLC_ERS;
		  else if (!strcmp(keyword,"ERS-2")) // ers
			m_readfilesinput.argvalue.sensor_id=SLC_ERS;
		  else if (!strcmp(keyword,"ERS2")) // ers
			m_readfilesinput.argvalue.sensor_id=SLC_ERS;
		  else if (!strcmp(keyword,"N1")) // envisat
			m_readfilesinput.argvalue.sensor_id=SLC_ASAR;
		  else if (!strcmp(keyword,"ASAR")) // envisat
			m_readfilesinput.argvalue.sensor_id=SLC_ASAR;
		  else if (!strcmp(keyword,"ENVISAT")) // envisat
			m_readfilesinput.argvalue.sensor_id=SLC_ASAR;
		  else if (!strcmp(keyword,"ATLANTIS")) // radarsat ceos reader
			m_readfilesinput.argvalue.sensor_id=SLC_RSAT;
		  else if (!strcmp(keyword,"RSAT")) // radarsat
			m_readfilesinput.argvalue.sensor_id=SLC_RSAT;
		  else if (!strcmp(keyword,"RADARSAT")) // radarsat
			m_readfilesinput.argvalue.sensor_id=SLC_RSAT;
		  else if (!strcmp(keyword,"RADARSAT-1")) // radarsat
			m_readfilesinput.argvalue.sensor_id=SLC_RSAT;
		  else if (!strcmp(keyword,"RADARSAT-2")) // radarsat
			m_readfilesinput.argvalue.sensor_id=SLC_RSAT;
		  else if (!strcmp(keyword,"JERS")) // jers
			m_readfilesinput.argvalue.sensor_id=SLC_JERS;
		  else if (!strcmp(keyword,"J-ERS")) // jers
			m_readfilesinput.argvalue.sensor_id=SLC_JERS;
		  else if (!strcmp(keyword,"ALOS")) // alos FBS [DON]
			m_readfilesinput.argvalue.sensor_id=SLC_ALOS;
		  else if (!strcmp(keyword,"TSX")) // TSX [PM]
			m_readfilesinput.argvalue.sensor_id=SLC_TSX;
		  else if (!strcmp(keyword,"TERRASARX")) // TSX
			m_readfilesinput.argvalue.sensor_id=SLC_TSX;
		  else if (!strcmp(keyword,"TERRASAR-X")) // TSX
			m_readfilesinput.argvalue.sensor_id=SLC_TSX;
		  else
			{
			ERROR << "M_IN_METHOD: method " << keyword << " not known for reading input files on line " << linecnt << ".";
				PRINT_ERROR(ERROR.get_str())
				throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"M_IN_VOL")) // volumefile filename
		  {
		  //m_readfilesinput.volfile =  word[1] ;   // pass keyword
		  m_readfilesinput.argvalue.volfile = word.charAt(1); // pass keyword
		  writearg(m_readfilesinput.argvalue.volfile);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"M_IN_LEA")) // leaderfile filename
		  {
		  m_readfilesinput.argvalue.leaderfile = word.charAt(1); // pass keyword
		  writearg(m_readfilesinput.argvalue.leaderfile);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"M_IN_NULL")) // nullfile filename
		  {
		  m_readfilesinput.argvalue.nullfile = word.charAt(1); // pass keyword
		  writearg(m_readfilesinput.argvalue.nullfile);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"M_IN_DAT")) // datafile filename
		  {
		  m_readfilesinput.argvalue.datfile = word.charAt(1); // pass keyword
		  writearg(m_readfilesinput.argvalue.datfile);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"M_RG_T_ERROR")) // master timing error
		  {
		   keyword = word.charAt(1);
		   String pLast = null;
		   m_readfilesinput.argvalue.rg_timing_error= strtod(keyword, pLast);
		   if (pLast.equals(keyword)) // fail
			{
			  ERROR << "M_RG_T_ERROR argument: " << keyword << " is not valid.";
			  PRINT_ERROR(ERROR.get_str())
			  throw(keyword_error);
			}
		  writearg(m_readfilesinput.argvalue.rg_timing_error);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"M_AZ_T_ERROR")) // master timing error
		  {
		   keyword = word.charAt(1);
		   String pLast = null;
		   m_readfilesinput.argvalue.az_timing_error = strtod(keyword, pLast);
		   if (pLast.equals(keyword)) // fail
			{
			  ERROR << "M_AZ_T_ERROR argument: " << keyword << " is not valid.";
			  PRINT_ERROR(ERROR.get_str())
			  throw(keyword_error);
			}
		  writearg(m_readfilesinput.argvalue.az_timing_error);
		  }

	// *******************************************************************
		else if (!strcmp(keyword,"S_IN_METHOD")) // ERS or ASAR ENVISAT
		  {
		  keyword = word.charAt(1);
		  writearg(keyword);
		  Character.toUpperCase(keyword);

		  if (!strcmp(keyword,"ERS"))
			s_readfilesinput.argvalue.sensor_id=SLC_ERS; // default
		  else if (!strcmp(keyword,"ERS-1")) // ers
			s_readfilesinput.argvalue.sensor_id=SLC_ERS;
		  else if (!strcmp(keyword,"ERS1")) // ers
			s_readfilesinput.argvalue.sensor_id=SLC_ERS;
		  else if (!strcmp(keyword,"ERS-2")) // ers
			s_readfilesinput.argvalue.sensor_id=SLC_ERS;
		  else if (!strcmp(keyword,"ERS2")) // ers
			s_readfilesinput.argvalue.sensor_id=SLC_ERS;
		  else if (!strcmp(keyword,"N1")) // envisat
			s_readfilesinput.argvalue.sensor_id=SLC_ASAR;
		  else if (!strcmp(keyword,"ASAR")) // envisat
			s_readfilesinput.argvalue.sensor_id=SLC_ASAR;
		  else if (!strcmp(keyword,"ENVISAT")) // envisat
			s_readfilesinput.argvalue.sensor_id=SLC_ASAR;
		  else if (!strcmp(keyword,"ATLANTIS")) // radarsat ceos reader
			s_readfilesinput.argvalue.sensor_id=SLC_RSAT;
		  else if (!strcmp(keyword,"RSAT")) // radarsat
			s_readfilesinput.argvalue.sensor_id=SLC_RSAT;
		  else if (!strcmp(keyword,"RADARSAT")) // radarsat
			s_readfilesinput.argvalue.sensor_id=SLC_RSAT;
		  else if (!strcmp(keyword,"RADARSAT-1")) // radarsat
			s_readfilesinput.argvalue.sensor_id=SLC_RSAT;
		  else if (!strcmp(keyword,"RADARSAT-2")) // radarsat
			s_readfilesinput.argvalue.sensor_id=SLC_RSAT;
		  else if (!strcmp(keyword,"JERS")) // jers
			s_readfilesinput.argvalue.sensor_id=SLC_JERS;
		  else if (!strcmp(keyword,"J-ERS")) // jers
			s_readfilesinput.argvalue.sensor_id=SLC_JERS;
		  else if (!strcmp(keyword,"ALOS")) // [DON]
			s_readfilesinput.argvalue.sensor_id=SLC_ALOS;
		  else if (!strcmp(keyword,"TSX")) // TSX [PM]
			s_readfilesinput.argvalue.sensor_id=SLC_TSX;
		  else if (!strcmp(keyword,"TERRASARX")) // TSX
			s_readfilesinput.argvalue.sensor_id=SLC_TSX;
		  else if (!strcmp(keyword,"TERRASAR-X")) // TSX
			s_readfilesinput.argvalue.sensor_id=SLC_TSX;
		  else
			{
			ERROR << "S_IN_METHOD: method " << keyword << " not known for reading input files on line " << linecnt << ".";
				PRINT_ERROR(ERROR.get_str())
				throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S_IN_VOL")) // volumefile filename
		  {
		  s_readfilesinput.argvalue.volfile = word.charAt(1); // pass keyword
		  writearg(s_readfilesinput.argvalue.volfile);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S_IN_LEA")) // leaderfile filename
		  {
		  s_readfilesinput.argvalue.leaderfile = word.charAt(1); // pass keyword
		  writearg(s_readfilesinput.argvalue.leaderfile);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S_IN_DAT")) // nullfile filename
		  {
		  s_readfilesinput.argvalue.datfile = word.charAt(1); // pass keyword
		  writearg(s_readfilesinput.argvalue.datfile);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S_IN_NULL")) // nullfile filename
		  {
		  s_readfilesinput.argvalue.nullfile = word.charAt(1); // pass keyword
		  writearg(s_readfilesinput.argvalue.nullfile);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S_RG_T_ERROR")) // slave timing error
		  {
		   keyword = word.charAt(1);
		   String pLast = null;
		   s_readfilesinput.argvalue.rg_timing_error = strtod(keyword, pLast);
		   if (pLast.equals(keyword)) // fail
			{
			  ERROR << "S_RG_T_ERROR argument: " << keyword << " is not valid.";
			  PRINT_ERROR(ERROR.get_str())
			  throw(keyword_error);
			}
		  writearg(s_readfilesinput.argvalue.rg_timing_error);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S_AZ_T_ERROR")) // slave timing error
		  {
		   keyword = word.charAt(1);
		   String pLast = null;
		   s_readfilesinput.argvalue.az_timing_error = strtod(keyword, pLast);
		   if (pLast.equals(keyword)) // fail
			{
			  ERROR << "S_AZ_T_ERROR argument: " << keyword << " is not valid.";
			  PRINT_ERROR(ERROR.get_str())
			  throw(keyword_error);
			}
		  writearg(s_readfilesinput.argvalue.az_timing_error);
		  }


	// **********************************************************************
	// *** ?_PORBITS
	// **********************************************************************
		else if (!strcmp(keyword,"M_ORBDIR")) // orbitfile filename
		  {
		  if (specified(porbitsinput.argvalue.m_orbdir))
				   WARNING.print("Prior occurence of M_ORBDIR ignored.");
		  porbitsinput.argvalue.m_orbdir = word.charAt(1); // pass keyword
		  writearg(porbitsinput.argvalue.m_orbdir);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"M_ORB_INTERVAL") || !strcmp(keyword,"S_ORB_INTERVAL"))
		  {
		   keyword = word.charAt(1); // pass keyword
		   String pLast = null;
		   porbitsinput.argvalue.timeinterval = strtol(keyword, pLast, BASE10); // int32
		   if (pLast.equals(keyword)) // fail
			{
			  ERROR << "[M|S]_ORB_INTERVAL argument: " << keyword << " is not valid.";
			  PRINT_ERROR(ERROR.get_str())
			  throw(keyword_error);
			}
		  writearg(porbitsinput.argvalue.timeinterval);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"M_ORB_EXTRATIME") || !strcmp(keyword,"S_ORB_EXTRATIME"))
		  {
		   keyword = word.charAt(1); // pass keyword
		   String pLast = null;
		   porbitsinput.argvalue.timebefore = strtol(keyword, pLast, BASE10); // int32
		   if (pLast.equals(keyword)) // fail
			{
			  ERROR << "[M|S]_ORB_EXTRATIME argument: " << keyword << " is not valid.";
			  PRINT_ERROR(ERROR.get_str())
			  throw(keyword_error);
			}
		  writearg(porbitsinput.argvalue.timebefore);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"M_ORB_DUMP")) // dump orbit with dt
		  {
		   keyword = word.charAt(1); // pass keyword
		   String pLast = null;
		   porbitsinput.argvalue.dumpmasterorbit = strtod(keyword, pLast);
		   if (pLast.equals(keyword)) // fail
			{
			  ERROR << "M_ORB_DUMP argument: " << keyword << " is not valid.";
			  PRINT_ERROR(ERROR.get_str())
			  throw(keyword_error);
			}
		  writearg(porbitsinput.argvalue.dumpmasterorbit);
		  INFO.print("dumping master orbit to ascii file: masterorbit.dat");
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S_ORBDIR")) // orbitfile filename
		  {
		  if (specified(porbitsinput.argvalue.s_orbdir))
				  WARNING.print("Prior occurence of S_ORBDIR ignored.");
		  porbitsinput.argvalue.s_orbdir = word.charAt(1); // pass keyword
		  writearg(porbitsinput.argvalue.s_orbdir);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S_ORB_DUMP")) // dump orbit with dt
		  {
		   keyword = word.charAt(1); // pass keyword
		   String pLast = null;
		   porbitsinput.argvalue.dumpslaveorbit = strtod(keyword, pLast);
		   if (pLast.equals(keyword)) // fail
			{
			  ERROR << "S_ORB_DUMP argument: " << keyword << " is not valid.";
			  PRINT_ERROR(ERROR.get_str())
			  throw(keyword_error);
			}
		  writearg(porbitsinput.argvalue.dumpslaveorbit);
		  INFO.print("dumping slave orbit to ascii file: slaveorbit.dat");
		  }


	// *******************************************************************
	// *** ?_CROP
	// *******************************************************************
		else if (!strcmp(keyword,"M_CROP_ID")) // identifier of run
		  {
		  m_cropinput.argvalue.idcrop = word.charAt(1); // pass keyword
		  writearg(m_cropinput.argvalue.idcrop);
		  }

	// *******************************************************************   CROP
		else if (!strcmp(keyword,"S_CROP_ID")) // identifier of run
		  {
		  s_cropinput.argvalue.idcrop = word.charAt(1); // pass keyword
		  writearg(s_cropinput.argvalue.idcrop);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"M_CROP_IN")) // SLC input image filename
		  {
		  m_cropinput.argvalue.filein1 = word.charAt(1); // pass keyword
		  writearg(m_cropinput.argvalue.filein1);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S_CROP_IN")) // SLC input image filename
		  {
		  s_cropinput.argvalue.filein1 = word.charAt(1); // pass keyword
		  writearg(s_cropinput.argvalue.filein1);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"M_CROP_OUT")) // SLC output image filename
		  {
		  m_cropinput.argvalue.fileout1 = word.charAt(1); // pass keyword
		  writearg(m_cropinput.argvalue.fileout1);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S_CROP_OUT")) // SLC output image filename
		  {
		  s_cropinput.argvalue.fileout1 = word.charAt(1); // pass keyword
		  writearg(s_cropinput.argvalue.fileout1);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"M_DBOW")) // min. line coord.
		  {
		   m_cropinput.argvalue.dbow.linelo = Integer.parseInt(word.charAt(1)); // pass keywords
		   m_cropinput.argvalue.dbow.linehi = Integer.parseInt(word.charAt(2));
		   m_cropinput.argvalue.dbow.pixlo = Integer.parseInt(word.charAt(3));
		   m_cropinput.argvalue.dbow.pixhi = Integer.parseInt(word.charAt(4));

		  writearg(m_cropinput.argvalue.dbow.linelo);
		  writearg(m_cropinput.argvalue.dbow.linehi);
		  writearg(m_cropinput.argvalue.dbow.pixlo);
		  writearg(m_cropinput.argvalue.dbow.pixhi);

	// ______initial check, later checked to #lines of image______
		  if (m_cropinput.argvalue.dbow.linelo <= 0 || m_cropinput.argvalue.dbow.pixlo <= 0 || m_cropinput.argvalue.dbow.linelo > m_cropinput.argvalue.dbow.linehi || m_cropinput.argvalue.dbow.pixlo > m_cropinput.argvalue.dbow.pixhi)
			{
			ERROR << "code 300: Arguments of M_DBOW card on line " << linecnt << " missing. [DBOW  min_line  max_line  min_pixel  max_pixel].";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S_DBOW")) // min. line coord.
		  {
		   s_cropinput.argvalue.dbow.linelo = Integer.parseInt(word.charAt(1)); // pass keywords
		   s_cropinput.argvalue.dbow.linehi = Integer.parseInt(word.charAt(2));
		   s_cropinput.argvalue.dbow.pixlo = Integer.parseInt(word.charAt(3));
		   s_cropinput.argvalue.dbow.pixhi = Integer.parseInt(word.charAt(4));

		  writearg(s_cropinput.argvalue.dbow.linelo);
		  writearg(s_cropinput.argvalue.dbow.linehi);
		  writearg(s_cropinput.argvalue.dbow.pixlo);
		  writearg(s_cropinput.argvalue.dbow.pixhi);
		  // ______initial check, later checked to #lines of image______
		  if (s_cropinput.argvalue.dbow.linelo <= 0 || s_cropinput.argvalue.dbow.pixlo <= 0 || s_cropinput.argvalue.dbow.linelo > s_cropinput.argvalue.dbow.linehi || s_cropinput.argvalue.dbow.pixlo > s_cropinput.argvalue.dbow.pixhi)
			{
			ERROR << "code 300: Arguments of S_DBOW card on line " << linecnt << " missing. [DBOW  min_line  max_line  min_pixel  max_pixel].";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		// BK 15-Dec-2003
		else if (!strcmp(keyword,"M_DBOW_GEO")) // 1e6*lat_0 [deg], lon_0, height, width [pix]
		  {
		  // use dbow to store tmp, compute later in crop if not zero
		  real8 tmp_lat_0;
		  real8 tmp_lon_0;
		  real8 tmp_height;
		  real8 tmp_width;
		  tmp_lat_0 = Double.parseDouble(word.charAt(1)); // pass keyword
		  tmp_lon_0 = Double.parseDouble(word.charAt(2)); // pass keyword
		  tmp_height = Double.parseDouble(word.charAt(3)); // pass keyword
		   tmp_width = Double.parseDouble(word.charAt(4)); // pass keyword
		  writearg(tmp_lat_0);
		  writearg(tmp_lon_0);
		  writearg(tmp_height);
		  writearg(tmp_width);
		  m_cropinput.argvalue.dbow_geo.linelo = uint((360.0+tmp_lat_0)*1e6);
		  m_cropinput.argvalue.dbow_geo.linehi = uint((360.0+tmp_lon_0)*1e6);
		  m_cropinput.argvalue.dbow_geo.pixlo = uint(tmp_height);
		  m_cropinput.argvalue.dbow_geo.pixhi = uint(tmp_width);
		  }

	// **********************************************************************
		// BK 15-Dec-2003
		else if (!strcmp(keyword,"S_DBOW_GEO")) // 1e6*lat_0 [deg], lon_0, height, width [pix]
		  {
		  // use dbow to store tmp, compute later in crop if not zero
		  real8 tmp_lat_0;
		  real8 tmp_lon_0;
		  real8 tmp_height;
		  real8 tmp_width;
		  tmp_lat_0 = Double.parseDouble(word.charAt(1)); // pass keyword
		  tmp_lon_0 = Double.parseDouble(word.charAt(2)); // pass keyword
		  tmp_height = Double.parseDouble(word.charAt(3)); // pass keyword
		   tmp_width = Double.parseDouble(word.charAt(4)); // pass keyword
		  writearg(tmp_lat_0);
		  writearg(tmp_lon_0);
		  writearg(tmp_height);
		  writearg(tmp_width);
		  s_cropinput.argvalue.dbow_geo.linelo = uint((360.0+tmp_lat_0)*1e6);
		  s_cropinput.argvalue.dbow_geo.linehi = uint((360.0+tmp_lon_0)*1e6);
		  s_cropinput.argvalue.dbow_geo.pixlo = uint(tmp_height);
		  s_cropinput.argvalue.dbow_geo.pixhi = uint(tmp_width);
		  }

	//____RaffaeleNutricato START MODIFICATION SECTION 9
	// *******************************************************************
	// *** ?_OVS
	// **********************************************************************
		else if (!strcmp(keyword,"M_OVS_OUT")) // oversampled SLC output image filename
		  {
		  m_oversample.argvalue.fileoutovs = word.charAt(1); // pass keyword
		  writearg(m_oversample.argvalue.fileoutovs);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S_OVS_OUT")) // oversampled SLC output image filename
		  {
		  s_oversample.argvalue.fileoutovs = word.charAt(1); // pass keyword
		  writearg(s_oversample.argvalue.fileoutovs);
		  }

	// *******************************************************************
		else if (!strcmp(keyword,"M_OVS_FACT_RNG")) // Oversampling ratio in the master range direction.
		  {
		  m_oversample.argvalue.OsrRange = strtol(word.charAt(1),null, BASE10); // pass keyword
		  writearg(m_oversample.argvalue.OsrRange);
		  if (m_oversample.argvalue.OsrRange < 1)
				  {
			PRINT_ERROR("(M_OVS_FACT_RNG < 1) Oversampling ratio must be at least 1.");
				  throw(keyword_error);
				  }
		  }

	// ********************************************************************** 
		else if (!strcmp(keyword,"S_OVS_FACT_RNG")) // Oversampling ratio in the slave range direction.
		  {
		  s_oversample.argvalue.OsrRange = strtol(word.charAt(1),null, BASE10); // pass keyword
		  writearg(s_oversample.argvalue.OsrRange);
		  if (s_oversample.argvalue.OsrRange < 1)
			{
			PRINT_ERROR("(S_OVS_FACT_RNG < 1) Oversampling ratio must be at least 1.");
			throw(keyword_error);
			}
		  }

	// *******************************************************************
		else if (!strcmp(keyword,"M_OVS_FACT_AZI")) // Oversampling ratio in the master azimuth direction.
		  {
		  m_oversample.argvalue.OsrAzimuth = strtol(word.charAt(1),null, BASE10); // pass keyword
		  writearg(m_oversample.argvalue.OsrAzimuth);
		  if (m_oversample.argvalue.OsrAzimuth < 1)
			{
			PRINT_ERROR("(M_OVS_FACT_AZI < 1) Oversampling ratio must be at least 1.");
			throw(keyword_error);
			}
		  if (m_oversample.argvalue.OsrAzimuth > 2)
			{
			PRINT_ERROR("(M_OVS_FACT_AZI > 2) Not implemented!");
			throw(keyword_error);
			}
		  }

	// ********************************************************************** 
		else if (!strcmp(keyword,"S_OVS_FACT_AZI")) // Oversampling ratio in the slave azimuth direction.
		  {
		  s_oversample.argvalue.OsrAzimuth = strtol(word.charAt(1),null, BASE10); // pass keyword
		  writearg(s_oversample.argvalue.OsrAzimuth);
		  if (s_oversample.argvalue.OsrAzimuth < 1)
			{
			PRINT_ERROR("(S_OVS_FACT_AZI < 1) Oversampling ratio must be at least 1.");
			throw(keyword_error);
			}
		  if (s_oversample.argvalue.OsrAzimuth > 2)
			{
			PRINT_ERROR("(S_OVS_FACT_AZI > 2) Not implemented!");
			throw(keyword_error);
			}
		  }

	// ********************************************************************** 
		else if (!strcmp(keyword,"M_OVS_KERNELSIZE")) // Length of the interpolation kernel.
		  {
		  m_oversample.argvalue.FilterSize = strtol(word.charAt(1),null, BASE10); // pass keyword
		  writearg(m_oversample.argvalue.FilterSize);
		  if (m_oversample.argvalue.FilterSize < 2)
			{
			 PRINT_ERROR("(M_OVS_KERNELSIZE < 2) Interpolation kernel length must be > 1.");
			throw(keyword_error);
			}
		  if (m_oversample.argvalue.FilterSize % 2)
			{
			PRINT_ERROR("(M_OVS_KERNELSIZE not even) Range Interpolation kernel length must be an even number.");
			throw(keyword_error);
			}
		  }

	// ********************************************************************** 
		else if (!strcmp(keyword,"S_OVS_KERNELSIZE")) // Length of the interpolation kernel.
		  {
		  s_oversample.argvalue.FilterSize = strtol(word.charAt(1),null, BASE10); // pass keyword
		  writearg(s_oversample.argvalue.FilterSize);
		  if (s_oversample.argvalue.FilterSize < 2)
			{
			PRINT_ERROR("code ???: (S_OVS_KERNELSIZE < 2) Interpolation kernel length must be > 1.");
			throw(keyword_error);
			}
		  if (s_oversample.argvalue.FilterSize % 2)
			{
			PRINT_ERROR("code ???: (S_OVS_KERNELSIZE not even) Interpolation kernel length must be an even number.");
			throw(keyword_error);
			}
		  }

	// ********************************************************************** 
		else if (!strcmp(keyword,"M_OVS_OUT_FORMAT")) // Output format [cr4] ci16, I suggest [cr4].
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"CR4"))
			m_oversample.argvalue.oformatflag = FORMATCR4; // default
		  else if (!strcmp(keyword,"CI2"))
			m_oversample.argvalue.oformatflag = FORMATCI2;
		  else
			{
			ERROR << "M_OVS_OUT_FORMAT: output format " << keyword << " not known for master range oversampling. line " << linecnt << ".";
			  PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }

	// ********************************************************************** 
		else if (!strcmp(keyword,"S_OVS_OUT_FORMAT")) // Output format [cr4] ci16, I suggest [cr4].
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"CR4"))
			s_oversample.argvalue.oformatflag = FORMATCR4; // default
		  else if (!strcmp(keyword,"CI2"))
			s_oversample.argvalue.oformatflag = FORMATCI2;
		  else
			{
			ERROR << "S_OVS_OUT_FORMAT: output format " << keyword << " not known for slave range oversampling. line " << linecnt << ".";
			  PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }
	//____RaffaeleNutricato END MODIFICATION SECTION 9

	// ____ start added by MA ____


	// **********************************************************************
	// *** SIMULATE AMPLITUDE FOR MASTER
	// **********************************************************************

		else if (!strcmp(keyword,"SAM_IN_DEM")) // input file
		  {
		  simampinput.argvalue.firefdem = word.charAt(1); // pass keyword
		  writearg(simampinput.argvalue.firefdem);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"SAM_IN_FORMAT")) // format input file
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"R4") || !strcmp(keyword,"REAL4"))
			simampinput.argvalue.iformatflag = FORMATR4;
		  else if (!strcmp(keyword,"I2") || !strcmp(keyword,"SHORT"))
			simampinput.argvalue.iformatflag = FORMATI2; // default
		  else if (!strcmp(keyword,"I2_BIGENDIAN") || !strcmp(keyword,"SHORT_BIGENDIAN"))
			simampinput.argvalue.iformatflag = FORMATI2_BIGENDIAN; // default
		  else if (!strcmp(keyword,"R8") || !strcmp(keyword,"REAL8"))
			simampinput.argvalue.iformatflag = FORMATR8;
		  else
			{
			ERROR << "SAM_IN_FORMAT: input format " << keyword << " not known (R4 R8 I2 (native) SHORT_BIGENDIAN); line " << linecnt << ".";
				  PRINT_ERROR(ERROR.get_str())
				  throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"SAM_IN_SIZE")) // nrow ncols (lat lon)
		  {
		  String pLast1;
		  String pLast2 = null;
		  simampinput.argvalue.demrows = strtoul(word.charAt(1), pLast1, BASE10);
		  simampinput.argvalue.demcols = strtoul(word.charAt(2), pLast2, BASE10);
		  if (pLast1.equals(word.charAt(1)) || pLast2.equals(word.charAt(2))) // fails to convert one of them to double.
		   {
			ERROR << "SAM_IN_SIZE: " << word.charAt(1) << " : " << word.charAt(2) << " are not valid.";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
		   }
		  writearg(simampinput.argvalue.demrows);
		  writearg(simampinput.argvalue.demcols);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"SAM_IN_DELTA")) // degrees delta lat lon
		  {
		  simampinput.argvalue.demdeltalat = Double.parseDouble(word.charAt(1)); // pass keyword
		  keyword = word.charAt(2); // pass keyword
		  writearg(simampinput.argvalue.demdeltalat);
		  writearg(keyword);
		  if (Character.isDigit(keyword.charAt(0)) || keyword.charAt(0)=='.') // likely to be 2 numbers
			simampinput.argvalue.demdeltalon = Double.parseDouble(keyword);
		  else // default same gridsize
			simampinput.argvalue.demdeltalon = simampinput.argvalue.demdeltalat;

		  // ______ Store as radians ______
		  simampinput.argvalue.demdeltalat = deg2rad(simampinput.argvalue.demdeltalat);
		  simampinput.argvalue.demdeltalon = deg2rad(simampinput.argvalue.demdeltalon);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"SAM_IN_UL")) // upperleft coordinates
		  {
		  String pLast1;
		  String pLast2 = null;
		  simampinput.argvalue.demlatleftupper = strtod(word.charAt(1), pLast1);
		  simampinput.argvalue.demlonleftupper = strtod(word.charAt(2), pLast2);
		  if (pLast1.equals(word.charAt(1)) || pLast2.equals(word.charAt(2))) // fails to convert one of them to double.
		   {
			ERROR << "SAM_IN_UL: " << word.charAt(1) << " : " << word.charAt(2) << " are not valid.";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
		   }
		  writearg(simampinput.argvalue.demlatleftupper);
		  writearg(simampinput.argvalue.demlonleftupper);
		  simampinput.argvalue.demlatleftupper = deg2rad(simampinput.argvalue.demlatleftupper);
		  simampinput.argvalue.demlonleftupper = deg2rad(simampinput.argvalue.demlonleftupper);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"SAM_IN_NODATA")) // flag for no data
		  {
		  String pLast = null;
		  simampinput.argvalue.demnodata = strtod(word.charAt(1), pLast);
		  if (pLast.equals(word.charAt(1))) // fails to convert to double.
		   {
			ERROR << "SAM_IN_NODATA: " << word.charAt(1) << " is not valid.";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
		   }
		  writearg(simampinput.argvalue.demnodata);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"SAM_OUT_DEM")) // name of output file
		  {
		  simampinput.argvalue.fodem = word.charAt(1); // pass keyword
		  writearg(simampinput.argvalue.fodem);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"SAM_OUT_FILE")) // name of output file
		  {
		  simampinput.argvalue.fosimamp = word.charAt(1); // pass keyword
		  writearg(simampinput.argvalue.fosimamp);
		  }

	// // **********************************************************************
	//     else if (!strcmp(keyword,"SAM_OUT_DEMI"))          // name of output file
	//       {
	//       strcpy(demassistinput.fodemi,  word[1] );      // pass keyword
	//       writearg(demassistinput.fodemi);
	//       }
	// 
	// // **********************************************************************
	//     else if (!strcmp(keyword,"SAM_OUT_DEM_LP"))        // name of output file
	//       {
	//       strcpy(demassistinput.forefdemhei,  word[1] );         // pass keyword
	//       writearg(demassistinput.forefdemhei);
	//       }


	// **********************************************************************
	// *** MASTER TIMING ERROR ESTIMATION using COREGISTRATION                             
	// **********************************************************************
		else if (!strcmp(keyword,"MTE_METHOD")) // method selector simamp coreg
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"MAGFFT"))
			{
			 mtiminginput.argvalue.method=cc_magfft; // default MTE_magfft
			}
		  else if (!strcmp(keyword,"MAGSPACE"))
			{
			 mtiminginput.argvalue.method=cc_magspace; // MTE_magspace
			}
		  else
			{
			 ERROR << "MTE_METHOD: method " << keyword << " not known for simamp correlation coregistration on line " << linecnt << ".";
			 PRINT_ERROR(ERROR.get_str())
				   throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"MTE_NWIN")) // #windows for simamp correlation
		  {
		  String pLast = null;
		  mtiminginput.argvalue.Nwin = strtoul(word.charAt(1), pLast, BASE10);
		  if (pLast.equals(word.charAt(1))) // fails to convert to double.
		   {
			ERROR << "MTE_NWIN: " << word.charAt(1) << " is not valid.";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
		   }
		  writearg(mtiminginput.argvalue.Nwin);
		  if (mtiminginput.argvalue.Nwin > 10000)
				 {
			PRINT_ERROR("Too many windows requested (MTE_NWIN > 10000).")
			throw(keyword_error);
				 }
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"MTE_IN_POS")) // file with #windows positions
		  {
		  mtiminginput.argvalue.ifpositions = word.charAt(1); // pass keyword
		  writearg(mtiminginput.argvalue.ifpositions);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"MTE_WINSIZE")) // windowsize for simamp correlation
		  {
		  String pLast1;
		  String pLast2 = null;
		  mtiminginput.argvalue.MasksizeL = strtoul(word.charAt(1), pLast1, BASE10);
		  mtiminginput.argvalue.MasksizeP = strtoul(word.charAt(2), pLast2, BASE10);
		  if (pLast1.equals(word.charAt(1)) || pLast2.equals(word.charAt(2))) // fails to convert one of them to double.
		   {
			ERROR << "MTE_WINSIZE: " << word.charAt(1) << " : " << word.charAt(2) << " are not valid.";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
		   }
		  writearg(mtiminginput.argvalue.MasksizeL);
		  writearg(mtiminginput.argvalue.MasksizeP);
		  if (mtiminginput.argvalue.MasksizeL > 4096 || mtiminginput.argvalue.MasksizeP > 4096)
				  {
			PRINT_ERROR("Too large correlation window (MTE_WINSIZE > 4096).");
			throw(keyword_error);
				  }
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"MTE_ACC")) // Searchdomain for correlation
		  { // simamp correlation
		  String pLast1;
		  String pLast2 = null;
		  mtiminginput.argvalue.AccL = strtoul(word.charAt(1), pLast1, BASE10);
		  mtiminginput.argvalue.AccP = strtoul(word.charAt(2), pLast2, BASE10);
		  if (pLast1.equals(word.charAt(1)) || pLast2.equals(word.charAt(2))) // fails to convert one of them to double.
		   {
			ERROR << "MTE_ACC: " << word.charAt(1) << " : " << word.charAt(2) << " are not valid.";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
		   }
		  writearg(mtiminginput.argvalue.AccL);
		  writearg(mtiminginput.argvalue.AccP);
		  if (mtiminginput.argvalue.AccL > 1000 || mtiminginput.argvalue.AccP > 1000)
				  {
			PRINT_ERROR("Too large searchwindow (MTE_ACC > 1000).");
			throw(keyword_error);
				  }
		  if (mtiminginput.argvalue.AccL == 0 || mtiminginput.argvalue.AccP == 0)
				  {
			PRINT_ERROR("Acc = 0 ?(MTE_ACC).");
			throw(keyword_error);
				  }
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"MTE_INITOFF")) // Initial offset
		  {
		  keyword = word.charAt(1); // pass keyword
		  String keyword2 = word.charAt(2); // pass keyword
		  writearg(keyword);
		  writearg(keyword2);
		  if (Character.isDigit(keyword2.charAt(0)) || keyword2.charAt(0)=='-') // likely to be 2 numbers
														// BK19/1/00 thanks to m.goos for '-'
			{
			mtiminginput.argvalue.initoffsetL = Integer.parseInt(keyword);
			mtiminginput.argvalue.initoffsetP = Integer.parseInt(keyword2);
			}
		  else
			{
			ERROR << "MTE_INITOFF: unknown input: " << keyword << ", " << keyword2 << " on line " << linecnt << ".";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }


	// ____ end added by MA ____


	// **********************************************************************
	// *** AZIMUTH FILTERING
	// **********************************************************************
		else if (!strcmp(keyword,"AF_BLOCKSIZE"))
		  {
		  String pLast = null;
		  filtaziinput.argvalue.fftlength = strtol(word.charAt(1), pLast, BASE10); // int32
		  if (pLast.equals(word.charAt(1))) // fails to convert to double.
		   {
			ERROR << "AF_BLOCKSIZE: " << word.charAt(1) << " is not valid.";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
		   }
		  writearg(filtaziinput.argvalue.fftlength);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"AF_HAMMING"))
		  {
		  String pLast = null;
		  filtaziinput.argvalue.hammingalpha = strtod(word.charAt(1), pLast);
		  if (pLast.equals(word.charAt(1))) // fails to convert to double.
		   {
			ERROR << "AF_HAMMING: " << word.charAt(1) << " is not valid.";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
		   }
		  writearg(filtaziinput.argvalue.hammingalpha);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"AF_OVERLAP"))
		  {
		  String pLast = null;
		  filtaziinput.argvalue.overlap = strtol(word.charAt(1), pLast, BASE10); // int32
		  if (pLast.equals(word.charAt(1))) // fails to convert to double.
		   {
			ERROR << "AF_OVERLAP: " << word.charAt(1) << " is not valid.";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
		   }
		  writearg(filtaziinput.argvalue.overlap);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"AF_OUT_MASTER"))
		  {
		  filtaziinput.argvalue.fomaster = word.charAt(1); // pass keyword
		  writearg(filtaziinput.argvalue.fomaster);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"AF_OUT_SLAVE"))
		  {
		  filtaziinput.argvalue.foslave = word.charAt(1); // pass keyword
		  writearg(filtaziinput.argvalue.foslave);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"AF_OUT_FORMAT")) // output format
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"CR4"))
			filtaziinput.argvalue.oformatflag = FORMATCR4; // default
		  else if (!strcmp(keyword,"CI2"))
			filtaziinput.argvalue.oformatflag = FORMATCI2;
		  else
			{
			ERROR << "AF_OUT_FORMAT: output format " << keyword << " not known for azimuth filtering. line " << linecnt << ".";
				  PRINT_ERROR(ERROR.get_str())
				  throw(keyword_error);
			}
		  }



	// **********************************************************************
	// *** COARSE CORR COREGISTRATION
	// **********************************************************************
		else if (!strcmp(keyword,"CC_METHOD")) // method selector coarse coreg
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"MAGFFT"))
			coarsecorrinput.argvalue.method=cc_magfft; // default
		  else if (!strcmp(keyword,"MAGSPACE"))
			coarsecorrinput.argvalue.method=cc_magspace;
		  else
			{
			ERROR << "CC_METHOD: method " << keyword << " not known for coarse correlation coregistration on line " << linecnt << ".";
				  PRINT_ERROR(ERROR.get_str())
				  throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CC_NWIN")) // #windows for coarse correlation
		  {
		  coarsecorrinput.argvalue.Nwin = Integer.parseInt(word.charAt(1)); // pass keyword
		  writearg(coarsecorrinput.argvalue.Nwin);
		  if (coarsecorrinput.argvalue.Nwin > 10000)
			{
			PRINT_ERROR("Too many windows requested (CC_NWIN > 10000).")
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CC_IN_POS")) // file with #windows positions
		  {
		  coarsecorrinput.argvalue.ifpositions = word.charAt(1); // pass keyword
		  writearg(coarsecorrinput.argvalue.ifpositions);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CC_WINSIZE")) // windowsize for coarse correlation
		  {
		  coarsecorrinput.argvalue.MasksizeL = Integer.parseInt(word.charAt(1)); // pass keyword
		  coarsecorrinput.argvalue.MasksizeP = Integer.parseInt(word.charAt(2)); // pass keyword
		  writearg(coarsecorrinput.argvalue.MasksizeL);
		  writearg(coarsecorrinput.argvalue.MasksizeP);
		  if (coarsecorrinput.argvalue.MasksizeL > 2048 || coarsecorrinput.argvalue.MasksizeP > 2048)
				  {
			PRINT_ERROR("Too large correlation window (CC_WINSIZE > 2048).");
			throw(keyword_error);
				  }
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CC_ACC")) // Searchdomain for correlation
		  { // coarse correlation
		  coarsecorrinput.argvalue.AccL = Integer.parseInt(word.charAt(1)); // pass keyword
		  coarsecorrinput.argvalue.AccP = Integer.parseInt(word.charAt(2)); // pass keyword
		  writearg(coarsecorrinput.argvalue.AccL);
		  writearg(coarsecorrinput.argvalue.AccP);
		  if (coarsecorrinput.argvalue.AccL > 1000 || coarsecorrinput.argvalue.AccP > 1000)
				  {
			PRINT_ERROR("Too large searchwindow (CC_ACC > 1000).");
			throw(keyword_error);
				  }
		  if (coarsecorrinput.argvalue.AccL == 0 || coarsecorrinput.argvalue.AccP == 0)
				  {
			PRINT_ERROR("Acc = 0 ?(CC_ACC).");
			throw(keyword_error);
				  }
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CC_INITOFF")) // Initial offset
		  {
		  keyword = word.charAt(1); // pass keyword
		  String keyword2 = word.charAt(2); // pass keyword
		  writearg(keyword);
		  writearg(keyword2);
		  if (!strcmp(keyword,"ORBIT") || !strcmp(keyword,"orbit"))
			{
			coarsecorrinput.argvalue.initoffsetL = NaN; // flag used in main, see processor.cc
			coarsecorrinput.argvalue.initoffsetP = NaN; // flag used in main
			INFO << "CC_INITOFF: \tInitial offsets from COARSEORB: " << generalinput.argvalue.i_resfile;
			INFO.print();
			}
		  else if (Character.isDigit(keyword2.charAt(0)) || keyword2.charAt(0)=='-') // likely to be 2 numbers
																			   // BK19/1/00 thanks to m.goos for '-'
			{
			//coarsecorrinput.initoffsetL = atof(keyword);
			//coarsecorrinput.initoffsetP = atof(keyword2);
			coarsecorrinput.argvalue.initoffsetL = Integer.parseInt(keyword);
			coarsecorrinput.argvalue.initoffsetP = Integer.parseInt(keyword2);
			}
		  else
			{
			ERROR << "CC_INITOFF: unknown input: " << keyword << ", " << keyword2 << " on line " << linecnt << ".";
				  PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }

	// **********************************************************************
	// *** FINE COREGISTRATION
	// **********************************************************************
		else if (!strcmp(keyword,"FC_NWIN")) // #windows for fine correlation
		  {
		  fineinput.argvalue.Nwin = Integer.parseInt(word.charAt(1)); // pass keyword
		  writearg(fineinput.argvalue.Nwin);
		  if (fineinput.argvalue.Nwin > 100000)
				  {
			PRINT_ERROR("Too many windows requested (FC_NWIN).")
			throw(keyword_error);
				  }
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"FC_IN_POS")) // file with #windows positions
		  {
		  fineinput.argvalue.ifpositions = word.charAt(1); // pass keyword
		  writearg(fineinput.argvalue.ifpositions);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"FC_WINSIZE")) // windowsize for fine correlation
		  {
		  fineinput.argvalue.MasksizeL = Integer.parseInt(word.charAt(1)); // pass keyword
		  fineinput.argvalue.MasksizeP = Integer.parseInt(word.charAt(2)); // pass keyword
		  writearg(fineinput.argvalue.MasksizeL);
		  writearg(fineinput.argvalue.MasksizeP);
		  if (fineinput.argvalue.MasksizeL > 1024 || fineinput.argvalue.MasksizeP > 1024)
			{
			PRINT_ERROR("Too large correlation window (FC_WINSIZE).")
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"FC_ACC")) // Searchdomain for correlation
		  { // fine correlation
		  fineinput.argvalue.AccL = Integer.parseInt(word.charAt(1)); // pass keyword
		  fineinput.argvalue.AccP = Integer.parseInt(word.charAt(2)); // pass keyword
		  writearg(fineinput.argvalue.AccL);
		  writearg(fineinput.argvalue.AccP);
		  if (fineinput.argvalue.AccL > 1000 || fineinput.argvalue.AccP > 1000)
			{
			PRINT_ERROR("Too large searchwindow (FC_ACC).")
			throw(keyword_error);
			}
		  if (fineinput.argvalue.AccL == 0 || fineinput.argvalue.AccP == 0)
			{
			PRINT_ERROR("Acc = 0 ?(FC_ACC).")
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"FC_INITOFF")) // Initial offset
		  {
		  keyword = word.charAt(1); // pass keyword
		  String keyword2 = word.charAt(2); // pass keyword
		  writearg(keyword);
		  writearg(keyword2);
		  if (!strcmp(keyword,"COARSECORR") || !strcmp(keyword,"coarsecorr")) // use value of coarse correlation
			{
			fineinput.argvalue.initoffsetL = NaN; // flag used in main
			fineinput.argvalue.initoffsetP = NaN; // flag used in main
			INFO << "FC_INITOFF: \tInitial offset from COARSECORR: " << generalinput.argvalue.i_resfile;
			INFO.print();
			}
		  else if (Character.isDigit(keyword2.charAt(0)) || keyword2.charAt(0)=='-') // likely to be 2 numbers
													// BK19/1/00 thanks to M.Goos for '-'
			{
			fineinput.argvalue.initoffsetL = Integer.parseInt(keyword);
			fineinput.argvalue.initoffsetP = Integer.parseInt(keyword2);
			}
		  else
			{
			ERROR << "FC_INITOFF: unknown input: " << keyword << ", " << keyword2 << " on line " << linecnt << ".";
				 PRINT_ERROR(ERROR.get_str())
				 throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"FC_METHOD")) // method selector fine coreg
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"CMPLXFFT"))
		   {
			//fineinput.method=fc_cmplxfft;                 // default
			PRINT_ERROR("CMPLXFFT not implemented in v1.0 of Doris.")
			throw(keyword_error);
		   }
		  else if (!strcmp(keyword,"CMPLXSPACE"))
		   {
			//fineinput.method=fc_cmplxspace;
		   PRINT_ERROR("CMPLXSPACE not implemented in v1.0 of Doris.")
		   throw(keyword_error);
		   }
		  else if (!strcmp(keyword,"MAGFFT"))
			fineinput.argvalue.method=fc_magfft;
		  else if (!strcmp(keyword,"MAGSPACE"))
			fineinput.argvalue.method=fc_magspace;
		  else if (!strcmp(keyword,"OVERSAMPLE"))
			fineinput.argvalue.method=fc_oversample;
		  else
			{
			ERROR << "FC_METHOD: method " << keyword << " not known for fine coregistration on line " << linecnt << ".";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"FC_OSFACTOR")) // oversampling factor
		  {
		  fineinput.argvalue.osfactor = Integer.parseInt(word.charAt(1)); // pass keyword
		  writearg(fineinput.argvalue.osfactor);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"FC_PLOT")) // plotting results
		  {
		  fineinput.argvalue.plotoffsets = true;
		  fineinput.argvalue.plotthreshold = Double.parseDouble(word.charAt(1)); // pass keyword
		  keyword = word.charAt(2); // pass keyword
		  writearg(fineinput.argvalue.plotthreshold);
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"BG"))
			fineinput.argvalue.plotmagbg = true;
		  else if (!strcmp(keyword,"NOBG"))
			fineinput.argvalue.plotmagbg = false;
		  else if (!strcmp(keyword,"ON")) // actually not allowed...
			fineinput.argvalue.plotoffsets = true;
		  else if (!strcmp(keyword,"OFF")) // actually not allowed...
			fineinput.argvalue.plotoffsets = false;
		  else
			WARNING.print("FC_PLOT: missing argument(s). (default: 0.4 NOBG)");
		  }

	// ____ start added by FvL ____

	// **********************************************************************
	// *** RELATIVE TIMING ERROR
	// **********************************************************************

		else if (!strcmp(keyword,"RTE_THRESHOLD")) // treshhold value
		  {
		  reltiminginput.argvalue.threshold = Double.parseDouble(word.charAt(1)); // pass keyword
		  writearg(reltiminginput.argvalue.threshold);
		  if (reltiminginput.argvalue.threshold > 1)
			{
			PRINT_ERROR("RTE_THRESHOLD: threshold > 1.")
				  throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"RTE_MAXITER")) // max number of offsets to reject
		  {
		  reltiminginput.argvalue.maxiter = Integer.parseInt(word.charAt(1)); // pass keyword
		  writearg(reltiminginput.argvalue.maxiter);
		  if (reltiminginput.argvalue.maxiter < 0)
			{
			WARNING.print("RTE_MAXITER: max. number of points to remove < 0? (using 0)");
			reltiminginput.argvalue.maxiter = 0;
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"RTE_K_ALPHA")) // critical value for outlier removal
		  {
		  reltiminginput.argvalue.k_alpha = Double.parseDouble(word.charAt(1)); // pass keyword
		  writearg(reltiminginput.argvalue.k_alpha);
		  if (reltiminginput.argvalue.k_alpha < 0)
			{
			WARNING.print("RTE_K_ALPHA: critical value < 0.0?");
			reltiminginput.argvalue.k_alpha = 1.97;
			}
		  }


	// **********************************************************************
	// *** DEM ASSISTED COREGISTRATION
	// **********************************************************************

		else if (!strcmp(keyword,"DAC_IN_DEM")) // input file
		  {
		  demassistinput.argvalue.firefdem = word.charAt(1); // pass keyword
		  writearg(demassistinput.argvalue.firefdem);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"DAC_IN_FORMAT")) // format input file
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"R4") || !strcmp(keyword,"REAL4"))
			demassistinput.argvalue.iformatflag = FORMATR4;
		  else if (!strcmp(keyword,"I2") || !strcmp(keyword,"SHORT"))
			demassistinput.argvalue.iformatflag = FORMATI2; // default
		  else if (!strcmp(keyword,"I2_BIGENDIAN") || !strcmp(keyword,"SHORT_BIGENDIAN"))
			demassistinput.argvalue.iformatflag = FORMATI2_BIGENDIAN; // default
		  else if (!strcmp(keyword,"R8") || !strcmp(keyword,"REAL8"))
			demassistinput.argvalue.iformatflag = FORMATR8;
		  else
			{
			ERROR << "DAC_IN_FORMAT: input format " << keyword << " not known (R4 R8 I2 (native) SHORT_BIGENDIAN); line " << linecnt << ".";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"DAC_IN_SIZE")) // nrow ncols (lat lon)
		  {
		  String pLast1;
		  String pLast2 = null;
		  demassistinput.argvalue.demrows = strtoul(word.charAt(1), pLast1, BASE10);
		  demassistinput.argvalue.demcols = strtoul(word.charAt(2), pLast2, BASE10);
		  if (pLast1.equals(word.charAt(1)) || pLast2.equals(word.charAt(2))) // fails to convert one of them to double.
		   {
			ERROR << "DAC_IN_SIZE: " << word.charAt(1) << " : " << word.charAt(2) << " are not valid.";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
		   }
		  writearg(demassistinput.argvalue.demrows);
		  writearg(demassistinput.argvalue.demcols);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"DAC_IN_DELTA")) // degrees delta lat lon
		  {
		  demassistinput.argvalue.demdeltalat = Double.parseDouble(word.charAt(1)); // pass keyword
		  keyword = word.charAt(2); // update keyword
		  writearg(demassistinput.argvalue.demdeltalat);
		  writearg(keyword); // lon
		  if (Character.isDigit(keyword.charAt(0)) || keyword.charAt(0)=='.') // likely to be 2 numbers
			demassistinput.argvalue.demdeltalon = Double.parseDouble(keyword);
		  else // default same gridsize
			demassistinput.argvalue.demdeltalon = demassistinput.argvalue.demdeltalat;

		  // ______ Store as radians ______
		  demassistinput.argvalue.demdeltalat = deg2rad(demassistinput.argvalue.demdeltalat);
		  demassistinput.argvalue.demdeltalon = deg2rad(demassistinput.argvalue.demdeltalon);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"DAC_IN_UL")) // upperleft coordinates
		  {
		  String pLast1;
		  String pLast2 = null;
		  demassistinput.argvalue.demlatleftupper = strtod(word.charAt(1), pLast1);
					  demassistinput.argvalue.demlonleftupper = strtod(word.charAt(2), pLast2);
		  if (pLast1.equals(word.charAt(1)) || pLast2.equals(word.charAt(2))) // fails to convert
		   {
			ERROR << "DAC_IN_UL: " << word.charAt(1) << " : " << word.charAt(2) << " are not valid.";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
		   }

		  writearg(demassistinput.argvalue.demlatleftupper);
		  writearg(demassistinput.argvalue.demlonleftupper);
		  demassistinput.argvalue.demlatleftupper = deg2rad(demassistinput.argvalue.demlatleftupper);
		  demassistinput.argvalue.demlonleftupper = deg2rad(demassistinput.argvalue.demlonleftupper);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"DAC_IN_NODATA")) // flag for no data
		  {
		  demassistinput.argvalue.demnodata = Double.parseDouble(word.charAt(1)); // pass keyword
		  writearg(demassistinput.argvalue.demnodata);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"DAC_OUT_DEM")) // name of output file
		  {
		  demassistinput.argvalue.fodem = word.charAt(1); // pass keyword
		  writearg(demassistinput.argvalue.fodem);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"DAC_OUT_DEMI")) // name of output file
		  {
		  demassistinput.argvalue.fodemi = word.charAt(1); // pass keyword
		  writearg(demassistinput.argvalue.fodemi);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"DAC_OUT_DEM_LP")) // name of output file
		  {
		  demassistinput.argvalue.forefdemhei = word.charAt(1); // pass keyword
		  writearg(demassistinput.argvalue.forefdemhei);
		  }

	// ____ end added by FvL ____


	// **********************************************************************
	// *** COMPUTATION OF COREGISTRATION PARAMETERS
	// **********************************************************************
		else if (!strcmp(keyword,"CPM_THRESHOLD")) // treshhold value
		  {
		  coregpminput.argvalue.threshold = Double.parseDouble(word.charAt(1)); // pass keyword
		  writearg(coregpminput.argvalue.threshold);
		  if (coregpminput.argvalue.threshold > 1)
			{
			PRINT_ERROR("CPM_THRESHOLD: threshold > 1.")
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CPM_DEGREE")) // degree of polynomial
		  {
		  coregpminput.argvalue.degree = Integer.parseInt(word.charAt(1)); // pass keyword
		  writearg(coregpminput.argvalue.degree);
		  if (coregpminput.argvalue.degree > 4)
			WARNING.print("CPM_DEGREE: degree > 4 dangerous at edges?");
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CPM_WEIGHT")) // weightmatrix
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"NONE"))
			coregpminput.argvalue.weightflag = 0; // default
		  else if (!strcmp(keyword,"LINEAR"))
			coregpminput.argvalue.weightflag = 1;
		  else if (!strcmp(keyword,"QUADRATIC"))
			coregpminput.argvalue.weightflag = 2;
		  else if (!strcmp(keyword,"BAMLER"))
			coregpminput.argvalue.weightflag = 3;
		  else
			{
			ERROR << "CPM_WEIGHT: data weighting option: " << keyword << " not known for computation of coregistration parameters on line " << linecnt << ".";
				PRINT_ERROR(ERROR.get_str())
				throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CPM_MAXITER")) // max number of offsets to reject
		  {
		  coregpminput.argvalue.maxiter = Integer.parseInt(word.charAt(1)); // pass keyword
		  writearg(coregpminput.argvalue.maxiter);
		  if (coregpminput.argvalue.maxiter < 0)
			{
			WARNING.print("CPM_MAXITER: max. number of points to remove < 0? (using 0)");
			coregpminput.argvalue.maxiter = 0;
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CPM_K_ALPHA")) // critical value for outlier removal
		  {
		  coregpminput.argvalue.k_alpha = Double.parseDouble(word.charAt(1)); // pass keyword
		  writearg(coregpminput.argvalue.k_alpha);
		  if (coregpminput.argvalue.k_alpha < 0)
			{
			WARNING.print("CPM_K_ALPHA: critical value < 0.0?");
			coregpminput.argvalue.k_alpha = 1.97;
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CPM_DUMP")) // boolean
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"OFF"))
			coregpminput.argvalue.dumpmodel = false;
		  else if (!strcmp(keyword,"ON") || !strncmp(keyword,"//",2) || !strncmp(keyword,"#",1) || !(keyword.charAt(0) == '\0')) // no keyword
			coregpminput.argvalue.dumpmodel = true;
		  else
			{
			WARNING << "CPM_DUMP: line " << linecnt << ": argument: " << keyword << " not recognized, no dumping to files.";
			WARNING.print();
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CPM_PLOT")) // plotting results
		  {
		  coregpminput.argvalue.plot = true;
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"BG"))
			coregpminput.argvalue.plotmagbg = true;
		  else if (!strcmp(keyword,"NOBG"))
			coregpminput.argvalue.plotmagbg = false;
		  else if (!strcmp(keyword,"ON")) // actually not allowed...
			coregpminput.argvalue.plot = true;
		  else if (!strcmp(keyword,"OFF")) // actually not allowed...
			coregpminput.argvalue.plot = false;
		  else
			WARNING.print("CPM_PLOT: missing argument. (default NOBG magnitude background)");
		  }


	// **********************************************************************
	// *** RANGE FILTERING
	// **********************************************************************
		else if (!strcmp(keyword,"RF_METHOD")) // method range filtering
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"ADAPTIVE"))
			filtrangeinput.argvalue.method = rf_adaptive;
		  else if (!strcmp(keyword,"PORBITS"))
			filtrangeinput.argvalue.method = rf_porbits;
		  else
			{
			ERROR << "RF_METHOD: method " << keyword << " not known for range filtering. line " << linecnt << ".";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"RF_SLOPE")) // terrain method porbits
		  {
		  filtrangeinput.argvalue.terrainslope = Double.parseDouble(word.charAt(1)); // pass keyword
		  writearg(filtrangeinput.argvalue.terrainslope);
		  deg2rad(filtrangeinput.argvalue.terrainslope);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"RF_THRESHOLD")) // threshhold value
		  {
		  filtrangeinput.argvalue.SNRthreshold = Double.parseDouble(word.charAt(1)); // pass keyword
		  writearg(filtrangeinput.argvalue.SNRthreshold);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"RF_HAMMING")) // alpha
		  {
		  filtrangeinput.argvalue.hammingalpha = Double.parseDouble(word.charAt(1)); // pass keyword
		  writearg(filtrangeinput.argvalue.hammingalpha);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"RF_WEIGHTCORR")) // boolean
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"OFF"))
			filtrangeinput.argvalue.doweightcorrel = false;
		  else if (!strcmp(keyword,"ON") || !strncmp(keyword,"//",2) || !strncmp(keyword,"#",1) || !(keyword.charAt(0) == '\0')) // no keyword
			filtrangeinput.argvalue.doweightcorrel = true;
		  else
			{
			filtrangeinput.argvalue.doweightcorrel = false; // default already...
			WARNING << "RF_WEIGHTCORREL: line " << linecnt << ": argument: " << keyword << " not recognized, weighting correlation set to OFF.";
			WARNING.print();
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"RF_OVERSAMPLE")) // int
		  {
		  filtrangeinput.argvalue.oversample = Integer.parseInt(word.charAt(1)); // pass keyword
		  writearg(filtrangeinput.argvalue.oversample);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"RF_NLMEAN")) // take mean over nlmean lines
		  {
		  filtrangeinput.argvalue.nlmean = Integer.parseInt(word.charAt(1)); // pass keyword
		  writearg(filtrangeinput.argvalue.nlmean);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"RF_FFTLENGTH")) // adaptive length
		  {
		  filtrangeinput.argvalue.fftlength = Integer.parseInt(word.charAt(1)); // pass keyword
		  writearg(filtrangeinput.argvalue.fftlength);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"RF_OVERLAP")) // overlap blocks
		  {
		  filtrangeinput.argvalue.overlap = Integer.parseInt(word.charAt(1)); // pass keyword
		  writearg(filtrangeinput.argvalue.overlap);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"RF_OUT_MASTER")) // filename
		  {
		  filtrangeinput.argvalue.fomaster = word.charAt(1); // pass keyword
		  writearg(filtrangeinput.argvalue.fomaster);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"RF_OUT_SLAVE")) // filename
		  {
		  filtrangeinput.argvalue.foslave = word.charAt(1); // pass keyword
		  writearg(filtrangeinput.argvalue.foslave);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"RF_OUT_FORMAT")) // output format
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"CR4"))
			filtrangeinput.argvalue.oformatflag = FORMATCR4; // default
		  else if (!strcmp(keyword,"CI2"))
			filtrangeinput.argvalue.oformatflag = FORMATCI2;
		  else
			{
			ERROR << "RF_OUT_FORMAT: output format " << keyword << " not known for range filtering. line " << linecnt << ".";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }


	// **********************************************************************
	// *** FLAT EARTH CORRECTION == compute reference phase since Feb-2000
	// **********************************************************************
		else if (!strcmp(keyword,"FE_METHOD")) // method selector flatearth
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"PORBITS"))
			comprefphainput.argvalue.method=fe_porbits;
		  else if (!strcmp(keyword,"METHOD2"))
			comprefphainput.argvalue.method=fe_method2;
		  else
			{
			ERROR << "FE_METHOD: argument: " << keyword << " not recognized.";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"FE_DEGREE")) // degree for flat earth correction
		  {
		  comprefphainput.argvalue.degree = Integer.parseInt(word.charAt(1)); // pass keyword
		  writearg(comprefphainput.argvalue.degree);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"FE_NPOINTS")) // number of points used
		  { // for estimation of polynomial
		  comprefphainput.argvalue.Npoints = Integer.parseInt(word.charAt(1)); // pass keyword // flat earth correction.
		  writearg(comprefphainput.argvalue.Npoints);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"FE_IN_POS")) // file with #windows positions
		  {
		  comprefphainput.argvalue.ifpositions = word.charAt(1); // pass keyword
		  writearg(comprefphainput.argvalue.ifpositions);
		  }


	// **********************************************************************
	// *** RESAMPLING (SLAVE)
	// **********************************************************************
		else if (!strcmp(keyword,"RS_METHOD")) // method selector resampling
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"CC4P"))
			resampleinput.argvalue.method = rs_cc4p; // default
		  else if (!strcmp(keyword,"CC6P"))
			resampleinput.argvalue.method = rs_cc6p;
		  else if (!strcmp(keyword,"TS6P"))
			resampleinput.argvalue.method = rs_ts6p;
		  else if (!strcmp(keyword,"TS8P"))
			resampleinput.argvalue.method = rs_ts8p;
		  else if (!strcmp(keyword,"TS16P"))
			resampleinput.argvalue.method = rs_ts16p;
		  else if (!strcmp(keyword,"KNAB4P"))
			resampleinput.argvalue.method = rs_knab4p;
		  else if (!strcmp(keyword,"KNAB6P"))
			resampleinput.argvalue.method = rs_knab6p;
		  else if (!strcmp(keyword,"KNAB8P"))
			resampleinput.argvalue.method = rs_knab8p;
		  else if (!strcmp(keyword,"KNAB10P"))
			resampleinput.argvalue.method = rs_knab10p;
		  else if (!strcmp(keyword,"KNAB16P"))
			resampleinput.argvalue.method = rs_knab16p;
		  else if (!strcmp(keyword,"RC6P"))
			resampleinput.argvalue.method = rs_rc6p;
		  else if (!strcmp(keyword,"RC12P"))
			resampleinput.argvalue.method = rs_rc12p;
		  else if (!strcmp(keyword,"RECT"))
			resampleinput.argvalue.method = rs_rect;
		  else if (!strcmp(keyword,"TRI"))
			resampleinput.argvalue.method = rs_tri;
		  else
			{
			ERROR << "RS_METHOD: method " << keyword << " not known for resampling. line " << linecnt << ".";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"RS_DBOW")) // database output window
		  {
		  resampleinput.argvalue.dbow.linelo = Integer.parseInt(word.charAt(1));
		  resampleinput.argvalue.dbow.linehi = Integer.parseInt(word.charAt(2));
		  resampleinput.argvalue.dbow.pixlo = Integer.parseInt(word.charAt(3));
		  resampleinput.argvalue.dbow.pixhi = Integer.parseInt(word.charAt(4));
		  writearg(resampleinput.argvalue.dbow.linelo);
		  writearg(resampleinput.argvalue.dbow.linehi);
		  writearg(resampleinput.argvalue.dbow.pixlo);
		  writearg(resampleinput.argvalue.dbow.pixhi);
		  if (resampleinput.argvalue.dbow.linelo <= 0 || resampleinput.argvalue.dbow.pixlo <= 0 || resampleinput.argvalue.dbow.linelo > resampleinput.argvalue.dbow.linehi || resampleinput.argvalue.dbow.pixlo > resampleinput.argvalue.dbow.pixhi)
			{
			ERROR << "code 300: Arguments of RS_DBOW card on line " << linecnt << " missing. [RS_DBOW  min_line  max_line  min_pixel  max_pixel].";
				  PRINT_ERROR(ERROR.get_str())
				  throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"RS_DBOW_GEO")) // database output
												  // window on GEO.LOCATION
		  {
		   // use dbow to store temp, compute later
		   real8 tmp_lat_0;
		   real8 tmp_lon_0;
		   real8 tmp_height; // total height, not half.window
		   real8 tmp_width; // total width, not half.windo

		   String pLast1;
		   String pLast2;
		   String pLast3;
		   String pLast4 = null;
		   tmp_lat_0 = strtod(word.charAt(1), pLast1);
		   tmp_lon_0 = strtod(word.charAt(2), pLast2);
		   tmp_height = strtod(word.charAt(3), pLast2);
		   tmp_width = strtod(word.charAt(4), pLast2);
		  if (pLast1.equals(word.charAt(1)) || pLast2.equals(word.charAt(2)) || pLast3.equals(word.charAt(3)) || pLast4.equals(word.charAt(4))) // fails to convert one of them to double.
		   {
			ERROR << "RS_DBOW_GEO: " << word.charAt(1) << " : " << word.charAt(2) << " : " << word.charAt(3) << " : " << word.charAt(4) << " are not valid.";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
		   }
			writearg(tmp_lat_0);
			writearg(tmp_lon_0);
			writearg(tmp_height);
			writearg(tmp_width);

			resampleinput.argvalue.dbow_geo.linelo = uint((360.0+tmp_lat_0)*1e6);
			resampleinput.argvalue.dbow_geo.linehi = uint((360.0+tmp_lon_0)*1e6);
			resampleinput.argvalue.dbow_geo.pixlo = uint(tmp_height);
			resampleinput.argvalue.dbow_geo.pixhi = uint(tmp_width);

		  }

	// **********************************************************************
		else if (!strcmp(keyword,"RS_OUT_FILE")) // name of output file
		  {
		  switch (priorrs_fileout)
			{
			case true:
			  WARNING << "RS_OUT_FILE: line: " << linecnt << ": " << "ignored due to prior occurence.";
			  WARNING.print();
			  break;
			default:
			  priorrs_fileout = true;
			  resampleinput.argvalue.fileout = word.charAt(1); // pass keyword (filename)
			  writearg(resampleinput.argvalue.fileout);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"RS_OUT_FORMAT")) // output format
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"CR4"))
			resampleinput.argvalue.oformatflag = FORMATCR4; // default
		  else if (!strcmp(keyword,"CI2"))
			resampleinput.argvalue.oformatflag = FORMATCI2;
		  else
			{
			ERROR << "RS_OUT_FORMAT: output format " << keyword << " not known for resampling. line " << linecnt << ".";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"RS_SHIFTAZI")) // true: shift before rs.
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"ON"))
			resampleinput.argvalue.shiftazi = true;
		  else if (!strcmp(keyword,"OFF") || !strncmp(keyword,"//",2) || !strncmp(keyword,"#",1) || !(keyword.charAt(0) == '\0')) // no keyword
			resampleinput.argvalue.shiftazi = false;
		  else
			{
			resampleinput.argvalue.shiftazi = true;
			WARNING << "RS_SHIFTAZI: line: " << linecnt << ": unknown argument: " << keyword << "; Set to ON (do shift azimuth spectrum).";
			WARNING.print();
			}
		  }


	// **********************************************************************
	// *** COMPUTATION OF INTERFEROGRAM
	// **********************************************************************
		else if (!strcmp(keyword,"INT_METHOD")) // method selector interfero
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"OLD"))
			interferoinput.argvalue.method = int_oldmethod;
		  else if (!strcmp(keyword,"OVERSAMPLE"))
			interferoinput.argvalue.method = int_oversample;
		  else
			{
			ERROR << "INT_METHOD: method " << keyword << " not known for interfero. line " << linecnt << ".";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"INT_OUT_INT")) // name of output file
		  {
		  interferoinput.argvalue.foint = word.charAt(1); // pass keyword
		  writearg(interferoinput.argvalue.foint);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"INT_OUT_CINT")) // name of output file
		  {
		  interferoinput.argvalue.focint = word.charAt(1); // pass keyword
		  writearg(interferoinput.argvalue.focint);
		  }

	// **********************************************************************
	//    else if (!strcmp(keyword,"INT_OUT_FE"))         // name of output file
	//      interferoinput.foflatearth =  word[1] ;         // pass keyword

	// **********************************************************************
		else if (!strcmp(keyword,"INT_MULTILOOK")) // multilookfactors
		  {
		  interferoinput.argvalue.multilookL = Integer.parseInt(word.charAt(1)); // pass keyword
		  interferoinput.argvalue.multilookP = Integer.parseInt(word.charAt(2)); // pass keyword
		  writearg(interferoinput.argvalue.multilookL);
		  writearg(interferoinput.argvalue.multilookP);
		  }

	// **********************************************************************
	// *** COMPUTATION OF COHERENCE
	// **********************************************************************
		else if (!strcmp(keyword,"COH_METHOD")) // method selector coherence
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"REFPHASE_ONLY")) // OLD
			coherenceinput.argvalue.method = coh_oldmethod;
		  else if (!strcmp(keyword,"INCLUDE_REFDEM")) // NEW
			coherenceinput.argvalue.method = coh_newmethod;
		  else
			{
			ERROR << "COH_METHOD: method " << keyword << " not known for coherence. line " << linecnt << ".";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"COH_OUT_COH")) // name of output file
		  {
		  coherenceinput.argvalue.focoh = word.charAt(1); // pass keyword
		  writearg(coherenceinput.argvalue.focoh);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"COH_OUT_CCOH")) // name of output file
		  {
		  coherenceinput.argvalue.foccoh = word.charAt(1); // pass keyword
		  writearg(coherenceinput.argvalue.foccoh);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"COH_MULTILOOK")) // multilookfactors
		  {
		  coherenceinput.argvalue.multilookL = Integer.parseInt(word.charAt(1)); // pass keyword
		  coherenceinput.argvalue.multilookP = Integer.parseInt(word.charAt(2)); // pass keyword
		  writearg(coherenceinput.argvalue.multilookL);
		  writearg(coherenceinput.argvalue.multilookP);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"COH_WINSIZE")) // estimator winsize
		  {
		  coherenceinput.argvalue.cohsizeL = Integer.parseInt(word.charAt(1)); // pass keyword
		  coherenceinput.argvalue.cohsizeP = Integer.parseInt(word.charAt(2)); // pass keyword
		  writearg(coherenceinput.argvalue.cohsizeL);
		  writearg(coherenceinput.argvalue.cohsizeP);
		  }

	// **********************************************************************
	// *** SUBTRACTION OF REFERENCE PHASE
	// **********************************************************************
		else if (!strcmp(keyword,"SRP_METHOD")) // method selector ref. phase
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"POLYNOMIAL"))
			subtrrefphainput.argvalue.method = srp_polynomial;
		  else if (!strcmp(keyword,"EXACT"))
			subtrrefphainput.argvalue.method = srp_exact;
		  else
			{
			ERROR << "SRP_METHOD: method " << keyword << " not known for subtraction of ref. phase. line " << linecnt << ".";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"SRP_OUT_CINT")) // name of output file
		  {
		  subtrrefphainput.argvalue.focint = word.charAt(1); // pass keyword
		  writearg(subtrrefphainput.argvalue.focint);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"SRP_OUT_REFPHA")) // name of output file
		  {
		  subtrrefphainput.argvalue.forefpha = word.charAt(1); // pass keyword
		  writearg(subtrrefphainput.argvalue.forefpha);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"SRP_DUMPREFPHA")) // true: dump ref.pha
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"ON"))
			subtrrefphainput.argvalue.dumponlyrefpha = true;
		  else if (!strcmp(keyword,"OFF") || !strncmp(keyword,"//",2) || !strncmp(keyword,"#",1) || !(keyword.charAt(0) == '\0')) // no keyword
			subtrrefphainput.argvalue.dumponlyrefpha = false;
		  else
			{
			subtrrefphainput.argvalue.dumponlyrefpha = false;
			WARNING << "SRP_DUMPREFPHA: line: " << linecnt << ": unknown argument: " << keyword << "; Set to OFF (no dump).";
			WARNING.print();
			}
		  }

	// **********************************************************************

	// ___________ added by FvL      
		else if (!strcmp(keyword,"SRP_OUT_H2PH")) // name of output file
		  {
		  subtrrefphainput.argvalue.foh2ph = word.charAt(1); // pass keyword
		  writearg(subtrrefphainput.argvalue.foh2ph);
		  }
	// ___________ end added by FvL


	// **********************************************************************
		else if (!strcmp(keyword,"SRP_MULTILOOK")) // multilookfactors
		  {
		  subtrrefphainput.argvalue.multilookL = Integer.parseInt(word.charAt(1)); // pass keyword
		  keyword = word.charAt(2); // pass keyword
		  writearg(subtrrefphainput.argvalue.multilookL);
		  writearg(keyword);
		  if (Character.isDigit(keyword.charAt(0)))
			subtrrefphainput.argvalue.multilookP = Integer.parseInt(keyword);
		  else // default same factor
			subtrrefphainput.argvalue.multilookP = subtrrefphainput.argvalue.multilookL;
		  }

	// **********************************************************************
	// *** PHASE FILTER
	// **********************************************************************
		else if (!strcmp(keyword,"PF_METHOD")) // method selector phase filtering
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"GOLDSTEIN"))
			filtphaseinput.argvalue.method = fp_goldstein;
		  else if (!strcmp(keyword,"SPATIALCONV"))
			filtphaseinput.argvalue.method = fp_spatialconv;
		  else if (!strcmp(keyword,"SPECTRAL"))
			filtphaseinput.argvalue.method = fp_spectral;
		  else
			{
			ERROR << "PF_METHOD: method " << keyword << " not known for phase filtering. line " << linecnt << ".";
				 PRINT_ERROR(ERROR.get_str())
				 throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"PF_OUT_FILE")) // filename
		  {
		  filtphaseinput.argvalue.fofiltphase = word.charAt(1); // pass keyword
		  writearg(filtphaseinput.argvalue.fofiltphase);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"PF_IN_FILE")) // filename
		  {
		  filtphaseinput.argvalue.fifiltphase = word.charAt(1); // pass keyword
		  String pLast = null;
		  filtphaseinput.argvalue.finumlines = strtoul(word.charAt(1), pLast, BASE10); // pass numoflines
		  if (pLast.equals(word.charAt(2))) // fails to convert one of them to double.
		   {
			ERROR << "PF_IN_FILE (numoflines): " << word.charAt(2) << " are not valid.";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
		   }
		  writearg(filtphaseinput.argvalue.fifiltphase);
		  writearg(filtphaseinput.argvalue.finumlines);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"PF_BLOCKSIZE")) // buffersize
		  {
		  filtphaseinput.argvalue.blocksize = Integer.parseInt(word.charAt(1)); // pass keyword
		  writearg(filtphaseinput.argvalue.blocksize);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"PF_ALPHA")) // alpha
		  {
		  filtphaseinput.argvalue.alpha = Double.parseDouble(word.charAt(1)); // pass keyword
		  writearg(filtphaseinput.argvalue.alpha);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"PF_OVERLAP")) // overlap
		  {
		  filtphaseinput.argvalue.overlap = Integer.parseInt(word.charAt(1)); // pass keyword
		  writearg(filtphaseinput.argvalue.overlap);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"PF_KERNEL")) // conv. kernel e.g. 3 1 1 1
		  {
		  keyword = word.charAt(1); // pass keyword
		  final int32 sizekernel = Integer.parseInt(keyword);
		  writearg(sizekernel);
		  if (!(isodd(sizekernel)))
			{
			 PRINT_ERROR("PF_KERNEL: size must be odd! (add 0, center around midpix)")
			 throw(keyword_error);
			}
		  filtphaseinput.argvalue.kernel.resize(1, sizekernel);
		  real4 sum =0.;
		  for (int32 argnum =0; argnum<sizekernel; ++argnum)
			{
			 keyword = word.charAt(1); // pass keyword
			 if (!(Character.isDigit(keyword.charAt(0))))
			 WARNING.print("kernel seems to be wrong?");
			 final real4 in = Double.parseDouble(keyword);
			 writearg(in);
			 filtphaseinput.argvalue.kernel(0,argnum) = in;
			 sum += Math.abs(in); // kernel -1 1
			}
			if (sum!=1)
			filtphaseinput.argvalue.kernel /= sum; // normalize
			INFO << "PF_KERNEL: Input kernel normalized by: " << sum;
			INFO.print();
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"PF_IN_KERNEL2D")) // filename for ascii 2d kernel
		  {
		  filtphaseinput.argvalue.fikernel2d = word.charAt(1); // pass keyword
		  writearg(filtphaseinput.argvalue.fikernel2d);
		  }

	// *******************************************************************
	// *** DINSAR
	// *******************************************************************
		else if (!strcmp(keyword,"DI_OUT_FILE")) // output filename
		  {
		  dinsarinput.argvalue.fodinsar = word.charAt(1); // pass keyword
		  writearg(dinsarinput.argvalue.fodinsar);
		  }

		else if (!strcmp(keyword,"DI_OUT_SCALED")) // output filename
		  {
		  dinsarinput.argvalue.foscaleduint = word.charAt(1); // pass keyword
		  writearg(dinsarinput.argvalue.foscaleduint);
		  }

		else if (!strcmp(keyword,"DI_IN_TOPOMASTER")) // input resultfilename
		  {
		  dinsarinput.argvalue.topomasterresfile = word.charAt(1); // pass keyword
		  writearg(dinsarinput.argvalue.fodinsar);
		  }

		else if (!strcmp(keyword,"DI_IN_TOPOSLAVE")) // input resultfilename
		  {
		  dinsarinput.argvalue.toposlaveresfile = word.charAt(1); // pass keyword
		  writearg(dinsarinput.argvalue.fodinsar);
		  }

		else if (!strcmp(keyword,"DI_IN_TOPOINT")) // input resultfilename
		  {
		  dinsarinput.argvalue.topointresfile = word.charAt(1); // pass keyword
		  writearg(dinsarinput.argvalue.fodinsar);
		  }


	// **********************************************************************
	// *** COMPUTATION OF REFERENCE DEM (phase)
	// **********************************************************************
	//    else if (!strcmp(keyword,"CRD_METHOD"))          // name of output file
	//      {
	//      filename =  word[1] ;   // pass keyword
	//      writearg(keyword2);
	//      toupper(keyword2);
	//      if (!strcmp(keyword2,"TRILINEAR") || !strcmp(keyword2,"TRI_LINEAR") ||
	//          !strcmp(keyword2,"TRILIN") || !strcmp(keyword2,"TRI_LIN"))
	//        comprefdeminput.method = crd_trilinear;
	//      else if (!strcmp(keyword2,"NN") || !strcmp(keyword2,"NEAREST") ||
	//               !strcmp(keyword2,"NEAREST_NEIGHBOR" ))
	//        comprefdeminput.method = crd_nearest;
	//      else
	//        {
	//        ERROR << "CRD_METHOD: "
	//             <<  filename
	//             << " not known (use TRILINEAR or NEAREST); line "
	//             << linecnt << ".";
	//      PRINT_ERROR(ERROR.get_str())
	//      throw(keyword_error);
	//        }
	//      }

	// **********************************************************************
		else if (!strcmp(keyword,"CRD_IN_DEM")) // input file
		  {
		  comprefdeminput.argvalue.firefdem = word.charAt(1); // pass keyword
		  writearg(comprefdeminput.argvalue.firefdem);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CRD_IN_FORMAT")) // format input file
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"R4") || !strcmp(keyword,"REAL4"))
			comprefdeminput.argvalue.iformatflag = FORMATR4;
		  else if (!strcmp(keyword,"I2") || !strcmp(keyword,"SHORT"))
			comprefdeminput.argvalue.iformatflag = FORMATI2; // default
		  else if (!strcmp(keyword,"I2_BIGENDIAN") || !strcmp(keyword,"SHORT_BIGENDIAN"))
			comprefdeminput.argvalue.iformatflag = FORMATI2_BIGENDIAN; // default
		  else if (!strcmp(keyword,"R8") || !strcmp(keyword,"REAL8"))
			comprefdeminput.argvalue.iformatflag = FORMATR8;
		  else
			{
			ERROR << "CRD_IN_FORMAT: input format " << keyword << " not known (R4 R8 I2 (native) SHORT_BIGENDIAN); line " << linecnt << ".";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CRD_IN_SIZE")) // nrow ncols (lat lon)
		  {
		  comprefdeminput.argvalue.demrows = strtoul(word.charAt(1), null, BASE10); // pass keyword [MA] atoi instead strtoul
		  comprefdeminput.argvalue.demcols = strtoul(word.charAt(2), null, BASE10); // pass keyword
		  writearg(comprefdeminput.argvalue.demrows);
		  writearg(comprefdeminput.argvalue.demcols);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CRD_IN_DELTA")) // degrees delta lat lon
		  {
		  comprefdeminput.argvalue.demdeltalat = Double.parseDouble(word.charAt(1)); // pass keyword
		  keyword = word.charAt(2); // update keyword
		  writearg(comprefdeminput.argvalue.demdeltalat);
		  writearg(keyword);
		  if (Character.isDigit(keyword.charAt(0)) || keyword.charAt(0)=='.') // likely to be 2 numbers
			comprefdeminput.argvalue.demdeltalon = Double.parseDouble(keyword);
		  else // default same gridsize
			comprefdeminput.argvalue.demdeltalon = comprefdeminput.argvalue.demdeltalat;

		  // ______ Store as radians ______
		  comprefdeminput.argvalue.demdeltalat = deg2rad(comprefdeminput.argvalue.demdeltalat);
		  comprefdeminput.argvalue.demdeltalon = deg2rad(comprefdeminput.argvalue.demdeltalon);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CRD_IN_UL")) // upperleft coordinates
		  {
		  String pLast1;
		  String pLast2 = null;
		  comprefdeminput.argvalue.demlatleftupper = strtod(word.charAt(1), pLast1);
					  comprefdeminput.argvalue.demlonleftupper = strtod(word.charAt(2), pLast2);
		  if (pLast1.equals(word.charAt(1)) || pLast2.equals(word.charAt(2))) // fails to convert one of them to double.
		   {
			ERROR << "CRD_IN_UL: " << word.charAt(1) << " : " << word.charAt(2) << " are not valid.";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
		   }
		  writearg(comprefdeminput.argvalue.demlatleftupper);
		  writearg(comprefdeminput.argvalue.demlonleftupper);
		  comprefdeminput.argvalue.demlatleftupper = deg2rad(comprefdeminput.argvalue.demlatleftupper);
		  comprefdeminput.argvalue.demlonleftupper = deg2rad(comprefdeminput.argvalue.demlonleftupper);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CRD_IN_NODATA")) // flag for no data
		  {
		  comprefdeminput.argvalue.demnodata = Double.parseDouble(word.charAt(1)); // pass keyword
		  writearg(comprefdeminput.argvalue.demnodata);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CRD_INCLUDE_FE")) // true: ref.pha incl. flat earth
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"OFF"))
			comprefdeminput.argvalue.includerefpha = false;
		  else if (!strcmp(keyword,"ON") || !strncmp(keyword,"//",2) || !strncmp(keyword,"#",1) || !(keyword.charAt(0) == '\0')) // no keyword
			comprefdeminput.argvalue.includerefpha = true;
		  else
			{
			comprefdeminput.argvalue.includerefpha = false;
			WARNING << "CRD_INCLUDE_FE: line: " << linecnt << ": unknown argument: " << keyword << "; Set to OFF (computing pure topo phase w.r.t. flat earth).";
		   WARNING.print();
		   }
		  }

	// **********************************************************************
	//    else if (!strcmp(keyword,"CRD_DENSE"))            // factor
	//      {
	//      comprefdeminput.extradense =  word[1] ;         // pass keyword
	//      writearg(comprefdeminput.extradense);
	//      }

	// **********************************************************************
		else if (!strcmp(keyword,"CRD_OUT_DEM")) // name of output file
		  {
		  comprefdeminput.argvalue.fodem = word.charAt(1); // pass keyword
		  writearg(comprefdeminput.argvalue.fodem);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CRD_OUT_DEMI")) // name of output file
		  {
		  comprefdeminput.argvalue.fodemi = word.charAt(1); // pass keyword
		  writearg(comprefdeminput.argvalue.fodemi);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"CRD_OUT_FILE")) // name of output file
		  {
		  comprefdeminput.argvalue.forefdem = word.charAt(1); // pass keyword
		  writearg(comprefdeminput.argvalue.forefdem);
		  }

	// ___________ added by FvL      
	// **********************************************************************
		else if (!strcmp(keyword,"CRD_OUT_H2PH")) // name of output file
		  {
		  comprefdeminput.argvalue.foh2ph = word.charAt(1); // pass keyword
		  writearg(comprefdeminput.argvalue.foh2ph);
		  }
	// ___________ end added by FvL

	// **********************************************************************
		else if (!strcmp(keyword,"CRD_OUT_DEM_LP")) // name of output file
		  {
		  comprefdeminput.argvalue.forefdemhei = word.charAt(1); // pass keyword
		  writearg(comprefdeminput.argvalue.forefdemhei);
		  }


	// **********************************************************************
	// *** SUBTRACTION OF REFERENCE DEM (phase)
	// **********************************************************************
		else if (!strcmp(keyword,"SRD_OUT_CINT")) // name of output file
		  {
		  subtrrefdeminput.argvalue.focint = word.charAt(1); // pass keyword
		  writearg(subtrrefdeminput.argvalue.focint);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"SRD_OFFSET")) // Line Pixel
		  {
		  subtrrefdeminput.argvalue.offsetL = Integer.parseInt(word.charAt(1)); // pass keyword
		  subtrrefdeminput.argvalue.offsetP = Integer.parseInt(word.charAt(2)); // pass keyword
		  writearg(subtrrefdeminput.argvalue.offsetL);
		  writearg(subtrrefdeminput.argvalue.offsetP);
		  }

	//  // **********************************************************************
	//  // *** ADDITION OF REFERENCE DEM (phase)
	//  // **********************************************************************
	//      else if (!strcmp(keyword,"ARD_OUT_CINT"))          // name of output file
	//        {
	//        addrefdeminput.focint =  word[1] ;    // pass keyword
	//        }

	// **********************************************************************
	// *** UNWRAPPING
	// **********************************************************************
		else if (!strcmp(keyword,"UW_METHOD")) // method selector unwrapping
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"RAMON"))
			unwrapinput.argvalue.method = uw_method1;
		  else if (!strcmp(keyword,"SNAPHU"))
			unwrapinput.argvalue.method = uw_method2; // default
		  else if (!strcmp(keyword,"MCF_DLR"))
			unwrapinput.argvalue.method = uw_method3;
		  else
			{
			ERROR << "UW_METHOD: method " << keyword << " not known for unwrapping on line " << linecnt << ".";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"UW_SEEDS")) // position of seeds
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  if (Character.isDigit(keyword.charAt(0)))
			{
			 unwrapinput.argvalue.deltaLseed = Integer.parseInt(keyword);
			 keyword = word.charAt(2); // update keyword
			 writearg(keyword);
			 if (Character.isDigit(keyword.charAt(0)))
				unwrapinput.argvalue.deltaPseed = Integer.parseInt(keyword);
			 else
				unwrapinput.argvalue.deltaPseed = unwrapinput.argvalue.deltaLseed;
			}
		  else // assume no numbers but filename with seeds
			{
			 unwrapinput.argvalue.seedfile = keyword;
			}
		   }

	// **********************************************************************
		else if (!strcmp(keyword,"UW_OUT_FILE")) // filename output unwrapped int.
		  {
		  unwrapinput.argvalue.fouint = word.charAt(1); // pass keyword
		  writearg(unwrapinput.argvalue.fouint);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"UW_OUT_FORMAT")) // output format
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"R4") || !strcmp(keyword,"REAL4"))
			unwrapinput.argvalue.oformatflag = FORMATR4;
		  else if (!strcmp(keyword,"CR4") || !strcmp(keyword,"COMPLEXR4"))
			{
			WARNING.print("UW_OUT_FORMAT = CR4 --> Using hgt format");
			unwrapinput.argvalue.oformatflag = FORMATHGT; // default
			}
		  else if (!strcmp(keyword,"HGT"))
			unwrapinput.argvalue.oformatflag = FORMATHGT; // default
		  else
			{
			unwrapinput.argvalue.oformatflag = FORMATHGT; // default
			WARNING << "UW_OUT_FORMAT: output format " << keyword << " not known (R4 or HGT, using HGT); line " << linecnt << ".";
			WARNING.print();
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"UW_OUT_REGIONS")) // filename output regions
		  {
		  unwrapinput.argvalue.foregions = word.charAt(1); // pass keyword
		  writearg(unwrapinput.argvalue.foregions);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"UW_SNAPHU_MODE")) // filename output regions
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"TOPO"))
			unwrapinput.argvalue.snaphu_mode = "TOPO"; // default TOPO
		  else if (!strcmp(keyword,"DEFO"))
			unwrapinput.argvalue.snaphu_mode = "DEFO";
		  else if (!strcmp(keyword,"SMOOTH"))
			unwrapinput.argvalue.snaphu_mode = "SMOOTH";
		  else if (!strcmp(keyword,"NOSTATCOSTS"))
			unwrapinput.argvalue.snaphu_mode = "NOSTATCOSTS";
		  else
			{
			ERROR << "UW_SNAPHU_MODE: " << keyword << " not known for unwrapping on line " << linecnt << ".";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"UW_SNAPHU_INIT"))
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"MST"))
			unwrapinput.argvalue.snaphu_init = "MST"; // default mst
		  else if (!strcmp(keyword,"MCF"))
			unwrapinput.argvalue.snaphu_init = "MCF";
		  else
			{
			WARNING << "UW_SNAPHU_INIT: " << keyword << " not known for unwrapping on line " << linecnt << " (using MST).";
			WARNING.print();
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"UW_SNAPHU_LOG"))
		  {
		  unwrapinput.argvalue.snaphu_log = word.charAt(1); // pass keyword
		  writearg(unwrapinput.argvalue.snaphu_log);
		  //strcpy(unwrapinput.snaphu_log,"-l ");
		  //strcat(unwrapinput.snaphu_log,keyword);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"UW_SNAPHU_COH"))
		  {
		  unwrapinput.argvalue.snaphu_coh = word.charAt(1); // pass keyword
		  writearg(unwrapinput.argvalue.snaphu_coh);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"UW_SNAPHU_VERBOSE"))
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"OFF"))
			unwrapinput.argvalue.snaphu_verbose = "FALSE";
		  else if (!strcmp(keyword,"ON") || !strncmp(keyword,"//",2) || !strncmp(keyword,"#",1) || !(keyword.charAt(0) == '\0')) // no keyword
			unwrapinput.argvalue.snaphu_verbose = "TRUE";
		  else
			{
			unwrapinput.argvalue.snaphu_verbose = "TRUE";
			WARNING << "UW_SNAPHU_VERBOSE: line: " << linecnt << ": unknown argument: " << keyword << "; Set to ON.";
		   WARNING.print();
			}
		  }

		else if (!strcmp(keyword,"UW_SNAPHU_NTILEROW"))
		  {
				 unwrapinput.argvalue.ntilerow = Integer.parseInt(word.charAt(1)); // pass keyword
				 writearg(unwrapinput.argvalue.ntilerow);
		   if (unwrapinput.argvalue.ntilerow > 50)
			 {
			  WARNING.print("UW_SNAPHU_NTILEROW > 100 tiles, is this okay? ");
			 }
		  INFO << "UW_SNAPHU_NTILEROW: \t " << unwrapinput.argvalue.ntilerow;
		  INFO.print();
		  }

		else if (!strcmp(keyword,"UW_SNAPHU_NTILECOL"))
		  {
		   unwrapinput.argvalue.ntilecol = Integer.parseInt(word.charAt(1)); // pass keyword
		   writearg(unwrapinput.argvalue.ntilecol);
		   if (unwrapinput.argvalue.ntilecol > 50)
		   {
			 WARNING.print("UW_SNAPHU_NTILECOL > 100 tiles, is this okay? ");
		   }

		   INFO << "UW_SNAPHU_NTILECOL: \t " << unwrapinput.argvalue.ntilecol;
		   INFO.print();
		  }

		else if (!strcmp(keyword,"UW_SNAPHU_ROWOVRLP"))
		  {
			unwrapinput.argvalue.rowovrlp = Integer.parseInt(word.charAt(1)); // pass keyword
			writearg(unwrapinput.argvalue.rowovrlp);
			INFO << "UW_SNAPHU_ROWOVRLP: \t " << unwrapinput.argvalue.rowovrlp;
			INFO.print();
		  }

		else if (!strcmp(keyword,"UW_SNAPHU_COLOVRLP"))
		  {
			unwrapinput.argvalue.colovrlp = Integer.parseInt(word.charAt(1)); // pass keyword
			writearg(unwrapinput.argvalue.colovrlp);
			INFO << "UW_SNAPHU_COLOVRLP: \t " << unwrapinput.argvalue.colovrlp;
			INFO.print();
		  }

		else if (!strcmp(keyword,"UW_SNAPHU_NPROC"))
		  {
			unwrapinput.argvalue.nproc = Integer.parseInt(word.charAt(1)); // pass keyword
			writearg(unwrapinput.argvalue.nproc);
			if (unwrapinput.argvalue.ntilecol > 2)
			{
			  WARNING.print("UW_SNAPHU_NPROC > 2CPUs, do you have a cluster?");
			}

			INFO << "UW_SNAPHU_NPROC: \t " << unwrapinput.argvalue.ntilecol;
			INFO.print();
		  }

		else if (!strcmp(keyword,"UW_SNAPHU_TILECOSTTHRESH"))
		  {
			unwrapinput.argvalue.tilecostthresh = Integer.parseInt(word.charAt(1)); // pass keyword
			writearg(unwrapinput.argvalue.tilecostthresh);
			if (unwrapinput.argvalue.ntilecol > 500)
			{
			  WARNING.print("UW_SNAPHU_TILECOSTTHRESH > 500, do you have a cluster?");
			}

			INFO << "UW_SNAPHU_TILECOSTTHRESH: \t " << unwrapinput.argvalue.tilecostthresh;
			INFO.print();
		  }

	// **********************************************************************
	// *** SLANT to HEIGHT CONVERSION
	// **********************************************************************
		else if (!strcmp(keyword,"S2H_METHOD")) // method selector slant2height
		  {
		  keyword = word.charAt(1); // pass keyword
		  writearg(keyword);
		  Character.toUpperCase(keyword);
		  if (!strcmp(keyword,"SCHWABISCH"))
			slant2hinput.argvalue.method = s2h_schwabisch;
		  else if (!strcmp(keyword,"AMBIGUITY"))
			slant2hinput.argvalue.method = s2h_ambiguity;
		  else if (!strcmp(keyword,"RODRIGUEZ"))
			slant2hinput.argvalue.method = s2h_rodriguez;
		  else
			{
			ERROR << "S2H_METHOD: method " << keyword << " not known for slant to height conversion, line " << linecnt << ".";
			PRINT_ERROR(ERROR.get_str())
			throw(keyword_error);
			}
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S2H_NPOINTS")) // number of points to use
		  {
		  slant2hinput.argvalue.Npoints = Integer.parseInt(word.charAt(1)); // pass keyword
		  writearg(slant2hinput.argvalue.Npoints);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S2H_DEGREE1D")) // degree of 1d polynomial
		  {
		  slant2hinput.argvalue.degree1d = Integer.parseInt(word.charAt(1)); // pass keyword
		  writearg(slant2hinput.argvalue.degree1d);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S2H_DEGREE2D")) // degree of 2d polynomial
		  {
		  slant2hinput.argvalue.degree2d = Integer.parseInt(word.charAt(1)); // pass keyword
		  writearg(slant2hinput.argvalue.degree2d);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S2H_NHEIGHTS")) // #heights to evaluate ref.pha
		  {
		  slant2hinput.argvalue.Nheights = Integer.parseInt(word.charAt(1)); // pass keyword
		  writearg(slant2hinput.argvalue.Nheights);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S2H_OUT_HEI")) // filename output height
		  {
		  slant2hinput.argvalue.fohei = word.charAt(1); // pass keyword
		  writearg(slant2hinput.argvalue.fohei);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S2H_OUT_PHI")) // filename output latitude
		  {
		  slant2hinput.argvalue.fophi = word.charAt(1); // pass keyword
		  writearg(slant2hinput.argvalue.fophi);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"S2H_OUT_LAM")) // filename output longitude
		  {
		  slant2hinput.argvalue.folam = word.charAt(1); // pass keyword
		  writearg(slant2hinput.argvalue.folam);
		  }

	// **********************************************************************
	// *** GEOCODING
	// **********************************************************************
	//      else if (!strcmp(keyword,"GEO_METHOD"))
	//        {
	//        geocodeinput.method =  word[1] ;      // pass keyword
	//        }

	// **********************************************************************
		else if (!strcmp(keyword,"GEO_OUT_PHI")) // output file latitude
		  {
		  geocodeinput.argvalue.fophi = word.charAt(1); // pass keyword
		  writearg(geocodeinput.argvalue.fophi);
		  }

	// **********************************************************************
		else if (!strcmp(keyword,"GEO_OUT_LAM")) // output file longitude
		  {
		  geocodeinput.argvalue.folam = word.charAt(1); // pass keyword
		  writearg(geocodeinput.argvalue.folam);
		  }


	// **********************************************************************
	// Assume wrong keyword, but if it starts with //, a c, or comment
	// then continue with a warning, no blank between key and arg
	// **********************************************************************
		else if (!strncmp(keyword,"COMMENT",7) || !strncmp(keyword,"C",1)) // but forgot a space
		  {
		  WARNING << "Obsolete or unknown keyword: \"" << keyword << "\" at line: " << linecnt << ". (Interpreted as comment, ignored)";
		  WARNING.print();
		  }
		// ______ Really cannot make anything from your input ______
		else
		  {
		  ERROR << "Unknown keyword: \"" << keyword << "\" at line: " << linecnt << ".";
		  PRINT_ERROR(ERROR.get_str())
		  throw(keyword_error);
		  }
		//optionsfile.getline(eachline,4*ONE27,'\n');               // goto next line. [MA] needs debug if the line is long loops forever.
		} // end while (file)






	// ====== Check input/ info to screen ======
	  checkgeneral(generalinput, onlyprocess); // also repair onlyprocess
	  INFO << "LISTINPUT: \tAppend input to logfile: \t" << listinput;
	  INFO.print();


	// ______ info ELLIPSOID card, compute and fill e2, e2b ______
	  ellipsinput.argvalue.showdata();


	// ====== Check input of step reading info from files ======
	  if (generalinput.argvalue.process[pr_m_readfiles])
		checkreadfiles(m_readfilesinput.argvalue, MASTERID); // check mandatory cards

	  if (generalinput.argvalue.process[pr_s_readfiles])
		checkreadfiles(s_readfilesinput.argvalue, SLAVEID); // check mandatory cards


	// ====== Check input of step conversion slc to raw ======
	  if (generalinput.argvalue.process[pr_m_crop])
		{
		checkcrop(m_cropinput.argvalue, MASTERID);
		if (!strcmp(s_cropinput.argvalue.fileout1,m_cropinput.argvalue.fileout1))
		  {
		  PRINT_ERROR("code 301: same name outputfile CROP_OUT for master and slave not allowed.")
		  throw(keyword_error);
		  }
		if (generalinput.argvalue.overwrit)
		  if (existed(m_cropinput.argvalue.fileout1))
			{
			INFO << "OVERWRIT: file " << m_cropinput.argvalue.fileout1 << " will be overwritten.";
			INFO.print();
			}
		}

	  if (generalinput.argvalue.process[pr_s_crop])
		{
		checkcrop(s_cropinput.argvalue, SLAVEID);
		if (!strcmp(s_cropinput.argvalue.fileout1,m_cropinput.argvalue.fileout1))
		  {
		  PRINT_ERROR("code 301: same name outputfile CROP_OUT for master and slave not allowed.")
		  throw(keyword_error);
		  }
		if (generalinput.argvalue.overwrit)
		  if (existed(s_cropinput.argvalue.fileout1))
			{
			INFO << "OVERWRIT: file " << s_cropinput.argvalue.fileout1 << " will be overwritten.";
			INFO.print();
			}
		}

	//____RaffaeleNutricato START MODIFICATION SECTION 10
	// ====== Check input of step oversample master======
	  if (generalinput.argvalue.process[pr_m_oversample])
		{
		checkoversample(m_oversample.argvalue, MASTERID);
		}
	// ====== Check input of step oversampling slave ======
	  if (generalinput.argvalue.process[pr_s_oversample])
		{
		checkoversample(s_oversample.argvalue, SLAVEID);
		}
	//____RaffaeleNutricato END MODIFICATION SECTION 10


	// ====== Check input of step porbits ======
	  if (generalinput.argvalue.process[pr_m_porbits])
		checkporbits(porbitsinput.argvalue, MASTERID);
	  if (generalinput.argvalue.process[pr_s_porbits])
		checkporbits(porbitsinput.argvalue, SLAVEID);

	// ====== Check input of SIMAMP step master====== [MA]
	if (generalinput.argvalue.process[pr_m_simamp])
	  {
	  checksimamp(simampinput.argvalue);
	  }

	// ====== Check input of MTIMING step master====== [MA]
	if (generalinput.argvalue.process[pr_m_mtiming])
	  {
	  checkmtiming(mtiminginput.argvalue);
	  }

	// ====== Check input of step azimuth filtering ======
	  if (generalinput.argvalue.process[pr_m_filtazi] || generalinput.argvalue.process[pr_s_filtazi])
		{
		// ______ Set defaults ______
		if (filtaziinput.argvalue.overlap == -1)
		  filtaziinput.argvalue.overlap = filtaziinput.argvalue.fftlength/8; // default
		// ______ Check input ______
		if (generalinput.argvalue.process[pr_m_filtazi] && generalinput.argvalue.process[pr_s_filtazi])
		  checkfiltazi(filtaziinput.argvalue, MASTERID+SLAVEID);
		else if (generalinput.argvalue.process[pr_m_filtazi])
		  checkfiltazi(filtaziinput.argvalue, MASTERID);
		else if (generalinput.argvalue.process[pr_s_filtazi])
		  checkfiltazi(filtaziinput.argvalue, SLAVEID);
		else
		  {
		  PRINT_ERROR("PANIC, this cannot be")
		  throw(keyword_error);
		  }
		}

	// ====== Check input of step range filtering ======
	  if (generalinput.argvalue.process[pr_m_filtrange])
		{
		// ______ Check defaults ______
		if (filtrangeinput.argvalue.fftlength==-999) // not specified
		  {
		  if (filtrangeinput.argvalue.method==rf_adaptive)
			filtrangeinput.argvalue.fftlength=64; // default
		  else if (filtrangeinput.argvalue.method==rf_porbits)
			filtrangeinput.argvalue.fftlength=1024; // default
		  }
		checkfiltrange(filtrangeinput.argvalue);
		}

	// ====== Check input of reserved EXTRA step master/slave ======
	  if (generalinput.argvalue.process[pr_m_EXTRA])
		{
		PRINT_ERROR("extra step master not implemented.")
		throw(keyword_error);
		}
	  if (generalinput.argvalue.process[pr_s_EXTRA])
		{
		PRINT_ERROR("extra step slave not implemented.")
		throw(keyword_error);
		}


	// ====== Check input of step coarse (orbits) ======
	// ______ no checks required. (no cards as well)
	  if (generalinput.argvalue.process[pr_i_coarse])
		{
		; // do nothing
		}


	// ====== Check coarse coregistration based on correlation ======
	// ______ Check + repair method selector coarse correlation ______
	  if (generalinput.argvalue.process[pr_i_coarse2]) // correlation
		{
		if (coarsecorrinput.argvalue.method == def_cc_method-999) // see at defaults, no method card used
		  {
		  coarsecorrinput.argvalue.method = def_cc_method; // default method
		  INFO.print("Default method will be used for coarse coregistration.");
		  }
		if (specified(coarsecorrinput.argvalue.ifpositions)) // filename specified
		  {
		  if (coarsecorrinput.argvalue.Nwin != def_cc_nwin+999) // something is specified
			{
			WARNING << "CC_NWIN: \t" << coarsecorrinput.argvalue.Nwin << " ignored due to existence of input file (CC_IN_POS) " << coarsecorrinput.argvalue.ifpositions;
			WARNING.print();
			}
		  coarsecorrinput.argvalue.Nwin = filelines(coarsecorrinput.argvalue.ifpositions);
		  }
		else if (coarsecorrinput.argvalue.Nwin == def_cc_nwin+999) // no inputfile, default
		  {
		  coarsecorrinput.argvalue.Nwin = def_cc_nwin;
		  INFO.print("Default number of windows will be used for coarse coregistration.");
		  }
		checkcoarsecorr(coarsecorrinput.argvalue);
		}


	// ====== Check input of step fine ======
	  if (generalinput.argvalue.process[pr_i_fine])
		{
		if (fineinput.argvalue.method == def_fc_method-999) // see at defaults, no method card used
		  {
		  fineinput.argvalue.method = def_fc_method; // default method
		  INFO.print("Default method will be used for fine coregistration.");
		  }
		if (specified(fineinput.argvalue.ifpositions)) // filename specified
		  {
		  if (fineinput.argvalue.Nwin != def_fc_nwin+999) // something is specified
			{
			WARNING << "FC_NWIN: \t" << fineinput.argvalue.Nwin << " ignored due to existence of input file (FC_IN_POS) " << fineinput.argvalue.ifpositions;
			WARNING.print();
			}
		  fineinput.argvalue.Nwin = filelines(fineinput.argvalue.ifpositions);
		  }
		else if (fineinput.argvalue.Nwin == def_fc_nwin+999) // no inputfile, default
		  {
		  fineinput.argvalue.Nwin = def_fc_nwin;
		  INFO.print("Default number of windows will be used for fine coregistration.");
		  }
		checkfine(fineinput.argvalue);
		}

	  // ====== Check input of TIMING step interferogram ====== [FvL]
	  if (generalinput.argvalue.process[pr_i_timing])
		{
		checkreltiming(reltiminginput.argvalue);
		}

	  // ====== Check input of DEMASSIST step interferogram ====== [FvL]
	  if (generalinput.argvalue.process[pr_i_demassist])
		{
		checkdemassist(demassistinput.argvalue);
		}

	  // ====== Check input step coregpm ======
	  if (generalinput.argvalue.process[pr_i_coregpm])
		{
		checkcoregpm(coregpminput.argvalue);
		}

	  // ====== Check + repair method selector resampling ======
	  if (generalinput.argvalue.process[pr_s_resample]) // request for process resample
		{
		if (resampleinput.argvalue.method == def_rs_method-999) // see at defaults, no method card used
		  {
		  resampleinput.argvalue.method = def_rs_method; // default method
		  INFO.print("RS_METHOD: Using default.");
		  }
		if (!priorrs_fileout)
		  INFO.print("RS_OUT_FILE: Using default.");
		checkresample(resampleinput.argvalue);
		}

	  // ====== Check + repair method selector flatearth correction ======
	  if (generalinput.argvalue.process[pr_i_comprefpha])
		{
		if (comprefphainput.argvalue.method == def_fe_method-999) // see at defaults, no method card used
		  {
		  comprefphainput.argvalue.method = def_fe_method; // default method
		  INFO.print("FE_METHOD: Using default.");
		  }
		if (comprefphainput.argvalue.Npoints == def_fe_Npoints-999) // default applies
		  {
		  comprefphainput.argvalue.Npoints = def_fe_Npoints; // default value
		  // INFO.print("FE_NPOINTS: Using default.");
		  }
		// ______ if read from file, count number of points ______
		if (specified(comprefphainput.argvalue.ifpositions)) // filename specified
		  {
		  INFO.print("FE_IN_POS: Using file to read positions.");
		  comprefphainput.argvalue.Npoints = filelines(comprefphainput.argvalue.ifpositions);
		  }
		if (comprefphainput.argvalue.degree == def_fe_degree-999) // default
		  {
		  comprefphainput.argvalue.degree = def_fe_degree; // default
		  INFO.print("FE_DEGREE: Using default.");
		  }
		checkcomprefpha(comprefphainput.argvalue);
		}

	  // ====== Check input of step interfero ======
	  if (generalinput.argvalue.process[pr_i_interfero])
		{
		checkinterfero(interferoinput.argvalue);
		}

	  // ====== Check input of SUBTRREFPHA step interferogram ======
	  if (generalinput.argvalue.process[pr_i_subtrrefpha])
		checksubtrrefpha(subtrrefphainput.argvalue);

	  // ====== Check input of COMPREFDEM step interferogram ======
	  if (generalinput.argvalue.process[pr_i_comprefdem])
		{
		checkcomprefdem(comprefdeminput.argvalue);
		}

	  // ====== Check input of SUBTRDEM step interferogram ======
	  if (generalinput.argvalue.process[pr_i_subtrrefdem])
		{
		checksubtrrefdem(subtrrefdeminput.argvalue);
		}

	  // ====== Check input of step coherence ======
	  if (generalinput.argvalue.process[pr_i_coherence])
		{
		checkcoherence(coherenceinput.argvalue);
		}

	  // ====== Check input of step phase filtering ======
	  if (generalinput.argvalue.process[pr_i_filtphase])
		{
		if (!specified(filtphaseinput.argvalue.fofiltphase)) // not specified, use default
		  {
		  if (filtphaseinput.argvalue.method==fp_goldstein)
			{
			String dummy127 = new String(new char[ONE27]);
			ostrstream omemfo = new ostrstream(dummy127,ONE27);
			omemfo << "cint." << filtphaseinput.argvalue.alpha << "filtered" << ends;
			filtphaseinput.argvalue.fofiltphase = dummy127;
			if (filtphaseinput.argvalue.kernel.size()==0)
			  {
			  INFO.print("PF_KERNEL: Using default kernel [1 2 3 2 1].");
			  filtphaseinput.argvalue.kernel.resize(1, 5);
			  filtphaseinput.argvalue.kernel(0,0) = 1./9.;
			  filtphaseinput.argvalue.kernel(0,1) = 2./9.;
			  filtphaseinput.argvalue.kernel(0,2) = 3./9.;
			  filtphaseinput.argvalue.kernel(0,3) = 2./9.;
			  filtphaseinput.argvalue.kernel(0,4) = 1./9.;
			  }
			}
		  else if (filtphaseinput.argvalue.method==fp_spatialconv)
			{
			filtphaseinput.argvalue.fofiltphase = "cint.filtered";
			// ______ Set default kernel ______
			if (!specified(filtphaseinput.argvalue.fikernel2d)) // no input file
			  {
			  if (filtphaseinput.argvalue.kernel.size()==0)
				{
				INFO.print("PF_KERNEL: Using default kernel [1 1 1].");
				filtphaseinput.argvalue.kernel.resize(1, 3);
				filtphaseinput.argvalue.kernel(0,0) = 1./3.;
				filtphaseinput.argvalue.kernel(0,1) = 1./3.;
				filtphaseinput.argvalue.kernel(0,2) = 1./3.;
				}
			  }
			}
		  else if (filtphaseinput.argvalue.method==fp_spectral)
			{
			filtphaseinput.argvalue.fofiltphase = "cint.filtered";
			}
		  else
			{
			PRINT_ERROR("Method filtphase not known")
			throw(keyword_error);
			}
		  }
		checkfiltphase(filtphaseinput.argvalue);
		}


	  // ====== Check input of step unwrapping ======
	  if (generalinput.argvalue.process[pr_i_unwrap]) // request for process unwrapping
		{
		if (unwrapinput.argvalue.method == def_uw_method-999) // see at defaults, no method card used
		  {
		  unwrapinput.argvalue.method = def_uw_method; // default method
		  INFO.print("Default method will be used for unwrapping.");
		  }
		checkunwrap(unwrapinput.argvalue);
		}


	  // ====== Check input of step slant2height ======
	  if (generalinput.argvalue.process[pr_i_slant2h])
		{
		if (slant2hinput.argvalue.method == def_s2h_method-999)
		  {
		  slant2hinput.argvalue.method = def_s2h_method; // default method
		  INFO.print("Default method will be used for slant2h.");
		  }
		if (slant2hinput.argvalue.Nheights <= slant2hinput.argvalue.degree1d)
		  {
		  WARNING << "S2H_NHEIGHTS: \t" << slant2hinput.argvalue.Nheights << " is too large because S2H_DEGREE1D=" << slant2hinput.argvalue.degree1d << "; set S2H_NHEIGHTS = " << slant2hinput.argvalue.degree1d+1 << " (minimum)";
		  WARNING.print();
		  slant2hinput.argvalue.Nheights = slant2hinput.argvalue.degree1d + 1;
		  }
		checkslant2h(slant2hinput.argvalue);
		}


	  // ====== Check input of step geocoding ======
	  if (generalinput.argvalue.process[pr_i_geocoding])
		{
		checkgeocode(geocodeinput.argvalue);
		}

	  // ====== Check input of dinsar step interferogram ======
	  if (generalinput.argvalue.process[pr_i_dinsar])
		{
		if (!strcmp(dinsarinput.argvalue.topomasterresfile,generalinput.argvalue.m_resfile))
		  setunspecified(dinsarinput.argvalue.topomasterresfile); // used in check...
		checkdinsar(dinsarinput.argvalue);
		if (!specified(dinsarinput.argvalue.topomasterresfile))
		  dinsarinput.argvalue.topomasterresfile = generalinput.argvalue.m_resfile;
		}

	  // ====== Check input of reserved EXTRA step interferogram ======
	  if (generalinput.argvalue.process[pr_i_EXTRA2])
		{
		PRINT_ERROR("extra step2 interferogram not implemented.")
		throw(keyword_error);
		}



	// ====== Copy input to logfile (via scratchfile, update in main) ======
	//   to avoid problems with not existing file it is always opened here
	  ofstream scratchreadinput = new ofstream("scratchreadinput", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchreadinput, "readinput: scratchreadinput", __FILE__, __LINE__);

	  if (listinput)
		{
		// ______Start writing______
		linecnt =0;
		scratchreadinput << "\n----------------------------------------------------\n" << "      BEGIN: copy of input (non-interpretive).\n" << "-line---keyword----argument-------------comment-----\n";
		optionsfile.seekg(0,ios.beg);
		optionsfile.clear(); // clear eofbit
		optionsfile.getline(eachline,4 *ONE27,'\n');
		while (optionsfile != null)
		  {
		  linecnt++;
		  scratchreadinput << setw(3) << linecnt << ": " << eachline << "\n";
		  optionsfile.getline(eachline,4 *ONE27,'\n');
		  }
		scratchreadinput << "\n----------------------------------------------------\n" << "      END: copy of input.\n" << "----------------------------------------------------\n\n";
		DEBUG.print("Finished writing to scratchreadinput.");
		}

	// ______Tidy up______
	  PROGRESS.print("Interpretation inputoptionsfile finished.");
	  scratchreadinput.close();
	  optionsfile.close();
	  } // END READINPUT
}



//----------------------------------------------------------------------------------------
//	Copyright © 2006 - 2008 Tangible Software Solutions Inc.
//
//	This class provides miscellaneous helper methods for strings.
//----------------------------------------------------------------------------------------
final class StringHelper
{
	//------------------------------------------------------------------------------------
	//	This method allows replacing a single character in a string, to help convert
	//	C++ code where a single character in a character array is replaced.
	//------------------------------------------------------------------------------------
	static String changeCharacter(String sourcestring, int charindex, char changechar)
	{
		return (charindex > 0 ? sourcestring.substring(0, charindex) : "")
			+ Character.toString(changechar) + (charindex < sourcestring.length() - 1 ? sourcestring.substring(charindex + 1) : "");
	}
}
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
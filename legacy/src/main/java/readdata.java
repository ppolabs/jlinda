public class GlobalMembersReaddata
{
//C++ TO JAVA CONVERTER TODO TASK: The implementation of the following method could not be found:
	//String strptime(String s, String format, RefObject<tm> tm);



	//***************************************************************
	// *    julday()                                                  *
	// * mm JAN=1, FEB=2, etc; id: day=1,2,3.., iyyy=1996             *
	// *                                                              *
	// #%// Bert Kampes, 07-Apr-2005
	// ***************************************************************
	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//#define IGREG (15+31L*(10+12L*1582))
	public static int32 julday(int32 id, int32 mm, int32 iyyy)
	  {
	  int32 jul;
	  int32 ja;
	  int32 jy;
	  int32 jm;
	  if (iyyy ==0)
		{
		//PRINT_ERROR("julday: error")
		//throw(some_error);
		WARNING.print("julday: year=0; impossible; continuing (only for Btemp computation)");
		return 0;
		}
	  if (iyyy<0)
		  ++iyyy;
	  if (mm>2)
		{
		jy =iyyy;
		jm =mm+1;
		}
	  else
		{
		jy =iyyy-1;
		jm =mm+13;
		}
	  jul = int32(Math.floor(365.25 *jy)+Math.floor(30.6001 *jm)+id+1720995);
	  if (id+31L*(mm+12L *iyyy) >= (15+31L*(10+12L *1582)))
		{
		ja = int32(0.01 *jy);
		jul += 2-ja+int32(0.25 *ja);
		}
	  return jul;
	  }
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#undef IGREG



	//***************************************************************
	// *    readvolume                                                *
	// *                                                              *
	// * reads volumefile                                             *
	// *  and writes to scratchlogfile and resfile                    *
	// *  checks if master is same as slave by id                     *
	// * see: annex C ERS SAR.SLC CCTand EXABYTE                      *
	// *      doc:er-is-epo-gs-5902.3                                 *
	// http://earth.esa.int/rootcollection/sysutil/01008.html
	// *                                                              *
	// * input:                                                       *
	// *  - struct: input options for readfiles                       *
	// *  - 3 checks with other (slave or master) volumefile          *
	// * output:                                                      *
	// *  - void, scratchfiles: ?                                     *
	// *                                                              *
	// *    Bert Kampes, 11-Dec-1998                                  *
	// * if radarsat, skip record number 4 of 360 bytes.              *
	// * according to specs                                           *
	// #%// Bert Kampes, 02-Aug-2004                                  *
	// #%// Davide Nitti (Don), 11-Nov-2008 Reader update for ALOS    *
	// ***************************************************************
	public static void readvolume(RefObject<input_readfiles> readfiles_arg, String checkvol1, String checkvol2, String checkvol3)
	  {
	  TRACE_FUNCTION("readvolume (BK 11-Dec-1998)")
	  final int16 sizeb1 = 1;
	  final int16 sizeb4 = 4;
	  final int16 sizei4 = 4;
	  final int16 sizei8 = 8;
	  final int16 sizea8 = 8;
	  final int16 sizea12 = 12;
	  final int16 sizea16 = 16;
	  final int16 sizea28 = 28;
	  final int16 sizea40 = 40;
	  final int16 sizea60 = 60;
	  uint lenrec1; // volfile...
	  uint lenrec2;
	  uint lenrec3;
	  uint lenrec4;
	  uint numpointrec;
	  uint numrec;
	  uint numvol;
	  String c4dummy = new String(new char[5]);
	  String c8date = new String(new char[9]);
	  String c8time = new String(new char[9]);
	  String c8agency = new String(new char[9]);
	  String c8nlins = new String(new char[9]);
	  String c12logvol = new String(new char[13]);
	  String c12country = new String(new char[13]);
	  String c12facility = new String(new char[13]);
	  String c16physid = new String(new char[17]);
	  String c16logvolid = new String(new char[17]);
	  String c16setid = new String(new char[17]);
	  String c16checkfilename = new String(new char[17]);
	  String c16dataref = new String(new char[17]);
	  String c28leaderrefclass = new String(new char[29]);
	  String c28datarefclass = new String(new char[29]);
	  String c40typespec = new String(new char[41]);
	  String c40physvolid = new String(new char[41]);
	  String c40sceneid = new String(new char[41]);
	  String c40sceneloc = new String(new char[41]);
	  String c60product = new String(new char[61]);
	  // --- Check for Atlantis processor (RSAT?) #%// Bert Kampes, 02-Aug-2004 ---
	  uint rec_seq; // type B4
	  byte rec_sub1; // type B1
	  byte rec_type;
	  byte rec_sub2;
	  byte rec_sub3;

	// ======Open files====== 
	  // ___ check if opened correctly, if not, try to use uppercase
	  // ___ from SCENE1/lea_01.001 #%// BK 27-Nov-2003
	  ifstream volumefile;
	  RefObject<ifstream> TempRefObject = new RefObject<ifstream>(volumefile);
	  openfstream(TempRefObject, readfiles_arg.argvalue.volfile);
	  volumefile = TempRefObject.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(volumefile, readfiles_arg.argvalue.volfile, __FILE__, __LINE__);

	// ______First check control arguments, allow check to be wrong ______
	  DEBUG.print("readvol: check on filename consistency.");
	  for (register int32 i =0; i<999; i++) // try some times
		{
		volumefile.seekg(44,ios.beg);
		volumefile.read((char)&c16physid, sizea16); // physical logical volume ID
		c16physid.charAt(16)='\0';
		volumefile.read((char)&c16logvolid, sizea16); // logical volume ID
		c16logvolid.charAt(16)='\0';
		volumefile.read((char)&c16setid, sizea16); // volume set ID
		c16setid.charAt(16)='\0';

		if (readfiles_arg.argvalue.sensor_id==SLC_ALOS) // added by don
		  {
		  String c8logvoltime = new String(new char[9]);
		  volumefile.seekg(120,ios.beg);
		  volumefile.read((char)&c8logvoltime, sizea8); // Logical volume creation time
		  c8logvoltime.charAt(8)='\0';
		  for (register int32 idx =0; idx<8; idx++)
			c16physid.charAt(idx+8)=c8logvoltime.charAt(idx);
		  }

		if ((!strcmp(checkvol1,c16physid)) && (!strcmp(checkvol2,c16logvolid)) && (!strcmp(checkvol3,c16setid)))
		  {
		  if (i ==0) // first time
			{
			WARNING.print("Volume file of master and slave seem to be the same, please change cd.");
			getanswer(); // wait
			}

		  else if (i ==1) // second time
			{
			WARNING << "ID of volume file 1: " << c16physid << ": " << c16logvolid << ": " << c16setid;
			WARNING.print();
			WARNING << "ID of volume file 2: " << checkvol1 << ":" << checkvol2 << ":" << checkvol3;
			WARNING.print();
			WARNING.print("Next time I will assume they are different.");
			getanswer(); // wait
			}

		  else // if i == 2, enough tries
			{
			break;
			}
		  }
		else // check passed
		  {
		  break;
		  }
		}


	  // ======Read volumefile======
	  // --- RECORD 1 ---
	  DEBUG.print("record 1 of volume file.");
	  volumefile.seekg(0,ios.beg); // place pointer at beginning of file
	  volumefile.read((char) rec_seq, sizeb4); // record number
	  rec_seq = ntohl(rec_seq); // Bert Kampes, 07-Apr-2005
	  volumefile.read((char)&rec_sub1, sizeb1); // first record sub type code
	  volumefile.read((char)&rec_type, sizeb1); // record type code
	  volumefile.read((char)&rec_sub2, sizeb1); // second record sub type code
	  volumefile.read((char)&rec_sub3, sizeb1); // third record sub type code
	  DEBUG.print("Expecting record 1 with code {192,192,18,18}");
	  DEBUG << "rec_seq: " << rec_seq << "; rec_sub1: " << (int)rec_sub1 << "; rec_type: " << (int)rec_type << "; rec_sub2: " << (int)rec_sub2 << "; rec_sub3: " << (int)rec_sub3;
	  DEBUG.print();

	  // --- Read important fields ---
	  volumefile.seekg(8,ios.beg);
	  volumefile.read((char) lenrec1, sizeb4); // length of record1
	  lenrec1 = ntohl(lenrec1); // bk 6 jul 2000, byteorder x86 machines.
	  if (lenrec1 != 360)
		{
		WARNING.print("probably something wrong here, byte order x86 (?).");
		WARNING << "readvolume: length of record 1 = \"" << lenrec1 << "\"; expected \"360\" for ESA SLC (full scene).";
		WARNING.print();
		}
	  volumefile.seekg(32,ios.beg);
	  volumefile.read((char)&c12logvol, sizea12); // logical volume etc.
	  c12logvol.charAt(12)='\0';
	  volumefile.seekg(112,ios.beg);
	  volumefile.read((char)&c8date, sizea8); // generating date
	  c8date.charAt(8)='\0';
	  volumefile.read((char)&c8time, sizea8); // generating time
	  c8time.charAt(8)='\0';
	  volumefile.read((char)&c12country, sizea12); // generating country
	  c12country.charAt(12)='\0';
	  volumefile.read((char)&c8agency, sizea8); // generating agency
	  c8agency.charAt(8)='\0';
	  volumefile.read((char)&c12facility, sizea12); // generating facility
	  c12facility.charAt(12)='\0';
	  volumefile.read((char)&c4dummy, sizei4); // #pointer records in vol.
	  c4dummy.charAt(4)='\0';
	  numpointrec = Integer.parseInt(c4dummy);

	  // if (readfiles_arg.sensor_id!=SLC_RSAT)
	  // Modified by LG for reading ALOS Fine
	  if ((readfiles_arg.argvalue.sensor_id!=SLC_RSAT) && (readfiles_arg.argvalue.sensor_id!=SLC_ALOS))
		{
		if (numpointrec!=2)
		  {
		  ERROR << "readvolume: number of pointer records = \"" << c4dummy << "\"; expected \"2\" for ESA SLC (full scene)." << ends;
		  WARNING.print(ERROR.get_str()); // Jan Kianicka report atlantis processor
		  WARNING.print("but for Atlantis processor this may be correct?");
		  ERROR.reset();
		  }
		}
	  else
		if (numpointrec!=3)
		  WARNING.print("expected 3 pointer records for Radarsat?");

	  volumefile.read((char)&c4dummy, sizei4); // #records in vol.
	  c4dummy.charAt(4)='\0';
	  numrec = Integer.parseInt(c4dummy);
	  // Modified by LG for reading ALOS Fine
	  if (readfiles_arg.argvalue.sensor_id == SLC_ALOS)
	  {
			  if (numrec!=1)
		  {
		  ERROR << "readvolume: number of records = \"" << c4dummy << "\"; expected \"1\" for JAXA SLC (full scene)." << ends;
		  WARNING.print(ERROR.get_str());
		  ERROR.reset();
		  }
	  }
	  else if (readfiles_arg.argvalue.sensor_id!=SLC_RSAT)
		{
		if (numrec!=4)
		  {
		  ERROR << "readvolume: number of records = \"" << c4dummy << "\"; expected \"4\" for ESA SLC (full scene)." << ends;
		  WARNING.print(ERROR.get_str());
		  ERROR.reset();
		  }
		}
	  else
		if (numrec!=5)
		  WARNING.print("expected 5 records for Radarsat?");
	  volumefile.read((char)&c4dummy, sizei4); // #vol.
	  c4dummy.charAt(4)='\0';
	  numvol = Integer.parseInt(c4dummy);
	  if (readfiles_arg.argvalue.sensor_id!=SLC_RSAT)
		if (numvol!=1)
		  {
		  ERROR << "readvolume: number of volumes = \"" << c4dummy << "\"; expected \"1\" for ESA SLC (full scene)." << ends;
		  WARNING.print(ERROR.get_str());
		  ERROR.reset();
		  }


	  // --- RECORD 2 ---
	  final uint startrec2 = lenrec1;
	  volumefile.seekg(startrec2,ios.beg);
	  volumefile.read((char) rec_seq, sizeb4); // record number
	  rec_seq = ntohl(rec_seq); // Bert Kampes, 07-Apr-2005
	  volumefile.read((char)&rec_sub1, sizeb1); // first record sub type code
	  volumefile.read((char)&rec_type, sizeb1); // record type code
	  volumefile.read((char)&rec_sub2, sizeb1); // second record sub type code
	  volumefile.read((char)&rec_sub3, sizeb1); // third record sub type code
	  DEBUG.print("Expecting record 2 with code {219,192,18,18}");
	  DEBUG << "rec_seq: " << rec_seq << "; rec_sub1: " << (int)rec_sub1 << "; rec_type: " << (int)rec_type << "; rec_sub2: " << (int)rec_sub2 << "; rec_sub3: " << (int)rec_sub3;
	  DEBUG.print();

	  // --- Read important fields ---
	  volumefile.seekg(startrec2+8,ios.beg); // leader file pointer rec.
	  DEBUG.print("record 2 of volume file.");
	  volumefile.read((char) lenrec2, sizeb4); // length of record2
	  lenrec2 = ntohl(lenrec2); // bk 6 jul 2000, byteorder x86 machines.
	  if (lenrec2 != 360)
		{
		WARNING.print("probably something wrong here, byte order x86. (?)");
		WARNING << "readvolume: length of record 2 = \"" << lenrec2 << "\"; expected \"360\" for ESA SLC (full scene).";
		WARNING.print();
		}
	  volumefile.seekg(startrec2+20,ios.beg);
	  volumefile.read((char)&c16checkfilename, sizea16); // referenced file name
	  c16checkfilename.charAt(16)='\0';
	  volumefile.read((char)&c28leaderrefclass, sizea28); // referenced file class
	  c28leaderrefclass.charAt(28)='\0';


	  // --- RECORD 3 ---
	  DEBUG.print("record 3 of volume file.");
	  final uint startrec3 = lenrec1 + lenrec2;
	  volumefile.seekg(startrec3,ios.beg);
	  volumefile.read((char) rec_seq, sizeb4); // record number
	  rec_seq = ntohl(rec_seq); // Bert Kampes, 07-Apr-2005
	  volumefile.read((char)&rec_sub1, sizeb1); // first record sub type code
	  volumefile.read((char)&rec_type, sizeb1); // record type code
	  volumefile.read((char)&rec_sub2, sizeb1); // second record sub type code
	  volumefile.read((char)&rec_sub3, sizeb1); // third record sub type code
	  DEBUG.print("Expecting record 3 with code {219,192,18,18}");
	  DEBUG << "rec_seq: " << rec_seq << "; rec_sub1: " << (int)rec_sub1 << "; rec_type: " << (int)rec_type << "; rec_sub2: " << (int)rec_sub2 << "; rec_sub3: " << (int)rec_sub3;
	  DEBUG.print();

	  // --- Read important fields ---
	  volumefile.seekg(startrec3+8,ios.beg); // data file pointer rec.
	  volumefile.read((char) lenrec3, sizeb4); // length of record3
	  lenrec3 = ntohl(lenrec3); // bk 6 jul 2000, byteorder x86 machines.
	  if (lenrec3 != 360)
		{
		WARNING.print("probably something wrong here, byte order x86. (?)");
		WARNING << "readvolume: length of record 3 = \"" << lenrec3 << "\"; expected \"360\" for ESA SLC (full scene).";
		WARNING.print();
		}
	  volumefile.seekg(startrec3+20,ios.beg);
	  volumefile.read((char)&c16dataref, sizea16); // referenced file name
	  c16dataref.charAt(16)='\0';
	  volumefile.read((char)&c28datarefclass, sizea28); // referenced file class
	  c28datarefclass.charAt(28)='\0';
	  volumefile.seekg(startrec3+100,ios.beg); // numlines for checking
	  volumefile.read((char)&c8nlins, sizei8);
	  c8nlins.charAt(8)='\0';


	  // --- RECORD 4 ---
	  uint startrec4 = lenrec1 + lenrec2 + lenrec3;
	  //if (readfiles_arg.sensor_id==SLC_RSAT)
	  //   startrec4=startrec4+360;// skip trailer record for RSAT
	  DEBUG.print("record 4 of volume file.");
	  DEBUG << "readvolume::rec4: start at byte: " << startrec4;
	  DEBUG.print();
	  // --- Check this record, RSAT has for us useless trailer info as rec.4 ---
	  // --- while next record of 360 is the text record of ERS ---
	  // ---start test for ATLANTIS (RSAT) ------------------------------------------
	  DEBUG.print("VMP: Expecting record 4 with code {18,63,18,18}");
	  volumefile.seekg(startrec4,ios.beg);
	  volumefile.read((char) rec_seq, sizeb4); // record number
	  rec_seq = ntohl(rec_seq); // Bert Kampes, 07-Apr-2005
	  volumefile.read((char)&rec_sub1, sizeb1); // first record sub type code
	  volumefile.read((char)&rec_type, sizeb1); // record type code
	  volumefile.read((char)&rec_sub2, sizeb1); // second record sub type code
	  volumefile.read((char)&rec_sub3, sizeb1); // third record sub type code
	  DEBUG << "rec_seq: " << rec_seq << "; rec_sub1: " << (int)rec_sub1 << "; rec_type: " << (int)rec_type << "; rec_sub2: " << (int)rec_sub2 << "; rec_sub3: " << (int)rec_sub3;
	  DEBUG.print();
	  if ((int)rec_sub1 ==18 && (int)rec_type ==63 && (int)rec_sub2 ==18 && (int)rec_sub3 ==18)
		{
		DEBUG.print("This is the expected text record with code {18,63,18,18}");
		readfiles_arg.argvalue.sar_processor = SARPR_VMP; // set determined sar processor/format
		}
	  else
		{
		//WARNING.print("This is not the expected text record, trying next one.")
		readfiles_arg.argvalue.sar_processor = SARPR_ATL; // set determined sar processor/format

			// Modified by LG for reading ALOS Fine
			if ((readfiles_arg.argvalue.sensor_id == SLC_ALOS) && ((int)rec_type ==192))
							 readfiles_arg.argvalue.sar_processor = SARPR_JAX;
			else if (readfiles_arg.argvalue.sensor_id!=SLC_RSAT)
		  {
		  DEBUG.print("This is NOT the expected text record with code {18,63,18,18}");
		  DEBUG.print("I assume in following Atlantis processed data (ERS or RSAT).");
		  DEBUG.print("If this is not the case, please contact doris_users@tudelft.nl");
		  }
		volumefile.seekg(startrec4+8,ios.beg); // text record
		volumefile.read((char) lenrec4, sizeb4); // length of record4
		lenrec4 = ntohl(lenrec4); // bk 6 jul 2000, byteorder x86 machines.
		startrec4 = startrec4 + lenrec4; // hopefully there is a next record
		}
	  // ---end test for ATLANTIS (RSAT) -------------------------------------------

	  // --- Read important fields ---
	  volumefile.seekg(startrec4+8,ios.beg); // text record
	  volumefile.read((char) lenrec4, sizeb4); // length of record4
	  lenrec4 = ntohl(lenrec4); // bk 6 jul 2000, byteorder x86 machines.
	  if (lenrec4 != 360)
		{
		WARNING.print("probably something wrong here, byte order x86. (?)");
		WARNING << "readvolume: length of record 3 = \"" << lenrec4 << "\"; expected \"360\" for ESA SLC (full scene).";
		WARNING.print();
		}
	  volumefile.seekg(startrec4+16,ios.beg);
	  volumefile.read((char)&c40typespec, sizea40); // product type specifier
	  c40typespec.charAt(40)='\0';
	  volumefile.read((char)&c60product, sizea60); // loc&date product gen.
	  c60product.charAt(60)='\0';
	  volumefile.read((char)&c40physvolid, sizea40); // physical vol id
	  c40physvolid.charAt(40)='\0';
	  volumefile.read((char)&c40sceneid, sizea40); // scene id
	  c40sceneid.charAt(40)='\0';
	  volumefile.read((char)&c40sceneloc, sizea40); // scene loc
	  c40sceneloc.charAt(40)='\0';

	  volumefile.close();


	// ====== Write results to scratch files ======
	  ofstream scratchlogfile = new ofstream("scratchlogvol", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchlogfile, "scratchlogvol", __FILE__, __LINE__);

	// ______Write information to scratchfiles______
	  scratchlogfile << "\n*******************************************************************" << "\n* EXTRACTED DATA FROM VOLUME FILE: " << readfiles_arg.argvalue.volfile << " *" << "\n*******************************************************************" << "\n\nVolume descriptor record" << "\n------------------------" << "\nLogical volume generating facility software " << "\n +release and revision level: \t\t\t" << c12logvol << "\nID of physical volume containing " << "\n +this volume descriptor: \t\t\t" << c16physid << "\nLogical volume identifier: \t\t\t" << c16logvolid << "\nVolume set identifier: \t\t\t\t" << c16setid << "\nLogical volume creation date (YYYYMMDD): \t" << c8date << "\nLogical volume creation time (HHMMSSDD): \t" << c8time << "\nLogical volume generation country: \t\t" << c12country << "\nLogical volume generating agency: \t\t" << c8agency << "\nLogical volume generation facility: \t\t" << c12facility << "\n\nLeader file pointer record" << "\n--------------------------" << "\nReferenced file name: \t\t\t\t" << c16checkfilename << "\nReferenced file class: \t\t\t\t" << c28leaderrefclass << "\n\nData file pointer record" << "\n------------------------" << "\nReferenced file name: \t\t\t\t" << c16dataref << "\nReferenced file class: \t\t\t\t" << c28datarefclass << "\nNumber of records in referenced file: \t\t" << c8nlins << "\n\nText record" << "\n-----------" << "\nProduct type specifier: \t\t\t" << c40typespec << "\nLocation and date/time of product creation: \t" << c60product << "\nPhysical volume identification: \t\t" << c40physvolid << "\nScene identification: \t\t\t\t" << c40sceneid << "\nScene location: \t\t\t\t" << c40sceneloc << "\nEND VOLUME FILE " << readfiles_arg.argvalue.volfile << "\n" << "\n";
	  scratchlogfile.close();

	  ofstream scratchresfile = new ofstream("scratchresvol", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile, "scratchresvol", __FILE__, __LINE__);
	  scratchresfile << "\n*******************************************************************";
	  if (readfiles_arg.argvalue.fileid==MASTERID)
		scratchresfile << "\n*_Start_" << processcontrol[pr_m_readfiles];
	  if (readfiles_arg.argvalue.fileid==SLAVEID)
		scratchresfile << "\n*_Start_" << processcontrol[pr_m_readfiles];
	  scratchresfile << "\n*******************************************************************" << "\nVolume file: \t\t\t\t\t" << readfiles_arg.argvalue.volfile << "\nVolume_ID: \t\t\t\t\t" << c16physid << "\nVolume_identifier: \t\t\t\t" << c16logvolid << "\nVolume_set_identifier: \t\t\t\t" << c16setid << "\n(Check)Number of records in ref. file: \t\t" << Integer.parseInt(c8nlins);
	  // --- write new string to determine sar processor ---
	  // --- this is used in cropping the data ---
	  if (readfiles_arg.argvalue.sar_processor==SARPR_VMP)
		scratchresfile << "\nSAR_PROCESSOR:  \t\t\t\tVMP";
	  if (readfiles_arg.argvalue.sar_processor==SARPR_ATL)
		scratchresfile << "\nSAR_PROCESSOR:  \t\t\t\tATLANTIS";
	  if (readfiles_arg.argvalue.sar_processor==SARPR_TUD)
		scratchresfile << "\nSAR_PROCESSOR:  \t\t\t\tTUD";
	  // start_added_by_don
	  if (readfiles_arg.argvalue.sar_processor==SARPR_JAX)
		scratchresfile << "\nSAR_PROCESSOR:  \t\t\t\tJAX";
	  // end_added_by_don
	  if (readfiles_arg.argvalue.sensor_id==SLC_RSAT)
		scratchresfile << "\nProduct type specifier: \t\t\t" << "RSAT"; // else it is just "PRODUCT:", but it is used in CROP step
	  // start_added_by_don
	  else if (readfiles_arg.argvalue.sensor_id==SLC_ALOS)
		scratchresfile << "\nProduct type specifier: \t\t\t" << "PRODUCT ALOS"; // else it is just "PRODUCT:", but it is used in CROP step
	  // end_added_by_don
	  else
		scratchresfile << "\nProduct type specifier: \t\t\t" << c40typespec;
	  scratchresfile << "\nLogical volume generating facility: \t\t\t" << c12facility << "\nLogical volume creation date: \t\t\t" << c8date << "\nLocation and date/time of product creation: \t" << c60product << "\nScene identification: \t\t\t\t" << c40sceneid << "\nScene location: \t\t\t\t" << c40sceneloc << "\n";
	  scratchresfile.close();

	// ______Tidy up______
	  PROGRESS.print("readvolume finished.");
	  } // END READVOLUME



	//***************************************************************
	// *    readleader                                                *
	// *                                                              *
	// * reads leaderfile                                             *
	// *  and writes to scratchlogfile, scratchresfile                *
	// * checks with volumefile #lines                                *
	// *                                                              *
	// * input:                                                       *
	// *  - struct with arguments for step0                           *
	// *  - check with volumefile                                     *
	// * output:                                                      *
	// *  - scratchlogfile                                            *
	// *  - scratchresfile                                            *
	// *  - scratchdatapoints                                         *
	// * see: annex C ERS SAR.SLC CCTand EXABYTE                      *
	// *      doc:er-is-epo-gs-5902.3                                 *
	// *                                                              *
	// *    Bert Kampes, 11-Dec-1998                                  *
	// * Included RSAT format based on document of ASF                *
	// #%// Bert Kampes, 03-Aug-2004                                  *
	// #%// Davide Nitti (Don), 11-Nov-2008  fixes for doppler        *
	// #     coefficient unit for Radarsat1 and ALOS                  * 
	// ***************************************************************
	public static void readleader(RefObject<input_readfiles> readfiles_arg, int32 checklines)
	  {
	  TRACE_FUNCTION("readleader (BK 11-Dec-1998)")
	  final int16 sizea2 = 2;
	  final int16 sizeb1 = 1;
	  final int16 sizeb4 = 4;
	  final int16 sizei4 = 4;
	  final int16 sizea4 = 4;
	  final int16 sizei8 = 8;
	  final int16 sizea8 = 8;
	  final int16 sizea12 = 12;
	  final int16 sizea16 = 16;
	  final int16 sizef16 = 16;
	  final int16 sizei16 = 16;
	  final int16 sized22 = 22;
	  final int16 sizea24 = 24;
	  final int16 sizea32 = 32;
	  final int16 sizea64 = 64;
	  uint lenrec1; // bk rsat record
	  uint lenrec2;
	  uint lenrec3;
	  uint lenrec4;
	  uint lenrec5;
	  uint lenrec6;
	  uint lenrec7;
	 String c2motioncomp = new String(new char[3]);
	 String c4dummy = new String(new char[5]);
	 String c4year = new String(new char[5]);
	 String c4month = new String(new char[5]);
	 String c4day = new String(new char[5]);
	 String c4dayofyear = new String(new char[5]);
	 String c4conversion = new String(new char[5]);
	 String c4compression = new String(new char[5]);
	 String c4clutterlock = new String(new char[5]);
	 String c4autofocus = new String(new char[5]);
	 String c8dummy = new String(new char[9]);
	 String c8orbitnr = new String(new char[9]);
	 String c8platformlat = new String(new char[9]);
	 String c8platformlon = new String(new char[9]);
	 String c8platformheading = new String(new char[9]);
	 String c8clockangle = new String(new char[9]);
	 String c8incidence = new String(new char[9]);
	 String c8freq = new String(new char[9]);
	 String c8systemid = new String(new char[9]);
	 String c8versionid = new String(new char[9]);
	 String c8satclockstep = new String(new char[9]);
	 String c8extindex = new String(new char[9]);
	 String c8qperch = new String(new char[9]);
	 String c8timeline = new String(new char[9]);
	 String c8timepix = new String(new char[9]);
	 String c8linecontent = new String(new char[9]);
	 String c12qdesc = new String(new char[13]);
	 String c16dummy = new String(new char[17]);
	 String c16lat11 = new String(new char[17]);
	 String c16lon11 = new String(new char[17]);
	 String c16lat1N = new String(new char[17]);
	 String c16lon1N = new String(new char[17]);
	 String c16latNN = new String(new char[17]);
	 String c16lonNN = new String(new char[17]);
	 String c16latN1 = new String(new char[17]);
	 String c16lonN1 = new String(new char[17]);
	 String c16leafilename = new String(new char[17]);
	 String c16centerlat = new String(new char[17]);
	 String c16centerlon = new String(new char[17]);
	 String c16centerheading = new String(new char[17]);
	 String c16ellipsoid = new String(new char[17]);
	 String c16semimajor = new String(new char[17]);
	 String c16semiminor = new String(new char[17]);
	 String c16GM = new String(new char[17]);
	 String c16J2 = new String(new char[17]);
	 String c16J3 = new String(new char[17]);
	 String c16J4 = new String(new char[17]);
	 String c16scenelength = new String(new char[17]);
	 String c16scenewidth = new String(new char[17]);
	 String c16platformid = new String(new char[17]);
	 String c16wavelength = new String(new char[17]);
	 String c16facilityid = new String(new char[17]);
	 String c16numpix = new String(new char[17]);
	 String c16numlin = new String(new char[17]);
	 String c16interpix = new String(new char[17]);
	 String c16interlin = new String(new char[17]);
	 String c16pulse = new String(new char[17]);
	 String c16ampconst = new String(new char[17]);
	 String c16amplinear = new String(new char[17]);
	 String c16ampquadratic = new String(new char[17]);
	 String c16ampcubic = new String(new char[17]);
	 String c16ampquartic = new String(new char[17]);
	 String c16phaseconst = new String(new char[17]);
	 String c16phaselinear = new String(new char[17]);
	 String c16phasequadratic = new String(new char[17]);
	 String c16phasecubic = new String(new char[17]);
	 String c16phasequartic = new String(new char[17]);
	 String c16sattimecode = new String(new char[17]);
	 String c16samplingrate = new String(new char[17]);
	 String c16rangedelay = new String(new char[17]);
	 String c16ranpulselen = new String(new char[17]);
	 String c16dci = new String(new char[17]);
	 String c16dcq = new String(new char[17]);
	 String c16boresight = new String(new char[17]);
	 String c16imbalance = new String(new char[17]);
	 String c16prf = new String(new char[17]);
	 String c16looksazi = new String(new char[17]);
	 String c16looksrange = new String(new char[17]);
	 String c16bandazi = new String(new char[17]);
	 String c16bandrange = new String(new char[17]);
	 String c16bandazitot = new String(new char[17]);
	 String c16bandrangetot = new String(new char[17]);
	 String c16inputsource = new String(new char[17]);
	 String c16resrange = new String(new char[17]);
	 String c16resazi = new String(new char[17]);
	 String c16linespace = new String(new char[17]);
	 String c16pixspace = new String(new char[17]);
	 String c16atdoppcconst = new String(new char[17]);
	 String c16atdoppclinear = new String(new char[17]);
	 String c16atdoppcquadratic = new String(new char[17]);
	 String c16xtdoppcconst = new String(new char[17]);
	 String c16xtdoppclinear = new String(new char[17]);
	 String c16xtdoppcquadratic = new String(new char[17]);
	 String c16atdopprconst = new String(new char[17]);
	 String c16atdopprlinear = new String(new char[17]);
	 String c16atdopprquadratic = new String(new char[17]);
	 String c16xtdopprconst = new String(new char[17]);
	 String c16xtdopprlinear = new String(new char[17]);
	 String c16xtdopprquadratic = new String(new char[17]);
	 String c16rcompdes = new String(new char[17]);
	 String c16zd1strange = new String(new char[17]);
	 String c16zdcenrange = new String(new char[17]);
	 String c16zdlstrange = new String(new char[17]);
	 String c16orien = new String(new char[17]);
	 String c16platincl = new String(new char[17]);
	 String c16platascn = new String(new char[17]);
	 String c16geocenter = new String(new char[17]);
	 String c16platalt = new String(new char[17]);
	 String c16plathead = new String(new char[17]);
	 String c16platgs = new String(new char[17]);
	 String c16refmajor = new String(new char[17]);
	 String c16refminor = new String(new char[17]);
	 String c16ltposerr = new String(new char[17]);
	 String c16ctposerr = new String(new char[17]);
	 String c16rposerr = new String(new char[17]);
	 String c22dummy = new String(new char[23]);
	 String c22seconds = new String(new char[23]);
	 String c22interval = new String(new char[23]);
	 String c22gmha = new String(new char[23]);
	 String c24zd1stazitime = new String(new char[25]);
	 String c24zdcenazitime = new String(new char[25]);
	 String c24zdlstazitime = new String(new char[25]);
	 String c32sceneref = new String(new char[33]);
	 String c32scenetime = new String(new char[33]);
	 String c32sensorid = new String(new char[33]);
	 String c32typespec = new String(new char[33]);
	 String c32algid = new String(new char[33]);
	 String c32projection = new String(new char[33]);
	 String c32sattime = new String(new char[33]);
	 String c32weightrange = new String(new char[33]);
	 String c32weightazi = new String(new char[33]);
	 String c32refellips = new String(new char[33]);
	 String c64rcs = new String(new char[65]);

	  // --- Check for RSAT #%// Bert Kampes, 02-Aug-2004 ---
	  uint rec_seq; // type B4
	  byte rec_sub1; // type B1
	  byte rec_type;
	  byte rec_sub2;
	  byte rec_sub3;



	// ======Open files====== 
	  ifstream leaderfile;
	  RefObject<ifstream> TempRefObject = new RefObject<ifstream>(leaderfile);
	  openfstream(TempRefObject, readfiles_arg.argvalue.leaderfile);
	  leaderfile = TempRefObject.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(leaderfile, readfiles_arg.argvalue.leaderfile, __FILE__, __LINE__);

	  // ======Read leaderfile======
	  // --- RECORD 1 ---
	  leaderfile.seekg(0,ios.beg); // place pointer at beginning of file
	  leaderfile.read((char) rec_seq, sizeb4); // record number
	  rec_seq = ntohl(rec_seq); // Bert Kampes, 07-Apr-2005
	  leaderfile.read((char)&rec_sub1, sizeb1); // first record sub type code
	  leaderfile.read((char)&rec_type, sizeb1); // record type code
	  leaderfile.read((char)&rec_sub2, sizeb1); // second record sub type code
	  leaderfile.read((char)&rec_sub3, sizeb1); // third record sub type code
	  DEBUG.print("Expecting record 1 with code {63,192,18,18}");
	  DEBUG << "rec_seq: " << rec_seq << "; rec_sub1: " << (int)rec_sub1 << "; rec_type: " << (int)rec_type << "; rec_sub2: " << (int)rec_sub2 << "; rec_sub3: " << (int)rec_sub3;
	  DEBUG.print();

	  // --- Read important data ---
	  leaderfile.seekg(8,ios.beg); // file descriptor record
	  leaderfile.read((char) lenrec1, sizeb4); // length of record1
	  lenrec1 = ntohl(lenrec1); // bk 6 jul 2000, byteorder x86 machines.
	  DEBUG << "readleader::record 1 of leader file: length: " << lenrec1;
	  DEBUG.print();
	  if (lenrec1 != 720)
		{
		WARNING.print("probably something wrong here, byte order x86.");
		WARNING << "readleader: length of record 1 = \"" << lenrec1 << "\"; expected \"720\" for ESA SLC (full scene).";
		WARNING.print();
		}
	  leaderfile.seekg(48,ios.beg);
	  leaderfile.read((char)&c16leafilename, sizea16); // file name
	  c16leafilename.charAt(16)='\0';


	  // --- RECORD 2 ---
	  DEBUG.print("readleader::reading record 2 of leader file.");
	  uint startrec2 = lenrec1;
	  leaderfile.seekg(startrec2,ios.beg); // place pointer at beginning of record
	  leaderfile.read((char) rec_seq, sizeb4); // record number
	  rec_seq = ntohl(rec_seq); // Bert Kampes, 07-Apr-2005
	  leaderfile.read((char)&rec_sub1, sizeb1); // first record sub type code
	  leaderfile.read((char)&rec_type, sizeb1); // record type code
	  leaderfile.read((char)&rec_sub2, sizeb1); // second record sub type code
	  leaderfile.read((char)&rec_sub3, sizeb1); // third record sub type code
	  DEBUG.print("ERS:  Expecting record 2 with code {10,10,31,20}");
	  DEBUG.print("RSAT: Expecting record 2 with code {18,10,18,20}");
	  DEBUG.print("RSAT record length is 4096, ERS 1886, but");
	  DEBUG.print("ERS contains more info on zero doppler times, etc.");
	  DEBUG.print("RSAT seems to have that info in the data file.");
	  DEBUG << "rec_seq: " << rec_seq << "; rec_sub1: " << (int)rec_sub1 << "; rec_type: " << (int)rec_type << "; rec_sub2: " << (int)rec_sub2 << "; rec_sub3: " << (int)rec_sub3;
	  DEBUG.print();

	  // ___ Read important parameters ___
	  leaderfile.seekg(startrec2+8,ios.beg); // slc data set summary record
	  leaderfile.read((char) lenrec2, sizeb4); // length of record2
	  lenrec2 = ntohl(lenrec2); // bk 6 jul 2000, byteorder x86 machines.
	  DEBUG << "readleader::record 2: start at: " << lenrec1 << "; length: " << lenrec2;
	  DEBUG.print();
	  //if (readfiles_arg.sensor_id==SLC_RSAT)
	  if (readfiles_arg.argvalue.sar_processor==SARPR_ATL)
		{
		if (lenrec2 != 4096)
		  {
		  WARNING.print("SARPR_ATL (RSAT) has 4096 record length, but other value found?");
		  }
		}
	  // start_added_by_don
	  else if (readfiles_arg.argvalue.sar_processor==SARPR_JAX)
		{
		if (lenrec2 != 4096)
		  {
		  WARNING.print("SARPR_JAX (ALOS) has 4096 record length, but other value found?");
		  }
		}
	  // end_added_by_don
	  else
		{
		if (lenrec2 != 1886)
		  {
		  WARNING.print("SARPR_ATL (RSAT) has 4096 record length");
		  WARNING.print("probably something wrong here, byte order on x86?");
		  WARNING << "readleader: length of record 2 = \"" << lenrec2 << "\"; expected \"1886\" for ESA SLC (full scene).  continuing";
		  WARNING.print();
		  }
		}


	  // ______Scene parameters______
	 // Modified by LG for reading ALOS Fine
	 if (readfiles_arg.argvalue.sensor_id == SLC_ALOS)
	 {
			 leaderfile.seekg(startrec2+20,ios.beg);
			 leaderfile.read((char)&c32sceneref, sizea32);
			 leaderfile.seekg(startrec2+20+sizea32+1,ios.beg);
	 }
	 else
	 {
			 leaderfile.seekg(startrec2+36,ios.beg);
			 leaderfile.read((char)&c32sceneref, sizea32); // scene ref. number
	 }
	  c32sceneref.charAt(32)='\0';
	  leaderfile.read((char)&c32scenetime, sizea32); // scene center time
	  c32scenetime.charAt(32)='\0';
	  leaderfile.seekg(startrec2+116,ios.beg);
	  leaderfile.read((char)&c16centerlat, sizef16); // centre latitude
	  c16centerlat.charAt(16)='\0';
	  leaderfile.read((char)&c16centerlon, sizef16); // centre longitude
	  c16centerlon.charAt(16)='\0';
	  leaderfile.read((char)&c16centerheading, sizef16); // center true heading
	  c16centerheading.charAt(16)='\0';
	  leaderfile.read((char)&c16ellipsoid, sizea16); // ell. designated
	  c16ellipsoid.charAt(16)='\0';
	  leaderfile.read((char)&c16semimajor, sizef16); // ell. semi major
	  c16semimajor.charAt(16)='\0';
	  leaderfile.read((char)&c16semiminor, sizef16); // ell. semi minor
	  c16semiminor.charAt(16)='\0';
	  leaderfile.read((char)&c16GM, sizef16); // GM
	  c16GM.charAt(16)='\0';
	  leaderfile.read((char)&c16dummy, sizef16); // dummy
	  leaderfile.read((char)&c16J2, sizef16); // J2
	  c16J2.charAt(16)='\0';
	  leaderfile.read((char)&c16J3, sizef16); // J3
	  c16J3.charAt(16)='\0';
	  leaderfile.read((char)&c16J4, sizef16); // J4
	  c16J4.charAt(16)='\0';
	  leaderfile.read((char)&c16dummy, sizef16); // dummy
	  leaderfile.read((char)&c16dummy, sizef16); // dummy
	  leaderfile.read((char)&c8dummy, sizei8); // center line#
	  c8dummy.charAt(8)='\0';
	  uint scenecenterline = Integer.parseInt(c8dummy);
	  leaderfile.read((char)&c8dummy, sizei8); // center pixel#
	  c8dummy.charAt(8)='\0';
	  uint scenecenterpixel = Integer.parseInt(c8dummy);
	  leaderfile.read((char)&c16scenelength, sizef16); // scene length
	  c16scenelength.charAt(16)='\0';
	  leaderfile.read((char)&c16scenewidth, sizef16); // scene width
	  c16scenewidth.charAt(16)='\0';

	// ______General mission / sensor parameters______
	  leaderfile.seekg(startrec2+396,ios.beg);
	  leaderfile.read((char)&c16platformid, sizea16); // platform mission id
	  c16platformid.charAt(16)='\0';
	  leaderfile.read((char)&c32sensorid, sizea32); // sensor id
	  c32sensorid.charAt(32)='\0';
	  leaderfile.read((char)&c8orbitnr, sizea8); // orbit number
	  c8orbitnr.charAt(8)='\0';
	  leaderfile.read((char)&c8platformlat, sizea8); // platform latitude
	  c8platformlat.charAt(8)='\0';
	  leaderfile.read((char)&c8platformlon, sizea8); // platform longitude
	  c8platformlon.charAt(8)='\0';
	  leaderfile.read((char)&c8platformheading, sizea8); // platform heading
	  c8platformheading.charAt(8)='\0';
	  leaderfile.read((char)&c8clockangle, sizea8); // sensor clock angle
	  c8clockangle.charAt(8)='\0';
	  leaderfile.read((char)&c8incidence, sizea8); // incidence angle
	  c8incidence.charAt(8)='\0';
	  leaderfile.read((char)&c8freq, sizea8); // radar frequency
	  c8freq.charAt(8)='\0';
	  leaderfile.read((char)&c16wavelength, sizea16); // radar wavelength
	  c16wavelength.charAt(16)='\0';
	  leaderfile.read((char)&c2motioncomp, sizea2); // indicator for compensation
	  c2motioncomp.charAt(2)='\0';
	  leaderfile.read((char)&c16pulse, sizea16); // range pulse code specifier
	  c16pulse.charAt(16)='\0';
	  leaderfile.read((char)&c16ampconst, sizef16); // amplitude constant term
	  c16ampconst.charAt(16)='\0';
	  leaderfile.read((char)&c16amplinear, sizef16); // amplitude linear term
	  c16amplinear.charAt(16)='\0';
	  leaderfile.read((char)&c16ampquadratic, sizef16); // amplitude quadrati term
	  c16ampquadratic.charAt(16)='\0';
	  leaderfile.read((char)&c16ampcubic, sizef16); // amplitude cubic term
	  c16ampcubic.charAt(16)='\0';
	  leaderfile.read((char)&c16ampquartic, sizef16); // amplitude quartic term
	  c16ampquartic.charAt(16)='\0';
	  leaderfile.read((char)&c16phaseconst, sizef16); // phase constant term
	  c16phaseconst.charAt(16)='\0';
	  leaderfile.read((char)&c16phaselinear, sizef16); // phase linear term
	  c16phaselinear.charAt(16)='\0';
	  leaderfile.read((char)&c16phasequadratic, sizef16); // phase quadratic term
	  c16phasequadratic.charAt(16)='\0';
	  leaderfile.read((char)&c16phasecubic, sizef16); // phase cubicterm
	  c16phasecubic.charAt(16)='\0';
	  leaderfile.read((char)&c16phasequartic, sizef16); // phase quartic term
	  c16phasequartic.charAt(16)='\0';
	  leaderfile.read((char)&c8extindex, sizei8); // chirp extraction
	  c8extindex.charAt(8)='\0';
	  leaderfile.read((char)&c8dummy, sizei8); // spare
	  leaderfile.read((char)&c16samplingrate, sizef16); // range sampling rate
	  c16samplingrate.charAt(16)='\0';
	  leaderfile.read((char)&c16rangedelay, sizef16); // delay
	  c16rangedelay.charAt(16)='\0';
	  leaderfile.read((char)&c16ranpulselen, sizef16); // range pulselength
	  c16ranpulselen.charAt(16)='\0';
	  leaderfile.read((char)&c4conversion, sizea4); // flag
	  c4conversion.charAt(4)='\0';
	  leaderfile.read((char)&c4compression, sizea4); // flag
	  c4compression.charAt(4)='\0';
	  leaderfile.read((char)&c16dummy, sizef16); // reserved
	  leaderfile.read((char)&c16dummy, sizef16); // reserved
	  leaderfile.read((char)&c8qperch, sizei8); // quantization
	  c8qperch.charAt(8)='\0';
	  leaderfile.read((char)&c12qdesc, sizea12); // quantization description
	  c12qdesc.charAt(12)='\0';
	  leaderfile.read((char)&c16dci, sizef16); // bias for i comp.
	  c16dci.charAt(16)='\0';
	  leaderfile.read((char)&c16dcq, sizef16); // bias for q comp.
	  c16dcq.charAt(16)='\0';
	  leaderfile.read((char)&c16imbalance, sizef16); // gain imbalance i&q
	  c16imbalance.charAt(16)='\0';
	  leaderfile.read((char)&c16dummy, sizef16); // spare
	  leaderfile.read((char)&c16dummy, sizef16); // spare
	  leaderfile.read((char)&c16dummy, sizef16); // reserved
	  leaderfile.read((char)&c16boresight, sizef16); // antenna
	  c16boresight.charAt(16)='\0';
	  leaderfile.read((char)&c4dummy, sizea4); // reserved
	  leaderfile.read((char)&c16prf, sizef16); // pulse repetition frequency
	  c16prf.charAt(16)='\0';

	// ______Sensor specific parameters______
	  leaderfile.seekg(startrec2+982,ios.beg);
	  leaderfile.read((char)&c16sattimecode, sizei16); // sat time code
	  c16sattimecode.charAt(16)='\0';
	  leaderfile.read((char)&c32sattime, sizea32); // sat time
	  c32sattime.charAt(32)='\0';
	  leaderfile.read((char)&c8satclockstep, sizei8); // sat clock step length
	  c8satclockstep.charAt(8)='\0';

	// ______General processing parameters______
	  leaderfile.seekg(startrec2+1046,ios.beg);
	  leaderfile.read((char)&c16facilityid, sizea16); // proc. facility id
	  c16facilityid.charAt(16)='\0';
	  leaderfile.read((char)&c8systemid, sizea8); // proc. system id
	  c8systemid.charAt(8)='\0';
	  leaderfile.read((char)&c8versionid, sizea8); // proc. version id
	  c8versionid.charAt(8)='\0';
	  leaderfile.read((char)&c16dummy, sizef16); // dummy
	  leaderfile.read((char)&c16dummy, sizef16); // dummy
	  leaderfile.read((char)&c32typespec, sizea32); // produkt type spec.
	  c32typespec.charAt(32)='\0';
	  leaderfile.read((char)&c32algid, sizea32); // proc. alg. id
	  c32algid.charAt(32)='\0';
	  leaderfile.read((char)&c16looksazi, sizef16); // number of looks
	  c16looksazi.charAt(16)='\0';
	  leaderfile.read((char)&c16looksrange, sizef16); // number of looks
	  c16looksrange.charAt(16)='\0';
	  leaderfile.read((char)&c16bandazi, sizef16); // bandwidth
	  c16bandazi.charAt(16)='\0';
	  leaderfile.read((char)&c16bandrange, sizef16); // bandwidth
	  c16bandrange.charAt(16)='\0';
	  leaderfile.read((char)&c16bandazitot, sizef16); // bandwidth
	  c16bandazitot.charAt(16)='\0';
	  leaderfile.read((char)&c16bandrangetot, sizef16); // bandwidth
	  c16bandrangetot.charAt(16)='\0';
	  leaderfile.read((char)&c32weightazi, sizea32); // weighting function
	  c32weightazi.charAt(32)='\0';
	  leaderfile.read((char)&c32weightrange, sizea32); // weighting function
	  c32weightrange.charAt(32)='\0';
	  leaderfile.read((char)&c16inputsource, sizea16); // data input
	  c16inputsource.charAt(16)='\0';
	  leaderfile.read((char)&c16resrange, sizef16); // resolution
	  c16resrange.charAt(16)='\0';
	  leaderfile.read((char)&c16resazi, sizef16); // resolution
	  c16resazi.charAt(16)='\0';
	  leaderfile.read((char)&c16dummy, sizef16); // reserved
	  leaderfile.read((char)&c16dummy, sizef16); // reserved
	  leaderfile.read((char)&c16atdoppcconst, sizef16); // along track centroid
	  c16atdoppcconst.charAt(16)='\0';
	  leaderfile.read((char)&c16atdoppclinear, sizef16); // along track centroid
	  c16atdoppclinear.charAt(16)='\0';
	  leaderfile.read((char)&c16atdoppcquadratic, sizef16); // along track centroid
	  c16atdoppcquadratic.charAt(16)='\0';
	  leaderfile.read((char)&c16dummy, sizef16); // spare
	  leaderfile.read((char)&c16xtdoppcconst, sizef16); // cross track centroid
	  c16xtdoppcconst.charAt(16)='\0';
	  leaderfile.read((char)&c16xtdoppclinear, sizef16); // cross track centroid
	  c16xtdoppclinear.charAt(16)='\0';
	  leaderfile.read((char)&c16xtdoppcquadratic, sizef16); // cross track centroid
	  c16xtdoppcquadratic.charAt(16)='\0';
	  leaderfile.read((char)&c8timepix, sizea8); // time direction
	  c8timepix.charAt(8)='\0';
	  leaderfile.read((char)&c8timeline, sizea8); // time direction
	  c8timeline.charAt(8)='\0';
	  leaderfile.read((char)&c16atdopprconst, sizef16); // along track rate
	  c16atdopprconst.charAt(16)='\0';
	  leaderfile.read((char)&c16atdopprlinear, sizef16); // along track rate
	  c16atdopprlinear.charAt(16)='\0';
	  leaderfile.read((char)&c16atdopprquadratic, sizef16); // along track rate
	  c16atdopprquadratic.charAt(16)='\0';
	  leaderfile.read((char)&c16dummy, sizef16); // spare
	  leaderfile.read((char)&c16xtdopprconst, sizef16); // cross track rate
	  c16xtdopprconst.charAt(16)='\0';
	  leaderfile.read((char)&c16xtdopprlinear, sizef16); // cross track rate
	  c16xtdopprlinear.charAt(16)='\0';
	  leaderfile.read((char)&c16xtdopprquadratic, sizef16); // cross track rate
	  c16xtdopprquadratic.charAt(16)='\0';
	  leaderfile.read((char)&c16dummy, sizef16); // spare
	  leaderfile.read((char)&c8linecontent, sizea8); // indicator
	  c8linecontent.charAt(8)='\0';
	  leaderfile.read((char)&c4clutterlock, sizea4); // flag
	  c4clutterlock.charAt(4)='\0';
	  leaderfile.read((char)&c4autofocus, sizea4); // flag
	  c4autofocus.charAt(4)='\0';
	  leaderfile.read((char)&c16linespace, sizef16);
	  c16linespace.charAt(16)='\0';
	  leaderfile.read((char)&c16pixspace, sizef16);
	  c16pixspace.charAt(16)='\0';
	  leaderfile.read((char)&c16rcompdes, sizea16); // range compression designator
	  c16rcompdes.charAt(16)='\0';



	  // --- RSAT does not fill this spare part of the record (blanks) ---
	  DEBUG.print("Assuming user has specified method RSAT if it is RSAT.");
	  DEBUG.print("although we could also easily detect it before");
	  boolean skipmapprojrecord = false; // some problem s with old IPAF?
	  uint numdatapoints =99999;
	  String c16incangle1strange = new String(new char[17]); //gk
	  String c16incanglecenrange = new String(new char[17]);
	  String c16incanglelstrange = new String(new char[17]);
	  String calK = new String(new char[17]);
	  String repplspwr = new String(new char[17]);
	  String c4numvalid = new String(new char[5]); //bk
	  String c4numinvalid = new String(new char[5]); //bk
	  c16incangle1strange = "skipped";
	  c16incanglecenrange = "skipped";
	  c16incanglelstrange = "skipped";
	  calK = "skipped";
	  repplspwr = "skipped";
	  c4numvalid = "999";
	  c4numinvalid = "999";
	  matrix<real8> STATE; // has to be declared at large scope, used later.
	  matrix<real8> STATE_INERTIAL; // has to be declared at large scope, used later.


	//if (readfiles_arg.sar_processor!=SARPR_ATL)
	// seems ERS data are same, also with Atlantis, so only check for RSAT here.
	if (readfiles_arg.argvalue.sensor_id!=SLC_RSAT)
	  {
	  DEBUG.print("Reading rest of this record, which for RSAT is empty.");

	  // ______Sensor specific local use segment______
	  leaderfile.seekg(startrec2+1766,ios.beg);
	  leaderfile.read((char)&c16zd1strange, sizef16); // zero doppler 1st pixel
	  c16zd1strange.charAt(16)='\0';
	  leaderfile.read((char)&c16zdcenrange, sizef16); // zero doppler centre pixel
	  c16zdcenrange.charAt(16)='\0';
	  leaderfile.read((char)&c16zdlstrange, sizef16); // zero doppler last pixel 2way
	  c16zdlstrange.charAt(16)='\0';
	  leaderfile.read((char)&c24zd1stazitime, sizea24); // zero doppler 1st pixel
	  c24zd1stazitime.charAt(24)='\0';
	  leaderfile.read((char)&c24zdcenazitime, sizea24); // zero doppler 1st pixel
	  c24zdcenazitime.charAt(24)='\0';
	  leaderfile.read((char)&c24zdlstazitime, sizea24); // zero doppler 1st pixel
	  c24zdlstazitime.charAt(24)='\0';


	  // --- RECORD 3 ---
	  DEBUG.print("record 3 of leader file (ERS).");
	  uint startrec3 = lenrec1 + lenrec2;
	  leaderfile.seekg(startrec3,ios.beg); // map projection data record
	  leaderfile.read((char) rec_seq, sizeb4); // record number
	  rec_seq = ntohl(rec_seq); // Bert Kampes, 07-Apr-2005
	  leaderfile.read((char)&rec_sub1, sizeb1); // first record sub type code
	  leaderfile.read((char)&rec_type, sizeb1); // record type code
	  leaderfile.read((char)&rec_sub2, sizeb1); // second record sub type code
	  leaderfile.read((char)&rec_sub3, sizeb1); // third record sub type code
	  DEBUG.print("ERS:  Expecting record 3 with code {10,20,31,20}");
	  DEBUG.print("RSAT: Expecting record 3 with code {18,60,18,20}");
	  DEBUG.print("RSAT record length is 4096, ERS 1886, but");
	  DEBUG.print("ERS contains more info on zero doppler times, etc.");
	  DEBUG.print("RSAT seems to have that in data file");
	  DEBUG << "rec_seq: " << rec_seq << "; rec_sub1: " << (int)rec_sub1 << "; rec_type: " << (int)rec_type << "; rec_sub2: " << (int)rec_sub2 << "; rec_sub3: " << (int)rec_sub3;
	  DEBUG.print();

	  // ___ Read important information ___
	  leaderfile.seekg(startrec3+8,ios.beg); // map projection data record
	  leaderfile.read((char) lenrec3, sizeb4); // length of record3
	  lenrec3 = ntohl(lenrec3); // bk 6 jul 2000, byteorder x86 machines.
	  DEBUG << "readleader::record 3: start at: " << startrec3 << "; length: " << lenrec3;
	  DEBUG.print();

	  //if (lenrec3 != 1620 ) // quarter scene 1046 (?)
	  if (lenrec3 < 1620) // quarter scene 1046 (?)
		{
		WARNING.print("Probably something wrong here, byte order on x86?");
		WARNING.print("you may also be trying to read leader of ESA.RAW not ESA.SLC.");
		WARNING.print("format (or quarter scene).");
		WARNING.print("We try to do something useful still, but be careful.");
		WARNING.print("Skipping map projection record, not in ESA.RAW leaderfile format.");
		WARNING << "readleader: length of record 3 = \"" << lenrec3 << "\"; expected \"1620\" for ESA SLC (full scene).";
		WARNING.print();
		lenrec3 = 0;
		skipmapprojrecord = true;
		}
	  // --- Trying to read anyway if it is longer ---
	  if (lenrec3 > 1620) // atlantis? rec2 is 4096 not 1886 (?)
		{
		//skipmapprojrecord = true;// no: want to read it sometimes...
		WARNING.print("record length longer than expected.  Atlantis?.  reading map.");
		WARNING.print("seems not ok. assuming mapprojection record not present");
		WARNING.print("trying to read further with next record");
		startrec3 = lenrec1+lenrec2+lenrec3; // try to continue with rec4??? rec3 seems platform data
		lenrec3 = 0; // so next time this record is read, but with platform data format...
		}

	  // skip map projection record if not present ...
	  // BK 18-Jul-2000
	  if (skipmapprojrecord) // not present, do not read.
		{
		WARNING.print("Skipping map projection record. Not found.");
		WARNING.print("Continuing with platform position data.");
		}
	  else // ESA.SLC format, do read map proj. record.
		{
		// ______ Map projection general information ______
		leaderfile.seekg(startrec3+28,ios.beg);
		leaderfile.read((char)&c32projection, sizea32); // map proj. descr.
		c32projection.charAt(32)='\0';
		leaderfile.read((char)&c16numpix, sizei16); // numpixels
		c16numpix.charAt(16)='\0';
		leaderfile.read((char)&c16numlin, sizei16); // numlines
		c16numlin.charAt(16)='\0';


		// ______ Continue ______
		leaderfile.read((char)&c16interpix, sizef16); // dist inter-pixel
		c16interpix.charAt(16)='\0';
		leaderfile.read((char)&c16interlin, sizef16); // dist inter-lines
		c16interlin.charAt(16)='\0';
		leaderfile.read((char)&c16orien, sizef16); // orientation at output
		c16orien.charAt(16)='\0';
		leaderfile.read((char)&c16platincl, sizef16); // actual platform inclination
		c16platincl.charAt(16)='\0';
		leaderfile.read((char)&c16platascn, sizef16); // actual ascending node
		c16platascn.charAt(16)='\0';
		leaderfile.read((char)&c16geocenter, sizef16);
		c16geocenter.charAt(16)='\0';
		leaderfile.read((char)&c16platalt, sizef16); // altitude
		c16platalt.charAt(16)='\0';
		leaderfile.read((char)&c16platgs, sizef16); // ground speed
		c16platgs.charAt(16)='\0';
		leaderfile.read((char)&c16plathead, sizef16); // heading
		c16plathead.charAt(16)='\0';
		leaderfile.read((char)&c32refellips, sizea32); // ellipsoid
		c32refellips.charAt(32)='\0';
		leaderfile.read((char)&c16refmajor, sizef16); // semi major
		c16refmajor.charAt(16)='\0';
		leaderfile.read((char)&c16refminor, sizef16); // semi minor
		c16refminor.charAt(16)='\0';

		// ______ Coordinates of four corner points ______
		leaderfile.seekg(startrec3+1072,ios.beg);
		leaderfile.read((char)&c16lat11, sizef16); // lat. 1st line 1st pix.
		c16lat11.charAt(16)='\0';
		leaderfile.read((char)&c16lon11, sizef16); // lon. 1st line 1st pix.
		c16lon11.charAt(16)='\0';
		leaderfile.read((char)&c16lat1N, sizef16); // lat. 1st line last pix.
		c16lat1N.charAt(16)='\0';
		leaderfile.read((char)&c16lon1N, sizef16); // lon. 1st line last pix.
		c16lon1N.charAt(16)='\0';
		leaderfile.read((char)&c16latNN, sizef16); // lat. last line last pix.
		c16latNN.charAt(16)='\0';
		leaderfile.read((char)&c16lonNN, sizef16); // lon. last line last pix.
		c16lonNN.charAt(16)='\0';
		leaderfile.read((char)&c16latN1, sizef16); // lat. last line 1st pix.
		c16latN1.charAt(16)='\0';
		leaderfile.read((char)&c16lonN1, sizef16); // lon. last line 1st pix.
		c16lonN1.charAt(16)='\0';
	  } // end skip map proj.


	  // --- RECORD 4 ---
	  DEBUG.print("record 4 of leader file.");
	  uint startrec4 =lenrec1+lenrec2+lenrec3;
	  leaderfile.seekg(startrec4+8,ios.beg); // slc platform position data record
	  leaderfile.read((char) lenrec4, sizeb4); // length of record4
	  lenrec4 = ntohl(lenrec4); // bk 6 jul 2000, byteorder x86 machines.
	  DEBUG << "readleader::record 4: start at: " << startrec4 << "; length (variable): " << lenrec4;
	  DEBUG.print();


	// ______ Positional data points ______
	  leaderfile.seekg(startrec4+140,ios.beg);
	  leaderfile.read((char)&c4dummy, sizei4); // number of data points
	  c4dummy.charAt(4)='\0';
	  numdatapoints = Integer.parseInt(c4dummy);
	  leaderfile.read((char)&c4year, sizei4); // year
	  c4year.charAt(4)='\0';
	  leaderfile.read((char)&c4month, sizei4); // month
	  c4month.charAt(4)='\0';
	  leaderfile.read((char)&c4day, sizei4); // day
	  c4day.charAt(4)='\0';
	  leaderfile.read((char)&c4dayofyear, sizei4); // day of year
	  c4dayofyear.charAt(4)='\0';
	  leaderfile.read((char)&c22seconds, sized22); // sec
	  c22seconds.charAt(22)='\0';
	  leaderfile.read((char)&c22interval, sized22); // interval time
	  c22interval.charAt(22)='\0';
	  leaderfile.read((char)&c64rcs, sizea64); // ref. coord. system
	  c64rcs.charAt(64)='\0';
	  leaderfile.read((char)&c22gmha, sized22); // greenwich mean hour angle
	  c22gmha.charAt(22)='\0';
	  leaderfile.read((char)&c16ltposerr, sizef16); // along track pos. error
	  c16ltposerr.charAt(16)='\0';
	  leaderfile.read((char)&c16ctposerr, sizef16); // across track pos. error
	  c16ctposerr.charAt(16)='\0';
	  leaderfile.read((char)&c16rposerr, sizef16); // radial pos. error
	  c16rposerr.charAt(16)='\0';

	// ______ Read statevector of data points ______
	  if (numdatapoints ==0) // test, sometimes on linux problems
		{
		WARNING.print("numdatapoints=0: probalbly something wrong with filepointer");
		WARNING.print("But continuing, be very careful!");
		numdatapoints =1; // arbitrary...
		STATE.resize(6, numdatapoints); // storage statevector
		}
	  else
		{
		STATE.resize(6, numdatapoints); // storage statevector
		for (register int32 i =0;i<numdatapoints;i++) // number of data points
		  {
		  leaderfile.seekg(startrec4+386+i *132,ios.beg); // start ith datarecord
		  for (register int32 j =0;j<6;j++) // read x,y,z,xdot,ydot,zdot
			{
			leaderfile.read((char)&c22dummy, sized22);
			c22dummy.charAt(22)='\0';
			STATE(j,i)=Double.parseDouble(c22dummy);
			}
		  }
		}

	  // --- RECORD 5 (Bianca Cassee) slc facility related data record [general type]) ---
	  DEBUG.print("record 5 of leader file."); //bc 18 dec 2003
	  uint startrec5 =lenrec1+lenrec2+lenrec3+lenrec4; //bc
	  leaderfile.seekg(startrec5+8,ios.beg); //slc facility related
	  leaderfile.read((char) lenrec5, sizeb4); //bc length of record4
	  lenrec5 = ntohl(lenrec5); //byteorder x86 machines.
	  DEBUG << "readleader::record 5: start at: " << startrec5 << "; length: " << lenrec5;
	  DEBUG.print();


	  // ______ Calibration information  ______              //bc
	  leaderfile.seekg(startrec5+582,ios.beg); //bc
	  leaderfile.read((char)&c16incangle1strange, sizef16); //bc
	  c16incangle1strange.charAt(16)='\0'; //bc
	  leaderfile.read((char)&c16incanglecenrange, sizef16); //bc
	  c16incanglecenrange.charAt(16)='\0'; // gk bc
	  leaderfile.read((char)&c16incanglelstrange, sizef16); //bc
	  c16incanglelstrange.charAt(16)='\0'; //bc
	  leaderfile.seekg(startrec5+662,ios.beg); //gk
	  leaderfile.read((char)&calK, sizef16); //gk
	  calK.charAt(16)='\0'; //gk
	  leaderfile.seekg(startrec5+566,ios.beg); //gk
	  leaderfile.read((char)&repplspwr, sizef16); //gk
	  repplspwr.charAt(16)='\0'; //gk
	  // ___ BK: number of invalid samples at end (DPAF/UKPAF problem) ___
	  String c4SWSTflag = new String(new char[5]); //bk
	  String c4SWSTchange = new String(new char[5]); //bk
	  String c4missingrawlines = new String(new char[5]); //bk
	  String c4validperline = new String(new char[5]); //bk
	  leaderfile.seekg(startrec5+98,ios.beg);
	  leaderfile.read((char)&c4SWSTflag, sizei4); // numsamples
	  c4SWSTflag.charAt(4)='\0';
	  if (Integer.parseInt(c4SWSTflag) != 0)
		WARNING.print("SWST was not constant");
	  leaderfile.seekg(startrec5+138,ios.beg);
	  leaderfile.read((char)&c4SWSTchange, sizei4); // numsamples
	  c4SWSTchange.charAt(4)='\0';
	  if (Integer.parseInt(c4SWSTchange) != 0)
		{
		WARNING << "c4SWSTchange: " << c4SWSTchange;
		WARNING.print();
		}
	  leaderfile.seekg(startrec5+142,ios.beg);
	  leaderfile.read((char)&c4missingrawlines, sizei4); // numsamples
	  c4missingrawlines.charAt(4)='\0';
	  if (Integer.parseInt(c4missingrawlines) != 0)
		{
		WARNING << "c4missingrawlines: " << c4missingrawlines;
		WARNING.print();
		}
	  leaderfile.seekg(startrec5+1722,ios.beg);
	  leaderfile.read((char)&c4validperline, sizei4); // numsamples
	  c4validperline.charAt(4)='\0';
	  DEBUG << "c4validperline: " << c4validperline;
	  DEBUG.print();

	//  
	//  leaderfile.seekg(startrec5+XXX,ios::beg);
	//  leaderfile.read((char*)&c4numvalid,sizei4);           // numsamples
	//  c4numvalid[4]='\0';
	//  leaderfile.read((char*)&c4numinvalid,sizei4);         // zero padded at end
	//  c4numinvalid[4]='\0';
	//  if (atoi(c4numinvalid) != 0)
	//    {
	//    WARNING << "Number of invalid samples " << c4numinvalid 
	//            << " not 0: rsr may be wrongly computed.";
	//    WARNING.print();
	//    }
	//  

	  //  // --- RECORD 6 (slc facility related data record [pcs type]) ---
	  //  DEBUG.print("record 6 of leader file.");                 //gk 28 jan 2004
	  //  uint startrec6=lenrec1+lenrec2+lenrec3+lenrec4+lenrec5;  //gk
	  //  leaderfile.seekg(startrec6+8,ios::beg);       //slc facility related 
	  //  leaderfile.read((char*)&lenrec6,sizeb4);      //gk length of record5
	  //  lenrec6 = ntohl(lenrec6);                     //byteorder x86 machines.
	  //  DEBUG << "readleader::record 6: start at: " << startrec6 
	  //        << "; length: " << lenrec6;
	  //  DEBUG.print();

	  } // endif ERS/RSAT
	else //RSAT method specified
	  {
	  skipmapprojrecord = true; // RSAT does not have this here
	  // --- RECORD 3 (skip) ---
	  DEBUG.print("record 3 of RSAT leader file (data quality).");
	  uint startrec3 = lenrec1 + lenrec2;
	  leaderfile.seekg(startrec3,ios.beg);
	  leaderfile.read((char) rec_seq, sizeb4); // record number
	  rec_seq = ntohl(rec_seq); // Bert Kampes, 07-Apr-2005
	  leaderfile.read((char)&rec_sub1, sizeb1); // first record sub type code
	  leaderfile.read((char)&rec_type, sizeb1); // record type code
	  leaderfile.read((char)&rec_sub2, sizeb1); // second record sub type code
	  leaderfile.read((char)&rec_sub3, sizeb1); // third record sub type code
	  DEBUG.print("ERS:  Expecting record 3 with code {10,20,31,20}");
	  DEBUG.print("RSAT: Expecting record 3 with code {18,60,18,20}");
	  DEBUG.print("RSAT record length should be 1620, ERS 1620");
	  DEBUG << "rec_seq: " << rec_seq << "; rec_sub1: " << (int)rec_sub1 << "; rec_type: " << (int)rec_type << "; rec_sub2: " << (int)rec_sub2 << "; rec_sub3: " << (int)rec_sub3;
	  DEBUG.print();
	  leaderfile.seekg(startrec3+8,ios.beg);
	  leaderfile.read((char) lenrec3, sizeb4); // length of record3
	  lenrec3 = ntohl(lenrec3); // bk 6 jul 2000, byteorder x86 machines.
	  DEBUG << "readleader::record 3: start at: " << startrec3 << "; length: " << lenrec3;
	  DEBUG.print();
	  if ((int)rec_sub1 ==18 && (int)rec_type ==60 && (int)rec_sub2 ==18 && (int)rec_sub3 ==20)
		DEBUG.print("This is the expected RSAT record with code {18,60,18,20}");
	  else
		WARNING.print("This is NOT the expected 3rd RSAT record with code {18,60,18,20}");


	  // --- RECORD 4 (skip) ---
	  DEBUG.print("record 4 of RSAT leader file (data hist).");
	  uint startrec4 = lenrec1 + lenrec2 + lenrec3;
	  leaderfile.seekg(startrec4,ios.beg);
	  leaderfile.read((char) rec_seq, sizeb4); // record number
	  rec_seq = ntohl(rec_seq); // Bert Kampes, 07-Apr-2005
	  leaderfile.read((char)&rec_sub1, sizeb1); // first record sub type code
	  leaderfile.read((char)&rec_type, sizeb1); // record type code
	  leaderfile.read((char)&rec_sub2, sizeb1); // second record sub type code
	  leaderfile.read((char)&rec_sub3, sizeb1); // third record sub type code
	  DEBUG.print("ERS:  Expecting record 4 with code {10,30,31,20}");
	  DEBUG.print("RSAT: Expecting record 4 with code {18,70,18,20}");
	  DEBUG.print("RSAT record length should be 16920");
	  DEBUG << "rec_seq: " << rec_seq << "; rec_sub1: " << (int)rec_sub1 << "; rec_type: " << (int)rec_type << "; rec_sub2: " << (int)rec_sub2 << "; rec_sub3: " << (int)rec_sub3;
	  DEBUG.print();
	  leaderfile.seekg(startrec4+8,ios.beg);
	  leaderfile.read((char) lenrec4, sizeb4); // length of record4
	  lenrec4 = ntohl(lenrec4); // bk 6 jul 2000, byteorder x86 machines.
	  DEBUG << "readleader::record 4: start at: " << startrec4 << "; length: " << lenrec4;
	  DEBUG.print();
	  if ((int)rec_sub1 ==18 && (int)rec_type ==70 && (int)rec_sub2 ==18 && (int)rec_sub3 ==20)
		DEBUG.print("This is the expected RSAT record with code {18,70,18,20}");
	  else
		WARNING.print("This is NOT the expected 4th RSAT record with code {18,70,18,20}");


	  // --- RECORD 5 (skip) ---
	  DEBUG.print("record 5 of RSAT leader file (data hist).");
	  uint startrec5 = lenrec1 + lenrec2 + lenrec3 + lenrec4;
	  leaderfile.seekg(startrec5,ios.beg);
	  leaderfile.read((char) rec_seq, sizeb4); // record number
	  rec_seq = ntohl(rec_seq); // Bert Kampes, 07-Apr-2005
	  leaderfile.read((char)&rec_sub1, sizeb1); // first record sub type code
	  leaderfile.read((char)&rec_type, sizeb1); // record type code
	  leaderfile.read((char)&rec_sub2, sizeb1); // second record sub type code
	  leaderfile.read((char)&rec_sub3, sizeb1); // third record sub type code
	  DEBUG.print("ERS:  Expecting record 5 with code {10,200,31,50}");
	  DEBUG.print("RSAT: Expecting record 5 with code {18,70,18,20}");
	  DEBUG.print("RSAT record length should be 16920");
	  DEBUG << "rec_seq: " << rec_seq << "; rec_sub1: " << (int)rec_sub1 << "; rec_type: " << (int)rec_type << "; rec_sub2: " << (int)rec_sub2 << "; rec_sub3: " << (int)rec_sub3;
	  DEBUG.print();
	  leaderfile.seekg(startrec5+8,ios.beg);
	  leaderfile.read((char) lenrec5, sizeb4); // length of record5
	  lenrec5 = ntohl(lenrec5); // bk 6 jul 2000, byteorder x86 machines.
	  DEBUG << "readleader::record 5: start at: " << startrec5 << "; length: " << lenrec5;
	  DEBUG.print();
	  if ((int)rec_sub1 ==18 && (int)rec_type ==70 && (int)rec_sub2 ==18 && (int)rec_sub3 ==20)
		DEBUG.print("This is the expected RSAT record with code {18,70,18,20}");
	  else
		WARNING.print("This is NOT the expected 5th RSAT record with code {18,70,18,20}");


	  // --- RECORD 6 (read a bit) ---
	  DEBUG.print("record 6 of RSAT leader file (detailed proc. param.).");
	  uint startrec6 = lenrec1 + lenrec2 + lenrec3 + lenrec4 + lenrec5;
	  leaderfile.seekg(startrec6,ios.beg);
	  leaderfile.read((char) rec_seq, sizeb4); // record number
	  rec_seq = ntohl(rec_seq); // Bert Kampes, 07-Apr-2005
	  leaderfile.read((char)&rec_sub1, sizeb1); // first record sub type code
	  leaderfile.read((char)&rec_type, sizeb1); // record type code
	  leaderfile.read((char)&rec_sub2, sizeb1); // second record sub type code
	  leaderfile.read((char)&rec_sub3, sizeb1); // third record sub type code
	  DEBUG.print("ERS:  Expecting record 6 with code {10,200,31,50}");
	  DEBUG.print("RSAT: Expecting record 6 with code {18,120,18,20}");
	  DEBUG.print("RSAT record length should be 7726");
	  DEBUG << "rec_seq: " << rec_seq << "; rec_sub1: " << (int)rec_sub1 << "; rec_type: " << (int)rec_type << "; rec_sub2: " << (int)rec_sub2 << "; rec_sub3: " << (int)rec_sub3;
	  DEBUG.print();
	  leaderfile.seekg(startrec6+8,ios.beg);
	  leaderfile.read((char) lenrec6, sizeb4); // length of record5
	  lenrec6 = ntohl(lenrec6); // bk 6 jul 2000, byteorder x86 machines.
	  DEBUG << "readleader::record 6: start at: " << startrec6 << "; length: " << lenrec6;
	  DEBUG.print();
	  if ((int)rec_sub1 ==18 && (int)rec_type ==120 && (int)rec_sub2 ==18 && (int)rec_sub3 ==20)
		DEBUG.print("This is the expected RSAT record with code {18,120,18,20}");
	  else
		WARNING.print("This is NOT the expected 6th RSAT record with code {18,120,18,20}");

	  DEBUG.print("todo: read asc/desc, BEAM, etc.");
	  // TODO: Annex B-11: LEA record 6
	  // if center lat not known, it seems here too, ignore them


	  // --- RECORD 7 (read platform data) ---
	  DEBUG.print("record 7 of RSAT leader file (platform pos. data).");
	  DEBUG.print("this is exact same as for ERS, but in later record.");
	  uint startrec7 = lenrec1 + lenrec2 + lenrec3 + lenrec4 + lenrec5 + lenrec6;
	  leaderfile.seekg(startrec7,ios.beg);
	  leaderfile.read((char) rec_seq, sizeb4); // record number
	  rec_seq = ntohl(rec_seq); // Bert Kampes, 07-Apr-2005
	  leaderfile.read((char)&rec_sub1, sizeb1); // first record sub type code
	  leaderfile.read((char)&rec_type, sizeb1); // record type code
	  leaderfile.read((char)&rec_sub2, sizeb1); // second record sub type code
	  leaderfile.read((char)&rec_sub3, sizeb1); // third record sub type code
	  DEBUG.print("RSAT: Expecting record 7 with code {18,30,18,20}");
	  DEBUG.print("RSAT record length should be 8960");
	  DEBUG << "rec_seq: " << rec_seq << "; rec_sub1: " << (int)rec_sub1 << "; rec_type: " << (int)rec_type << "; rec_sub2: " << (int)rec_sub2 << "; rec_sub3: " << (int)rec_sub3;
	  DEBUG.print();
	  leaderfile.seekg(startrec7+8,ios.beg);
	  leaderfile.read((char) lenrec7, sizeb4); // length of record5
	  lenrec7 = ntohl(lenrec7); // bk 6 jul 2000, byteorder x86 machines.
	  DEBUG << "readleader::record 7: start at: " << startrec7 << "; length: " << lenrec7;
	  DEBUG.print();
	  if ((int)rec_sub1 ==18 && (int)rec_type ==30 && (int)rec_sub2 ==18 && (int)rec_sub3 ==20)
		DEBUG.print("This is the expected RSAT record with code {18,30,18,20}");
	  else
		WARNING.print("This is NOT the expected 7th RSAT record with code {18,30,18,20}");
	  // ______ Keppler elements may be useful to get some sort of orbit data points ______
	  // ______ by evaluating them myself, in earth fixed system ______
	  String c16semimajor = new String(new char[17]);
	  String c16inclination = new String(new char[17]);
	  String c16eccentricity = new String(new char[17]);
	  String c16argofperi = new String(new char[17]);
	  String c16lonofnode = new String(new char[17]);
	  String c16meananomaly = new String(new char[17]);
	  // Modified by LG for reading JERS Fine
	  if (readfiles_arg.argvalue.sensor_id == SLC_ALOS)
					startrec7 = startrec3;
	  leaderfile.seekg(startrec7+44,ios.beg);
	  leaderfile.read((char)&c16semimajor, sizef16);
	  c16semimajor.charAt(16)='\0';
	  leaderfile.read((char)&c16inclination, sizef16);
	  c16inclination.charAt(16)='\0';
	  leaderfile.read((char)&c16eccentricity, sizef16);
	  c16eccentricity.charAt(16)='\0';
	  leaderfile.read((char)&c16argofperi, sizef16);
	  c16argofperi.charAt(16)='\0';
	  leaderfile.read((char)&c16lonofnode, sizef16);
	  c16lonofnode.charAt(16)='\0';
	  leaderfile.read((char)&c16meananomaly, sizef16);
	  c16meananomaly.charAt(16)='\0';
	  // Modified by LG for reading JERS Fine
	  if (readfiles_arg.argvalue.sensor_id != SLC_ALOS)
	  {
		INFO.print("Orbit Kepplerian elements follow:");
		INFO << "c16semimajor [km]:     " << c16semimajor;
		INFO.print();
		INFO << "c16inclination [rad]:  " << c16inclination;
		INFO.print();
		INFO << "c16eccentricity [-]:   " << c16eccentricity;
		INFO.print();
		INFO << "c16argofperi [rad]:    " << c16argofperi;
		INFO.print();
		INFO << "c16lonofnode [rad]:    " << c16lonofnode;
		INFO.print();
		INFO << "c16meananomaly [rad]:  " << c16meananomaly;
		INFO.print();
	  }
	  //real8 semimajor = atof(c16semimajor);

	  // ______ Positional data points ______
	  leaderfile.seekg(startrec7+140,ios.beg);
	  leaderfile.read((char)&c4dummy, sizei4); // number of data points
	  c4dummy.charAt(4)='\0';
	  numdatapoints = Integer.parseInt(c4dummy);
	  leaderfile.read((char)&c4year, sizei4); // year of first data point
	  c4year.charAt(4)='\0';
	  leaderfile.read((char)&c4month, sizei4); // month of first data point
	  c4month.charAt(4)='\0';
	  leaderfile.read((char)&c4day, sizei4); // day of first data point
	  c4day.charAt(4)='\0';
	  leaderfile.read((char)&c4dayofyear, sizei4); // daynumber of year
	  c4dayofyear.charAt(4)='\0';
	  leaderfile.read((char)&c22seconds, sized22); // sec of day of first point
	  c22seconds.charAt(22)='\0';
	  leaderfile.read((char)&c22interval, sized22); // interval time
	  c22interval.charAt(22)='\0';
	  leaderfile.read((char)&c64rcs, sizea64); // ref. coord. system
	  c64rcs.charAt(64)='\0';
	  leaderfile.read((char)&c22gmha, sized22); // greenwich mean hour angle
	  c22gmha.charAt(22)='\0';
	  leaderfile.read((char)&c16ltposerr, sizef16); // along track pos. error
	  c16ltposerr.charAt(16)='\0';
	  leaderfile.read((char)&c16ctposerr, sizef16); // across track pos. error
	  c16ctposerr.charAt(16)='\0';
	  leaderfile.read((char)&c16rposerr, sizef16); // radial pos. error
	  c16rposerr.charAt(16)='\0';

	  // --- Maybe Julian Day is useful for conversion ---
	  INFO << "Seconds of day of first data point: " << c22seconds;
	  INFO.print();
	  INFO << "Day of data:        " << c4day;
	  INFO.print();
	  INFO << "Month of data:      " << c4month;
	  INFO.print();
	  INFO << "Year of data:       " << c4year;
	  INFO.print();
	  int32 jd_statevector = julday(Integer.parseInt(c4day), Integer.parseInt(c4month), Integer.parseInt(c4year));
	  INFO << "Julian Day of data: " << jd_statevector;
	  INFO.print();


	  // ______ Read statevector of data points ______
	  STATE.resize(6, numdatapoints); // storage statevector
	  for (register int32 i =0;i<numdatapoints;i++) // number of data points
		{
		leaderfile.seekg(startrec7+386+i *132,ios.beg); // start ith datarecord
		for (register int32 j =0;j<6;j++) // read x,y,z,xdot,ydot,zdot
		  {
		  leaderfile.read((char)&c22dummy, sized22);
		  c22dummy.charAt(22)='\0';
		  STATE(j,i)=Double.parseDouble(c22dummy);
		  }
		}
	  // --- if this is an inertial system (rsat) we need to convert to earth fixed ---
	  // --- Moreover: we want to reduce the number of data points to local arc ------- 
	  WARNING.print("Convert orbit data to earth fixed (please check this).");
	  if (!(strcmp(c64rcs,"INERTIAL")))
		INFO.print("Inertial system for orbit: transforming to Earth fixed.");
	  else
		WARNING.print("system for orbit: transforming to Earth fixed (?).");
	  // This angle refers to ascending node of equator crossing (time), ie. where Z=0
	  // not to first state vectors time!  so following is wrong... ?
	  // note that we have also the Keppler elements, etc.
	  real8 csi2cse = -1.0; // 1 for initial to earth fixed, -1 for back transform
	  real8 GMHA = deg2rad(Double.parseDouble(c22gmha)); // Greenwich Mean Hour Angle in [deg]
	  real8 earthrot = 7.292115856E-5; // [rad/s] 2pi rad in 23h56'04.091" (siderial day)
	  // can we simply use the given Greenwich mean angle to rotate around Z?
	  // or do we need to take precesion.nutation polarwobble into account?
	  // to what (time) does the annotated GMHA refer?
	  DEBUG << "GMHA [rad]: " << GMHA;
	  DEBUG.print();
	  DEBUG << "Convertion from inertial to earth fixed [1/-1]: " << csi2cse;
	  DEBUG.print();
	  DEBUG << "earthrot [rad/s]: " << earthrot;
	  DEBUG.print();
	  // --- Create a new state vector matrix ---
	  // --- these computation could be checked by using TIEPOINT, the center point -----
	  // --- Limit to only 2 points before, 2 after scene, and use polyfit(3) interp. ---
	  // --- Still this seems to have a large error... ---
	  STATE_INERTIAL = STATE; // copy the Z values
	  for (register int32 i =0;i<numdatapoints;i++) // number of data points
		{
		real8 dt = new real8(i); // since annotated GMHA???
		real8 angle = csi2cse*(GMHA+earthrot *dt); // current angle of Greenwich
		DEBUG << "current angle for this data point [deg]: " << rad2deg(angle);
		DEBUG.print();
		STATE(0,i) = Math.cos(angle)*STATE_INERTIAL(0,i)-Math.sin(angle)*STATE_INERTIAL(1,i); // x
		STATE(1,i) = Math.sin(angle)*STATE_INERTIAL(0,i)+Math.cos(angle)*STATE_INERTIAL(1,i); // y
		//STATE(2,i)     = STATE_INERTIAL(2,i);// Z is same
		STATE(3,i) = -1.0 *csi2cse *earthrot* (Math.sin(angle)*STATE_INERTIAL(0,i)+Math.cos(angle)*STATE_INERTIAL(1,i))+ Math.cos(angle)*STATE_INERTIAL(3,i)-Math.sin(angle)*STATE_INERTIAL(4,i); // xdot
		STATE(4,i) = csi2cse *earthrot* (Math.cos(angle)*STATE_INERTIAL(0,i)-Math.sin(angle)*STATE_INERTIAL(1,i))+ Math.sin(angle)*STATE_INERTIAL(3,i)+Math.cos(angle)*STATE_INERTIAL(4,i); // ydot
		//STATE(5,i)     = STATE_INERTIAL(5,i);// Zdot is same ?
		}
	  INFO.print("info on fast/slow time is obtained in readdat, not here.");
	  WARNING.print("todo: reduce the number of orbital data to arc around area only");
	  // --- Limit the data points to 2 points before/after: 32 minutes ---
	  // --- better yet, interpolate using keppler to 5 points 5 seconds dt over take ---
	}

	  leaderfile.close();



	// ====== Compute prf and rsr and use these (based on data is better) ======
	  // Modified by LG for reading ALOS Fine
	  if ((readfiles_arg.argvalue.sar_processor==SARPR_ATL) || (readfiles_arg.argvalue.sar_processor == SARPR_JAX))
		skipmapprojrecord = true;
	  real8 prfcomputed; // [Hz] written to result file
	  real8 rsrcomputed; // [MHz] written to result file
	  if (skipmapprojrecord ==true) // not present in old? ESA.SLC versions?
		{
		c32projection = "skipped";
		c16numpix = "skipped";
		c16numlin = "skipped";
		//strcpy(c16numpix,     c16scenewidth);
		//strcpy(c16numlin,     c16scenelength);
		c16interpix = "skipped";
		c16interlin = "skipped";
		c16orien = "skipped";
		c16platincl = "skipped";
		c16platascn = "skipped";
		c16geocenter = "skipped";
		c16platalt = "skipped";
		c16platgs = "skipped";
		c16plathead = "skipped";
		c32refellips = "skipped";
		c16refmajor = "skipped";
		c16refminor = "skipped";
		c16lat11 = "skipped";
		c16lon11 = "skipped";
		c16lat1N = "skipped";
		c16lon1N = "skipped";
		c16latNN = "skipped";
		c16lonNN = "skipped";
		c16latN1 = "skipped";
		c16lonN1 = "skipped";
		WARNING.print("Using nominal values for PRF, RSR, not computed.");
		WARNING.print(" + because map projection record not present(?).");
		prfcomputed = Double.parseDouble(c16prf); // Hz
		// start_added_by_don
		if (readfiles_arg.argvalue.sensor_id==SLC_ALOS)
		  prfcomputed /= 1e3;
		// end_added_by_don
		rsrcomputed = Double.parseDouble(c16samplingrate); // MHz
		}
	  else
		{
		// ______ Check prf ______
		tm tijdstart;
		tm tijdend;
		RefObject<tm> TempRefObject2 = new RefObject<tm>(tijdstart);
		strptime(c24zd1stazitime, "%d-%b-%Y %T", TempRefObject2);
		tijdstart = TempRefObject2.argvalue;
		RefObject<tm> TempRefObject3 = new RefObject<tm>(tijdend);
		strptime(c24zdlstazitime, "%d-%b-%Y %T", TempRefObject3);
		tijdend = TempRefObject3.argvalue;
		int32 index =0;
		while (c24zd1stazitime.charAt(index) != '.')
		  ++index;
		String fracsec1st = new String(new char[25]);
		fracsec1st = c24zd1stazitime.charAt(index+1);
		index =0;
		while (c24zdlstazitime.charAt(index) != '.')
		  ++index;
		String fracseclst = new String(new char[25]);
		fracseclst = c24zdlstazitime.charAt(index+1);
		final real8 ta1 = Double.parseDouble(fracsec1st)/1000. + real8(tijdstart.tm_sec + 60 *tijdstart.tm_min + 3600 *tijdstart.tm_hour);
		final real8 taN = Double.parseDouble(fracseclst)/1000. + real8(tijdend.tm_sec + 60 *tijdend.tm_min + 3600 *tijdend.tm_hour);
		// BK 28-Sep-2000: numlin -1 !
		prfcomputed = (Double.parseDouble(c16numlin) - 1.) / (taN-ta1);

		// ______ compute rsr ______
		//  const real8 tr1     = atof(c16zd1strange);              // 2way ms
		//  const real8 trN     = atof(c16zdlstrange);              // 2way ms
		// BK 28-Sep-2000: numpix -1 !
		rsrcomputed = 0.001 * (Double.parseDouble(c16numpix)-1.0) / (Double.parseDouble(c16zdlstrange)-Double.parseDouble(c16zd1strange)); // MHz
		} // else skipped map projection record

	  // BK 28-Oct-2003, for Atlantis processor
	  // ___ Check if rsr is in MHz, assume about 20 MHz ---
	  //if (rsrcomputed>10000000.0 && rsrcomputed<30000000.0)
	  if (rsrcomputed>10000.0)
		{
		WARNING.print("Adapting RSR from Hz to MHz.");
		WARNING << "Old value rsrcomputed: " << rsrcomputed << "; new value rsrcomputed: " << rsrcomputed/1e6;
		WARNING.print();
		rsrcomputed/=1e6;
		}
	  real8 rangebandwidth = Double.parseDouble(c16bandrangetot);
	  // start_added_by_don
	  if (readfiles_arg.argvalue.sensor_id==SLC_ALOS)
		rangebandwidth /= 1e3;
	  // end_added_by_don
	  if (rangebandwidth < 1.0) // e.g. JERS not annotated.
		{
		WARNING.print("computed rangebandwidth < 1.0 ignoring this.");
		rangebandwidth = rsrcomputed; // set this to no oversampling
		}
	  //if (rangebandwidth>10000000.0 && rangebandwidth<30000000.0)
	  if (rangebandwidth>10000.0)
		{
		WARNING.print("Adapting RBW from Hz to MHz.");
		WARNING << "Old value rangebandwidth: " << rangebandwidth << "; new value rangebandwidth: " << rangebandwidth/1e6;
		WARNING.print();
		rangebandwidth/=1e6;
		}


	// ====== Write information to scratchfiles ======
	  ofstream scratchlogfile = new ofstream("scratchloglea", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchlogfile, "readleader: scratchloglea", __FILE__, __LINE__);

	  scratchlogfile << "\n\n*******************************************************************" << "\n* EXTRACTED DATA FROM LEADER FILE: " << readfiles_arg.argvalue.leaderfile << " *" << "\n*******************************************************************" << "\n\nFile descriptor record" << "\n----------------------" << "\nFile name: " << c16leafilename << "\n\nSLC data set summary record: scene parameters" << "\n---------------------------------------------" << "\nScene reference number: \t\t\t\t" << c32sceneref << "\nScene centre time (UTC) <YYYYMMDDhhmmssttt>: \t\t" << c32scenetime << "\nProcessed scene centre geodetic latitude" << "\n +(positive for North latitude) (degrees): \t\t" << c16centerlat << "\nProcessed scene centre longitude" << "\n +(negative for West longitude) (degrees): \t\t" << c16centerlon << "\nProcessed scene centre true heading (degrees): \t\t" << c16centerheading << "\nEllipsoid designator: \t\t\t\t\t" << c16ellipsoid << "\nEllipsoid semimajor axis (km): \t\t\t\t" << c16semimajor << "\nEllipsoid semiminor axis (km): \t\t\t\t" << c16semiminor << "\nEarth mass times grav. const. (M.G) (kg.m/s^2): \t" << c16GM << "\nEllipsoid J2 parameter: \t\t\t\t" << c16J2 << "\nEllipsoid J3 parameter: \t\t\t\t" << c16J3 << "\nEllipsoid J4 parameter: \t\t\t\t" << c16J4 << "\nScene center line number (the line number at " << "\n +the scene centre including zero fill): \t\t" << scenecenterline << "\nScene center pixel number (the pixel number at" << "\n +the scene centre including zero fill): \t\t" << scenecenterpixel << "\nProcessed scene length incl. zero fill (km): \t\t" << c16scenelength << "\nProcessed scene width incl. zero fill (km): \t\t" << c16scenewidth << "\n\nSLC data set summary record: general mission/sensor parameters" << "\n--------------------------------------------------------------" << "\nSensor platform mission identifier: \t\t\t" << c16platformid << "\nSensor ID, mode of op. for this channel: \t\t" << c32sensorid << "\nOrbit number: \t\t\t\t\t\t\t" << c8orbitnr << "\nSensor platform geodetic latitude at nadir" << "\n +corresponding to scene centre" << "\n +(positive for North latitude) (degrees): \t\t\t" << c8platformlat << "\nSensor platform geodetic longitude at nadir" << "\n +corresponding to scene centre" << "\n +(negative for West longitude) (degrees): \t\t\t" << c8platformlon << "\nSensor platform heading at nadir" << "\n +corresponding to scene centre" << "\n +(clockwise positive from North) (degrees): \t\t\t" << c8platformheading << "\nSensor clock angle as measured relative to" << "\n +sensor platform flight direction (degrees): \t" << c8clockangle << "\nIncidence angle at scene centre as derived" << "\n +from sensor platform orientation (degrees): \t" << c8incidence << "\nRadar frequency (GHz): \t\t\t\t\t\t" << c8freq << "\nRadar wavelength (meters): \t\t\t" << c16wavelength << "\nMotion compensator identifier" << "\n +(00=no, 01=on board, 10=in processor, 11=both): \t" << c2motioncomp << "\nRange pulse specifier: \t\t\t\t\t" << c16pulse << "\nNominal range pulse (chirp) amplitude coefficient," << "\n +constant term: \t\t\t\t\t" << c16ampconst << "\n +linear term (sec-1): \t\t\t\t\t" << c16amplinear << "\n +quadratic term (sec-2): \t\t\t\t" << c16ampquadratic << "\n +cubic term (sec-3): \t\t\t\t\t" << c16ampcubic << "\n +quartic term (sec-4): \t\t\t\t" << c16ampquartic << "\nNominal range pulse (chirp) phase coefficient," << "\n +constant term (cycles): \t\t\t\t" << c16phaseconst << "\n +linear term (Hz): \t\t\t\t\t" << c16phaselinear << "\n +quadratic term (Hz/sec): \t\t\t\t" << c16phasequadratic << "\n +cubic term (Hz/sec2): \t\t\t\t" << c16phasecubic << "\n +quartic term (Hz/sec3): \t\t\t\t" << c16phasequartic << "\nDown linked chirp extraction index (samples): \t\t\t" << c8extindex << "\nRange sampling rate (MHz): \t\t\t\t" << c16samplingrate << "\nRange gate delay at early ege (in time)" << "\n +at the start of the image (microsec): \t" << c16rangedelay << "\nRange pulse length (microsec): \t\t\t\t" << c16ranpulselen << "\nBase band conversion flag: \t\t\t" << c4conversion << "\nRange compressed flag (YES=range compressed data): \t" << c4compression << "\nQuantization per channel I & Q" << "\n +(5I 5Q/6I 6Q for OGRB/OBRC) (bits): \t\t\t\t" << c8qperch << "\nQuantizer descriptor: \t\t\t\t\t" << c12qdesc << "\nDC Bias for I-component (actual value): \t\t" << c16dci << "\nDC Bias for Q-component (actual value): \t\t" << c16dcq << "\nGain imbalance for I & Q (actual value): \t\t" << c16imbalance << "\nAntenna mechanical boresight angle" << "\n +relative to platform vertical axis: \t\t\t" << c16boresight << "\nPulse Repetition Frequency (PRF) (actual value): \t" << c16prf << "\n\nSLC data set summary record: sensor specific parameters" << "\n-------------------------------------------------------" << "\nSatellite encoded binary time code: \t\t\t" << c16sattimecode << "\nSatellite clock time (UTC) (YYYYMMDDhhmmssttt): \t" << c32sattime << "\nSatellite clock step length (nanosec): \t\t\t\t" << c8satclockstep << "\n\nSLC data set summary record: general processing parameters" << "\n----------------------------------------------------------" << "\nProcessing facility identifier: \t\t\t" << c16facilityid << "\nProcessing system identifier: \t\t\t\t" << c8systemid << "\nProcessing version identifier: \t\t\t\t" << c8versionid << "\nProduct type specifier: \t\t\t\t" << c32typespec << "\nProcessing algorithm identifier: \t\t\t" << c32algid << "\nNominal number of looks processed in azimuth (looks): \t" << c16looksazi << "\nNominal number of looks processed in range (looks): \t" << c16looksrange << "\nBandwidth per look in azimuth (null-to-null) (Hz): \t" << c16bandazi << "\nBandwidth per look in range (MHz): \t\t\t" << c16bandrange << "\nTotal processor bandwidth in azimuth (Hz): \t\t" << c16bandazitot << "\nTotal processor bandwidth in range (MHz): \t\t" << c16bandrangetot << "\nWeighting function designator in azimuth: \t\t" << c32weightazi << "\nWeighting function designator in range: \t\t" << c32weightrange << "\nData input source: \t\t\t\t\t" << c16inputsource << "\nNominal resolution in range (3-dB width) (m): \t\t" << c16resrange << "\nNominal resolution in azimuth (3-dB width) (m): \t" << c16resazi << "\nAlong track Doppler frequency centroid" << "\n +at early edge of image" << "\n +constant term (Hz): \t\t\t\t\t" << c16atdoppcconst << "\n +linear term (Hz/sec): \t\t\t\t" << c16atdoppclinear << "\n +quadratic term (Hz/sec/sec): \t\t\t\t" << c16atdoppcquadratic << "\nCross track Doppler frequency centroid" << "\n +at early edge of image" << "\n +constant term (Doppler centroid) (Hz): \t\t" << c16xtdoppcconst << "\n +linear term (Slope of Doppler centroid) (Hz/sec): \t" << c16xtdoppclinear << "\n +quadratic term (Hz/sec/sec): \t\t\t\t" << c16xtdoppcquadratic << "\nTime direction indicator along pixel direction: \t" << c8timepix << "\nTime direction indicator along line direction: \t\t" << c8timeline << "\nAlong track Doppler frequency rate" << "\n +at early edge of image" << "\n +constant term (Hz/sec): \t\t\t\t" << c16atdopprconst << "\n +linear term (Hz/sec/sec): \t\t\t\t" << c16atdopprlinear << "\n +quadratic term (Hz/sec/sec/sec): \t\t\t" << c16atdopprquadratic << "\nCross track Doppler frequency rate" << "\n +at early edge of image" << "\n +constant term (Azimuth FM rate) (Hz/sec): \t\t" << c16xtdopprconst << "\n +linear term (Slope of Azimuth FM rate) (Hz/sec/sec): \t" << c16xtdopprlinear << "\n +quadratic term (Hz/sec/sec/sec): \t\t\t" << c16xtdopprquadratic << "\nLine content indicator: \t\t\t\t" << c8linecontent << "\nClutterlock applied flag: \t\t\t\t" << c4clutterlock << "\nAutofocussing applied flag: \t\t\t\t" << c4autofocus << "\nLine spacing (m): \t\t\t\t\t" << c16linespace << "\nPixel spacing (in range) (m): \t\t\t\t" << c16pixspace << "\nProcessor range compression designator: \t\t" << c16rcompdes << "\n\nSLC data set summary record: sensor specific local use segment" << "\n--------------------------------------------------------------" << "\nZero-doppler range time (two-way)" << "\n +of first range pixel (millisec): \t\t\t" << c16zd1strange << "\n +of centre range pixel (millisec): \t\t\t" << c16zdcenrange << "\n +of last range pixel (millisec): \t\t\t" << c16zdlstrange << "\nZero-doppler azimuth time" << "\n +of first azimuth pixel (UTC): \t\t" << c24zd1stazitime << "\n +of centre azimuth pixel (UTC): \t\t" << c24zdcenazitime << "\n +of last azimuth pixel (UTC): \t\t\t" << c24zdlstazitime << "\n\nMap projection data record: general information" << "\n-----------------------------------------------" << "\nMap projection descriptor: \t\t\t\t" << c32projection << "\nNumber of pixels per line of image: \t\t\t" << c16numpix << "\nNumber of lines: \t\t\t\t\t" << c16numlin << "\nNominal inter-pixel distance in output scene (m): \t" << c16interpix << "\nNominal inter-line distance in output scene (m): \t" << c16interlin << "\nOrientation at output scene centre" << "\n +[for geodetic products this is simply the" << "\n +convergence of the meridians, ie: the angle" << "\n +between geographic north and map grid north" << "\n +(angle of projection axis from true North)] (degrees): " << c16orien << "\nActual platform orbital inclination (degrees): \t\t" << c16platincl << "\nActual ascending node (longitude at equator) (degrees): " << c16platascn << "\nGeocentre to platform distance at input scene centre: \t" << c16geocenter << "\nPlatform geodetic altitude over the ellipsoid: \t\t" << c16platalt << "\nGround speed at nadir at input scene centre time: \t" << c16platgs << "\nPlatform heading at nadir" << "\n +corresponding to scene centre (degrees): \t\t" << c16plathead << "\nName of reference ellipsoid: \t\t\t\t" << c32refellips << "\nSemimajor axis of ref.ellipsoid (km): \t\t\t" << c16refmajor << "\nSemiminor axis of ref.ellipsoid (km): \t\t\t" << c16refminor << "\n\nMap projection data record: coordinates of four corner points" << "\n-------------------------------------------------------------" << "\n1st line 1st pixel geodetic latitude" << "\n +(positive for North latitude) (degrees): \t\t" << c16lat11 << "\n1st line 1st pixel geodetic longitude" << "\n +(negative for West longitude) (degrees): \t\t" << c16lon11 << "\n1st line last pixel geodetic latitude (degrees): \t" << c16lat1N << "\n1st line last pixel geodetic longitude (degrees): \t" << c16lon1N << "\nlast line last pixel geodetic latitude (degrees): \t" << c16latNN << "\nlast line last pixel geodetic longitude (degrees): \t" << c16lonNN << "\nlast line 1st pixel geodetic latitude (degrees): \t" << c16latN1 << "\nlast line 1st pixel geodetic longitude (degrees): \t" << c16lonN1 << "\n\nSLC platform position data record: positional data points" << "\n---------------------------------------------------------" << "\nNumber of data points: \t\t\t\t" << numdatapoints << "\nYear of data point <YYYY>: \t\t\t" << c4year << "\nMonth of data point <$$MM>: \t\t\t" << c4month << "\nDay of data point <$$DD>: \t\t\t" << c4day << "\nDay in the year <GMT> (jan 1st=day 1): \t\t" << c4dayofyear << "\nSeconds of day of data: \t\t\t" << c22seconds << "\nTime interval between data points: \t\t" << c22interval << "\nReference coordinate system: \t\t\t" << c64rcs << "\nGreenwich mean hour angle (degrees): \t\t" << c22gmha << "\nAlong track position error (meters): \t\t" << c16ltposerr << "\nAcross track position error (meters): \t\t" << c16ctposerr << "\nRadial position error (meters): \t\t" << c16rposerr << "\n\nSLC facility related data record [general type]: calibration information" << "\n------------------------------------------------------------------------" << "\nIncidence angle at first range pixel (at mid-azimuth): \t" << c16incangle1strange << "\nIncidence angle at centre range pixel (at mid-azimuth): " << c16incanglecenrange << "\nIncidence angle at last range pixel (at mid-azimuth): \t" << c16incanglelstrange << "\nAbsolute calibration constant K: \t" << calK << "\nLenrec6: \t" << lenrec6 << "\nReplica Pulse Power: \t" << repplspwr << "\nNumber of valid samples:   \t\t\t" << c4numvalid << "\nNumber of invalid samples: \t\t\t" << c4numinvalid << "\n";
	  int32 point;

	  // --- STATE VECTORS in INERTIAL SYSTEM ---
	  // only convert these if RSAT and ATLANTIS?
	//if (readfiles_arg.sensor_id==SLC_RSAT)
	if (readfiles_arg.argvalue.sensor_id==SLC_RSAT && readfiles_arg.argvalue.sar_processor==SARPR_ATL)
	  {
	  scratchlogfile << "\nRSAT: INERTIAL system orbital data as in leaderfile";
	  for (register int32 k =0;k<numdatapoints;k++) // number of data points
		{
		point = k+1;
		real8 secofday = Double.parseDouble(c22seconds) + real8(k)*Double.parseDouble(c22interval);
	  scratchlogfile << "\nSLC platform position data record: data point: " << point << "\n------------------------------------------------\n" << point << " data point - Seconds of day (s):    \t\t" << setprecision(13) << secofday << "\n" << point << " data point - Position vector X (m): \t\t" << setprecision(13) << STATE_INERTIAL(0,k) << "\n" << point << " data point - Position vector Y (m): \t\t" << setprecision(13) << STATE_INERTIAL(1,k) << "\n" << point << " data point - Position vector Z (m): \t\t" << setprecision(13) << STATE_INERTIAL(2,k) << "\n" << point << " data point - Velocity vector X (mm/s): \t" << setprecision(13) << STATE_INERTIAL(3,k) << "\n" << point << " data point - Velocity vector Y (mm/s): \t" << setprecision(13) << STATE_INERTIAL(4,k) << "\n" << point << " data point - Velocity vector Z (mm/s): \t" << setprecision(13) << STATE_INERTIAL(5,k) << "\n";
		}
	  scratchlogfile << "\nRSAT: and converted to earth fixed using GMHA:";
	  }

	  // --- STATE VECTORS in EARTH FIXED SYSTEM ---
	  for (register int32 k =0;k<numdatapoints;k++) // number of data points
		{
		point = k+1;
		real8 secofday = Double.parseDouble(c22seconds) + real8(k)*Double.parseDouble(c22interval);
	  scratchlogfile << "\nSLC platform position data record: data point: " << point << "\n------------------------------------------------\n" << point << " data point - Seconds of day (s):    \t\t" << setprecision(13) << secofday << "\n" << point << " data point - Position vector X (m): \t\t" << setprecision(13) << STATE(0,k) << "\n" << point << " data point - Position vector Y (m): \t\t" << setprecision(13) << STATE(1,k) << "\n" << point << " data point - Position vector Z (m): \t\t" << setprecision(13) << STATE(2,k) << "\n" << point << " data point - Velocity vector X (mm/s): \t" << setprecision(13) << STATE(3,k) << "\n" << point << " data point - Velocity vector Y (mm/s): \t" << setprecision(13) << STATE(4,k) << "\n" << point << " data point - Velocity vector Z (mm/s): \t" << setprecision(13) << STATE(5,k) << "\n";
		}

	  scratchlogfile << "\nEND LEADER FILE " << readfiles_arg.argvalue.leaderfile << "\n" << "\n";
	  scratchlogfile.close();



	  // ___ "repair" bandwidth if not given, e.g., in JERS? ___
	  if (Double.parseDouble(c16bandazitot) < 1.0)
		{
		WARNING.print("could not find azimuth band width, using prf...");
		c16bandazitot = c16prf;
		}
	  if (Double.parseDouble(c16bandrangetot) < 1.0)
		{
		WARNING.print("could not find range band width, using rsr...");
		c16bandrangetot = c16samplingrate;
		}

	  real8 wavelength_computed = Double.parseDouble(c16wavelength); // seems freq. not ok for RSAT
	  // --- seems not all information there for Atlantis CEOS ---
	  //if (readfiles_arg.sensor_id!=SLC_RSAT)
	  if (readfiles_arg.argvalue.sar_processor!=SARPR_ATL)
		{
		INFO << "Computing wavelength from radar freq. " << c8freq << "GHz"; // GHz
		INFO.print();
		wavelength_computed = (0.000000001 *SOL/Double.parseDouble(c8freq)); // seems more reliable, BK 03/04
		INFO << "Annotated wavelength: " << c16wavelength << " re-computed: " << wavelength_computed;
		INFO.print();
		}

	  // Email KKMohanty: time indicator is different then line counter, i.e, image
	  // can be flipped up/down left/right.  The utilities converting line number to
	  // az_time can thus be confused.  Repair this.
	  //  Time direction indicator along pixel direction:         INCREASE or DECREASE
	  //  Time direction indicator along line direction:          INCREASE or DECREASE
	  // BK 20-Apr-2004
	  INFO.print("Checking time indicator along line direction.");
	  // ___ check indicator in line direction ___
	  if (!strcmp(c8timeline,"INCREASE"))
		{
		INFO.print("Time indicator along line direction is INCREASE, this is OK.");
		}
	  else if (!strcmp(c8timeline,"DECREASE"))
		{
		WARNING.print("Time indicator along line direction is DECREASE: should adapt times.");
		WARNING.print("TODO: done for RSAT, not for ERS.");
		// ___ Adapt the time of first line/last line, such that orbit is correctly ___
		// ___ interpolated and line2ta works OK ___
		}
	  else
		{
		WARNING.print("unexpected: time indicator along line direction not INCREASE nor DECREASE.");
		}
	  // ___ check indicator in pixel direction ___
	  if (!strcmp(c8timepix,"INCREASE"))
		{
		INFO.print("Time indicator along pixel direction is INCREASE, this is OK.");
		}
	  else if (!strcmp(c8timepix,"DECREASE"))
		{
		WARNING.print("Time indicator along pixel direction is DECREASE: should adapt times.");
		WARNING.print("TODO: done for RSAT, not for ERS.");
		}
	  else
		{
		WARNING.print("unexpected: time indicator along pixel direction not INCREASE nor DECREASE.");
		}



	  ofstream scratchresfile = new ofstream("scratchreslea", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile, "readleader: scratchreslea", __FILE__, __LINE__);
	  // Modified by LG for reading ALOS Fine
	  if(!strcmp(c16centerlat,"                "))
			  c16centerlat = "0";
	  if(!strcmp(c16centerlon,"                "))
			  c16centerlon = "0";

	  scratchresfile << "Leader file:                                 \t" << readfiles_arg.argvalue.leaderfile << "\nSensor platform mission identifer:         \t" << c32sensorid << "\nScene_centre_latitude:                     \t" << c16centerlat << "\nScene_centre_longitude:                    \t" << c16centerlon;
	  // start_added_by_don
	  if (readfiles_arg.argvalue.sensor_id == SLC_ALOS)
		{
		scratchresfile << "\nRadar_wavelength (m):                      \t" << c16wavelength;
		}
	  else
	  // end_added_by_don
		{
		scratchresfile << "\nRadar_wavelength (m):                      \t" << setprecision(16) << wavelength_computed;
		// <<  c16wavelength;
		}

		// ______ Azimuth ______
	  //if (readfiles_arg.sar_processor==SARPR_ATL)
	  // start_added_by_don
	  if ((readfiles_arg.argvalue.sensor_id == SLC_ALOS) || (readfiles_arg.argvalue.sar_processor==SARPR_ATL))
	  // end_added_by_don
		{
		WARNING.print("First_pixel_azimuth_time not in leader, but done in readdat.");
		scratchresfile << "\nRSAT_First_pixel_azimuth_time (UTC):            \t" << "skipped_see_datfile_below";
		  //<<  c24zd1stazitime;
		}
	  else
		{
		scratchresfile << "\nFirst_pixel_azimuth_time (UTC):            \t" << c24zd1stazitime;
		}
		// email to ESA, they said this value is most accurate and should be used.
		//<< "\nPulse_Repetition_Frequency (actual, Hz):   \t"
		//<<  c16prf
		// (but that seems strange to me...)
	  scratchresfile << "\nPulse_Repetition_Frequency (computed, Hz): \t" << setprecision(16) << prfcomputed << "\nTotal_azimuth_band_width (Hz):             \t" << c16bandazitot << "\nWeighting_azimuth:                         \t" << c32weightazi << "\nXtrack_f_DC_constant (Hz, early edge):     \t" << c16xtdoppcconst;

	  // start_added_by_don
	  if ((readfiles_arg.argvalue.sar_processor==SARPR_JAX) || (readfiles_arg.argvalue.sar_processor==SARPR_ATL))
		{
		real8 xtdoppclinear = Double.parseDouble(c16xtdoppclinear) * rsrcomputed * 1e6;
		real8 xtdoppcquadratic = Double.parseDouble(c16xtdoppcquadratic) * sqr(rsrcomputed *1e6);
		scratchresfile << "\nXtrack_f_DC_linear (Hz/s, early edge):     \t" << setprecision(16) << xtdoppclinear << "\nXtrack_f_DC_quadratic (Hz/s/s, early edge): \t" << xtdoppcquadratic;
		}
	  else
	  // start_added_by_don
		{
		scratchresfile << "\nXtrack_f_DC_linear (Hz/s, early edge):     \t" << c16xtdoppclinear << "\nXtrack_f_DC_quadratic (Hz/s/s, early edge): \t" << c16xtdoppcquadratic;
	  }

	  // ______ Range ______
	  //if (readfiles_arg.sar_processor==SARPR_ATL)
	  // Modified by LG for reading ALOS Fine
	  if ((readfiles_arg.argvalue.sensor_id == SLC_ALOS) || (readfiles_arg.argvalue.sar_processor==SARPR_ATL))
		{
		WARNING.print("Range_time_to_first_pixel not in leader, but done in readdat.");
		scratchresfile << "\nRSAT_Range_time_to_first_pixel (2way) (ms):     \t" << "skipped_see_datfile_below";
		  //<<  c16zd1strange;
		}
	  else
		{
		real8 range_t1 = Double.parseDouble(c16zd1strange);
		// Bert Kampes, 20-Sep-2005: analysis of Corner Reflectors showed about 4
		// samples offset: systematic error, correct for ERS1 and ERS2.
		// this means that the radarcoding is now better.
		// If you don;t like this, you can use card [MS]_T_RG_ERROR to correct
		// but it is not that important.
		//
		// No, we did not like this and turned it off. FvL, 13-6-2006
		//
		//real8 range_t1 = atof(c16zd1strange);
		//if (readfiles_arg.sar_processor==SARPR_VMP)
		//  {
		//  INFO.print("ERS range time correction: -0.00019235 msec (two-way); ");
		//  range_t1 += -0.00019235;// ERS CR computed by BK ~75m
		//  }
		scratchresfile << "\nRange_time_to_first_pixel (2way) (ms):     \t" << setprecision(16) << range_t1;
		}
	  scratchresfile << "\nRange_sampling_rate (computed, MHz):       \t" << setprecision(16) << rsrcomputed << "\nTotal_range_band_width (MHz):               \t" << rangebandwidth << "\nWeighting_range:                            \t" << c32weightrange << "\n";

	  scratchresfile.setf(ios.fixed | ios.floatfield);
	  scratchresfile << "\n\n*******************************************************************" << "\n*_Start_leader_datapoints" << "\n*******************************************************************" << "\nt(s)\t\tX(m)\t\tY(m)\t\tZ(m)" << "\nNUMBER_OF_DATAPOINTS: \t\t\t" << numdatapoints << "\n";

	  for (register int32 l =0;l<numdatapoints;l++) // number of data points
		{
		scratchresfile << setprecision(11) << Double.parseDouble(c22seconds)+real4(l)*Double.parseDouble(c22interval);
		for (register int32 m =0;m<3;m++) // no velocities
		  {
		  scratchresfile << " \t" << setprecision(13) << STATE(m,l);
		  }
		scratchresfile << "\n";
		}
	  scratchresfile << "\n*******************************************************************" << "\n* End_leader_datapoints:_NORMAL" << "\n*******************************************************************\n";
	  scratchresfile.close();


	// ______ Check with volumefile (but first write what is read) ______
	  if (Integer.parseInt(c16numlin) != checklines)
		{
		WARNING << "Number of lines of leader file seems " << "not to be consistent with volume file: " << c16numlin << " != " << checklines;
		WARNING.print();
		WARNING.print("This may be caused by wrong? format SLC of this PAF, check parameterpool.");
		}
	  if (100*((Double.parseDouble(c16prf)-prfcomputed)/prfcomputed) > 2.)
		WARNING.print("deviation PRF computed from PRF read in leaderfile > 2%");
	  if (100*((Double.parseDouble(c16samplingrate)-rsrcomputed)/rsrcomputed) > 2.)
		WARNING.print("deviation range sampling rate computed from RSR read in leaderfile > 2%");


	// ______ Write some info ______
	  INFO << "Number of lines, pixels:               " << c16numlin << " " << c16numpix;
	  INFO.print();
	  INFO << "Pulse repetition frequency (Hz):       " << c16prf << ends;
	  INFO.print();
	  INFO << "Pulse repetition frequency (computed): " << setprecision(16) << prfcomputed;
	  INFO.print();
	  INFO << "Range sampling rate (Mhz):             " << c16samplingrate << ends;
	  INFO.print();
	  INFO << "Range sampling rate (computed Mhz):    " << setprecision(16) << rsrcomputed;
	  INFO.print();
	  INFO << "UTC of first azimuth line:             " << c24zd1stazitime;
	  INFO.print();
	  INFO << "UTC of last azimuth line:              " << c24zdlstazitime;
	  INFO.print();
	  INFO << "Range time to first pixel (2way ms):   " << c16zd1strange;
	  INFO.print();
	  INFO << "Range time to last pixel (2way ms):    " << c16zdlstrange;
	  INFO.print();
	  INFO << "First corner of image (lat,lon):       " << c16lat11 << " " << c16lon11;
	  INFO.print();
	  INFO << "Second corner of image (lat,lon):      " << c16lat1N << " " << c16lon1N;
	  INFO.print();
	  INFO << "Third corner of image (lat,lon):       " << c16latNN << " " << c16lonNN;
	  INFO.print();
	  INFO << "Fourth corner of image (lat,lon):      " << c16latN1 << " " << c16lonN1;
	  INFO.print();

	  INFO << "Weighting function designator azimuth: " << c32weightazi;
	  INFO.print();
	  Character.toUpperCase(c32weightazi);
	  if (strncmp(c32weightazi,"HAMMING",7))
		{
		WARNING << INFO.get_str() << " not HAMMING.";
		WARNING.print();
		}
	  INFO << "Weighting function designator range:   " << c32weightrange;
	  INFO.print();
	  Character.toUpperCase(c32weightrange);
	  if (strncmp(c32weightrange,"HAMMING",7))
		{
		WARNING << INFO.get_str() << " not HAMMING.";
		WARNING.print();
		}

	  PROGRESS.print("readleader finished.");
	  } // END READLEADER



	//***************************************************************
	// *    readnull                                                  *
	// *                                                              *
	// * reads nullfile                                               *
	// *  and writes to scratchlogfile, scratchresfile                *
	// * see: annex C ERS SAR.SLC CCTand EXABYTE                      *
	// *      doc:er-is-epo-gs-5902.3                                 *
	// *                                                              *
	// *    Bert Kampes, 11-Dec-1998                                  *
	// ***************************************************************
	public static void readnull(input_readfiles readfiles_arg)
	  {
	// ======Open files====== 
	  ifstream nullfile;
	  RefObject<ifstream> TempRefObject = new RefObject<ifstream>(nullfile);
	  openfstream(TempRefObject, readfiles_arg.nullfile);
	  nullfile = TempRefObject.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(nullfile, readfiles_arg.nullfile, __FILE__, __LINE__);

	  ofstream scratchlogfile = new ofstream("scratchlognull", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchlogfile, "readnull: scratchlognull", __FILE__, __LINE__);

	  ofstream scratchresfile = new ofstream("scratchresnull", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile, "readnull: scratchresnull", __FILE__, __LINE__);

	// ======Read nullfile======
	  DEBUG.print("nullfile should be read?");

	// ______Tidy up______
	  nullfile.close();
	  scratchlogfile.close();
	  scratchresfile.close();
	  PROGRESS.print("readnull finished.");
	  } // END READNUL



	//***************************************************************
	// *    readdat                                                   *
	// *                                                              *
	// * Extract info from data file,                                 *
	// *  see also appendix C of CD-R distribution (ESA).             *
	// * checks with volumefile #lines                                *
	// * checks some things like numchannels                          *
	// *                                                              *
	// * info is written to scratchresfile?                           *
	// *                                                              *
	// *    Bert Kampes, 21-Dec-1998                                  *
	// ***************************************************************
	public static void readdat(RefObject<input_readfiles> readfiles_arg, int32 checklines)
	  {
	  TRACE_FUNCTION("readdat (BK 21-Dec-1998)")
	  final int16 sizeb4 = 4;
	  final int16 sizeb1 = 1;
	  final int16 sizei4 = 4;
	  final int16 sizei6 = 6;
	  final int16 sizei8 = 8;
	  uint numchannels; // offset (=0?)
	  uint leftborder;
	  uint rightborder;
	  uint bottomborder;
	  uint topborder;
	  uint numdatarec; // SAR DATA records length
	  uint numlines;
	  uint numpixels;
	  uint numbytesdata;
	  uint lendatarec2;
	  String c4 = new String(new char[5]); // correctly 9 for \0
	  String c6 = new String(new char[7]);
	  String c8 = new String(new char[9]);


	  // --- Check for RSAT #%// Bert Kampes, 02-Aug-2004 ---
	  uint rec_seq; // type B4
	  byte rec_sub1; // type B1
	  byte rec_type;
	  byte rec_sub2;
	  byte rec_sub3;



	// ______Open files______ 
	  ifstream datfile;
	  RefObject<ifstream> TempRefObject = new RefObject<ifstream>(datfile);
	  openfstream(TempRefObject, readfiles_arg.argvalue.datfile);
	  datfile = TempRefObject.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(datfile, readfiles_arg.argvalue.datfile, __FILE__, __LINE__);


	  // --- RECORD 1 ---
	  DEBUG.print("record 1 of data file (ERS and RSAT).");
	  datfile.seekg(0,ios.beg);
	  datfile.read((char) rec_seq, sizeb4); // record number
	  rec_seq = ntohl(rec_seq); // Bert Kampes, 07-Apr-2005
	  datfile.read((char)&rec_sub1, sizeb1); // first record sub type code
	  datfile.read((char)&rec_type, sizeb1); // record type code
	  datfile.read((char)&rec_sub2, sizeb1); // second record sub type code
	  datfile.read((char)&rec_sub3, sizeb1); // third record sub type code
	  DEBUG.print("ERS/RSAT: Expecting record 1 with code {63,192,18,18}");
	  DEBUG << "rec_seq: " << rec_seq << "; rec_sub1: " << (int)rec_sub1 << "; rec_type: " << (int)rec_type << "; rec_sub2: " << (int)rec_sub2 << "; rec_sub3: " << (int)rec_sub3;
	  DEBUG.print();
	  if ((int)rec_sub1 ==63 && (int)rec_type ==192 && (int)rec_sub2 ==18 && (int)rec_sub3 ==18)
		DEBUG.print("This is the expected record with code {63,192,18,18}");
	  else
		WARNING.print("This is NOT the expected record with code {63,192,18,18}");
	  uint lenrec1; // length of record1
	  datfile.seekg(8,ios.beg);
	  datfile.read((char) lenrec1, sizeb4); // length of record
	  lenrec1 = ntohl(lenrec1); // bk 6 jul 2000, byteorder x86 machines.
	  DEBUG.print("ERS:  Expecting record 1 with length 10012");
	  DEBUG.print("RSAT: Expecting record 1 with length 16252");
	  DEBUG << "readdat::record 1: start at: " << 0 << "; length: " << lenrec1;
	  DEBUG.print();


	  // ====== Get general info ======
	  datfile.seekg(180,ios.beg);
	  datfile.read((char)&c6, sizei6); // number of SAR DATA records (lines)
	  c6.charAt(6)='\0';
		numdatarec = Integer.parseInt(c6);
	  datfile.read((char)&c6, sizei6); // SAR DATA record length
	  c6.charAt(6)='\0';
		lendatarec2 = Integer.parseInt(c6);
	  datfile.seekg(232,ios.beg); // SAR Related data
	  datfile.read((char)&c4, 4);
	  c4.charAt(4)='\0';
	  numchannels = Integer.parseInt(c4);
	  datfile.read((char)&c8, sizei8);
	  c8.charAt(8)='\0';
	  numlines = Integer.parseInt(c8);

	  datfile.read((char)&c4, sizei4);
		c4.charAt(4)='\0';
		leftborder = Integer.parseInt(c4);
	  datfile.read((char)&c8, sizei8); // number of pixels
		c8.charAt(8)='\0';
		numpixels = Integer.parseInt(c8);
	  datfile.read((char)&c4, sizei4);
		c4.charAt(4)='\0';
		rightborder = Integer.parseInt(c4);
	  datfile.read((char)&c4, sizei4);
		c4.charAt(4)='\0';
		topborder = Integer.parseInt(c4);
	  datfile.read((char)&c4, sizei4);
		c4.charAt(4)='\0';
		bottomborder = Integer.parseInt(c4);

	  datfile.seekg(280,ios.beg); // Record data
	  datfile.read((char)&c8, sizei8);
	  c8.charAt(8)='\0';
	  numbytesdata = Integer.parseInt(c8);


	// ======Write to scratchfiles======
	  ofstream scratchlogdat = new ofstream("scratchlogdat", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchlogdat, "readdat: scratchlogdat", __FILE__, __LINE__);

	  scratchlogdat << "\n\n*******************************************************************" << "\n* EXTRACTED DATA FROM DATA FILE: " << readfiles_arg.argvalue.datfile << " *" << "\n*******************************************************************" << "\nNumber of SAR channels in file:         \t" << numchannels << "\nNumber of lines:                        \t" << numlines << "\nNumber of pixels:                       \t" << numpixels << "\nNumber of left border pixels per line:  \t" << leftborder << "\nNumber of right border pixels per line: \t" << rightborder << "\nNumber of topborder lines:              \t" << topborder << "\nNumber of bottom border lines:          \t" << bottomborder << "\nNumber of bytes per data group:         \t" << numbytesdata << "\n";
	  scratchlogdat.close();


	  // --- write RESULTFILE, add info for RSAT which was not in LEADER ---
	  ofstream scratchresdat = new ofstream("scratchresdat", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresdat, "readdat: scratchresdat", __FILE__, __LINE__);
	  scratchresdat << "Datafile: \t\t\t\t\t" << readfiles_arg.argvalue.datfile << "\nNumber_of_lines_original: \t\t\t" << numlines << "\nNumber_of_pixels_original: \t\t\t" << numpixels;

	  // --- if rsat, read datfile further, write az/rg time ----------------
	//if (readfiles_arg.sensor_id==SLC_RSAT)
	//if (readfiles_arg.sar_processor==SARPR_ATL)
	  if ((readfiles_arg.argvalue.sar_processor==SARPR_ATL) || (readfiles_arg.argvalue.sar_processor==SARPR_JAX)) //LG
	  {
	  // --- RECORD 2 ---
	  uint startrec2 = lenrec1;
	  DEBUG.print("record 2 of data file (RSAT).");
	  datfile.seekg(startrec2,ios.beg);
	  datfile.read((char) rec_seq, sizeb4); // record number
	  rec_seq = ntohl(rec_seq); // Bert Kampes, 07-Apr-2005
	  datfile.read((char)&rec_sub1, sizeb1); // first record sub type code
	  datfile.read((char)&rec_type, sizeb1); // record type code
	  datfile.read((char)&rec_sub2, sizeb1); // second record sub type code
	  datfile.read((char)&rec_sub3, sizeb1); // third record sub type code
	  DEBUG.print("RSAT: Expecting record 2 with code {50,11,18,20}");
	  DEBUG << "rec_seq: " << rec_seq << "; rec_sub1: " << (int)rec_sub1 << "; rec_type: " << (int)rec_type << "; rec_sub2: " << (int)rec_sub2 << "; rec_sub3: " << (int)rec_sub3;
	  DEBUG.print();
	  if ((int)rec_sub1 ==50 && (int)rec_type ==11 && (int)rec_sub2 ==18 && (int)rec_sub3 ==20)
		DEBUG.print("This is the expected record with code {50,11,18,20}");
	  else
		WARNING.print("This is NOT the expected record with code {50,11,18,20}");
	  uint lenrec2; // length of record2
	  datfile.seekg(startrec2+8,ios.beg);
	  datfile.read((char) lenrec2, sizeb4); // length of record
	  lenrec2 = ntohl(lenrec2); // bk 6 jul 2000, byteorder x86 machines.
	  DEBUG.print("RSAT: Expecting record 2 with length variable");
	  DEBUG << "readdat::record 2: start at: " << startrec2 << "; length: " << lenrec2;
	  DEBUG.print();
	  uint startrec3 = lenrec1+lenrec2;
	  uint startrecN = lenrec1+(numlines-1)*lenrec2; // start of last record
	  // --- azimuth time to first line (depends on decrease/increase): ---
	  uint zdmsecofday1 = 99999; // B4
	  uint zdmsecofday2 = 99999; // B4
	  uint zdmsecofdayN = 99999; // B4
	  datfile.seekg(startrec2+44,ios.beg);
	  datfile.read((char) zdmsecofday1, sizeb4); // range to first pix
	  zdmsecofday1 = ntohl(zdmsecofday1); // Bert Kampes, 07-Apr-2005
	  datfile.seekg(startrec3+44,ios.beg);
	  datfile.read((char) zdmsecofday2, sizeb4); // range to first pix
	  zdmsecofday2 = ntohl(zdmsecofday2); // Bert Kampes, 07-Apr-2005
	  datfile.seekg(startrecN+44,ios.beg);
	  datfile.read((char) zdmsecofdayN, sizeb4); // range to first pix
	  zdmsecofdayN = ntohl(zdmsecofdayN); // Bert Kampes, 07-Apr-2005
	  INFO << "zdmsecofday1: " << zdmsecofday1;
	  INFO.print();
	  DEBUG << "zdmsecofday2: " << zdmsecofday2;
	  DEBUG.print();
	  INFO << "zdmsecofdayN: " << zdmsecofdayN;
	  INFO.print();
	  // --- Check PRF/RSR ---
	  real8 prf_check = new real8(numlines-1);
	  INFO << "PRF check (computed [Hz]): " << prf_check;
	  INFO.print();


	  // format should be: "22-AUG-1997 18:22:10.246"
	  if (zdmsecofday1 < zdmsecofdayN) // increase, use ZD time of first line
		{
		INFO.print("INCREASING azimuth time detected RSAT");
		}
	  else // decrease: use first-numlines*prf? or read last record.
		{
		INFO.print("DECREASING azimuth time detected (flip image up-down) RSAT");
		zdmsecofday1 = zdmsecofdayN; // use last record (line) as first line (flip)
		startrec2 = startrecN; // use last record (line) as first line (flip)
		}
	  // --- Format the string UTC ---
	  uint acq_year = 1999;
	  uint acq_day = 19;
	  datfile.seekg(startrec2+36,ios.beg);
	  datfile.read((char) acq_year, sizeb4); // range to first pix
	  acq_year = ntohl(acq_year); // Bert Kampes, 07-Apr-2005
	  datfile.read((char) acq_day, sizeb4); // range to first pix
	  acq_day = ntohl(acq_day); // Bert Kampes, 07-Apr-2005
	  INFO << "acq_year: " << acq_year;
	  INFO.print();
	  INFO << "acq_day: " << acq_day;
	  INFO.print();
	  // --- fill tm struct using daynumber with strptime, re-format it using strftime ---
	  String datestring = new String(new char[13]); // e.g., "25-Jan-1999"
	  tm tm_tmp;
	  String buf = new String(new char[9]); // e.g., "1999 191";
	  String.format(buf,"%4d %03d", acq_year,acq_day);
	  RefObject<tm> TempRefObject2 = new RefObject<tm>(tm_tmp);
	  strptime(buf, "%Y %j", TempRefObject2);
	  tm_tmp = TempRefObject2.argvalue;
	  int q = strftime(datestring, 12, "%d-%b-%Y", tm_tmp); // q: number of bytes written
	  INFO << "numbytes in datestring: (12?) " << q;
	  INFO.print();

	  // --- Fill important part: (sec of day) ---
	  String c24zd1stazitime = new String(new char[25]);
	  int32 hour = real8(zdmsecofday1)/1000.0/60.0/60.0; // floor
	  int32 min = (real8(zdmsecofday1)-hour *1000.0 *60.0 *60.0)/1000.0/60.0; // floor
	  int32 sec = (real8(zdmsecofday1)-hour *1000.0 *60.0 *60.0-min *1000.0 *60.0)/1000.0; // floor
	  int32 msec = (real8(zdmsecofday1)-hour *1000.0 *60.0 *60.0-min *1000.0 *60.0-1000.0 *sec); // floor
	  //sprintf(c24zd1stazitime,"01-JAN-1990 %02d:%02d:%02d.%03d", hour,min,sec,msec);
	  String.format(c24zd1stazitime,"%11s %02d:%02d:%02d.%03d", datestring, hour,min,sec,msec);
	  INFO << "c24zd1stazitime: " << c24zd1stazitime;
	  INFO.print();

	  // --- range time to first pixel is distance -------------------------
	  uint range1st = 99999;
	  uint rangelst = 99999;
	  // start_added_by_don
	  real8 zd1strange;
	  if (readfiles_arg.argvalue.sar_processor==SARPR_JAX)
		{
		datfile.seekg(startrec2+116,ios.beg);
		datfile.read((char) range1st, sizeb4); // range to first pix
		range1st = ntohl(range1st); // Bert Kampes, 07-Apr-2005
		zd1strange = 2000.0 *range1st/SOL;
		INFO << "range1st: " << range1st;
		INFO.print();
		 }
	  if (readfiles_arg.argvalue.sar_processor==SARPR_ATL)
	  // end_added_by_don
		{
	  datfile.seekg(startrec2+64,ios.beg);
	  datfile.read((char) range1st, sizeb4); // range to first pix
	  range1st = ntohl(range1st); // Bert Kampes, 07-Apr-2005
	  datfile.seekg(startrec2+72,ios.beg);
	  datfile.read((char) rangelst, sizeb4); // range to last pix
	  rangelst = ntohl(rangelst); // Bert Kampes, 07-Apr-2005
	  zd1strange = (range1st<rangelst) ? 2000.0 *range1st/SOL : 2000.0 *rangelst/SOL;
	  INFO << "range1st: " << range1st;
	  INFO.print();
	  INFO << "rangelst: " << rangelst;
	  INFO.print();
	  // --- Check PRF/RSR ---
	  real8 rsr_check = new real8(numpixels-1);
	  INFO << "RSR check (computed [MHz]): " << rsr_check;
	  INFO.print();
	  }
	  //char dummydate[] = "01-JAN-1990 ";// not used except maybe getorb later
	  //WARNING.print("RSAT: using a dummy date for orbit, only secofday important.");
		//    << "22-AUG-1997 18:22:10.246"
	  scratchresdat << "\nFirst_pixel_azimuth_time (UTC):            \t" << c24zd1stazitime << "\nRange_time_to_first_pixel (2way) (ms):     \t" << setprecision(16) << zd1strange;
	  }

	  scratchresdat << "\n*******************************************************************";
	  if (readfiles_arg.argvalue.fileid == MASTERID)
		scratchresdat << "\n* End_" << processcontrol[pr_m_readfiles] << "_NORMAL";
	  if (readfiles_arg.argvalue.fileid == SLAVEID)
		scratchresdat << "\n* End_" << processcontrol[pr_s_readfiles] << "_NORMAL";
	  scratchresdat << "\n*******************************************************************" << "\n";
	  scratchresdat.close();
	  datfile.close();

	// ______Tidy up______
	  if (numchannels != 1) // ??
		{
		WARNING << "code 904: Number of channels in file: " << readfiles_arg.argvalue.datfile << " = " << numchannels << " != 1 ";
		WARNING.print();
		WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
		}

	// ______ Check with volumefile ______
	  if (numlines != checklines)
		{
		WARNING << "code 902: data file: " << readfiles_arg.argvalue.datfile << " numlin=" << numlines << " vs.  volume file: " << readfiles_arg.argvalue.volfile << " numlin=" << checklines;
		WARNING.print();
		WARNING.print(" +this means data and volume file seem not to correspond.");
		}

	// ______ Check with previous section ______
	  if (numlines != numdatarec)
		{
		WARNING << "code 904: Number of lines seems not to be consistent in file: " << readfiles_arg.argvalue.datfile << " : " << numlines << " != " << numdatarec;
		WARNING.print();
		WARNING.print(" +this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
		}
	  if (bottomborder != topborder != leftborder != rightborder != 0)
		{
		WARNING << "code 904: Not implemented: offset border: left,right,bottom,top: " << leftborder << "," << rightborder << "," << bottomborder << "," << topborder << " in file: " << readfiles_arg.argvalue.datfile;
		WARNING.print();
		WARNING.print(" +this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
		}

	  // start_added_by_don
	  if (readfiles_arg.argvalue.sar_processor==SARPR_JAX)
		 {
		 if ((numbytesdata / 8) != numpixels)
		   {
		   WARNING << "code 904AAA: Number of pixels seems to be inconsistent in file: " << readfiles_arg.argvalue.datfile << ": " << numpixels << " != " << (numbytesdata / 8);
		   WARNING.print();
		   WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
		   }
		 }
	  else
		 {
	  if ((numbytesdata / 4) != numpixels)
		{
		WARNING << "code 904: Number of pixels seems to be inconsistent in file: " << readfiles_arg.argvalue.datfile << ": " << numpixels << " != " << (numbytesdata / 4);
		WARNING.print();
		WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
		}
	  }
	  PROGRESS.print("readdat (header info) finished.");
	  } // END readdat



	//***************************************************************
	// *    writeslc                                                  *
	// *                                                              *
	// * Inputfile in ceos slc format is converted to                 *
	// *  raw format outputfile.                                      *
	// * Some info is extracted from header,                          *
	// *  see also appendix C of CD-R distribution (ESA).             *
	// * checks with volumefile #lines                                *
	// *                                                              *
	// * A buffer (BUFFERMEMSIZE) is used to increase speed.          *
	// * Input is the inputfile and outputfile,                       *
	// *  an identifier for this run,                                 *
	// *  the area which should be converted (default is total image) *
	// * info is written to scratchresfile                            *
	// http://earth.esa.int/rootcollection/sysutil/01008.html
	// *                                                              *
	// *    Bert Kampes, 11-Dec-1998                                  *
	// ***************************************************************
	public static void writeslc(input_gen generalinput, input_crop crop_arg, int checklines)
	  {
	  final int16 sizeb4 = 4;
	  final int16 sizei4 = 4;
	  final int16 sizei6 = 6;
	  final int16 sizei8 = 8;
	  uint lenrec1; // length of general record1
	  uint lenrec2; // (nominal) length of data records
	  String c4 = new String(new char[5]); // correctly 9 for \0
	  String c6 = new String(new char[7]);
	  String c8 = new String(new char[9]);

	  // ______ Write some info ______ 
	  TRACE_FUNCTION("writeslc (BK 11-Dec-1998)")
	  PROGRESS.print("Start cropping slc data.");
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __X86PROCESSOR__
	  INFO.print("Swapping Big Endian (CEOS input) to Little Endian (your platform).");
	//  #else
	  INFO.print("NO byte swapping performed, you must be on Big Endian platform.");

	//  #endif
	  // ______ Open files ______ 
	  ifstream datfile;
	  RefObject<ifstream> TempRefObject = new RefObject<ifstream>(datfile);
	  openfstream(TempRefObject, crop_arg.filein1);
	  datfile = TempRefObject.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(datfile, crop_arg.filein1, __FILE__, __LINE__);

	  // ====== Get data such as recordlength ======
	  datfile.seekg(8,ios.beg);
	  datfile.read((char) lenrec1, sizeb4); // length of record1
	  lenrec1 = ntohl(lenrec1); // bk 6 jul 2000, byteorder x86 machines.
	  DEBUG.print("record 1 of data file.");
	  if (lenrec1 != 19612) // quarter scene 10012? JERS 22196?
		{
		WARNING << "writeslc: length of record 1 = \"" << lenrec1 << "\"; expected \"19612\" for ERS SLC (CEOS, full scene).";
		WARNING.print();
		DEBUG.print("10012 seems ERS quarter scene?");
		DEBUG.print("22196 seems JERS scene?");
		}

	  datfile.seekg(180,ios.beg);
	  datfile.read((char)&c6, sizei6); // number of SAR DATA records (lines)
		c6.charAt(6)='\0';
		final uint numdatarec = Integer.parseInt(c6);
	  DEBUG << "numdatarec: " << numdatarec;
	  DEBUG.print();
	  //datfile.seekg(186,ios::beg);
	  datfile.read((char)&c6, sizei6); // SAR DATA record length
		c6.charAt(6)='\0';
		final uint lendatarec2 = Integer.parseInt(c6);
	  DEBUG << "lendatarec2: " << lendatarec2;
	  DEBUG.print();
	  datfile.seekg(232,ios.beg); // SAR Related data
	  datfile.read((char)&c4, 4);
		c4.charAt(4)='\0';
		final uint numchannels = Integer.parseInt(c4);
	  DEBUG << "numchannels: " << numchannels;
	  DEBUG.print();

	  datfile.read((char)&c8, sizei8);
		c8.charAt(8)='\0';
		uint numlines = Integer.parseInt(c8);
	  DEBUG << "numlines: " << numlines;
	  DEBUG.print();

	  datfile.read((char)&c4, sizei4);
		c4.charAt(4)='\0';
		final uint leftborder = Integer.parseInt(c4);
	  DEBUG << "leftborder: " << leftborder;
	  DEBUG.print();
	  datfile.read((char)&c8, sizei8); // number of pixels
		c8.charAt(8)='\0';
		uint numpixels = Integer.parseInt(c8);
	  DEBUG << "numpixels: " << numpixels;
	  DEBUG.print();
	  datfile.read((char)&c4, sizei4);
		c4.charAt(4)='\0';
		final uint rightborder = Integer.parseInt(c4);
	  DEBUG << "rightborder: " << rightborder;
	  DEBUG.print();
	  datfile.read((char)&c4, sizei4);
		c4.charAt(4)='\0';
		final uint topborder = Integer.parseInt(c4);
	  DEBUG << "topborder: " << topborder;
	  DEBUG.print();
	  datfile.read((char)&c4, sizei4);
		c4.charAt(4)='\0';
		final uint bottomborder = Integer.parseInt(c4);
	  DEBUG << "bottomborder: " << bottomborder;
	  DEBUG.print();

	  datfile.seekg(280,ios.beg); // Record data
	  datfile.read((char)&c8, sizei8);
		c8.charAt(8)='\0';
		final uint numbytesdata = Integer.parseInt(c8);
	  DEBUG << "numbytesdata: " << numbytesdata;
	  DEBUG.print();


	// ====== Check with volumefile / internal ======
	  if (numlines != checklines)
		{
		WARNING << "code 902: data file: " << crop_arg.filein1 << " numlin=" << numlines << " vs. volume file: " << crop_arg.filein1 << " numlin=" << checklines;
		WARNING.print();
		WARNING.print(" +this means data and volume file seem not to correspond.");
		}

	// ______ Check with previous section ______
	  if (numlines != numdatarec)
		{
		WARNING << "code 904: Number of lines seems not to be consistent in file: " << crop_arg.filein1 << " : " << numlines << " != " << numdatarec;
		WARNING.print();
		WARNING.print(" +this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
		}
	  if ((numbytesdata / 4) != numpixels)
		{
		WARNING << "code 904: Number of pixels seems to be inconsistent in file: " << crop_arg.filein1 << ": " << numpixels << " != " << (numbytesdata / 4);
		WARNING.print();
		WARNING.print(" +this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
		}


	// ====== Start copy input to output (raw) format with buffer======
	// ______ Check and process optional offset parameters______
	// ______ Lcnlow is corner line, lcnhi is other corner, pcnlow, pixel coord. low etc.
	  uint linestart = 1; // counters for loops
	  uint lineend = numlines;
	  uint pixelstart = 1;
	  uint pixelend = numpixels; // only for resultfile

	  if (crop_arg.dbow.linehi!=0 && crop_arg.dbow.linelo!=0 && crop_arg.dbow.pixhi!=0 && crop_arg.dbow.pixlo!=0)
		{
		window tempdbow = new window(crop_arg.dbow.linelo, crop_arg.dbow.linehi, crop_arg.dbow.pixlo, crop_arg.dbow.pixhi);
		if (crop_arg.dbow.linehi>numlines)
		  {
		  WARNING << "Specified input DBOW linehi > numlines: " << crop_arg.dbow.linehi << " > " << numlines << ". I set linehi = " << numlines;
		  WARNING.print();
		  tempdbow.linehi=numlines;
		  }
		if (crop_arg.dbow.pixhi>numpixels)
		  {
		  WARNING << "Specified input DBOW pixhi > numpixels: " << crop_arg.dbow.pixhi << " > " << numpixels << ". I set pixhi = " << numpixels;
		  WARNING.print();
		  tempdbow.pixhi=numpixels;
		  }
	// ______ Only hi values are possibly adapted, low is a constant ______
		numlines = tempdbow.linehi - crop_arg.dbow.linelo + 1;
		numpixels = tempdbow.pixhi - crop_arg.dbow.pixlo + 1;

		linestart = crop_arg.dbow.linelo;
		lineend = tempdbow.linehi;
		pixelstart = crop_arg.dbow.pixlo;
		pixelend = tempdbow.pixhi; // only for resultfile
		}

	// ______ Note complex<short> not in ANSI c ______
	  matrix LINE = new matrix(1, 2 *numpixels); // size of short

	// ====== Process requested lines ======
	  ofstream datoutfile;
	  RefObject<ofstream> TempRefObject2 = new RefObject<ofstream>(datoutfile);
	  openfstream(TempRefObject2, crop_arg.fileout1, generalinput.overwrit);
	  datoutfile = TempRefObject2.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(datoutfile, crop_arg.fileout1, __FILE__, __LINE__);

	  // ______ info on data, to avoid X86 problems ______
	  // ______ according to CEOS specs, byte 13-16 is first complex pixel, etc. ______
	  datfile.seekg(lenrec1+12,ios.beg);
	  matrix TMPSHORT = new matrix(1, 2);
	  datfile >> TMPSHORT; // read in first complex pixel for test
	  real8 tmpmag = Math.sqrt(real8(int16(ntohs(TMPSHORT(0,0)))*int16(ntohs(TMPSHORT(0,0)))) + real8(int16(ntohs(TMPSHORT(0,1)))*int16(ntohs(TMPSHORT(0,1)))));
	  DEBUG << "First complex element in datafile: (" << int16(ntohs(TMPSHORT(0,0))) << "," << int16(ntohs(TMPSHORT(0,1))) << "); mag = " << tmpmag;
	  DEBUG.print();
	  if (tmpmag > 10000.)
		{
		WARNING.print(DEBUG.get_str());
		WARNING.print("this is a byteorder problem on X86? (use ntohs)");
		}
	  DEBUG << "TEST: (realpart): " << TMPSHORT(0,0) << ", (imagpart): " << TMPSHORT(0,1);
	  DEBUG.print();
	  DEBUG << "TEST: htons(realpart): " << htons(TMPSHORT(0,0)) << ", htons(imagpart): " << htons(TMPSHORT(0,1));
	  DEBUG.print();
	  DEBUG << "TEST: ntohs(realpart): " << ntohs(TMPSHORT(0,0)) << ", ntohs(imagpart): " << ntohs(TMPSHORT(0,1));
	  DEBUG.print();
	  DEBUG << "TEST: short int(ntohs(realpart)): " << int16(ntohs(TMPSHORT(0,0))) << ", (imagpart): " << int16(ntohs(TMPSHORT(0,1)));
	  DEBUG.print();


	  // ====== perline is faster than perbuffer, less memory etc. BK1998 ======
	  datfile.seekg(lenrec1+(linestart-1)*lendatarec2+8,ios.beg);
	  datfile.read((char) lenrec2, sizeb4); // length of first record
	  lenrec2 = ntohl(lenrec2); // bk 6 jul 2000, byteorder x86 machines.

	  if (lenrec2 != lendatarec2)
		{
		ERROR << "code 904: Length of datarecords seems to be inconsistent in file: " << crop_arg.filein1 << ": " << lenrec2 << " != " << lendatarec2;
		WARNING.print(ERROR.get_str());
		ERROR.reset();
		}

	  final int32 TEN = 10;
	  final int32 TENPERCENT = (.5 *TEN+numlines)/TEN; // number of lines
	  int32 percentage = 0; // initialization
	  final int32 tmpstart = lenrec1+12-lendatarec2+(pixelstart-1)*4; // sizeof=4
	  for (register int32 linecnt =linestart; linecnt<=lineend; linecnt++)
		{
		if (!((linecnt-linestart)%TENPERCENT))
		  {
		  PROGRESS << "WRITESLC: " << setw(3) << percentage << "%";
		  PROGRESS.print();
		  percentage += TEN;
		  }
		datfile.seekg(tmpstart+linecnt *lendatarec2,ios.beg);
		datfile >> LINE;
		// ______ BK 13 July 2000: swapbytes for X86 (intel) linux cpus ______
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//	#if __X86PROCESSOR__
		for (int ii =0; ii<LINE.pixels(); ++ii)
		  LINE(0,ii) = int16(ntohs(LINE(0,ii))); // changed from htons 171100 BK datoutfile << LINE;
	//	#else
		datoutfile << LINE;
	//	#endif
		}
	  datfile.close(); // close files
	  datoutfile.close();



	// ====== Write results to scratchfile ======
	  ofstream scratchresfile = new ofstream("scratchres2raw", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile, "writeslc: scratchres2raw", __FILE__, __LINE__);
	  scratchresfile << "\n\n*******************************************************************\n";
		//<< "\n*_Start_crop:\t\t\t"
		//<<  crop_arg.idcrop
	  if (crop_arg.fileid == MASTERID)
		scratchresfile << "*_Start_" << processcontrol[pr_m_crop];
	  if (crop_arg.fileid == SLAVEID)
		scratchresfile << "*_Start_" << processcontrol[pr_s_crop];
	// ______ updateslcimage greps these ______
	  scratchresfile << "\t\t\t" << crop_arg.idcrop << "\n*******************************************************************" << "\nData_output_file: \t\t\t\t" << crop_arg.fileout1 << "\nData_output_format: \t\t\t\t" << "complex_short" << "\nFirst_line (w.r.t. original_image): \t\t" << linestart << "\nLast_line (w.r.t. original_image): \t\t" << lineend << "\nFirst_pixel (w.r.t. original_image): \t\t" << pixelstart << "\nLast_pixel (w.r.t. original_image): \t\t" << pixelend << "\n*******************************************************************";
	  if (crop_arg.fileid == MASTERID)
		scratchresfile << "\n* End_" << processcontrol[pr_m_crop] << "_NORMAL";
	  if (crop_arg.fileid == SLAVEID)
		scratchresfile << "\n* End_" << processcontrol[pr_s_crop] << "_NORMAL";
	  scratchresfile << "\n*******************************************************************" << "\n";
	  scratchresfile.close();

	// ______ Tidy up do checks here ______
	  if (numchannels != 1) // ??
		{
		WARNING << "code 904: Number of channels in file: " << crop_arg.filein1 << " = " << numchannels << " != 1 ";
		WARNING.print();
		WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
		}
	  if (bottomborder != topborder != leftborder != rightborder != 0)
		{
		WARNING << "code 904: Not implemented: offset border: left,right,bottom,top: " << leftborder << "," << rightborder << "," << bottomborder << "," << topborder << " in file: " << crop_arg.filein1;
		WARNING.print();
		WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
		}
	  PROGRESS.print("WRITESLC: 100%");
	  } // END writeslc



	//***************************************************************
	// *    envisatdump_data                                          *
	// *                                                              *
	// * Via a system call to the c-program envisat_dump_data, the    *
	// * SLC data is wrtten to file.  The resfile is here created.    *
	// * envisatdumpdata writes SLC data out in host order.           *
	// * it is important that crop_arg.dbow is correctly filled.      *
	// *    Bert Kampes, 16-JUN-2003                                  *
	// ***************************************************************
	public static void envisat_dump_data(input_crop crop_arg)
	  {
	  // ______ Write some info ______ 
	  TRACE_FUNCTION("envisat_dump_data (BK 16-Jun-2003)")
	  // ______ Build command ______
	  // ______ make sure l0 etc. are correctly defined ______
	  // ____ assume these are filled correctly ___
	  INFO.reset();
	  if (crop_arg.dbow.linehi!=0 && crop_arg.dbow.linelo!=0 && crop_arg.dbow.pixhi!=0 && crop_arg.dbow.pixlo!=0)
		INFO << "envisat_dump_data " << crop_arg.filein1 << " " << crop_arg.fileout1 << " " << crop_arg.dbow.linelo << " " << crop_arg.dbow.linehi << " " << crop_arg.dbow.pixlo << " " << crop_arg.dbow.pixhi << ends;
	  else
		INFO << "envisat_dump_data " << crop_arg.filein1 << " " << crop_arg.fileout1 << ends;
	  String cmd = new String(new char[512]); // command string
	  cmd = INFO.get_str();
	  INFO.print("With following command the envisat data was cropped.");
	  INFO.print(cmd);
	  PROGRESS.print("system call may take some time...");
	  system(cmd); // this does the work
	  INFO.reset();
	  INFO.print();

	  // ====== Write results to scratchfile ======
	  ofstream scratchresfile = new ofstream("scratchres2raw", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile, "writeslc: scratchres2raw", __FILE__, __LINE__);
	  scratchresfile << "\n\n*******************************************************************\n";
	  if (crop_arg.fileid == MASTERID)
		scratchresfile << "*_Start_" << processcontrol[pr_m_crop];
	  if (crop_arg.fileid == SLAVEID)
		scratchresfile << "*_Start_" << processcontrol[pr_s_crop];
		// ______ updateslcimage greps these ______
	  scratchresfile << "\t\t\t" << crop_arg.idcrop << "\n*******************************************************************" << "\nData_output_file: \t\t\t\t" << crop_arg.fileout1 << "\nData_output_format: \t\t\t\t" << "complex_short" << "\nFirst_line (w.r.t. original_image): \t\t" << crop_arg.dbow.linelo << "\nLast_line (w.r.t. original_image): \t\t" << crop_arg.dbow.linehi << "\nFirst_pixel (w.r.t. original_image): \t\t" << crop_arg.dbow.pixlo << "\nLast_pixel (w.r.t. original_image): \t\t" << crop_arg.dbow.pixhi << "\n*******************************************************************";
	  if (crop_arg.fileid == MASTERID)
		scratchresfile << "\n* End_" << processcontrol[pr_m_crop] << "_NORMAL";
	  if (crop_arg.fileid == SLAVEID)
		scratchresfile << "\n* End_" << processcontrol[pr_s_crop] << "_NORMAL";
	  scratchresfile << "\n*******************************************************************" << "\n";
	  scratchresfile.close();
	  } // END envisat_dump_data


	//***************************************************************
	// *    tsxdump_data                                              *
	// *                                                              *
	// * Via a system call to the python tsx_dump_data, the           *
	// * SLC data is written to file.  The resfile is here created.   *
	// * tsxdumpdata writes SLC data out in host order.               *
	// * it is important that crop_arg.dbow is correctly filled.      *
	// *                                                              *
	// * Dependencies: GDAL                                           *
	// ***************************************************************
	public static void tsx_dump_data(input_crop crop_arg)
	{
	  // ______ Write some info ______ 
	  TRACE_FUNCTION("tsx_dump_data (PM 06-Apr-2009)")
		// ______ Build command ______
		// ______ make sure l0 etc. are correctly defined ______
		// ____ assume these are filled correctly ___
		int16 status = 0; // [MA] check exit status of system calls for proper error handling
		INFO.reset();
	  if (crop_arg.dbow.linehi!=0 && crop_arg.dbow.linelo!=0 && crop_arg.dbow.pixhi!=0 && crop_arg.dbow.pixlo!=0)
			 //<< " " << crop_arg.dbow.linelo - 1
			 //<< " " << crop_arg.dbow.pixlo - 1
		INFO << "tsx_dump_data.py " << crop_arg.filein1 << " " << crop_arg.fileout1 << " " << crop_arg.dbow.linelo << " " << crop_arg.dbow.linehi << " " << crop_arg.dbow.pixlo << " " << crop_arg.dbow.pixhi << ends;
	  else
		INFO << "tsx_dump_data.py " << crop_arg.filein1 << " " << crop_arg.fileout1 << ends;
	  String cmd = new String(new char[512]); // command string
	  cmd = INFO.get_str();
	  INFO.print("With following command TSX cosar data was cropped.");
	  INFO.print(cmd);
	  PROGRESS.print("system call may take some time...");
	  status =system(cmd); // this does the work
	  if (status != 0) // [MA] TODO make it a function
		{
		ERROR << "tsx_dump_data.py: failed with exit code: " << status;
		PRINT_ERROR(ERROR.get_str())
		throw(some_error);
		}
	  INFO.reset();
	  INFO.print();

	  // ====== Write results to scratchfile ======
	  ofstream scratchresfile = new ofstream("scratchres2raw", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile, "writeslc: scratchres2raw", __FILE__, __LINE__);
	  scratchresfile << "\n\n*******************************************************************\n";
	  if (crop_arg.fileid == MASTERID)
		scratchresfile << "*_Start_" << processcontrol[pr_m_crop];
	  if (crop_arg.fileid == SLAVEID)
		scratchresfile << "*_Start_" << processcontrol[pr_s_crop];
		// ______ updateslcimage greps these ______
	  scratchresfile << "\t\t\t" << crop_arg.idcrop << "\n*******************************************************************" << "\nData_output_file: \t\t\t\t" << crop_arg.fileout1 << "\nData_output_format: \t\t\t\t" << "complex_short" << "\nFirst_line (w.r.t. original_image): \t\t" << crop_arg.dbow.linelo << "\nLast_line (w.r.t. original_image): \t\t" << crop_arg.dbow.linehi << "\nFirst_pixel (w.r.t. original_image): \t\t" << crop_arg.dbow.pixlo << "\nLast_pixel (w.r.t. original_image): \t\t" << crop_arg.dbow.pixhi << "\n*******************************************************************";
	  if (crop_arg.fileid == MASTERID)
		scratchresfile << "\n* End_" << processcontrol[pr_m_crop] << "_NORMAL";
	  if (crop_arg.fileid == SLAVEID)
		scratchresfile << "\n* End_" << processcontrol[pr_s_crop] << "_NORMAL";
	  scratchresfile << "\n*******************************************************************" << "\n";
	  scratchresfile.close();
	} // END tsx_dump_data


	//***************************************************************
	// *    radarsat_dump_data                                        *
	// *                                                              *
	// * Inputfile in ceos slc format is converted to                 *
	// *  raw format outputfile.                                      *
	// * image is flipped if time direction is not INCREASE,INCREASE  *
	// * crop is thus done either geometric correct, or if specified  *
	// #%// Bert Kampes, 04-Aug-2004
	// * DBOW a bit strange ... ?
	// #%// Fix for DBOW for cropping radarsat images - Andy  Mar,2009
	// ***************************************************************
	public static void radarsat_dump_data(input_gen generalinput, input_crop crop_arg)
	  {
	  final int16 sizeb4 = 4;
	  final int16 sizeb1 = 1;
	  final int16 sizei4 = 4;
	  final int16 sizei6 = 6;
	  final int16 sizei8 = 8;
	  uint lenrec1; // length of general record1
	  uint lenrec2; // (nominal) length of data records
	  String c4 = new String(new char[5]); // correctly 9 for \0
	  String c6 = new String(new char[7]);
	  String c8 = new String(new char[9]);
	  // --- Check for RSAT #%// Bert Kampes, 02-Aug-2004 ---
	  uint rec_seq; // type B4
	  byte rec_sub1; // type B1
	  byte rec_type;
	  byte rec_sub2;
	  byte rec_sub3;

	  // ______ Write some info ______ 
	  TRACE_FUNCTION("radarsat_dump_data (Bert Kampes 04-Aug-2004)")
	  PROGRESS.print("Start cropping slc data for RADARSAT.");
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __X86PROCESSOR__
	  INFO.print("Swapping Big Endian (CEOS input) to Little Endian (your platform).");
	//  #else
	  INFO.print("NO byte swapping performed, you must be on Big Endian platform.");

	//  #endif
	  // ______ Open files ______ 
	  ifstream datfile;
	  RefObject<ifstream> TempRefObject = new RefObject<ifstream>(datfile);
	  openfstream(TempRefObject, crop_arg.filein1);
	  datfile = TempRefObject.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(datfile, crop_arg.filein1, __FILE__, __LINE__);

	  // ====== Get data such as recordlength ======
	  // --- RECORD 1 ---
	  DEBUG.print("record 1 of data file (ERS and RSAT).");
	  datfile.seekg(0,ios.beg);
	  datfile.read((char) rec_seq, sizeb4); // record number
	  rec_seq = ntohl(rec_seq); // bk 6 jul 2000, byteorder x86 machines.
	  datfile.read((char)&rec_sub1, sizeb1); // first record sub type code
	  datfile.read((char)&rec_type, sizeb1); // record type code
	  datfile.read((char)&rec_sub2, sizeb1); // second record sub type code
	  datfile.read((char)&rec_sub3, sizeb1); // third record sub type code
	  DEBUG.print("RSAT: Expecting record 1 with code {63,192,18,18}");
	  DEBUG << "rec_seq: " << rec_seq << "; rec_sub1: " << (int)rec_sub1 << "; rec_type: " << (int)rec_type << "; rec_sub2: " << (int)rec_sub2 << "; rec_sub3: " << (int)rec_sub3;
	  DEBUG.print();
	  if ((int)rec_sub1 ==63 && (int)rec_type ==192 && (int)rec_sub2 ==18 && (int)rec_sub3 ==18)
		DEBUG.print("This is the expected record with code {63,192,18,18}");
	  else
		WARNING.print("This is NOT the expected record with code {63,192,18,18}");
	  datfile.seekg(8,ios.beg);
	  datfile.read((char) lenrec1, sizeb4); // length of record
	  lenrec1 = ntohl(lenrec1); // bk 6 jul 2000, byteorder x86 machines.
	  DEBUG.print("RSAT: Expecting record 1 with length 16252");
	  DEBUG << "radarsat_dump_data::record 1: start at: " << 0 << "; length: " << lenrec1;
	  DEBUG.print();

	  // --- Get some info in RECORD 1 ---
	  datfile.seekg(180,ios.beg);
	  datfile.read((char)&c6, sizei6); // number of SAR DATA records (lines)
		c6.charAt(6)='\0';
		final uint numdatarec = Integer.parseInt(c6);
	  DEBUG << "numdatarec: " << numdatarec;
	  DEBUG.print();
	  //datfile.seekg(186,ios::beg);
	  datfile.read((char)&c6, sizei6); // SAR DATA record length
		c6.charAt(6)='\0';
		final uint lendatarec2 = Integer.parseInt(c6);
	  DEBUG << "lendatarec2: " << lendatarec2;
	  DEBUG.print();
	  datfile.seekg(232,ios.beg); // SAR Related data
	  datfile.read((char)&c4, 4);
		c4.charAt(4)='\0';
		final uint numchannels = Integer.parseInt(c4);
	  DEBUG << "numchannels: " << numchannels;
	  DEBUG.print();

	  datfile.read((char)&c8, sizei8);
		c8.charAt(8)='\0';
		uint numlines = Integer.parseInt(c8);
	  DEBUG << "numlines: " << numlines;
	  DEBUG.print();

	  datfile.read((char)&c4, sizei4);
		c4.charAt(4)='\0';
		final uint leftborder = Integer.parseInt(c4);
	  DEBUG << "leftborder: " << leftborder;
	  DEBUG.print();
	  datfile.read((char)&c8, sizei8); // number of pixels
		c8.charAt(8)='\0';
		uint numpixels = Integer.parseInt(c8);
	  DEBUG << "numpixels: " << numpixels;
	  DEBUG.print();
	  datfile.read((char)&c4, sizei4);
		c4.charAt(4)='\0';
		final uint rightborder = Integer.parseInt(c4);
	  DEBUG << "rightborder: " << rightborder;
	  DEBUG.print();
	  datfile.read((char)&c4, sizei4);
		c4.charAt(4)='\0';
		final uint topborder = Integer.parseInt(c4);
	  DEBUG << "topborder: " << topborder;
	  DEBUG.print();
	  datfile.read((char)&c4, sizei4);
		c4.charAt(4)='\0';
		final uint bottomborder = Integer.parseInt(c4);
	  DEBUG << "bottomborder: " << bottomborder;
	  DEBUG.print();

	  datfile.seekg(280,ios.beg); // Record data
	  datfile.read((char)&c8, sizei8);
		c8.charAt(8)='\0';
		final uint numbytesdata = Integer.parseInt(c8);
	  DEBUG << "numbytesdata: " << numbytesdata;
	  DEBUG.print();

	// ______ Check with previous section ______
	  if (numlines != numdatarec)
		{
		WARNING << "code 904: Number of lines seems not to be consistent in file: " << crop_arg.filein1 << " : " << numlines << " != " << numdatarec;
		WARNING.print();
		WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
		}
	  if ((numbytesdata / 4) != numpixels)
		{
		WARNING << "code 904: Number of pixels seems to be inconsistent in file: " << crop_arg.filein1 << ": " << numpixels << " != " << (numbytesdata / 4);
		WARNING.print();
		WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
		}


	// ====== Start copy input to output (raw) format with buffer======
	// ______ Check and process optional offset parameters______
	// ______ Lcnlow is corner line, lcnhi is other corner, pcnlow, pixel coord. low etc.
	  uint linestart = 1; // counters for loops, first=1;
	  uint lineend = numlines;
	  uint pixelstart = 1; // first pix is 1
	  uint pixelend = numpixels; // only for resultfile
	  uint orig_numlines = numlines;
	  uint orig_numpixels = numpixels;

	  if (crop_arg.dbow.linehi!=0 && crop_arg.dbow.linelo!=0 && crop_arg.dbow.pixhi!=0 && crop_arg.dbow.pixlo!=0)
		{
		WARNING.print("cropping data may be difficult, due to INC/DECREASE storage");
		window tempdbow = new window(crop_arg.dbow.linelo, crop_arg.dbow.linehi, crop_arg.dbow.pixlo, crop_arg.dbow.pixhi);
		if (crop_arg.dbow.linehi>numlines)
		  {
		  WARNING << "Specified input DBOW linehi > numlines: " << crop_arg.dbow.linehi << " > " << numlines << ". I set linehi = " << numlines;
		  WARNING.print();
		  tempdbow.linehi=numlines;
		  }
		if (crop_arg.dbow.pixhi>numpixels)
		  {
		  WARNING << "Specified input DBOW pixhi > numpixels: " << crop_arg.dbow.pixhi << " > " << numpixels << ". I set pixhi = " << numpixels;
		  WARNING.print();
		  tempdbow.pixhi=numpixels;
		  }
		// ______ Only hi values are possibly adapted, low is a constant ______
		numlines = tempdbow.linehi - crop_arg.dbow.linelo + 1;
		numpixels = tempdbow.pixhi - crop_arg.dbow.pixlo + 1;

		linestart = crop_arg.dbow.linelo;
		lineend = tempdbow.linehi;
		pixelstart = crop_arg.dbow.pixlo;
		pixelend = tempdbow.pixhi; // only for resultfile
		}


	  // --- Find out if data is stored increasing in line/pix ---
	  // --- RECORD 2 ---
	  uint startrec2 = lenrec1;
	  DEBUG.print("record 2 of data file (RSAT).");
	  datfile.seekg(startrec2,ios.beg);
	  datfile.read((char) rec_seq, sizeb4); // record number
	  rec_seq = ntohl(rec_seq); // bk 6 jul 2000, byteorder x86 machines.
	  datfile.read((char)&rec_sub1, sizeb1); // first record sub type code
	  datfile.read((char)&rec_type, sizeb1); // record type code
	  datfile.read((char)&rec_sub2, sizeb1); // second record sub type code
	  datfile.read((char)&rec_sub3, sizeb1); // third record sub type code
	  DEBUG.print("RSAT: Expecting record 2 with code {50,11,18,20}");
	  DEBUG << "rec_seq: " << rec_seq << "; rec_sub1: " << (int)rec_sub1 << "; rec_type: " << (int)rec_type << "; rec_sub2: " << (int)rec_sub2 << "; rec_sub3: " << (int)rec_sub3;
	  DEBUG.print();
	  if ((int)rec_sub1 ==50 && (int)rec_type ==11 && (int)rec_sub2 ==18 && (int)rec_sub3 ==20)
		DEBUG.print("This is the expected record with code {50,11,18,20}");
	  else
		WARNING.print("This is NOT the expected record with code {50,11,18,20}");
	  datfile.seekg(startrec2+8,ios.beg);
	  datfile.read((char) lenrec2, sizeb4); // length of record
	  lenrec2 = ntohl(lenrec2); // bk 6 jul 2000, byteorder x86 machines.
	  DEBUG.print("RSAT: Expecting record 2 with length variable");
	  DEBUG << "radarsat_dump_data::record 2: start at: " << startrec2 << "; length: " << lenrec2;
	  DEBUG.print();
	  uint startrec3 = lenrec1+lenrec2;
	  uint startrecN = lenrec1+(numlines-1)*lenrec2; // start of last record
	  // --- azimuth time to first line (depends on decrease/increase): ---
	  uint zdmsecofday1 = 99999; // B4
	  uint zdmsecofday2 = 99999; // B4
	  uint zdmsecofdayN = 99999; // B4
	  datfile.seekg(startrec2+44,ios.beg);
	  datfile.read((char) zdmsecofday1, sizeb4); // range to first pix
	  zdmsecofday1 = ntohl(zdmsecofday1); // bk 6 jul 2000, byteorder x86 machines.
	  datfile.seekg(startrec3+44,ios.beg);
	  datfile.read((char) zdmsecofday2, sizeb4); // range to first pix
	  zdmsecofday2 = ntohl(zdmsecofday2); // bk 6 jul 2000, byteorder x86 machines.
	  datfile.seekg(startrecN+44,ios.beg);
	  datfile.read((char) zdmsecofdayN, sizeb4); // range to first pix
	  zdmsecofdayN = ntohl(zdmsecofdayN); // bk 6 jul 2000, byteorder x86 machines.
	  INFO << "zdmsecofday1: " << zdmsecofday1;
	  INFO.print();
	  DEBUG << "zdmsecofday2: " << zdmsecofday2;
	  DEBUG.print();
	  INFO << "zdmsecofdayN: " << zdmsecofdayN;
	  INFO.print();
	  boolean increasing_line = true; // assume this
	  // --- I assume linestart is smaller than end, so swap them if required ---
	  if (zdmsecofday1 < zdmsecofdayN) // increase, use ZD time of first line
		{
		INFO.print("INCREASE line direction detected (OK).");
		}
	  else // decreasing lines: flip up-down
		{
		WARNING.print("DECREASE line direction detected: I will flip up-down the RSAT data");
		increasing_line = false;
		}

	  // --- range time to first pixel is distance -------------------------
	  uint range1st = 99999;
	  uint rangelst = 99999;
	  datfile.seekg(startrec2+64,ios.beg);
	  datfile.read((char) range1st, sizeb4); // range to first pix
	  range1st = ntohl(range1st); // bk 6 jul 2000, byteorder x86 machines.
	  datfile.seekg(startrec2+72,ios.beg);
	  datfile.read((char) rangelst, sizeb4); // range to last pix
	  rangelst = ntohl(rangelst); // bk 6 jul 2000, byteorder x86 machines.
	  INFO << "range1st: " << range1st;
	  INFO.print();
	  INFO << "rangelst: " << rangelst;
	  INFO.print();
	  boolean increasing_pix = true; // assume this
	  if (range1st < rangelst) // increase
		{
		INFO.print("INCREASE pixel direction detected (OK).");
		}
	  else // decrease: flip data left-right
		{
		WARNING.print("DECREASE pixel direction detected: I will flip left-right the RSAT data");
		increasing_pix = false;
		}



	  // ====== Process requested lines ======
	  ofstream datoutfile;
	  RefObject<ofstream> TempRefObject2 = new RefObject<ofstream>(datoutfile);
	  openfstream(TempRefObject2, crop_arg.fileout1, generalinput.overwrit);
	  datoutfile = TempRefObject2.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(datoutfile, crop_arg.fileout1, __FILE__, __LINE__);

	  // ______ info on data, to avoid X86 problems ______
	  // ______ according to RSAT CEOS specs, byte 192-195 is first complex pixel, etc. ______
	  datfile.seekg(lenrec1+192,ios.beg);
	  matrix TMPSHORT = new matrix(1, 2);
	  datfile >> TMPSHORT; // read in first complex pixel for test
	  real8 tmpmag = Math.sqrt(real8(int16(ntohs(TMPSHORT(0,0)))*int16(ntohs(TMPSHORT(0,0)))) + real8(int16(ntohs(TMPSHORT(0,1)))*int16(ntohs(TMPSHORT(0,1)))));
	  DEBUG << "First complex element in datafile: (" << int16(ntohs(TMPSHORT(0,0))) << "," << int16(ntohs(TMPSHORT(0,1))) << "); mag = " << tmpmag;
	  DEBUG.print();
	  if (tmpmag > 10000.)
		{
		WARNING.print(DEBUG.get_str());
		WARNING.print("this is a byteorder problem on X86? (use ntohs)");
		}
	  DEBUG << "TEST: (realpart): " << TMPSHORT(0,0) << ", (imagpart): " << TMPSHORT(0,1);
	  DEBUG.print();
	  DEBUG << "TEST: htons(realpart): " << htons(TMPSHORT(0,0)) << ", htons(imagpart): " << htons(TMPSHORT(0,1));
	  DEBUG.print();
	  DEBUG << "TEST: ntohs(realpart): " << ntohs(TMPSHORT(0,0)) << ", ntohs(imagpart): " << ntohs(TMPSHORT(0,1));
	  DEBUG.print();
	  DEBUG << "TEST: short int(ntohs(realpart)): " << int16(ntohs(TMPSHORT(0,0))) << ", (imagpart): " << int16(ntohs(TMPSHORT(0,1)));
	  DEBUG.print();


	  // --- Simple way pix by pix reading so we can easily flip if req. ---
	  datfile.seekg(lenrec1+(linestart-1)*lendatarec2+8,ios.beg);
	  datfile.read((char) lenrec2, sizeb4); // length of first record
	  lenrec2 = ntohl(lenrec2); // bk 6 jul 2000, byteorder x86 machines.
	  if (lenrec2 != lendatarec2)
		{
		ERROR << "code 904: Length of datarecords seems to be inconsistent in file: " << crop_arg.filein1 << ": " << lenrec2 << " != " << lendatarec2;
		WARNING.print(ERROR.get_str());
		ERROR.reset();
		}


	  // ______ Note complex<short> not in ANSI c ______
	  int32 percentage = 0; // initialization
	  final int32 TENPERCENT = (5+numlines)/10; // number of lines
	  for (register int32 linecnt =linestart; linecnt<=lineend; linecnt++)
		{
		// 03/2009 AH
		//int32 line2read = (increasing_line==true) ? linecnt : lineend-(linecnt-linestart);
		int32 line2read = (increasing_line ==true) ? linecnt : orig_numlines-linecnt+1;
		for (register int32 pixcnt =pixelstart; pixcnt<=pixelend; pixcnt++)
		  {
		  // 03/2009 AH
		  //int32 pix2read = (increasing_pix==true) ? pixcnt : pixelend-(pixcnt-pixelstart);
		  int32 pix2read = (increasing_pix ==true) ? pixcnt : orig_numpixels-pixcnt+1;
		  // --- set pointer before complex pixel to read ---
		  uint tmpstart = lenrec1+192+(line2read-1)*lendatarec2+(pix2read-1)*4;
		  datfile.seekg(tmpstart,ios.beg);
		  datfile >> TMPSHORT;
		  // ______ BK 13 July 2000: swapbytes for X86 (intel) linux cpus ______
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//	  #if __X86PROCESSOR__
		  TMPSHORT(0,0) = int16(ntohs(TMPSHORT(0,0))); // byteswap
		  TMPSHORT(0,1) = int16(ntohs(TMPSHORT(0,1))); // byteswap datoutfile << TMPSHORT;
	//	  #else
		  datoutfile << TMPSHORT;
	//	  #endif
		  }
		if (!((linecnt-linestart)%TENPERCENT))
		  {
		  PROGRESS << "radarsat_dump_data: " << setw(3) << percentage << "%";
		  PROGRESS.print();
		  percentage += 10;
		  }
		}
	  datfile.close(); // close files
	  datoutfile.close();



	// ====== Write results to scratchfile ======
	  ofstream scratchresfile = new ofstream("scratchres2raw", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile, "writeslc: scratchres2raw", __FILE__, __LINE__);
	  scratchresfile << "\n\n*******************************************************************\n";
	  if (crop_arg.fileid == MASTERID)
		scratchresfile << "*_Start_" << processcontrol[pr_m_crop];
	  if (crop_arg.fileid == SLAVEID)
		scratchresfile << "*_Start_" << processcontrol[pr_s_crop];
		// ______ updateslcimage greps these ______
	  scratchresfile << "\t\t\t" << crop_arg.idcrop << "\n*******************************************************************" << "\nData_output_file: \t\t\t\t" << crop_arg.fileout1 << "\nData_output_format: \t\t\t\t" << "complex_short" << "\nFirst_line (w.r.t. original_image): \t\t" << linestart << "\nLast_line (w.r.t. original_image): \t\t" << lineend << "\nFirst_pixel (w.r.t. original_image): \t\t" << pixelstart << "\nLast_pixel (w.r.t. original_image): \t\t" << pixelend << "\n*******************************************************************";
	  if (crop_arg.fileid == MASTERID)
		scratchresfile << "\n* End_" << processcontrol[pr_m_crop] << "_NORMAL";
	  if (crop_arg.fileid == SLAVEID)
		scratchresfile << "\n* End_" << processcontrol[pr_s_crop] << "_NORMAL";
	  scratchresfile << "\n*******************************************************************" << "\n";
	  scratchresfile.close();

	// ______ Tidy up do checks here ______
	  if (numchannels != 1) // ??
		{
		WARNING << "code 904: Number of channels in file: " << crop_arg.filein1 << " = " << numchannels << " != 1 ";
		WARNING.print();
		WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
		}
	  if (bottomborder != topborder != leftborder != rightborder != 0)
		{
		WARNING << "code 904: Not implemented: offset border: left,right,bottom,top: " << leftborder << "," << rightborder << "," << bottomborder << "," << topborder << " in file: " << crop_arg.filein1;
		WARNING.print();
		WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
		}
	  PROGRESS.print("radarsat_dump_data: 100%");
	  } // END radarsat_dump_data




	//____RaffaeleNutricato START MODIFICATION SECTION 2
	//***************************************************************
	// *    OversampleSLC                                             *
	// *                                                              *
	// * Oversamples the SLC by an integer factor.                    *
	// * For now only range oversampling is performed.                *
	// *    Raffaele Nutricato, 12-Jan-2004                           *
	// * Azimuth oversampling with factor 2 added.                    *
	// *    Bert Kampes, 30-Jul-2005                                  *
	// **************************************************************
	public static void OversampleSLC(input_gen generalinput, slcimage imageinfo, input_oversample oversampleinput, int16 fileid)
	  {
	  TRACE_FUNCTION("OversampleSLC (Raffaele Nutricato 12-Jan-2004)")
	  String infile = new String(new char[EIGHTY]); // Input file which is master/slave.raw renamed as .old
	  String outfile = new String(new char[EIGHTY]); // Output file which is the oversampled version.
	  infile = imageinfo.file;
	  outfile = oversampleinput.fileoutovs;
	  final int32 OsrRange = oversampleinput.OsrRange; // Range oversampling ratio.
	  final int32 OsrAzimuth = oversampleinput.OsrAzimuth; // Azimuth oversampling ratio.
	  final int32 FilterSize = oversampleinput.FilterSize; // Length of the kernel for the oversampling in range.
	  if (OsrAzimuth!=1 && OsrAzimuth!=2)
		{
		ERROR.print("oversampling in azimuth: only factor 2");
		throw(some_error);
		}

	  // ______Open input file_____
	  ifstream ifile;
	  RefObject<ifstream> TempRefObject = new RefObject<ifstream>(ifile);
	  openfstream(TempRefObject, infile);
	  ifile = TempRefObject.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(ifile, infile, __FILE__, __LINE__);

	  // ______Open output file_____
	  ofstream ofile;
	  RefObject<ofstream> TempRefObject2 = new RefObject<ofstream>(ofile);
	  openfstream(TempRefObject2, outfile, generalinput.overwrit);
	  ofile = TempRefObject2.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(ofile, outfile, __FILE__, __LINE__);

	  // ______ Compute the size of original cropped image ______
	  final int32 numlines = imageinfo.currentwindow.lines();
	  final int32 numpixels = imageinfo.currentwindow.pixels();

	  // ______ Define the accuracy of the digital signal processings ______
	  // ______ in range (in azimuth I use float, BK ______

	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//  #define RN_DSP_ACCURACY double
	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//  #define RN_DSP_CPXACCURACY complr8
	  // ______ Interpolation kernel section ______
	  matrix LINE_IN = new matrix(1, OsrRange*(numpixels-1)+1); // zero alternating
	  matrix LINE_OUT = new matrix(1, OsrRange *numpixels); // ovs line
	  final int32 interp_size = FilterSize * OsrRange - 1;
	  matrix INTERP_KERNEL = new matrix(1, interp_size);

	  // ______ Generate the range interpolator impulse response ______
	  final _double invosr = 1.0/(_double)OsrRange;
	  _double interpsamplepos = invosr - (_double)FilterSize/2.0;
	  for (int32 i =0; i<interp_size; i++)
		{
		INTERP_KERNEL(0,i) = complr8(sinc(interpsamplepos),0);
		interpsamplepos += invosr;
		}

	  // ______ Normalize kernel (BK) ______
	  // ______ (i would use 6pn real4 rc kernel, but ok.)
	  _double qsum = 0.0;
	  for (int32 i =0; i<interp_size; ++i)
		qsum += real(INTERP_KERNEL(0,i));
	  INTERP_KERNEL /= qsum; // complex divided by float

	  // ====== Variables for oversampling in azimuth [BK] ======
	  // ______ Factor two oversampling is supported in azimuth ______
	  // ______ We use a buffer of 6 lines (oversampled in range) ______
	  // ______ these samples are interpolated in between using a ______
	  // ______ raised cosine kernel of 6 points (claimed 0.9999) ______
	  // ______ A big matrix is used that contains all shifted kernels ______
	  final int32 NP_kernel_az = 6;
	  matrix<complr4> BUFFER_AZ = new matrix(NP_kernel_az, LINE_OUT.pixels()); // rotating buffer
	  matrix<complr4> KERNEL_AZ = new matrix(NP_kernel_az, LINE_OUT.pixels()); // raised cosine
	  matrix<complr4> LINE_OUT2 = new matrix(1, LINE_OUT.pixels()); // interpolated in azi.
	  if (OsrAzimuth!=1)
		{
		INFO.print("Initializing azimuth kernel");
		final real4 CHI = imageinfo.prf/imageinfo.abw; // oversampling factor az
		matrix<real4> x_axis = new matrix(NP_kernel_az, 1);
		for (int32 i =0; i<NP_kernel_az; ++i)
		  x_axis(i,0) = -NP_kernel_az/2 + 0.5 + i; // [-2.5 -1.5 -0.5 0.5 1.5 2.5]
		matrix<complr4> tmp_kernel = mat2cr4rc_kernel(x_axis, CHI, NP_kernel_az);
		DEBUG.print("Normalizing kernel");
		real4 qsum = 0.0;
		for (int32 i =0; i<NP_kernel_az; ++i)
		  qsum += real(tmp_kernel(i,0));
		tmp_kernel /= qsum; // complr4
		DEBUG.print("Shifting kernels with Doppler");
		// ___ Doppler centroid is function of range only ____
		// ___ to shift spectrum of convolution kernel to fDC of data, multiply
		// ___ in the space domain with a phase trend of -2pi*t*fdc/prf
		// ___ (to shift back (no need) you would use +fdc), see manual;
		for (int32 x =0; x<LINE_OUT.pixels(); ++x)
		  {
		  final real4 pix = real4(imageinfo.currentwindow.pixlo)+real4(x)/2.0;
		  final real4 slope = 2.0 *PI *imageinfo.pix2fdc(pix)/imageinfo.prf;
		  for (int32 i =0; i<NP_kernel_az; ++i)
			{
			// ___ Modify kernel, shift spectrum to fDC ___
			final real4 t = x_axis(i,0)*slope; // the phase ramp
			KERNEL_AZ(i,x) = tmp_kernel(i,0)*complr4(Math.cos(t),-Math.sin(t)); // note '-' (see manual)
			}
		  }
		} //if azimuth ovs


	  // ====== Loop on the lines to oversample ======
	  final int32 TEN = 10;
	  final int32 TENPERCENT = Math.floor((0.5 *TEN+numlines)/TEN); // number of lines
	  int32 percentage = 0; // initialization
	  for (register int32 linecnt =0; linecnt<numlines; linecnt++)
		{
		if (!(linecnt%TENPERCENT))
		  {
		  PROGRESS << "OVERSAMPLESLC: " << setw(3) << percentage << "%";
		  PROGRESS.print();
		  percentage += TEN;
		  }
		// ______ Read input data in larger line ______
		switch (imageinfo.formatflag)
		  {
		  case FORMATCR4:
			{
			matrix<real4> bufferrreal4 = new matrix(1, 1);
			matrix<real4> bufferrimag4 = new matrix(1, 1);
			for (int32 ii =0; ii<numpixels; ++ii)
			  {
			  ifile >> bufferrreal4;
			  ifile >> bufferrimag4;
			  // ______Generate a zero filled copy of LINE______
			  // RN LINE_IN must be cleaned!!!
			  LINE_IN(0,OsrRange *ii) = complr8(bufferrreal4(0,0),bufferrimag4(0,0));
			  }
			break;
			}
		  // ______ Convert first to ci2 before writing to file ______
		  case FORMATCI2:
			{
			matrix<int16> bufferrealint16 = new matrix(1, 1);
			matrix<int16> bufferimagint16 = new matrix(1, 1);
			for (int32 ii =0; ii<numpixels; ++ii)
			  {
			  ifile >> bufferrealint16;
			  ifile >> bufferimagint16;
			  // ______Generate a zero filled copy of LINE______
			  // RN LINE_IN must be cleaned!!!
			  LINE_IN(0,OsrRange *ii) = complr8(bufferrealint16(0,0),bufferimagint16(0,0));
			  }
			break;
			}
		  default:
			PRINT_ERROR("Unknown input format for the cropped image.");
			throw(unhandled_case_error);
		  } // end switch reading input line

		// ______Spatial convolution between LINE_IN and INTERP_KERNEL______
		int jmin;
		int jmax;
		int RN_k;
		int minpos;
		int maxpos;
		RN_k = 0;
		minpos = (interp_size-1)/2;
		maxpos = (interp_size-1)/2 + (LINE_IN.pixels() - 1) + OsrRange - 1;
		for (int ii =minpos; ii<=maxpos; ii++)
		  {
		  LINE_OUT(0,RN_k) = 0;
		  jmin = max(int32(0), int32(ii-interp_size+1));
		  jmax = min(int32(ii),int32(LINE_IN.pixels()-1));
		  for (int j =jmin; j<=jmax; j++)
			LINE_OUT(0,RN_k) += LINE_IN(0,j) * INTERP_KERNEL(0,ii-j);
		  RN_k++;
		  }

		// ______ Oversample in azimuth ______
		// ______ e.g., kernel is 6 points, buffer is 6 lines ______
		// ______ if linecnt==5, then buffer is filled with ______
		// ______ first 6 lines of file.  LINE_OUT2 is the ______
		// ______ interpolated line at BUFFER[2.5,*], i.e., ______
		// ______ to write in correct order each time we write ______
		// ______ BUFFER[2,*] and LINE_OUT2[*];  not 5 and 5.5 ______
		// ______ To correct for this offset, we do: ______
		// ______ linecnt==0: do not write (add at end)
		// ______ linecnt==1: do not write (add at end)
		// ______ linecnt==2: do not write (add at end)
		// ______ linecnt==3: write 0 and 0.5
		// ______ linecnt==4: write 1 and 1.5
		// ______ linecnt==5: write 2 and 2.5
		// ______ linecnt==6: write etc.
		// ______ linecnt==last: write last + 6 extra lines.
		if (OsrAzimuth!=1) //i.e., 2
		  {
		  // ______ Rotating BUFFER_AZ ______
		  // ______ I dont trust shifting pointers so I copy lines ______
		  // ______ this is slower, but guarantees continuous in mem. ______
		  for (int32 x =0; x<LINE_OUT.pixels(); ++x)
			for (int32 i =0; i<NP_kernel_az-1; ++i)
			  BUFFER_AZ(i,x) = BUFFER_AZ(i+1,x); //rotate buffer by copy
		  for (int32 x =0; x<LINE_OUT.pixels(); ++x)
			BUFFER_AZ(NP_kernel_az-1,x) = complr4(LINE_OUT(0,x)); // add new line
		  // ______ Oversample in azimuth (interpolate) ______
		  LINE_OUT2 = sum(dotmult(BUFFER_AZ, KERNEL_AZ), 1); // at half+0.5
		  }

		// ______ Write LINE_OUT in the output file ______
		switch (oversampleinput.oformatflag)
		  {
		  // ______ Convert first to cr4 before writing to file ______
		  case FORMATCR4:
			{
			complr4 buffercr4;
			if (OsrAzimuth!=1)
			  {
			  if (linecnt>=NP_kernel_az/2)
				{
				// _____ write line 2.0 ______
				for (int ii =0; ii<LINE_OUT.pixels(); ii++)
				  {
				  buffercr4 = complr4(BUFFER_AZ(NP_kernel_az/2-1,ii));
				  ofile.write((char) buffercr4, sizeof(complr4));
				  }
				// _____ write line 2.5 ______
				for (int ii =0; ii<LINE_OUT.pixels(); ii++)
				  {
				  buffercr4 = complr4(LINE_OUT2(0,ii));
				  ofile.write((char) buffercr4, sizeof(complr4));
				  }
				// _____ write additional zero lines at end ______
				if (linecnt==numlines-1) //write extra lines at end
				  {
				  buffercr4 = complr4(0.0,0.0);
				  for (int i =0; i<NP_kernel_az; ++i)
					for (int ii =0; ii<LINE_OUT.pixels(); ii++)
					  ofile.write((char) buffercr4, sizeof(complr4));
				  }
				}
			  }
			else // no azimuth oversampling
			  {
			  for (int ii =0; ii<LINE_OUT.pixels(); ii++)
				{
				buffercr4 = complr4(LINE_OUT(0,ii));
				ofile.write((char) buffercr4, sizeof(complr4));
				}
			  }
			break; //switch
			}
		  // ______ Convert first to ci2 before writing to file ______
		  case FORMATCI2:
			{
			compli16 bufferci16;
			if (OsrAzimuth!=1)
			  {
			  if (linecnt>=NP_kernel_az/2)
				{
				// _____ write line 2.0 ______
				for (int ii =0; ii<LINE_OUT.pixels(); ii++)
				  {
				  bufferci16 = cr4toci2(complr4(BUFFER_AZ(NP_kernel_az/2-1,ii)));
				  ofile.write((char) bufferci16, sizeof(compli16));
				  }
				// _____ write line 2.5 ______
				for (int ii =0; ii<LINE_OUT.pixels(); ii++)
				  {
				  bufferci16 = cr4toci2(complr4(LINE_OUT2(0,ii)));
				  ofile.write((char) bufferci16, sizeof(compli16));
				  }
				// _____ write extra line at end ______
				if (linecnt==numlines-1) // write extra lines at end
				  {
				  bufferci16 = compli16(0,0);
				  for (int i =0; i<NP_kernel_az; ++i)
					for (int ii =0; ii<LINE_OUT.pixels(); ii++)
					  ofile.write((char) bufferci16, sizeof(compli16));
				  }
				}
			  }
			else // no azimuth oversampling
			  {
			  for (int ii =0; ii<LINE_OUT.pixels(); ii++)
				{
				bufferci16 = cr4toci2(complr4(LINE_OUT(0,ii)));
				ofile.write((char) bufferci16, sizeof(compli16));
				}
			  }
			break;
			}
			default:
			  PRINT_ERROR("Unknown output format for the oversampled image.");
			  throw(unhandled_case_error);
		  } // end switch writing output
		LINE_IN.clean();
		} // end for loop over lines


	  // ______ Close files ______
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #undef RN_DSP_ACCURACY
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #undef RN_DSP_CPXACCURACY
	  ifile.close();
	  ofile.close();


	  // ====== Write results to scratchfile ======
	  ofstream scratchresfile = new ofstream("scratchoversample", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile, "OversampleSLC: scratchoversample", __FILE__, __LINE__);
	  scratchresfile << "\n\n*******************************************************************\n";
	  if (fileid == MASTERID)
	  {
		scratchresfile << "*_Start_" << processcontrol[pr_m_oversample];
		scratchresfile << "\t\t\t" << "master";
	  }
	  if (fileid == SLAVEID)
		{
		scratchresfile << "*_Start_" << processcontrol[pr_s_oversample];
		scratchresfile << "\t\t\t" << "slave";
		}
	  // ______ updateslcimage greps these ______
	  // ______ [BK] write new file size, etc. here ______
	  scratchresfile << "\n*******************************************************************" << "\nData_output_file: \t\t\t\t" << outfile << "\nData_output_format: \t\t\t\t";
		if (oversampleinput.oformatflag==FORMATCR4)
		  scratchresfile << "complex_real4";
		if (oversampleinput.oformatflag==FORMATCI2)
		  scratchresfile << "complex_short";
	  // ______ updateslcimage does not grep these ______
	  scratchresfile << "\nFirst_line (w.r.t. ovs_image):       \t\t" << (imageinfo.currentwindow.linelo-1)*OsrAzimuth+1 << "\nLast_line (w.r.t. ovs_image):        \t\t" << imageinfo.currentwindow.linehi *OsrAzimuth << "\nFirst_pixel (w.r.t. ovs_image):      \t\t" << (imageinfo.currentwindow.pixlo-1)*OsrRange+1 << "\nLast_pixel (w.r.t. ovs_image):       \t\t" << imageinfo.currentwindow.pixhi *OsrRange << "\nMultilookfactor_azimuth_direction:   \t\t" << 1.0/OsrAzimuth << "\nMultilookfactor_range_direction:     \t\t" << 1.0/OsrRange << "\nNumber of lines (oversampled):       \t\t" << OsrAzimuth *numlines << "\nNumber of pixels (oversampled):      \t\t" << OsrRange *numpixels << "\n#First_line (w.r.t. original_image): \t\t" << imageinfo.currentwindow.linelo << "\n#Last_line (w.r.t. original_image):  \t\t" << imageinfo.currentwindow.linehi << "\n#First_pixel (w.r.t. original_image):\t\t" << imageinfo.currentwindow.pixlo << "\n#Last_pixel (w.r.t. original_image): \t\t" << imageinfo.currentwindow.pixhi;
	  scratchresfile << "\n*******************************************************************";
	  if (fileid == MASTERID)
		scratchresfile << "\n* End_" << processcontrol[pr_m_oversample] << "_NORMAL";
	  if (fileid == SLAVEID)
		scratchresfile << "\n* End_" << processcontrol[pr_s_oversample] << "_NORMAL";
	  scratchresfile << "\n*******************************************************************" << "\n";
	  scratchresfile.close();

	  PROGRESS.print("OVERSAMPLESLC: 100%");
	  } // END OversampleSLC
	//____RaffaeleNutricato END MODIFICATION SECTION 2

	// Copy code from writeslc
	// do some modification
	// Modified by LG for reading ALOS Fine
	// To read complex real4 data 
	public static void palsar_fine_dump_data(input_gen generalinput, input_crop crop_arg, int checklines)
	{

	  final int16 sizeb4 = 4;
	  final int16 sizei4 = 4;
	  final int16 sizei6 = 6;
	  final int16 sizei8 = 8;
	  uint lenrec1; // length of general record1
	  uint lenrec2; // (nominal) length of data records
	  String c4 = new String(new char[5]); // correctly 9 for \0
	  String c6 = new String(new char[7]);
	  String c8 = new String(new char[9]);

	  // ______ Write some info ______ 
	  TRACE_FUNCTION("palsar_fine_dump_data (LG 28-Dec-2005)")
	  PROGRESS.print("Start cropping slc data.");

	  // ______ Open files ______ 
	  ifstream datfile;
	  RefObject<ifstream> TempRefObject = new RefObject<ifstream>(datfile);
	  openfstream(TempRefObject, crop_arg.filein1);
	  datfile = TempRefObject.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(datfile, crop_arg.filein1, __FILE__, __LINE__);

	  // ====== Get data such as recordlength ======
	  datfile.seekg(8,ios.beg);
	  datfile.read((char) lenrec1, sizeb4); // length of record1
	  lenrec1 = ntohl(lenrec1); // bk 6 jul 2000, byteorder x86 machines.
	  DEBUG.print("record 1 of data file.");
	  if (lenrec1 != 720) // PALSAR data = 720 + Rec_bytes*lines
		{
		WARNING << "palsar_fine_dump_data : length of record 1 = \"" << lenrec1 << "\"; expected \"720\" for PALSAR FINE SLC (CEOS, full scene).";
		WARNING.print();

		}

	  datfile.seekg(180,ios.beg);
	  datfile.read((char)&c6, sizei6); // number of SAR DATA records (lines)
		c6.charAt(6)='\0';
		final uint numdatarec = Integer.parseInt(c6);
	  DEBUG << "numdatarec: " << numdatarec;
	  DEBUG.print();
	  //datfile.seekg(186,ios::beg);
	  datfile.read((char)&c6, sizei6); // SAR DATA record length
		c6.charAt(6)='\0';
		final uint lendatarec2 = Integer.parseInt(c6);
	  DEBUG << "lendatarec2: " << lendatarec2;
	  DEBUG.print();
	  datfile.seekg(232,ios.beg); // SAR Related data
	  datfile.read((char)&c4, 4);
		c4.charAt(4)='\0';
		final uint numchannels = Integer.parseInt(c4);
	  DEBUG << "numchannels: " << numchannels;
	  DEBUG.print();

	  datfile.read((char)&c8, sizei8);
		c8.charAt(8)='\0';
		uint numlines = Integer.parseInt(c8);
	  DEBUG << "numlines: " << numlines;
	  DEBUG.print();

	  datfile.read((char)&c4, sizei4);
		c4.charAt(4)='\0';
		final uint leftborder = Integer.parseInt(c4);
	  DEBUG << "leftborder: " << leftborder;
	  DEBUG.print();
	  datfile.read((char)&c8, sizei8); // number of pixels
		c8.charAt(8)='\0';
		uint numpixels = Integer.parseInt(c8);
	  DEBUG << "numpixels: " << numpixels;
	  DEBUG.print();
	  datfile.read((char)&c4, sizei4);
		c4.charAt(4)='\0';
		final uint rightborder = Integer.parseInt(c4);
	  DEBUG << "rightborder: " << rightborder;
	  DEBUG.print();
	  datfile.read((char)&c4, sizei4);
		c4.charAt(4)='\0';
		final uint topborder = Integer.parseInt(c4);
	  DEBUG << "topborder: " << topborder;
	  DEBUG.print();
	  datfile.read((char)&c4, sizei4);
		c4.charAt(4)='\0';
		final uint bottomborder = Integer.parseInt(c4);
	  DEBUG << "bottomborder: " << bottomborder;
	  DEBUG.print();

	  datfile.seekg(280,ios.beg); // Record data
	  datfile.read((char)&c8, sizei8);
		c8.charAt(8)='\0';
		final uint numbytesdata = Integer.parseInt(c8);
	  DEBUG << "numbytesdata: " << numbytesdata;
	  DEBUG.print();


	// ====== Check with volumefile / internal ======
	// It seems that the lines (N+1) get from volume file is wrong,
	  // it seems to be the pixle number in one line 
	  if (numlines != checklines)
		{
		WARNING << "code 902: data file: " << crop_arg.filein1 << " numlin=" << numlines << " vs. volume file: " << crop_arg.filein1 << " numlin=" << checklines;
		WARNING.print();
		WARNING.print(" +this means data and volume file seem not to correspond.");
		}

	// ______ Check with previous section ______
	  if (numlines != numdatarec)
		{
		WARNING << "code 904: Number of lines seems not to be consistent in file: " << crop_arg.filein1 << " : " << numlines << " != " << numdatarec;
		WARNING.print();
		WARNING.print(" +this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
		}
	  if ((numbytesdata / 8) != numpixels)
		{
		WARNING << "code 904: Number of pixels seems to be inconsistent in file: " << crop_arg.filein1 << ": " << numpixels << " != " << (numbytesdata / 8);
		WARNING.print();
		WARNING.print(" +this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
		}


	// ====== Start copy input to output (raw) format with buffer======
	// ______ Check and process optional offset parameters______
	// ______ Lcnlow is corner line, lcnhi is other corner, pcnlow, pixel coord. low etc.
	  uint linestart = 1; // counters for loops
	  uint lineend = numlines;
	  uint pixelstart = 1;
	  uint pixelend = numpixels; // only for resultfile

	  if (crop_arg.dbow.linehi!=0 && crop_arg.dbow.linelo!=0 && crop_arg.dbow.pixhi!=0 && crop_arg.dbow.pixlo!=0)
		{
		window tempdbow = new window(crop_arg.dbow.linelo, crop_arg.dbow.linehi, crop_arg.dbow.pixlo, crop_arg.dbow.pixhi);
		if (crop_arg.dbow.linehi>numlines)
		  {
		  WARNING << "Specified input DBOW linehi > numlines: " << crop_arg.dbow.linehi << " > " << numlines << ". I set linehi = " << numlines;
		  WARNING.print();
		  tempdbow.linehi=numlines;
		  }
		if (crop_arg.dbow.pixhi>numpixels)
		  {
		  WARNING << "Specified input DBOW pixhi > numpixels: " << crop_arg.dbow.pixhi << " > " << numpixels << ". I set pixhi = " << numpixels;
		  WARNING.print();
		  tempdbow.pixhi=numpixels;
		  }
	// ______ Only hi values are possibly adapted, low is a constant ______
		numlines = tempdbow.linehi - crop_arg.dbow.linelo + 1;
		numpixels = tempdbow.pixhi - crop_arg.dbow.pixlo + 1;

		linestart = crop_arg.dbow.linelo;
		lineend = tempdbow.linehi;
		pixelstart = crop_arg.dbow.pixlo;
		pixelend = tempdbow.pixhi; // only for resultfile
		}

	// ______ Note complex<short> not in ANSI c ______
	  // 
		  matrix LINE = new matrix(1, 2 *numpixels); // size of real4

	// ====== Process requested lines ======
	  ofstream datoutfile;
	  RefObject<ofstream> TempRefObject2 = new RefObject<ofstream>(datoutfile);
	  openfstream(TempRefObject2, crop_arg.fileout1, generalinput.overwrit);
	  datoutfile = TempRefObject2.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(datoutfile, crop_arg.fileout1, __FILE__, __LINE__);

	  // ______ info on data, to avoid X86 problems ______
	  // ______ according to CEOS specs, byte 413 is first complex pixel, etc. ______
	  // 720 + 84124*linenum + 412 
	  datfile.seekg(lenrec1 + 412,ios.beg);

	  matrix TMPREAL4 = new matrix(1, 2);
	  datfile >> TMPREAL4; // read in first complex pixel for test

	  real8 tmpmag = Math.sqrt(real8(real4(ntohl(TMPREAL4(0,0)))*real4(ntohl(TMPREAL4(0,0)))) + real8(real4(ntohl(TMPREAL4(0,1)))*real4(ntohl(TMPREAL4(0,1)))));

	  DEBUG << "First complex element in datafile: (" << real4(ntohl(TMPREAL4(0,0))) << "," << real4(ntohl(TMPREAL4(0,1))) << "); mag = " << tmpmag;
	  DEBUG.print();
	  if (tmpmag > 10000.)
		{
		WARNING.print(DEBUG.get_str());
		WARNING.print("this is a byteorder problem on X86? (use ntohs)");
		}
	  DEBUG << "TEST: (realpart): " << TMPREAL4(0,0) << ", (imagpart): " << TMPREAL4(0,1);
	  DEBUG.print();
	  DEBUG << "TEST: htons(realpart): " << ntohl(TMPREAL4(0,0)) << ", htons(imagpart): " << ntohl(TMPREAL4(0,1));
	  DEBUG.print();
	  DEBUG << "TEST: ntohs(realpart): " << ntohl(TMPREAL4(0,0)) << ", ntohs(imagpart): " << ntohl(TMPREAL4(0,1));
	  DEBUG.print();
	  DEBUG << "TEST: short int(ntohs(realpart)): " << real4(ntohl(TMPREAL4(0,0))) << ", (imagpart): " << real4(ntohl(TMPREAL4(0,1)));
	  DEBUG.print();


	  // ====== perline is faster than perbuffer, less memory etc. BK1998 ======

	  datfile.seekg(lenrec1+(linestart-1)*lendatarec2 + 8,ios.beg);
	  datfile.read((char) lenrec2, sizeb4); // length of first record
	  lenrec2 = ntohl(lenrec2); // bk 6 jul 2000, byteorder x86 machines.

	  if (lenrec2 != lendatarec2)
		{
		ERROR << "code 904: Length of datarecords seems to be inconsistent in file: " << crop_arg.filein1 << ": " << lenrec2 << " != " << lendatarec2;
		WARNING.print(ERROR.get_str());
		ERROR.reset();
		}

	  final int32 TEN = 10;
	  final int32 TENPERCENT = (.5 *TEN+numlines)/TEN; // number of lines
	  int32 percentage = 0; // initialization
	  final int32 tmpstart = lenrec1+412-lendatarec2+(pixelstart-1)*8; // sizeof=8
	  char pD;
	  String pc;
	  for (register int32 linecnt =linestart; linecnt<=lineend; linecnt++)
		{
		if (!((linecnt-linestart)%TENPERCENT))
		  {
		  PROGRESS << "WRITESLC: " << setw(3) << percentage << "%";
		  PROGRESS.print();
		  percentage += TEN;
		  }
		datfile.seekg(tmpstart+linecnt *lendatarec2,ios.beg);
		datfile >> LINE;
		// ______ LG 28 DEC 2005: swapbytes for X86 (intel) linux cpus ______
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//	#if __X86PROCESSOR__
		for (int ii =0; ii<LINE.pixels(); ++ii)
			{
					pc = (char) LINE(0,ii);
					pD = pc;
					pc = *(pc+3);
					*(pc+3) = pD;
					pD = *(pc+1);
					*(pc+1) = *(pc+2);
					*(pc+2) = pD;
		//  LINE(0,ii) = LINE(0,ii)));      // changed from htons 171100 BK
			}
	//	#endif
		datoutfile << LINE;
		}
	  datfile.close(); // close files
	  datoutfile.close();



	// ====== Write results to scratchfile ======
	  ofstream scratchresfile = new ofstream("scratchres2raw", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile, "writeslc: scratchres2raw", __FILE__, __LINE__);
	  scratchresfile << "\n\n*******************************************************************\n";
		//<< "\n*_Start_crop:\t\t\t"
		//<<  crop_arg.idcrop
	  if (crop_arg.fileid == MASTERID)
		scratchresfile << "*_Start_" << processcontrol[pr_m_crop];
	  if (crop_arg.fileid == SLAVEID)
		scratchresfile << "*_Start_" << processcontrol[pr_s_crop];
	// ______ updateslcimage greps these ______
	  scratchresfile << "\t\t\t" << crop_arg.idcrop << "\n*******************************************************************" << "\nData_output_file: \t\t\t\t" << crop_arg.fileout1 << "\nData_output_format: \t\t\t\t" << "complex_real4" << "\nFirst_line (w.r.t. original_image): \t\t" << linestart << "\nLast_line (w.r.t. original_image): \t\t" << lineend << "\nFirst_pixel (w.r.t. original_image): \t\t" << pixelstart << "\nLast_pixel (w.r.t. original_image): \t\t" << pixelend << "\n*******************************************************************";
	  if (crop_arg.fileid == MASTERID)
		scratchresfile << "\n* End_" << processcontrol[pr_m_crop] << "_NORMAL";
	  if (crop_arg.fileid == SLAVEID)
		scratchresfile << "\n* End_" << processcontrol[pr_s_crop] << "_NORMAL";
	  scratchresfile << "\n*******************************************************************" << "\n";
	  scratchresfile.close();

	// ______ Tidy up do checks here ______
	  if (numchannels != 1) // ??
		{
		WARNING << "code 904: Number of channels in file: " << crop_arg.filein1 << " = " << numchannels << " != 1 ";
		WARNING.print();
		WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
		}
	  if (bottomborder != topborder != leftborder != rightborder != 0)
		{
		WARNING << "code 904: Not implemented: offset border: left,right,bottom,top: " << leftborder << "," << rightborder << "," << bottomborder << "," << topborder << " in file: " << crop_arg.filein1;
		WARNING.print();
		WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
		}
	  PROGRESS.print("WRITESLC: 100%");
	} // end palsar_fine_dump_data
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
// * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/readdata.cc,v $
// * $Revision: 3.37 $
// * $Date: 2005/10/06 11:09:20 $
// * $Author: kampes $
// *
// * routines for initial work, reading SLC,
// * writing to internal format.
// * see e.g: http://earth.esa.int:81/sarslc
// ***************************************************************



//#ifdef WIN32
//  #include "winsock2.h"                 // Jia changed this
//#else
//#endif


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
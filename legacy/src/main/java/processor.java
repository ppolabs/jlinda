public class GlobalMembersProcessor
{
	//
	// @file   processor.cc Doris InSAR processor.
	// @brief  main routine calling all modules.
	//
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
	// * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/processor.cc,v $  *
	// * $Revision: 4.02 $                                            *
	// * $Date: 2009/05/19 15:17:21 $                                 *
	// * $Author: TUDelft $                                           *
	// *                                                              *
	// * - main()                                                     *
	// * - initmessages()                                             *
	// * - handleinput()                                              *
	// * - usage()                                                    *
	// * - copyright()                                                *
	// * - quote()                                                    *
	// * - preview(...)                                               *
	// ***************************************************************


	// ====== Include files used here =========================================





	// ====== Global variables, declared extern in constants.h ======
	// === Declared extern in constants.h, set in main to be able to use them ===
	// === in readinput.cc they are set via the input file ======================
	public static bk_messages TRACE;
	public static bk_messages DEBUG;
	public static bk_messages INFO;
	public static bk_messages PROGRESS;
	public static bk_messages WARNING;
	public static bk_messages ERROR;
	public static bk_messages matDEBUG;
	public static bk_messages matERROR;

	// ______ used in matrix_bk.cc to keep track of total allocated ______
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if __DEBUGMAT2
	public static uint totalallocated =0; // [B] matrixdebuging void initmessages();

//***************************************************************
// @brief  initmessages initializes logging messages
// @param  none no arguments
// * initmessages                                                 *
// * set defaults for global messages objects                     *
// * #%// BK 12-Apr-2003                                          *
// ***************************************************************
	//#else
	public static void initmessages()
	  {
	  TRACE.setidentifyer("TRACE");
	  TRACE.dostderr(0);
	  TRACE.doprint(0);
	  TRACE.bellrings(0);
	  DEBUG.setidentifyer("DEBUG");
	  DEBUG.dostderr(0);
	  DEBUG.doprint(0);
	  DEBUG.bellrings(0);
	  INFO.setidentifyer("INFO");
	  INFO.dostderr(0);
	  INFO.doprint(1);
	  INFO.bellrings(0);
	  PROGRESS.setidentifyer("PROGRESS");
	  PROGRESS.dostderr(1);
	  PROGRESS.doprint(1);
	  PROGRESS.bellrings(1);
	  WARNING.setidentifyer("WARNING");
	  WARNING.dostderr(1);
	  WARNING.doprint(1);
	  WARNING.bellrings(2);
	  ERROR.setidentifyer("ERROR");
	  ERROR.dostderr(1);
	  ERROR.doprint(1);
	  ERROR.bellrings(3);
	  ERROR.doexit(0);
	  matDEBUG.setidentifyer("mtxDEBUG");
	  matDEBUG.dostderr(0);
	  matDEBUG.doprint(1); // default print, define prevents too much
	  matERROR.setidentifyer("mtxERROR");
	  matERROR.dostderr(1);
	  matERROR.doprint(1);
	  matERROR.bellrings(10); // should not occur anyway
	  matERROR.doexit(1);
	  } // END initmessages
	//#endif

//***************************************************************
// @brief  usage echos synopsis
// @param  programname name of this program
// * usage                                                        *
// * BK 01-Dec-1999                                               *
// ***************************************************************



	// ====== Prototypes functions in this file ======
	public static void usage(RefObject<String> programname)
	  {
	  TRACE_FUNCTION("usage");
	  String libs = checklibs();
	  cerr << "\n\t ***     **     **    *    **" << "\n\t *  *   *  *   *  *   *   *" << "\n\t *  *   *  *   *  *   *    *" << "\n\t *  *   *  *   * *    *     *" << "\n\t ***     **    *  *   *   **\n\n\n";
	  //     << "\n\t    Created by: Bert Kampes.\n\n"
		   //<< "\n\t\t\t(default: \"inputoptionsfile\")\n"
	  cerr << "\n  Program: \"" << programname.argvalue << "\" " << SWVERSION << "\n\tInterferometric processor for SAR SLC data.\n" << "\n\t(c) Delft University of Technology, the Netherlands." << "\n\n" << "  SYNOPSIS:\n\t" << programname.argvalue << " infile | -h [searchterm] | -v | -c | -q\n\n" << "\t  infile:       input file for " << programname.argvalue << "\n" << "\t  -h [term]:    call \"helpdoris\" (script with searchable help)\n" << "\t  -c:           return copyright notice.\n" << "\t  -q:           return random quote (not so random).\n" << "\t  -v:           return version number.\n\n\n" << "  LIBRARIES (used): " << libs << "\n" << "\n";
	  throw(usage_error); // normal exit
	  } // END usage

//***************************************************************
// @brief  handleinput parses input
// @param  argc original command line input to program
// @param  argv original command line input to program
// @param  input_gen general input struct
// * handleinput                                                  *
// * returns inputfilename as input_gen.logfile                   *
// * BK 01-Dec-1999                                               *
// ***************************************************************
	public static void handleinput(int argc, RefObject<String[]> argv, RefObject<input_gen> input_general)
	  {
	  TRACE_FUNCTION("handleinput");
	  // ====== Handle command line arguments ======
	  switch (argc)
		{
		case 3: // help search term
		  ; // do nothing, assume argv[1] is -h ...
		  //--- fall through --- //
		case 2: // optionsfilename, -h, or -ver
		  if (!strcmp(argv.argvalue.charAt(1),"-ver") || !strcmp(argv.argvalue.charAt(1),"-v") || !strcmp(argv.argvalue.charAt(1),"-VER") || !strcmp(argv.argvalue.charAt(1),"-V")) // version number asked
			{
			cerr << "Software name:    " << SWNAME << "\nSoftware version: " << SWVERSION << "\n";
			exit(0);
			}
		  // ______ Help wanted ______
		  else if (!strcmp(argv.argvalue.charAt(1),"-h") || !strcmp(argv.argvalue.charAt(1),"-H") || !strcmp(argv.argvalue.charAt(1),"-help") || !strcmp(argv.argvalue.charAt(1),"-HELP"))
			{
			// g++ makes: helpdoris (null) if empty argv[2] 
			if (argc ==2)
			  {
			  INFO.print("Making system call:");
			  INFO.print("helpdoris");
			  system("helpdoris");
			  }
			else
			  {
			  INFO << "helpdoris " << argv.argvalue.charAt(2) << ends; // search pattern
			  String cmd = new String(new char[512]); // command string
			  cmd = INFO.get_str();
			  INFO.print("Making system call:");
			  INFO.print(cmd);
			  system(cmd);
			  }
			exit(0);
			}
		  // ______ Copyright notice ______
		  else if (!strcmp(argv.argvalue.charAt(1),"-c") || !strcmp(argv.argvalue.charAt(1),"-C"))
			{
			copyright();
			}
		  // ______ Random quote ______
		  else if (!strcmp(argv.argvalue.charAt(1),"-q") || !strcmp(argv.argvalue.charAt(1),"-Q"))
			{
			quote();
			exit(0);
			}
		  // ______ Assume input file ______
		  else
			input_general.argvalue.logfile = argv.argvalue.charAt(1); // store in logfile.
		  break; // --- stop falling, no more options --- //
		// default no name for inputfile
		//case 1:                                           // no arguments
		//  strcpy(input_general.logfile,"inputoptionsfile");       // default name
		//  break;
		default:
		  usage(argv.argvalue.charAt(0));
		} // switch handle arguments

	  // ______ Check command line arguments ______
	  if (existed(input_general.argvalue.logfile)==true)
		{
		INFO << "input file: \"" << input_general.argvalue.logfile << "\"";
		INFO.print();
		}
	  else
		{
		ERROR << argv.argvalue.charAt(0) << " input file: \"" << input_general.argvalue.logfile << "\" not found.";
		ERROR.print();
		throw(input_error); // check if exception is caught in try block
		//usage(argv[0]);
		}
	  } // END handleinput

//***************************************************************
// @brief  quote like fortune program
// * quote                                                        *
// * I copied some quotes from the internet etc.                  *
// * sorry not to have a source for them.                         *
// * BK 15-Nov-2000                                               *
// ***************************************************************
	public static void quote()
	  {
	  TRACE_FUNCTION("quote");
	  // ______ Initialize quotes (send me new ones) ______
	  final int32 largestquote = ONE27;
	  String[] QUOTES = { "tip: Check out \"http://enterprise.lr.tudelft.nl/doris/\"", "tip: Check out \"http://enterprise.lr.tudelft.nl/dig/\"", "tip: See the online manual at \"http://enterprise.lr.tudelft.nl/doris/\"", "tip: Type: \"helpdoris\" for a synopsis of all input cards.", "tip: The \"run\" script can be used to generate template input files.", "tip: Get ASAR/ERS[12] precise orbits via \"http://www.deos.tudelft.nl/ers/precorbs/\"", "tip: For RADARSAT, use polyfit orbit interpolation, not approximation.", "tip: cpxfiddle is a nifty little utility.", "tip: dateconv is a nifty little utility.", "GMT is a great package.", "The only thing we have to fear is fear itself.", "2b || !2b = 1", "Sent me better quotes. Quick.", "The whole is more than the sum of its parts.", "To Thales the primary question was not what do we know,\n\tbut how do we know it.", "For the things of this world cannot be made known without\n\ta knowledge of mathematics.", "Life is a school of probability.", "\"Obvious\" is the most dangerous word in mathematics.", "An expert is a man who has made all the mistakes, which can be made,\n\tin a very narrow field.", "Structures are the weapons of the mathematician.", "A witty statesman said, you might prove anything by figures.", "Men pass away, but their deeds abide.", "It is a good thing for an uneducated man to read books of quotations.", "Mathematics is written for mathematicians.", "Revolutions never occur in mathematics.", "It is easier to square the circle than to get round a mathematician.", "Cognito ergo sum. \"I think, therefore I am.\"", "It is not enough to have a good mind. The main thing is to use it well.", "From a drop of water a logician could predict an Atlantic or a Niagara.", "I don't believe in mathematics.", "Imagination is more important than knowledge.", "A Mathematician is a machine for turning coffee into theorems.", "Whenever you can, count.", "Mathematics is a language.", "One should always generalize.", "Statistics: the mathematical theory of ignorance.", "When we ask advice, we are usually looking for an accomplice.", "What we know is not much. What we do not know is immense.", "Medicine makes people ill, mathematics make them sad and theology makes them sinful.", "A great truth is a truth whose opposite is also a great truth.", "I feign no hypotheses.", "It is not certain that everything is uncertain.", "Though this be madness, yet there is method in it.", "I have no faith in political arithmetic.", "Fourier is a mathematical poem.", "We think in generalities, but we live in details.", "I think that God in creating man somewhat overestimated his ability.", "Add little to little and there will be a big pile.", "Computers in the future may weigh no more than 1.5 tons. (1949)", "There is no reason anyone would want a computer in their home. (1977)", "Heavier-than-air flying machines are impossible. (1895)", "Everything that can be invented has been invented. (1899)", "640K ought to be enough for anybody. (1981)", "Pentiums melt in your PC, not in your hand.", "Always remember you're unique, just like everyone else.", "Ever notice how fast Windows runs? Neither did I.", "Double your drive space - delete Windows!", "Circular definition: see Definition, circular.", "43.3% of statistics are meaningless.", "Very funny, Scotty. Now beam down my clothes.", "I've got a problem. I say everything twice", "Don't insult the alligator till after you cross the river.", "Black holes are where God divided by zero.", "% make fire\n\tMake: Don't know how to make fire. Stop.", "% why not?\n\tNo match.", "% gotta light?\n\tNo match.", "% !1984\n\t 1984: Event not found. # (on some systems)", "% How's my lovemaking?\n\t Unmatched '.", "% \"How would you rate his incompetence?\n\tUnmatched \".", "% [Where is Jimmy Hoffa?\n\tMissing ].", "% [Where is my brain?\n\tMissing ].", "% ^How did the sex change^ operation go?\n\tModifier failed.", "% ar x \"my love life\"\n\tar: my love life does not exist", "This time it will surely run.", "Bug? That's not a bug, that's a feature.", "It's redundant! It's redundant!", "cpxfiddle is a great program.", "The shortest path between two truths in the real domain\n\tpasses through the complex domain.", "You have a tendency to feel you are superior to most computers.", "The first 90% of a project takes 90% of the time.", "The last 10% of a project takes 90% of the time.", "Any given program, when running, is obsolete.", "Any given program costs more and takes longer.", "If a program is useful, it will have to be changed.", "If a program is useless, it will have to be documented.", "Any given program will expand to fill all available memory.", "The value of a program is porportional to the weight of its output.", "Program complexity grows until it exceeds the capability\n\tof the programmer who must maintain it.", "Make it possible for programmers to write programs in English\n\tand you will find that programmers cannot write in English.", "On a helmet mounted mirror used by US cyclists:\n\t\"REMEMBER, OBJECTS IN THE MIRROR ARE ACTUALLY BEHIND YOU.\"", "On a New Zealand insect spray\n\t\"THIS PRODUCT NOT TESTED ON ANIMALS.\"", "In some countries, on the bottom of Coke bottles:\n\t\"OPEN OTHER END.\"", "On a Sears hairdryer:\n\t\"DO NOT USE WHILE SLEEPING.\"", "On a bar of Dial soap:\n\t\"DIRECTIONS - USE LIKE REGULAR SOAP.\"", "On a Korean kitchen knife:\n\t\"WARNING KEEP OUT OF CHILDREN.\"", "On an American Airlines packet of nuts:\n\t\"INSTRUCTIONS - OPEN PACKET, EAT NUTS.\"", "On a child's superman costume:\n\t\"WEARING OF THIS GARMENT DOES NOT ENABLE YOU TO FLY.\"", "Looking at wrapped interferograms is like running in circles.", "Conversation equals conservation (proposition with thesis Ramon Hanssen).", "Unlikely children's book title:\n\t\"Curios George and the High-Voltage Fence\".", "Unlikely children's book title:\n\t\"Controlling the Playground: Respect through Fear\".", "Unlikely children's book title:\n\t\"Mr Fork and Ms Electrical Outlet Become Friends\".", "Unlikely children's book title:\n\t\"Strangers Have the Best Candy\".", "Unlikely children's book title:\n\t\"Daddy Drinks Because You Cry\".", "Stanley looked quite bored and somewhat detached, but then penguins often do.", "Trouble with Windows? Reboot. Trouble with Unix? Be root.", "The good thing about standards is that there are so many to choose from.", "You can always count on people to do the right thing,\n\tafter they have exhausted all the alternatives.", "Where there is matter, there is geometry.", "The simplest schoolboy is now familiar with facts for which Archimedes\n\twould have sacrificed his life.", "Get the fastest fourier transform in the west at http://www.fftw.org/", "See http://www.gnu.org/ for compiler updates, etc.", "You can only find truth with logic if you have already found truth without it.", "Everything should be made as simple as possible, but not simpler.", "Seek simplicity, and distrust it.", "Descartes commanded the future from his study more\n\tthan Napoleon from the throne.", "Say what you know, do what you must, come what may.", "To the devil with those who published before us.", "The words figure and fictitious both derive from\n\tthe same Latin root, fingere. Beware!", "The best material model of a cat is another, or preferably the same, cat.", "He who loves practice without theory is like the sailor who boards ship\n\twithout a rudder and compass.", "Nature uses as little as possible of anything.", "Mathematics is not yet ready for such problems.", "You can only find truth with logic if you have already found truth without it.", "Natural selection is a mechanism for generating an exceedingly high\n\tdegree of improbability.", "Similar figures exist only in pure geometry.", "Writing briefly takes far more time than writing at length.", "Like the crest of a peacock so is mathematics at the head of all knowledge.", "The real danger is not that computers will begin to think like men,\n\tbut that men will begin to think like computers.", "Why did the blond stare at the orange juice?\n\tit said concentrate.", "Hofstadter's Law: It always takes longer than you expect,\neven when you take into account Hofstadter's Law.", "It can be of no practical use to know that Pi is irrational, but if we\n\tcan know, it surely would be intolerable not to know.", "Life is a complex, it has real and imaginary parts.", "Beauty can be perceived but not explained.", "Programming is like sex: one mistake and you have to\n\tsupport it for the rest of your life. [Michael Sinz]", "I wish you were here and I was there [me].", "In mathematics you dont understand things, you just get used to them [Neumann].", "A man is incomplete until he is married. After that, he is finished [Zsa Zsa Gabor].", "Unlikely children's book title:\n\t\"The Kids' Guide to Hitchhiking\".", "Unlikely children's book title:\n\t\"Whining, Kicking and Crying to Get Your Way\".", "Unlikely children's book title:\n\t\"Dads New Wife Robert\".", "Unlikely children's book title:\n\t\"The Little Sissy Who Snitched\".", ""};

	  // ______ Get random seed and print quote +endl ______
	  final int32 NUMQUOTES = 141;
	  final int32 seed = time(null);
	  final int32 quotenumber = seed%NUMQUOTES;
	  cerr << "\n  ..." << QUOTES.charAt(quotenumber) << "\n\n";
	  } // END quote

//***************************************************************
// @brief  copyright echos copyriht info
// * copyright                                                    *
// * BK 15-Nov-2000                                               *
// ***************************************************************
	public static void copyright()
	  {
	  TRACE_FUNCTION("copyright");
	  cerr << "\n\t ***     **     **   *    **" << "\n\t *  *   *  *   *     *   *" << "\n\t *  *   *  *   *     *    *" << "\n\t *  *   *  *   *     *     *" << "\n\t ***     **    *     *   **\n" << "\n" << "\n" << "COPYRIGHT NOTICE:\n" << "-----------------\n" << " Doris program is free software; you can redistribute it and/or modify\n" << "it under the terms of the GNU General Public License as published by\n" << "the Free Software Foundation; either version 2 of the License, or\n" << "(at your option) any later version.\n" << "\n" << "Doris is distributed in the hope that it will be useful,\n" << "but WITHOUT ANY WARRANTY; without even the implied warranty of\n" << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n" << "GNU General Public License for more details.\n" << "\n" << "You should have received a copy of the GNU General Public License\n" << "along with this program; if not, write to the Free Software\n" << "Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA\n" << "--- (END COPYRIGHT NOTICE) ---\n" << "\n" << "Publications that contain results produced by the Doris software should \n" << "contain an acknowledgment. (For example: The interferometric processing \n" << "was performed using the freely available Doris software package developed \n" << "by the Delft Institute of Earth Observation and Space Systems (DEOS), Delft \n" << "University of Technology, or include a reference to: Bert Kampes and \n" << "Stefania Usai. \"Doris: The Delft Object-oriented Radar Interferometric \n" << "software.\" In: proceedings 2nd ITC ORS symposium, August 1999. (cdrom)).\n" << "\n[Written by Bert Kampes, (c)1998-2003].\n\n";
	//
	//       << "COPYRIGHT NOTICE:\n"
	//       << "Doris is a scientific-purpose software and cannot be commercialized, nor \n"
	//       << "can parts or products of it be commercialized. Parties interested in using \n"
	//       << "Doris or its products for any commercial purposes are requested to contact \n"
	//       << "Prof. Dr. Roland Klees of DEOS (r.klees@citg.tudelft.nl)\n"
	//       << "\n"
	//       << "Our version of the software is the only official one. Please do not \n"
	//       << "distribute the Doris software to third parties, instead refer to the Doris \n"
	//       << "home page. This in order to guarantee uniformity in the distribution of \n"
	//       << "updates and information.\n"
	//       << "\n"
	//       << "Delft University of Technology is not responsible for any damage caused by \n"
	//       << "errors in the software or the documentation.\n"
	//       << "\n"
	//       << "Users are very welcome to extend the capabilities of the Doris software by \n"
	//       << "implementing new algorithms or improving the existing ones. It is intended \n"
	//       << "that if new software is developed based on Doris, that this also is made \n"
	//       << "available for free to the other users (through us).\n"
	//       << "\n"
	//       << "We would appreciate if any addition or modification of the software would \n"
	//       << "be announced first to us, so that it can be included in the official \n"
	//       << "(next) version of the software.\n"
	//       << "\n"
	//       << "Publications that contain results produced by the Doris software should \n"
	//       << "contain an acknowledgment. (For example: The interferometric processing \n"
	//       << "was performed using the freely available Doris software package developed \n"
	//       << "by the Delft Institute of Earth Observation and Space Systems (DEOS), Delft \n"
	//       << "University of Technology, or include a reference to: Bert Kampes and \n"
	//       << "Stefania Usai. \"Doris: The Delft Object-oriented Radar Interferometric \n"
	//       << "software.\" In: proceedings ITC 2nd ORS symposium, August 1999. (cdrom)).\n"
	//       << "\n[Written by Bert Kampes, (c)1998-2003].\n\n";
	//
	  exit(0);
	  } // END copyright

//***************************************************************
// @brief  preview generates cpxfiddle scripts for quicklooks
// @param  previewflag etc.
// * preview                                                      *
// * create scritps to generate sunraster files.                  *
// * previewflag set by PREVIEW card:                             *
// * i.e.: 0: OFF; 1:ON (generate sh script); 2:xv (call shscript)*
// * BK 07-Feb-2002                                               *
// * MA 24 Dec 2008 updates due to g++ 4.3                        *
// ***************************************************************
	//void preview(int32 previewflag, int32 width, int32 informat, const char *infile, const char *outfile, const char *opts);
	public static void preview(int32 previewflag, int32 width, int32 formatflag, String infile, String outfile, String opts)
	  {
	  TRACE_FUNCTION("preview (BK 07-Feb-2002)");
	  DEBUG << "previewflag: " << previewflag;
	  DEBUG.print();
	  if (previewflag ==0)
		{
		DEBUG.print("NO preview requested.");
		return;
		}
	  String scriptname = new String(new char[127]);
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if WIN32 || WINDOWS
	  scriptname = outfile.c_str();
	  scriptname += ".bat";
	//  #else
	  scriptname = "./"; // prevent error if "." is not in path (unix/linux/cygwin slash)
	  //strcat(scriptname,outfile);
	  scriptname += outfile.c_str();
	  scriptname += ".sh";
	//  #endif
	  INFO.reset(); // be sure to get rid of old chars before system call
	  if (previewflag<3)
		{
		PROGRESS.print("Start PREVIEW generation.");
		// ______ find out format of inputfile ______
		switch (formatflag)
		  {
		  case FORMATCR4:
			{
			INFO << "cpxfiddle -w " << width << " " << opts;
			INFO << " -f cr4 -l1 -p1 -P" << width; // -L?
			break;
			}
		  case FORMATCI2:
			{
			INFO << "cpxfiddle -w " << width << " " << opts;
			INFO << " -f ci2 -l1 -p1 -P" << width; // -L?
			break;
			}
		  // ___ trick cpxfiddle, make it think this complex ___
		  // ___ width is in complex pixels, but hgt is band interleaved ___
		  // ___ read second band (phase) as part of complex file___
		  // ___ -q normal must be set to trick cpxfiddle ___
		  case FORMATHGT:
			{
			PROGRESS.print("Trying to trick cpxfiddle for HGT format (should work)");
			int32 start_band_2 = Math.floor(width/2.0)+2.0; // (start at 1)
			DEBUG << "start of second band for cpxfiddle: " << start_band_2;
			DEBUG.print();
			INFO << "cpxfiddle -w " << width << " " << opts;
			INFO << " -f cr4 -l1 -p" << start_band_2 << " -P" << width;
			break;
			}
		  // ___ trick cpxfiddle, make it think this complex ___
		  // ___ this may not always work (in case of odd width)___
		  // ___ -q normal must be set to trick cpxfiddle ___
		  case FORMATR4:
			{
			// [added real4 option to cpxfiddle BK 13-Apr-2003]
			PROGRESS.print("Trying to trick cpxfiddle for REAL4 format (should work)");
			INFO << "cpxfiddle -w " << width << " " << opts;
			INFO << " -f r4 -l1 -p1 -P" << width;
			break;
			}
		  default:
			PRINT_ERROR("PREVIEW generation: Unknown formatflag.")
			throw(unhandled_case_error); // exit
		  }
		INFO << " " << infile << " > " << outfile << ends;
		INFO.print("With following command the SUNraster image can be generated again.");
		INFO.print(INFO.get_str()); // buffer should not be touched by previous print
		// ______ Create the script ______ //
		ofstream scriptfile;
		scriptfile.open(scriptname, ios.out);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		bk_assert(scriptfile, scriptname, __FILE__, __LINE__);
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if WIN32
		scriptfile << "REM !/bin/sh\n" << "REM  Script generated by Doris software to preview binary file.\n" << "REM  The program cpxfiddle is called to create a SUNraster file.\n" << "REM  For unwrapped interferograms (hgt format) cpxfiddle is tricked.\n" << "REM  Send comments to: doris-users <doris_users@tudelft.nl>\n" << "REM  to run this script type at prompt:\n" << "REM  " << scriptname << "\n" << INFO.get_str() << "\n";
		if (previewflag ==2)
		  {
		  scriptfile << "cximage " << outfile << "\n"; // Jia uses this windows viewer
		  }
		scriptfile << "REM EOF\n";
		scriptfile.close();
	//#else
		scriptfile << "#!/bin/sh\n" << "# Script generated by Doris software to preview binary file.\n" << "# The program cpxfiddle is called to create a SUNraster file.\n" << "# For unwrapped interferograms (hgt format) cpxfiddle is tricked.\n" << "# Send comments to: doris-users <doris_users@tudelft.nl>\n" << "# to run, type at prompt:\n" << "# " << scriptname << "\n" << INFO.get_str() << "\n";
		if (previewflag ==2)
		  {
		  scriptfile << "xv " << outfile << "\n";
		  }
		scriptfile << "# EOF\n";
		scriptfile.close();
		//
		// ______ Make file executable via system call to chmod ______ //
		DEBUG << "chmod +x " << scriptname << "\n" << ends;
		DEBUG.print(); // rewinds buffer (for next call)
		system(DEBUG.get_str());
	//#endif
		// ______ Actual make the call ______ //
		INFO.reset(); // make sure buffer is empty before system call
		INFO << scriptname << "&" << "\n" << ends;
		system(INFO.get_str());
		INFO.print(); // resets buffer (for next call)
		PROGRESS << "SUNraster file created of: " << infile << " (see also file: " << scriptname << ")";
		PROGRESS.print();
		}
	  else
		{
		PRINT_ERROR("PREVIEW generation: Unknown option.")
		throw(unhandled_case_error); // exit
		}
	  } // END preview






	//***************************************************************
	// @brief  main main program
	// @param  argc argument count
	// @param  argv value of arguments
	// *    main                                                      *
	// * A big switch to select the modules the user requests.        *
	// * Handles some general information such as size of             *
	// * current file as well.                                        *
	// *    Bert Kampes, 11-Dec-1998                                  *
	// ***************************************************************
	static int Main(int argc, RefObject<String[]> argv)
	  {
	  // ______ Catch math errors, zero division, etc. ______
	  // ______ v3.15: seems g++ 4 generates SIGSEGV errors due to over-optimization? ___
	  CatchSignals(handle_signal); // floating point errors, segmentation faults, etc.
	  // try should work with current compilers
	  try
	  {
	  // ====== Set defaults for message class so we can use it here ======
	  initmessages(); // set DEBUG etc.
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if ! NO_FASTTRIG
	  init_fasttrig(); // lookup table for fast_sin() ... etc. char ident[] = "@(#)Doris InSAR software, $Revision: 4.02.6 $, $Author: TUDelft $";
	//  #else
	  String ident = "@(#)Doris InSAR software, $Revision: 4.02.6 $, $Author: TUDelft $";

	//  #endif
	  // ====== ident string for `what doris` ======
	  cerr << "\n";
	  INFO.print(ident); //use ident so it is not optimized away

	  // ====== Structs for storing input ======
	  input_gen input_general; // general inputparameters
													// misuse logfile for input file name
	  input_ell input_ellips; // ellipsoid pm.
	  input_pr_orbits input_porbits; // inputparam. readfiles slave
	  input_readfiles input_m_readfiles; // inputparam. readfiles master
	  input_crop input_m_crop; // inputparam. crop master
	  input_oversample input_m_oversample; // by Raffaele Nutricato
	  input_simamp input_m_simamp; // inputparam. simulated amplitude [FvL, MA]
	  input_mtiming input_m_mtiming; // inputparam. simu. amplitude correl. + master timing error e.i absolute timing [MA]

	  input_readfiles input_s_readfiles; // inputparam. readfiles slave
	  input_crop input_s_crop; // inputparam. crop slave
	  input_oversample input_s_oversample; // by Raffaele Nutricato
	  input_filtazi input_ms_filtazi; // inputparam. azimuth filtering
	  input_resample input_s_resample; // inputparam. resample

	  input_coarsecorr input_i_coarsecorr; // inputparam. coarse correlation
	  input_fine input_i_fine; // inputparam. fine correlation
	  input_reltiming input_i_timing; // inputparam. relative timing error [FvL]
	  input_demassist input_i_demassist; // inputparam. dem assisted coregistraton [FvL]
	  input_coregpm input_i_coregpm; // inputparam. coreg pm comp.
	  input_filtrange input_ms_filtrange; // inputparam. range filtering
	  input_interfero input_i_interfero; // inputparam. interfero.
	  input_coherence input_i_coherence; // inputparam. coherence
	  input_comprefpha input_i_comprefpha; // inputparam. flat earth corr.
	  input_subtrrefpha input_i_subtrrefpha; // inputparam. subtr earth corr.
	  input_comprefdem input_i_comprefdem; // inputparam. ref. dem
	  input_subtrrefdem input_i_subtrrefdem; // inputparam. subtr ref. dem
	  input_filtphase input_i_filtphase; // inputparam. phase filter
	  input_dinsar input_i_dinsar; // inputparam. 3pass diff. insar
	  input_unwrap input_i_unwrap; // inputparam. unwrapping
	  input_slant2h input_i_slant2h; // inputparam. slant2height
	  input_geocode input_i_geocode; // inputparam. geocode


	// ______ Info on images/products ______
	  slcimage master; // info on master image
	  productinfo simamp; // info on simulated master amplitude
	  slcimage slave; // info on slave image
	  productinfo interferogram; // interferogram
	  productinfo unwrappedinterf; // interferogram unwrapped
	  //productinfo         differentialinterf;     // diff. interferogram
	  productinfo radarcodedrefdem; // reference DEM
	  productinfo coherence; // coherence image
	  productinfo heightinradarsystem; // s2h, in radarcoord.


	// ______ Orbital data (interpolation) ______
	  orbit masterorbit; // class orbit
	  orbit slaveorbit; // class orbit


	// ______ Matrices for storing polynomial coefficients ______
	  matrix<real8> coeff_cpmL; // coregistration parameters
	  matrix<real8> coeff_cpmP; // coregistration parameters
	  matrix<real8> coeff_flat; // flatearth correction
	  matrix<real8> coeff_h2ph; // h2ph factors, added by FvL



	// ====== Initialize options ("inputoptionsfile"), files etc ======
	  RefObject<input_gen> TempRefObject = new RefObject<input_gen>(input_general);
	  handleinput(argc, argv, TempRefObject); // handle input
	  input_general = TempRefObject.argvalue;
	  printcpu(1); // initialize timers
	  readinput(input_general, input_ellips, input_porbits, input_m_readfiles, input_m_crop, input_m_oversample, input_m_simamp, input_m_mtiming, input_s_readfiles, input_s_crop, input_s_oversample, input_ms_filtazi, input_i_coarsecorr, input_i_fine, input_i_timing, input_i_demassist, input_i_coregpm, input_s_resample, input_ms_filtrange, input_i_interfero, input_i_coherence, input_i_comprefpha, input_i_subtrrefpha, input_i_comprefdem, input_i_subtrrefdem, input_i_filtphase, input_i_dinsar, input_i_unwrap, input_i_slant2h, input_i_geocode);
	  // ___ Perform some initial tests of platform after SCREEN is set ___
	  inittest(); // initial tests
	  // ___ transfer general input to orbit ___
	  masterorbit.set_interp_method(input_general.orb_interp); // transfer input
	  slaveorbit.set_interp_method(input_general.orb_interp); // transfer input
	  // ___ trace function after readinput (screen level known) ___
	  TRACE_FUNCTION("main");



	// ====== Write initial info if appropriate ======
	  initwrite(input_general.logfile, LOGID); // general info logfile(=1)
	  RefObject<input_gen> TempRefObject2 = new RefObject<input_gen>(input_general);
	  boolean processmaster = doinitwrite(TempRefObject2, MASTERID);
	  input_general = TempRefObject2.argvalue;
	  RefObject<input_gen> TempRefObject3 = new RefObject<input_gen>(input_general);
	  boolean processlave = doinitwrite(TempRefObject3, SLAVEID);
	  input_general = TempRefObject3.argvalue;
	  RefObject<input_gen> TempRefObject4 = new RefObject<input_gen>(input_general);
	  boolean processinterf = doinitwrite(TempRefObject4, INTERFID);
	  input_general = TempRefObject4.argvalue;
	  if (processmaster)
		initwrite(input_general.m_resfile, MASTERID); // general info master result=2
	  if (processlave)
		initwrite(input_general.s_resfile, SLAVEID); // general info slave result=3
	  if (processinterf)
		initwrite(input_general.i_resfile, INTERFID); // general info interf_out=4
	  updatefile("scratchreadinput", input_general.logfile); // copy input2log (readinput)



	// ====== Check requested processing and fill alreadyprocessed ======
	  int16[] alreadyprocessed = new int16[NUMPROCESSES];
	  checkprocessing(input_general, alreadyprocessed); // initialize alreadyprocessed

	  int16 status = 0; // [MA] check exit status of system calls for proper error handling

	  DEBUG.print("Time spent for initialization:");
	  printcpu();
	  PROGRESS.print("Finished initialization");




	// *********************************************************
	// **** start master                                    ****
	// *********************************************************

	// ====== PROCESS INFOFILES CDROM: MASTER ======
	  if (input_general.process[pr_m_readfiles])
		{
		PROGRESS.print("Start M_READFILES.");
		alreadyprocessed[pr_m_readfiles] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing readfiles for master.";
		  getanswer();
		  }

		// ______ Select method for READFILES______
		switch (input_m_readfiles.sensor_id)
		  {
		  case SLC_ALOS:
			INFO.print("ALOS: Under development: is ALOS really CEOS/Atlantis like ERS?");
		  // ______ JERS also ceos? (falls through to ers) ______
		  case SLC_JERS:
			INFO.print("JERS: Under development: is JERS really CEOS like ERS?");
		  // ______ RSAT also CEOS? (falls through to ERS, account in reader) ______
		  case SLC_RSAT:
			WARNING.print("RSAT: for orbit highest polyfit recommended.");
			WARNING.print("RSAT: Under development: CEOS reader seems to work.");
		  // ______ ERS-1/2 ______
		  case SLC_ERS:
			// ______ Get checks slave volumefile ______
			String c16checkvol1 = new String(new char[17]);
			String c16checkvol2 = new String(new char[17]);
			String c16checkvol3 = new String(new char[17]); // check id of volume file
			if (alreadyprocessed[pr_s_readfiles])
			  {
			  readres(c16checkvol1,sizeof(c16checkvol1), input_general.s_resfile,"Volume_ID:");
			  readres(c16checkvol2,sizeof(c16checkvol2), input_general.s_resfile,"Volume_identifier:");
			  readres(c16checkvol3,sizeof(c16checkvol3), input_general.s_resfile,"Volume_set_identifier:");
			  }
			// ______ Extract info ______
			readvolume(input_m_readfiles, c16checkvol1, c16checkvol2, c16checkvol3);
			// ___ written in readvolume ___
			updatefile("scratchresvol", input_general.m_resfile);
			// ___ written in readvolume ___
			updatefile("scratchlogvol", input_general.logfile);
			// ______Get check with master volumefile______
			String c8checkleadat = new String(new char[9]); // check vol met lea en dat
			readres(c8checkleadat,sizeof(c8checkleadat), input_general.m_resfile,"(Check)Number",5);
			// ______ Extract info ______
			// read info
			readleader(input_m_readfiles, Integer.parseInt(c8checkleadat)-1);
			// made in readleader
			updatefile("scratchreslea", input_general.m_resfile);
			// made in readleader
			updatefile("scratchloglea", input_general.logfile);
			// ______ DONT PROCESS NULLFILE ______
			//    readnull(input_m_readfiles);
			//    updatefile("scratchresnull",input_general.m_resfile);
			//    updatefile("scratchlognull",input_general.logfile);
			//
			// read data file info
			readdat(input_m_readfiles, Integer.parseInt(c8checkleadat)-1);
			// update resfile
			updatefile("scratchresdat", input_general.m_resfile);
			// update logfile
			updatefile("scratchlogdat", input_general.logfile);
			break;
		  // ______ ENVISAT ASAR ______
		  case SLC_ASAR:
			INFO.reset(); // make sure nothing in buffer
			INFO << "envisat_dump_header2doris.csh " << input_m_readfiles.datfile << " > scratchres_envisat" << "\n" << ends;
			String cmd = new String(new char[512]); // command string
			cmd = INFO.get_str();
			INFO.print("With following command the envisat header was read.");
			INFO.print(cmd);
			PROGRESS.print("Making system call to envisat_dump_header2doris.csh");
			PROGRESS.print("(also requires envisat_dump_header program)");
			system(cmd); // this does the work
			INFO.reset();
			PROGRESS.print("Finished system call to envisat_dump_header2doris.csh");
			// ___ update resfile ___
			updatefile("scratchres_envisat", input_general.m_resfile);
			// ___ update logfile ___
			updatefile("envisat_dump_header.log", input_general.logfile);
			break;
		// ______ TSX COSAR ______
		case SLC_TSX:
		  INFO.reset(); // make sure nothing in buffer
		  INFO << "tsx_dump_header2doris.py " << input_m_readfiles.leaderfile << " > scratchres_tsx" << "\n" << ends;
		  //char cmd[512];// command string
		  cmd = INFO.get_str();
		  INFO.print("With following command the TSX header was read.");
		  INFO.print(cmd);
		  PROGRESS.print("Making system call to tsx_dump_header2doris.py");
		  PROGRESS.print("(also requires python, gdal and XML libraries)");
		  status =system(cmd); // this does the work
		  if (status != 0) // [MA] TODO make it a function
			{
			ERROR << "tsx_dump_header2doris.py: failed with exit code: " << status;
			PRINT_ERROR(ERROR.get_str())
			throw(some_error);
			}
		  INFO.reset();
		  PROGRESS.print("Finished system call to tsx_dump_header2doris.py");
		  // ___ update resfile ___
		  updatefile("scratchres_tsx", input_general.m_resfile);
		  // ___ update logfile ___
		  // updatefile("tsx_dump_header.log",input_general.logfile);
		  break;
		  default:
			PRINT_ERROR("M_READFILES step is impossible, check input and data for the sensor parameters.")
			throw(input_error); // exit
		  }
		PROGRESS.print("Finished M_READFILES.");
		DEBUG.print("Time spent for reading files master:");
		printcpu();
		}

	// ______Fill slcimage struct______
	  if (existed(input_general.m_resfile))
		{
		master.fillslcimage(input_general.m_resfile);
		interferogram.win = master.currentwindow;
		}

	  // ______ Account for timing errors ______
	  if (input_m_readfiles.az_timing_error !=0.0 || input_m_readfiles.rg_timing_error != 0.0)
		{
		 PROGRESS.print("The correction values for master image timing errors are specified.");
		 PROGRESS << "m_az_timing_error: " << input_m_readfiles.az_timing_error;
		 PROGRESS.print();
		 PROGRESS << "m_rg_timing_error: " << input_m_readfiles.rg_timing_error;
		 PROGRESS.print();
		 master.add_rg_t_error(input_m_readfiles.rg_timing_error);
		 master.add_az_t_error(input_m_readfiles.az_timing_error);
		 master.timingerror_flag = true; // (true) then time is updated manually [MA]
		}


	// ====== GET PRECISE ORBITS (GETORB): MASTER ======
	// BK 15-Dec-2003 (moved before crop since required if dbow_geo
	  if (input_general.process[pr_m_porbits])
		{
		PROGRESS.print("Start M_PORBITS.");
		alreadyprocessed[pr_m_porbits] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing get precise orbits for master.";
		  getanswer();
		  }
		// ______ Getorb does the real work ______
		getorb(master,input_porbits,MASTERID);
		convertgetorbout(MASTERID);
		// ______ result of conv.getorb master ______
		updatefile("scratchdatapoints", input_general.m_resfile);
		// ______ remove data section by leader ______
		removedatleader(input_general.m_resfile);
		PROGRESS.print("Finished M_PORBITS.");
		DEBUG.print("Time spent for obtaining precise orbits master:");
		printcpu();
		}
	// ______ Get coefficients for orbit interpolation if data in file ______
	  //(alreadyprocessed[pr_m_readfiles] && isspace(input_porbits.m_orbdir[0])))
	  if (alreadyprocessed[pr_m_porbits] || alreadyprocessed[pr_m_readfiles])
		masterorbit.initialize(input_general.m_resfile);

	  // ______ Dump interpolated orbit if requested ______
	  // ______ also possible to dump orbit of readfiles with card M_ORB_DUMP ___
	  // ___ this returns if no orbit available ___
	  masterorbit.dumporbit(input_porbits, MASTERID);



	// ====== PROCESS DATAFILE CDROM: MASTER ======
	  if (input_general.process[pr_m_crop])
		{
		PROGRESS.print("Start M_CROP.");
		alreadyprocessed[pr_m_crop] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing crop for master.";
		  getanswer();
		  }

		// ______ Assume orbits loaded, use GEO card if present ______
		if (input_m_crop.dbow_geo.pixhi != 0)
		  {
		  PROGRESS.print("Computing m_crop based on center lat/lon/height/width.");
		  real8 centerline;
		  real8 centerpix;
		  real8 centerphi = deg2rad(input_m_crop.dbow_geo.linelo/1e6-360.0);
		  real8 centerlambda = deg2rad(input_m_crop.dbow_geo.linehi/1e6-360.0);
		  real8 centerheight = input_general.terrain_height;
		  DEBUG.print("Using height of HEIGHT card");
		  cn centerpos; // compute phi/lambda/h-->x,y,z
		  centerpos = input_ellips.ell2xyz(centerphi,centerlambda,centerheight);
		  DEBUG << "Converted M_DBOW_GEO phi,lambda,hei --> x,y,z: " << centerphi << ", " << centerlambda << ", " << centerheight << " --> " << centerpos.x << " " << centerpos.y << " " << centerpos.z;
		  DEBUG.print();
		  DEBUG << "Center master SLC: x,y,z: " << master.approxcentreoriginal.x << " " << master.approxcentreoriginal.y << " " << master.approxcentreoriginal.z;
		  DEBUG.print();
		  if (Math.abs(master.approxcentreoriginal.x-centerpos.x)>60000.0) // 50km
			WARNING.print("M_DBOW_GEO: coordinates seem to be outside SLC area? (X)");
		  if (Math.abs(master.approxcentreoriginal.y-centerpos.y)>60000.0) // 50km
			WARNING.print("M_DBOW_GEO: coordinates seem to be outside SLC area? (Y)");
		  if (Math.abs(master.approxcentreoriginal.z-centerpos.z)>60000.0) // 50km
			WARNING.print("M_DBOW_GEO: coordinates seem to be outside SLC area? (Z)");
		  RefObject<real8> TempRefObject5 = new RefObject<real8>(centerline);
		  RefObject<real8> TempRefObject6 = new RefObject<real8>(centerpix);
		  RefObject<orbit> TempRefObject7 = new RefObject<orbit>(masterorbit);
		  xyz2lp(TempRefObject5, TempRefObject6, master, TempRefObject7, centerpos, 10, 1e-3); // MAXITERATIONS, CONVERGENCE_TIME
		  centerline = TempRefObject5.argvalue;
		  centerpix = TempRefObject6.argvalue;
		  masterorbit = TempRefObject7.argvalue;
		  DEBUG << "CROP: center line/pix: " << centerline << " " << centerpix;
		  DEBUG.print();
		  int32 l0 = int32(centerline+0.5) - input_m_crop.dbow_geo.pixlo/2;
		  int32 lN = l0 + input_m_crop.dbow_geo.pixlo - 1;
		  int32 p0 = int32(centerpix+0.5) - input_m_crop.dbow_geo.pixhi/2;
		  int32 pN = p0 + input_m_crop.dbow_geo.pixhi - 1;
		  if (l0 < int32(master.currentwindow.linelo))
			  l0 = master.currentwindow.linelo;
		  if (lN > int32(master.currentwindow.linehi))
			  lN = master.currentwindow.linehi;
		  if (p0 < int32(master.currentwindow.pixlo))
			  p0 = master.currentwindow.pixlo;
		  if (pN > int32(master.currentwindow.pixhi))
			  pN = master.currentwindow.pixhi;
		  INFO << "DBOW from GEO: " << l0 << " " << lN << " " << p0 << " " << pN;
		  INFO.print();
		  // ___ Simply fill dbow with it and it will be done correctly! ___
		  input_m_crop.dbow.linelo = uint(l0);
		  input_m_crop.dbow.linehi = uint(lN);
		  input_m_crop.dbow.pixlo = uint(p0);
		  input_m_crop.dbow.pixhi = uint(pN);
		  }

		// ______ For each sensor, of course the formats differ... ______
		// ______ make sure radarsat_dump_data is called for ATLANTIS ______
		switch (master.sensor)
		  {
			case SLC_ALOS:
					String c8checkleadata = new String(new char[9]); // check vol met lea en dat
				readres(c8checkleadata,sizeof(c8checkleadata), input_general.m_resfile,"(Check)Number",5);
				palsar_fine_dump_data(input_general,input_m_crop, Integer.parseInt(c8checkleadata)-1);
				updatefile("scratchres2raw", input_general.m_resfile); // update resfile
				break;

		 // ______ JERS also CEOS? (falls through to ers) ______
		  case SLC_JERS:
			WARNING.print("JERS: Under development: is JERS really CEOS like ERS?");
		  // ______ ERS-1/2 ______
		  case SLC_ERS:
			switch (master.sar_processor)
			  {
			  case SARPR_VMP:
				// --- default, ESA SLCs ---
				// ______Get check with volumefile______
				String c8checkleadat = new String(new char[9]); // check vol met lea en dat
				readres(c8checkleadat,sizeof(c8checkleadat), input_general.m_resfile,"(Check)Number",5);
				writeslc(input_general,input_m_crop, Integer.parseInt(c8checkleadat)-1);
				updatefile("scratchres2raw", input_general.m_resfile); // update resfile
				break;
			  case SARPR_ATL:
				// --- self processed from RAW ---
				radarsat_dump_data(input_general,input_m_crop);
				updatefile("scratchres2raw", input_general.m_resfile); // update resfile
				WARNING.print("ATLANTIS data cropping not fully checked yet.");
				break;
			  case SARPR_TUD:
				//break;
			  default:
				PRINT_ERROR("unrecognized sar_processor.")
				throw(input_error); // exit
			  }
			break;
		  // ______ ENVISAT ASAR ______
		  case SLC_ASAR:
			PROGRESS.print("System call to get ASAR SLC data (requires envisat_dump_data program)");
			// ______ Data window must be set correctly ______
			if (input_m_crop.dbow.linelo == 0 || input_m_crop.dbow.linehi == 0 || input_m_crop.dbow.pixlo == 0 || input_m_crop.dbow.pixhi == 0)
			  {
			  input_m_crop.dbow = master.currentwindow; // needs to be set correctly after readfiles
			  }
			envisat_dump_data(input_m_crop);
			PROGRESS.print("Finished system call to envisat_dump_data");
			updatefile("scratchres2raw", input_general.m_resfile); // update resfile
			break;
		  // ______ TSX COSAR ______
		case SLC_TSX:
		  PROGRESS.print("System call to get TSX SLC/COSAR data (requires tsx_dump_data program)");
		  // ______ Data window must be set correctly ______
		  if (input_m_crop.dbow.linelo == 0 || input_m_crop.dbow.linehi == 0 || input_m_crop.dbow.pixlo == 0 || input_m_crop.dbow.pixhi == 0)
			{
			  input_m_crop.dbow = master.currentwindow; // needs to be set correctly after readfiles
			}
		  tsx_dump_data(input_m_crop);
		  PROGRESS.print("Finished system call to tsx_dump_data");
		  updatefile("scratchres2raw", input_general.m_resfile); // update resfile
		  break;
		  // ______ RADARSAT (detected by Product specifier) ______
		  case SLC_RSAT:
			radarsat_dump_data(input_general,input_m_crop);
			updatefile("scratchres2raw", input_general.m_resfile); // update resfile
			WARNING.print("RSAT (Atlantis) not fully checked yet.");
			break;
		  default:
			PRINT_ERROR("M_CROP step is impossible, check input and data for the sensor parameters.")
			//PRINT_ERROR("impossible error (sensor not ERS JERS RSAT nor ASAR).")
			throw(input_error); // exit
		  }
		PROGRESS.print("Finished M_CROP.");
		DEBUG.print("Time spent for writing raw format master:");
		printcpu();
		}

	// ______ Update info on image ______
	  if (alreadyprocessed[pr_m_crop])
		{
		INFO.print("");
		INFO.print("master: latest known processing stage: crop");
		String SECTIONID = new String(new char[ONE27]);
		SECTIONID = "*_Start_";
		SECTIONID += processcontrol[pr_m_crop];
		master.updateslcimage(input_general.m_resfile,SECTIONID);
		}

	// ______ Generate magnitude preview of master SLC if requested ______
	  if (input_general.process[pr_m_crop])
		{
		PROGRESS.print("calling preview for cropped master");
		if (master.sensor!=SLC_RSAT && master.sensor!=SLC_ALOS && master.sensor!=SLC_TSX) // master, slave , resampled slave TODO: convert to preview function
		  preview(input_general.preview, master.currentwindow.pixels(), master.formatflag, master.file, "master_mag.ras", "-e 0.5 -s 1.0 -q mag -o sunraster -b -c gray -M 2/10"); // general
	//    else if (master.sensor == SLC_ALOS)
	//      preview(input_general.preview,
	//              master.currentwindow.pixels(), master.formatflag,
	//              master.file, "master_mag.ras", 
	//              "-e 0.5 -s 1.0 -q mag -o sunraster -b -c gray -M 2/2");   // ALOS VV/HH az/ra 5/5 m
		else if (master.sensor == SLC_ALOS || master.sensor == SLC_TSX)
		  preview(input_general.preview, master.currentwindow.pixels(), master.formatflag, master.file, "master_mag.ras", "-e 0.5 -s 1.0 -q mag -o sunraster -b -c gray -M 10/10"); // [MA] higher mlook to reduce filesize
		else //multilooking depends on beam
		  preview(input_general.preview, master.currentwindow.pixels(), master.formatflag, master.file, "master_mag.ras", "-e 0.5 -s 1.0 -q mag -o sunraster -b -c gray -M 5/6");
		}



	//____ RaffaeleNutricato START MODIFICATION SECTION 5 ______
	// ====== OVERSAMPLE SLC: MASTER ======
	  if (input_general.process[pr_m_oversample])
		{
		PROGRESS.print("Start M_OVS.");
		alreadyprocessed[pr_m_oversample] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing oversample for master.";
		  getanswer();
		  }
		// ______Master and Slave oversampling ratios should be identical______
		if (alreadyprocessed[pr_s_oversample])
		  {
		  String checkmultilook = new String(new char[20]);
		  int32 checkOsrRange;
		  readres(checkmultilook,sizeof(checkmultilook), input_general.s_resfile,"Multilookfactor_range_direction:");
		  checkOsrRange = int32(1.0/Double.parseDouble(checkmultilook)+0.5);
		  INFO << "m_ovs: " << input_m_oversample.OsrRange << "checkOsrRange: " << checkOsrRange;
		  INFO.print();
		  if (input_m_oversample.OsrRange != checkOsrRange)
			WARNING.print("Master and Slave range oversampling factors should be identical!!!");
		  }
		if (alreadyprocessed[pr_s_oversample])
		  {
		  String checkmultilook = new String(new char[20]);
		  int32 checkOsrAzimuth;
		  readres(checkmultilook,sizeof(checkmultilook), input_general.s_resfile,"Multilookfactor_azimuth_direction:");
		  checkOsrAzimuth = int32(1.0/Double.parseDouble(checkmultilook)+0.5);
		  if (input_m_oversample.OsrAzimuth != checkOsrAzimuth)
			WARNING.print("Master and Slave azimuth oversampling factors should be identical!!!");
		  }

		// ______ Oversample master image ______
		OversampleSLC(input_general,master,input_m_oversample,MASTERID);
		updatefile("scratchoversample", input_general.m_resfile); // update resfile

		PROGRESS.print("Finished M_OVS.");
		DEBUG.print("Time spent for oversampling master:");
		printcpu();
		}

	  // ______ Update info on image ______
	  if (alreadyprocessed[pr_m_oversample])
		{
		INFO.print("");
		INFO.print("master: latest known processing stage: oversample");
		String SECTIONID = new String(new char[ONE27]);
		SECTIONID = "*_Start_";
		SECTIONID += processcontrol[pr_m_oversample];
		master.updateslcimage(input_general.m_resfile,SECTIONID); // Bert Kampes, 29-Jul-2005
		interferogram.win = master.currentwindow; // larger size
		}

	  // ______ Generate magnitude preview of master SLC if requested ______
	  // -M 2/10 may make no sense after oversampling (depends).
	  if (input_general.process[pr_m_oversample])
		{
		PROGRESS.print("calling preview for oversampled master");
		preview(input_general.preview, master.currentwindow.pixels(), master.formatflag, master.file, "master_ovs_mag.ras", "-e 0.5 -s 1.0 -q mag -o sunraster -b -c gray -M 4/4");
		}


	// ____ start added by MA ____

	//***************************************************************
	// *              SIMAMP + MASTER TIMING ERROR                    *
	// *                                                              *
	// * Estimate timing error (azimuth, range) of master image       *
	// * based on simulated amplitude image (DEM based)               *
	// *                                                              *
	// * 1. Simulate Amplitude                                        *
	// * 2. Correlate MASTER with DEM (SIM. Ampli.) for timing        *
	// *                                                              *
	// ***************************************************************

	// ====== SIMULATE AMPLITUDE (MASTER) ======
	  if (input_general.process[pr_m_simamp])
		{
		PROGRESS.print("Start SIMULATE AMPLITUDE.");
		alreadyprocessed[pr_m_simamp] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nSimulating master amplitude using DEM.\n";
		  getanswer();
		  }

		if (master.timingerror_flag == true) // [MA] time is already updated, this will shift dem
		  {
		   WARNING << "SAM: M_AZ_T_ERROR = " << input_m_readfiles.az_timing_error << " sec.";
		   WARNING.print();
		   WARNING << "SAM: M_RG_T_ERROR = " << input_m_readfiles.rg_timing_error << " sec.";
		   WARNING.print();
		   WARNING << "SAM: timing error cards are already set, please check. This will effect the simulation.";
		   WARNING.print();
		  }

	//    if ( master.timingerror_flag == false || master.timingerror_flag == true )   // this has no effect, just want to keep it for the future [MA]
	//      {
			sim_amplitude(input_general, input_ellips, input_m_simamp, master, masterorbit);
	//      }

		// ______ Update log files ______
		updatefile("scratchlogsimamp", input_general.logfile);
		updatefile("scratchressimamp", input_general.m_resfile);

		PROGRESS.print("Finished SIMULATE AMPLITUDE.");
		DEBUG.print("Time spent for simulating master amplitude:");
		printcpu();
		}

	// ______ Fill input struct for  ______
	  if (alreadyprocessed[pr_m_simamp]) // update simamp image
		{
		INFO.print("");
		INFO.print("master: latest known processing stage: simamplitude");
		String SECTIONID = new String(new char[ONE27]);
		SECTIONID = "*_Start_";
		SECTIONID += processcontrol[pr_m_simamp];
		simamp.fillproductinfo(input_general.m_resfile,SECTIONID);
		}


	// ______ Generate magnitude preview of simulated amlitude if requested ______
	  if (input_general.process[pr_m_simamp])
		{
		PROGRESS.print("calling preview for simulated amplitude");
		if (master.sensor!=SLC_RSAT)
		  preview(input_general.preview, int32((simamp.win.pixels()/simamp.multilookP)), simamp.formatflag, simamp.file, "master_simamp.ras", "-e 1.0 -s 1.0 -q normal -o sunraster -b -c gray -M 2/10");
		else //multilooking depends on beam
		  preview(input_general.preview, int32((simamp.win.pixels()/simamp.multilookP)), simamp.formatflag, master.file, "master_simamp.ras", "-e 1.0 -s 1.0 -q normal -o sunraster -b -c gray -M 5/6");
		}

	// ====== AMPLITUDE CORRELATION MASTER v.s. DEM (SIMULATED AMPLITUDE)
	// ====== in order to estimate master timing error

	  if (input_general.process[pr_m_mtiming])
		{
		PROGRESS.print("Start MASTER TIMING.");
		alreadyprocessed[pr_m_mtiming] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing for master timing error using correlation.";
		  getanswer();
		  }

	//    if ( master.timingerror_flag == true )  // [MA] manual timing error cards no direct effect on correlation but influences sim. amp. 
	//      {
	//       WARNING << "MTE: M_AZ_T_ERROR and M_RG_T_ERROR cards are already set, please check. This does effect amplitude simulation step." ;
	//       WARNING.print(); 
	//      }  // TODO check that this is not necassary addition

		// ______Select method______
		switch (input_m_mtiming.method)
		  {
		  case cc_magfft:
			mtiming_correlfft(input_m_mtiming, master, simamp); // infoimage
			break;
		  case cc_magspace:
			mtiming_correl(input_m_mtiming, master, simamp); // infoimage
			break;
		  default:
			PRINT_ERROR("impossible error for checked input.")
			throw(input_error); // exit
		  }

		// ______ Update log files ______
		updatefile("scratchlogmtiming", input_general.logfile);
		updatefile("scratchresmtiming", input_general.m_resfile);
		PROGRESS.print("Finished MASTER TIMING.");
		DEBUG.print("Time spent for master timing error estimation with correlation:");
		printcpu();
		}

	// ______ Read estimates to interferogram info struct ______
	  if (alreadyprocessed[pr_m_mtiming]) // [MA]
		{
		INFO.print("");
		INFO.print("master: latest known processing stage: master timing");

		if (master.timingerror_flag == true) // [MA] manual timing error cards no direct effect on correlation but influences sim. amp.
		  { // if sb put timing cards and forgot to delete mtiming results, give warning since time is updated twice.
		   WARNING << "MTE: M_AZ_T_ERROR and M_RG_T_ERROR cards are already set, please check. Make sure you don't update master time twice by accident.";
		   WARNING.print();
		  }

		   String c6aztimingerrorlines = new String(new char[7]); // estimate azimuth timing error (lines)
		   String c6rtimingerrorpixels = new String(new char[7]); // estimate range timing error (pixels)
		   String c16aztimingerrorsec = new String(new char[17]); // estimate azimuth timing error (sec)
		   String c16rtimingerrorsec = new String(new char[17]); // estimate range timing error (sec)
		   readres(c6aztimingerrorlines,sizeof(c6aztimingerrorlines), input_general.m_resfile,"Coarse_correlation_translation_lines", 1);
		   readres(c6rtimingerrorpixels,sizeof(c6rtimingerrorpixels), input_general.m_resfile,"Coarse_correlation_translation_pixels", 1);
		   readres(c16aztimingerrorsec,sizeof(c16aztimingerrorsec), input_general.m_resfile,"Master_azimuth_timing_error", 1);
		   readres(c16rtimingerrorsec,sizeof(c16rtimingerrorsec), input_general.m_resfile,"Master_range_timing_error", 1);

		   DEBUG << "original master azimuth time   : " << setprecision(16) << master.t_azi1 << "\t   range time\t\t  : " << master.t_range1;
		   DEBUG.print();
		   DEBUG << "master azimuth time error lines: " << c6aztimingerrorlines << "\t\t   range time error pixels: " << c6rtimingerrorpixels;
		   DEBUG.print();
		   DEBUG << "master azimuth time error [sec]: " << Double.parseDouble(c16aztimingerrorsec) << "\t   range time error [sec] : " << Double.parseDouble(c16rtimingerrorsec);
		   DEBUG.print();

		   // _____ Update Master Azimuth and Range time _____
		   // master.add_offsetL2az_t_error(-atoi(c6aztimingerrorlines)); // read offsetL (master-dem) and update time (more precise?)
		   // master.add_offsetP2rg_t_error(-atoi(c6rtimingerrorpixels)); // read offsetP (master-dem) and update time
		   master.add_az_t_error(Double.parseDouble(c16aztimingerrorsec)); // read offsetL in seconds from the res filei and update master timing.
		   master.add_rg_t_error(Double.parseDouble(c16rtimingerrorsec));

		   INFO << "updated master azimuth time    : " << master.t_azi1 << " range time\t\t  : " << master.t_range1;
		   INFO.print();

		   // master.az_timing_error = -atoi(c6aztimingerrorlines);   // (never used) [MA]
		   // master.r_timing_error  = -atoi(c6rtimingerrorpixels);


		  // NOTE: time update is not necessary before reading orbital data
		  //       using getorb since we have extra-time. However if there 
		  //       is shift in more than PRF (lines) than it is better to set M_AZ_T_ERROR cards to read orbits with updated timing.

		} // alreadyprocessed mtiming



	// ____ end added by MA ____


	// *********************************************************
	// **** start slave                                     ****
	// *********************************************************

	// ====== PROCESS INFOFILE CDROM: SLAVE ======
	  if (input_general.process[pr_s_readfiles])
		{
		PROGRESS.print("Start S_READFILES.");
		// ___ update alreadyprocessed for this step ___
		alreadyprocessed[pr_s_readfiles] =1;
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing readfiles for slave.";
		  getanswer();
		  }

		// ______ Select method for SLAVE READFILES______
		switch (input_s_readfiles.sensor_id)
		  {
		  case SLC_ALOS:
			INFO.print("ALOS: Under development: is ALOS really CEOS/Atlantis like ERS?");
		  // ______ JERS also ceos? (falls through to ers) ______
		  case SLC_JERS:
			INFO.print("JERS: Under development: is JERS really CEOS like ERS?");
		  // ______ RSAT also CEOS? (falls through to ERS, account in reader) ______
		  case SLC_RSAT:
			WARNING.print("RSAT: CEOS reader seems to work.");
		  // ______ ERS-1/2 ______
		  case SLC_ERS:
			// ______ Get checks with master volumefile ______
			String c16checkvol1 = new String(new char[17]);
			String c16checkvol2 = new String(new char[17]);
			String c16checkvol3 = new String(new char[17]); // check id of volume file
			if (alreadyprocessed[pr_m_readfiles])
			  {
			  readres(c16checkvol1,sizeof(c16checkvol1), input_general.m_resfile,"Volume_ID:");
			  readres(c16checkvol2,sizeof(c16checkvol2), input_general.m_resfile,"Volume_identifier:");
			  readres(c16checkvol3,sizeof(c16checkvol3), input_general.m_resfile,"Volume_set_identifier:");
			  }

			// ______ Extract info ______
			readvolume(input_s_readfiles, c16checkvol1, c16checkvol2, c16checkvol3);
			// created in readvolume
			updatefile("scratchresvol", input_general.s_resfile);
			// created in readvolume
			updatefile("scratchlogvol", input_general.logfile);

			// ______ Get check with slave volumefile ______
			String c8checkleadat = new String(new char[9]); // check vol met lea en dat
			readres(c8checkleadat,sizeof(c8checkleadat), input_general.s_resfile,"(Check)Number",5);

			// ______ Extract info ______
			readleader(input_s_readfiles, Integer.parseInt(c8checkleadat)-1);
			updatefile("scratchreslea", input_general.s_resfile); // made in readleader
			updatefile("scratchloglea", input_general.logfile); // made in readleader

			//  // DONT PROCESS NULLFILE
			//readnull(input_s_readfiles);
			//updatefile("scratchresnull",input_general.s_resfile);
			//updatefile("scratchlognull",input_general.logfile);
			// ___ read data file info ___
			readdat(input_s_readfiles, Integer.parseInt(c8checkleadat)-1);
			updatefile("scratchresdat", input_general.s_resfile); // update resfile
			updatefile("scratchlogdat", input_general.logfile); // update logfile
			break; // (this break was missing in v3.4, thanks Raffaele Nutricato)
		  // ______ ENVISAT ASAR ______
		  case SLC_ASAR:
			INFO.reset(); // make sure nothing in buffer
			INFO << "envisat_dump_header2doris.csh " << input_s_readfiles.datfile << " > scratchres_envisat" << "\n" << ends;
			String cmd = new String(new char[512]); // command string
			cmd = INFO.get_str();
			INFO.print("With following command the envisat header was read.");
			INFO.print(cmd);
			PROGRESS.print("Making system call to envisat_dump_header2doris.csh");
			PROGRESS.print("(also requires envisat_dump_header program)");
			system(cmd); // This does the actual work
			INFO.reset();
			PROGRESS.print("Finished system call to envisat_dump_header2doris.csh");
			updatefile("scratchres_envisat", input_general.s_resfile); // update resfile
			updatefile("envisat_dump_header.log", input_general.logfile); // update logfile
			break;
		// ______ TSX COSAR ______
		case SLC_TSX:
		  INFO.reset(); // make sure nothing in buffer
		  INFO << "tsx_dump_header2doris.py " << input_s_readfiles.leaderfile << " > scratchres_tsx" << "\n" << ends;
		  //      char cmd[512];// command string
		  cmd = INFO.get_str();
		  INFO.print("With following command the TSX header was read.");
		  INFO.print(cmd);
		  PROGRESS.print("Making system call to tsx_dump_header2doris.py");
		  PROGRESS.print("(also requires python, gdal and XML libraries)");
		  status =system(cmd); // this does the work
		  if (status != 0) // [MA] TODO make it a function
			{
			ERROR << "tsx_dump_header2doris.py: failed with exit code: " << status;
			PRINT_ERROR(ERROR.get_str())
			throw(some_error);
			}
		  INFO.reset();
		  PROGRESS.print("Finished system call to tsx_dump_header2doris.py");
		  // ___ update resfile ___
		  updatefile("scratchres_tsx", input_general.s_resfile);
		  // ___ update logfile ___
		  // updatefile("tsx_dump_header.log",input_general.logfile);
		  break;
		  default:
			PRINT_ERROR("S_READFILES step is impossible, check input and data for the sensor parameters.")
			throw(input_error); // exit
		  }
		PROGRESS.print("Finished S_READFILES.");
		DEBUG.print("Time spent for reading files slave:");
		printcpu();
		}

	  // ______ Fill slcimage struct ______
	  if (existed(input_general.s_resfile))
		slave.fillslcimage(input_general.s_resfile);

	  // ______ Check sensor correspondence ______
	  if (alreadyprocessed[pr_m_readfiles] && alreadyprocessed[pr_s_readfiles])
		if (master.sensor != slave.sensor)
		  WARNING.print("master.sensor not same as slave.sensor (ERS-ERS, ASAR-ASAR only)");


	  // ______ Account for slave timing errors ______
	  if (input_s_readfiles.az_timing_error !=0.0 || input_s_readfiles.rg_timing_error != 0.0)
		{
		 PROGRESS.print("The correction values for slave image timing errors are specified.");
		 PROGRESS << "s_az_timing_error: " << input_s_readfiles.az_timing_error;
		 PROGRESS.print();
		 PROGRESS << "s_rg_timing_error: " << input_s_readfiles.rg_timing_error;
		 PROGRESS.print();
		 slave.add_rg_t_error(input_s_readfiles.rg_timing_error);
		 slave.add_az_t_error(input_s_readfiles.az_timing_error);
		 slave.timingerror_flag = true; // (true) then time is updated manually [MA]
		}



	// BK 15-Dec-2003 (moved before crop since required if dbow_geo
	// ====== GET PRECISE ORBITS (GETORB): SLAVE ======
	  if (input_general.process[pr_s_porbits])
		{
		PROGRESS.print("Start S_PORBITS.");
		alreadyprocessed[pr_s_porbits] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing get precise orbits for slave.";
		  getanswer();
		  }
		// ______getorb does the real work______
		getorb(slave,input_porbits,SLAVEID);
		convertgetorbout(SLAVEID);
		updatefile("scratchdatapoints", input_general.s_resfile); // result of conv.getorb slave
		removedatleader(input_general.s_resfile); // remove data section by leader
		PROGRESS.print("Finished S_PORBITS.");
		DEBUG.print("Time spent for obtaining precise orbits slave:");
		printcpu();
		}
	  // ______ Get coefficients for orbit interpolation if data in file ______
	  // (alreadyprocessed[pr_s_readfiles] && isspace(input_porbits.s_orbdir[0])))
	  if (alreadyprocessed[pr_s_porbits] || alreadyprocessed[pr_s_readfiles])
		slaveorbit.initialize(input_general.s_resfile);

	  // ______ Dump interpolated orbit if requested ______
	  //if (alreadyprocessed[pr_s_porbits])
	  // ___ compbaseline return if no orbit available ___
	  // ______ also possible to dump orbit of readfiles with card M_ORB_DUMP ___
	  slaveorbit.dumporbit(input_porbits, SLAVEID);


	  // --- initialize BASELINE class if orbits available ---
	  // Bert Kampes, 04-Apr-2005
	  BASELINE baseline;
	  baseline.model_parameters(master,slave,masterorbit,slaveorbit,input_ellips);


	  // ______ Stdout baseline parametrizations as INFO ______
	  //if (alreadyprocessed[pr_m_porbits] && alreadyprocessed[pr_s_porbits])
	  // ___ compbaseline return if no orbit available ___
	  //compbaseline(input_general, master, slave, masterorbit, slaveorbit);
	  for (register int32 i =0; i<input_general.dumpbaselineL; ++i) // azimuthlines
		{
		final real8 line = (input_general.dumpbaselineL==1) ? master.currentwindow.linelo + master.currentwindow.lines()/2.0 : master.currentwindow.linelo + i *master.currentwindow.lines()/ (input_general.dumpbaselineL-1.0);
		for (register int32 j =0; j<input_general.dumpbaselineP; ++j) // rangepixels
		  {
		  final real8 pixel = (input_general.dumpbaselineP==1) ? master.currentwindow.pixlo + master.currentwindow.pixels()/2.0 : master.currentwindow.pixlo + j *master.currentwindow.pixels()/ (input_general.dumpbaselineP-1.0);
		  for (register int32 height =0; height<1001; height=height+1000) // height
			baseline.dump(line,pixel,height); // eval parameters on grid
		  }
		}


	  // ====== Get tiepoint coordinates etc., this returns if nothing to do ======
	  DEBUG.print("Computating integration constant based on tiepoint");
	  tiepoint(input_general, master, slave, masterorbit, slaveorbit, input_ellips);



	// ====== PROCESS DATAFILE CDROM: SLAVE ======
	  if (input_general.process[pr_s_crop])
		{
		PROGRESS.print("Start S_CROP.");
		// ___ update alreadypr. ___
		alreadyprocessed[pr_s_crop] =1;
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing crop for slave.";
		  getanswer();
		  }

		// ______ Assume orbits loaded, use GEO card if present ______
		if (input_s_crop.dbow_geo.pixhi != 0)
		  {
		  PROGRESS.print("Computing s_crop based on center lat/lon/height/width.");
		  real8 centerline;
		  real8 centerpix;
		  real8 centerphi = deg2rad(input_s_crop.dbow_geo.linelo/1e6-360.0);
		  real8 centerlambda = deg2rad(input_s_crop.dbow_geo.linehi/1e6-360.0);
		  real8 centerheight = input_general.terrain_height;
		  DEBUG.print("Using height of HEIGHT card");
		  cn centerpos; // compute phi/lambda/h-->x,y,z
		  //ell2xyz(input_ellips,centerpos,centerphi,centerlambda,centerheight);
		  centerpos = input_ellips.ell2xyz(centerphi,centerlambda,centerheight);
		  DEBUG << "Converted S_DBOW_GEO phi,lambda,hei --> x,y,z: " << centerphi << ", " << centerlambda << ", " << centerheight << " --> " << centerpos.x << " " << centerpos.y << " " << centerpos.z;
		  DEBUG.print();
		  DEBUG << "Center slave SLC: x,y,z: " << slave.approxcentreoriginal.x << " " << slave.approxcentreoriginal.y << " " << slave.approxcentreoriginal.z;
		  DEBUG.print();
		  if (Math.abs(slave.approxcentreoriginal.x-centerpos.x)>60000.0) // 60km
			WARNING.print("S_DBOW_GEO: coordinates seem to be outside SLC area? (X)");
		  if (Math.abs(slave.approxcentreoriginal.y-centerpos.y)>60000.0) // 60km
			WARNING.print("S_DBOW_GEO: coordinates seem to be outside SLC area? (Y)");
		  if (Math.abs(slave.approxcentreoriginal.z-centerpos.z)>60000.0) // 60km
			WARNING.print("S_DBOW_GEO: coordinates seem to be outside SLC area? (Z)");
		  RefObject<real8> TempRefObject8 = new RefObject<real8>(centerline);
		  RefObject<real8> TempRefObject9 = new RefObject<real8>(centerpix);
		  RefObject<orbit> TempRefObject10 = new RefObject<orbit>(slaveorbit);
		  xyz2lp(TempRefObject8, TempRefObject9, slave, TempRefObject10, centerpos, 10, 1e-3); // MAXITERATIONS, CONVERGENCE_TIME
		  centerline = TempRefObject8.argvalue;
		  centerpix = TempRefObject9.argvalue;
		  slaveorbit = TempRefObject10.argvalue;
		  DEBUG << "CROP: center line/pix: " << centerline << " " << centerpix;
		  DEBUG.print();
		  int32 l0 = int32(centerline+0.5) - input_s_crop.dbow_geo.pixlo/2;
		  int32 lN = l0 + input_s_crop.dbow_geo.pixlo - 1;
		  int32 p0 = int32(centerpix+0.5) - input_s_crop.dbow_geo.pixhi/2;
		  int32 pN = p0 + input_s_crop.dbow_geo.pixhi - 1;
		  if (l0 < int32(slave.currentwindow.linelo))
			  l0 = slave.currentwindow.linelo;
		  if (lN > int32(slave.currentwindow.linehi))
			  lN = slave.currentwindow.linehi;
		  if (p0 < int32(slave.currentwindow.pixlo))
			  p0 = slave.currentwindow.pixlo;
		  if (pN > int32(slave.currentwindow.pixhi))
			  pN = slave.currentwindow.pixhi;
		  INFO << "DBOW from GEO: " << l0 << " " << lN << " " << p0 << " " << pN;
		  INFO.print();
		  // ___ Simply fill dbow with it and it will be done correctly! ___
		  input_s_crop.dbow.linelo = uint(l0);
		  input_s_crop.dbow.linehi = uint(lN);
		  input_s_crop.dbow.pixlo = uint(p0);
		  input_s_crop.dbow.pixhi = uint(pN);
		  }

		// ______ For each sensor, of course the formats differ... ______
		switch (slave.sensor)
		  {
		  // _____ start added by don       
		  case SLC_ALOS:
				 String c8checkleadata = new String(new char[9]); // check vol met lea en dat
			readres(c8checkleadata,sizeof(c8checkleadata), input_general.s_resfile,"(Check)Number",5);
			palsar_fine_dump_data(input_general,input_s_crop, Integer.parseInt(c8checkleadata)-1);
			updatefile("scratchres2raw", input_general.s_resfile); // update resfile
			   break;
		  // _____ end added by don 
		  // ______ JERS also CEOS? (falls through to ers) ______
		  case SLC_JERS:
			WARNING.print("JERS: Under development: is JERS really CEOS like ERS?");
		  // ______ ERS-1/2 ______
		  case SLC_ERS:
			switch (slave.sar_processor)
			  {
			  case SARPR_VMP:
				// --- default, ESA SLCs ---
				// ______ Get check with volumefile ______
				String c8checkleadat = new String(new char[9]); // check vol met lea en dat
				readres(c8checkleadat,sizeof(c8checkleadat), input_general.s_resfile,"(Check)Number",5);
				writeslc(input_general,input_s_crop, Integer.parseInt(c8checkleadat)-1);
				updatefile("scratchres2raw", input_general.s_resfile); // update resfile
				break;
			  case SARPR_ATL:
				// --- self processed from RAW using Atlantis ---
				radarsat_dump_data(input_general,input_s_crop);
				updatefile("scratchres2raw", input_general.s_resfile); // update resfile
				WARNING.print("ATLANTIS data cropping not fully checked yet.");
				break;
			  case SARPR_TUD:
				//break;
			  default:
				PRINT_ERROR("unrecognized sar_processor.")
				throw(input_error); // exit
			  }
			break;
		  // ______ ENVISAT ASAR ______
		  case SLC_ASAR:
			PROGRESS.print("System call to get ASAR SLC data (requires envisat_dump_data program)");
			// ______ Data window must be set correctly ______
			if (input_s_crop.dbow.linelo == 0 || input_s_crop.dbow.linehi == 0 || input_s_crop.dbow.pixlo == 0 || input_s_crop.dbow.pixhi == 0)
			  input_s_crop.dbow = slave.currentwindow;
			envisat_dump_data(input_s_crop);
			PROGRESS.print("Finished system call to envisat_dump_data");
			updatefile("scratchres2raw", input_general.s_resfile); // update resfile
			break;
		  // ______ TSX COSAR ______
		case SLC_TSX:
		  PROGRESS.print("System call to get TSX SLC/COSAR data (requires tsx_dump_data program)");
		  // ______ Data window must be set correctly ______
		  if (input_s_crop.dbow.linelo == 0 || input_s_crop.dbow.linehi == 0 || input_s_crop.dbow.pixlo == 0 || input_s_crop.dbow.pixhi == 0)
			{
			  input_s_crop.dbow = slave.currentwindow; // needs to be set correctly after readfiles
			}
		  tsx_dump_data(input_s_crop);
		  PROGRESS.print("Finished system call to tsx_dump_data");
		  updatefile("scratchres2raw", input_general.s_resfile); // update resfile
		  break;
		  // ______ RADARSAT ______
		  case SLC_RSAT:
			radarsat_dump_data(input_general,input_s_crop);
			updatefile("scratchres2raw", input_general.s_resfile); // update resfile
			WARNING.print("RSAT (Atlantis) not fully checked yet.");
			break;
		  default:
			PRINT_ERROR("S_CROP step is impossible, check input and data for the sensor parameters.")
			//PRINT_ERROR("impossible error (sensor not (J)ERS or ASAR).")
			throw(input_error); // exit
		  }
		PROGRESS.print("Finished S_CROP.");
		DEBUG.print("Time spent for writing raw format slave:");
		printcpu();
		}

	  // ______ Update info on image ______
	  if (alreadyprocessed[pr_s_crop])
		{
		INFO.print("");
		INFO.print("slave: latest known processing stage: crop");
		String SECTIONID = new String(new char[ONE27]);
		SECTIONID = "*_Start_";
		SECTIONID += processcontrol[pr_s_crop];
		slave.updateslcimage(input_general.s_resfile,SECTIONID);
		}

	  // ______ Check wavelength ______
	  if (alreadyprocessed[pr_m_crop] && alreadyprocessed[pr_s_crop])
		{
		if (Math.abs(master.wavelength-slave.wavelength) > EPS)
		  WARNING.print("wavelength master not equal to wavelength slave.");
		}


	  // ______ Generate magnitude preview if requested ______
	  if (input_general.process[pr_s_crop])
		{
		PROGRESS.print("calling preview for cropped slave");
		if (slave.sensor!=SLC_RSAT && slave.sensor!=SLC_ALOS && slave.sensor!=SLC_TSX)
		  preview(input_general.preview, slave.currentwindow.pixels(), slave.formatflag, slave.file, "slave_mag.ras", "-e 0.5 -s 1.0 -q mag -o sunraster -b -c gray -M 2/10");
		else if (slave.sensor == SLC_ALOS || slave.sensor == SLC_TSX)
		  preview(input_general.preview, slave.currentwindow.pixels(), slave.formatflag, slave.file, "slave_mag.ras", "-e 0.5 -s 1.0 -q mag -o sunraster -b -c gray -M 10/10"); // [MA] higher mlook to reduce filesize
		else // multilooking depends on beam
		  preview(input_general.preview, slave.currentwindow.pixels(), slave.formatflag, slave.file, "slave_mag.ras", "-e 0.5 -s 1.0 -q mag -o sunraster -b -c gray -M 5/6");
		}



	//____RaffaeleNutricato START MODIFICATION SECTION 6
	// ====== OVERSAMPLE SLC: SLAVE ======
	// For the moment only along the range direction
	  if (input_general.process[pr_s_oversample])
		{
		PROGRESS.print("Start S_OVS.");
		alreadyprocessed[pr_s_oversample] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing oversample for slave.";
		  getanswer();
		  }
		// ______Master and Slave oversampling ratioes should be identical______
		if (alreadyprocessed[pr_m_oversample])
		  {
		  String checkmultilook = new String(new char[20]);
		  int32 checkOsrRange;
		  RefObject<String> TempRefObject11 = new RefObject<String>(checkmultilook);
		  readres(TempRefObject11, sizeof(checkmultilook), input_general.m_resfile, "Multilookfactor_range_direction:", 0, 1);
		  checkmultilook = TempRefObject11.argvalue;
		  checkOsrRange = int32(1.0/Double.parseDouble(checkmultilook)+0.5);
		  if (input_s_oversample.OsrRange != checkOsrRange)
			WARNING.print("Master and Slave range oversampling factors should be identical!!!");
		  }
		if (alreadyprocessed[pr_m_oversample])
		  {
		  String checkmultilook = new String(new char[20]);
		  int32 checkOsrAzimuth;
		  RefObject<String> TempRefObject12 = new RefObject<String>(checkmultilook);
		  readres(TempRefObject12, sizeof(checkmultilook), input_general.m_resfile, "Multilookfactor_azimuth_direction:", 0, 1);
		  checkmultilook = TempRefObject12.argvalue;
		  checkOsrAzimuth = int32(1.0/Double.parseDouble(checkmultilook)+0.5);
		  if (input_s_oversample.OsrAzimuth != checkOsrAzimuth)
			WARNING.print("Master and Slave azimuth oversampling factors should be identical!!!");
		  }

		// ______ Oversample slave image ______
		OversampleSLC(input_general,slave,input_s_oversample,SLAVEID);
		updatefile("scratchoversample", input_general.s_resfile); // update resfile
		PROGRESS.print("Finished S_OVS.");
		DEBUG.print("Time spent for oversampling slave:");
		printcpu();
		}

	  // ______ Update info on image ______
	  if (alreadyprocessed[pr_s_oversample])
		{
		INFO.print("");
		INFO.print("slave: latest known processing stage: oversample");
		String SECTIONID = new String(new char[ONE27]);
		SECTIONID = "*_Start_";
		SECTIONID += processcontrol[pr_s_oversample];
		slave.updateslcimage(input_general.s_resfile,SECTIONID); // Bert Kampes, 29-Jul-2005
		}

	  // ______ Generate magnitude preview of slave SLC if requested ______
	  // -M 2/10 make no sense after oversampling!! (depends on factors..)
	  if (input_general.process[pr_s_oversample])
		{
		PROGRESS.print("calling preview for oversampled slave");
		preview(input_general.preview, slave.currentwindow.pixels(), slave.formatflag, slave.file, "slave_ovs_mag.ras", "-e 0.5 -s 1.0 -q mag -o sunraster -b -c gray -M 4/4");
		}
	//____RaffaeleNutricato END MODIFICATION SECTION 6




	// *********************************************************
	// **** start interferogram                             ****
	// *********************************************************

	// ====== COARSE CO-REGISTRATION BASED ON ORBITS ======
	  if (input_general.process[pr_i_coarse])
		{
		PROGRESS.print("Start COARSE_ORB.");
		alreadyprocessed[pr_i_coarse] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing coarse co-registration based on orbits.";
		  getanswer();
		  }
		RefObject<orbit> TempRefObject13 = new RefObject<orbit>(masterorbit);
		RefObject<orbit> TempRefObject14 = new RefObject<orbit>(slaveorbit);
		coarseporbit(input_ellips, master, slave, TempRefObject13, TempRefObject14, baseline);
		masterorbit = TempRefObject13.argvalue;
		slaveorbit = TempRefObject14.argvalue;

		// ______ Update log files ______
		updatefile("scratchlogcoarse", input_general.logfile);
		updatefile("scratchrescoarse", input_general.i_resfile);
		PROGRESS.print("Finished COARSE_ORB.");
		DEBUG.print("Time spent for coarse coregistration with orbits:");
		printcpu();
		}

	// ______ Read estimates to interferogram info struct ______
	  if (alreadyprocessed[pr_i_coarse])
		{
		INFO.print("");
		INFO.print("ifg: latest known processing stage: coarse (orbits)");
		String c6initoffL = new String(new char[7]); // initial offset lines
		String c6initoffP = new String(new char[7]); // initial offset pixels
		readres(c6initoffL,sizeof(c6initoffL),input_general.i_resfile, "Coarse_orbits_translation_lines:", 0);
		readres(c6initoffP,sizeof(c6initoffP),input_general.i_resfile, "Coarse_orbits_translation_pixels:", 0);
		slave.coarseoffsetL = Integer.parseInt(c6initoffL); // used as initial value
		slave.coarseoffsetP = Integer.parseInt(c6initoffP); // used as initial value
		slave.coarseorbitoffsetL = Integer.parseInt(c6initoffL); // to estimate timing error[FvL]
		slave.coarseorbitoffsetP = Integer.parseInt(c6initoffP); // to estimate timing error[FvL]
		//courseoffset is updated by coarse correlation, therefore extra
		//parameter courseorbitoffset [FvL]
		master.coarseoffsetL = -slave.coarseoffsetL; // (never used)
		master.coarseoffsetP = -slave.coarseoffsetP; // azifilt

		// ______ Corners of current slave in master coordinate system ______
		final int32 sL0 = slave.currentwindow.linelo - slave.coarseoffsetL;
		final int32 sLN = slave.currentwindow.linehi - slave.coarseoffsetL;
		final int32 sP0 = slave.currentwindow.pixlo - slave.coarseoffsetP;
		final int32 sPN = slave.currentwindow.pixhi - slave.coarseoffsetP;
		interferogram.win.linelo = max(sL0,int32(master.currentwindow.linelo));
		interferogram.win.linehi = min(sLN,int32(master.currentwindow.linehi));
		interferogram.win.pixlo = max(sP0,int32(master.currentwindow.pixlo));
		interferogram.win.pixhi = min(sPN,int32(master.currentwindow.pixhi));
		}



	// ====== COARSE CO-REGISTRATION BASED ON CORRELATION ======
	  if (input_general.process[pr_i_coarse2])
		{
		PROGRESS.print("Start COARSE_CORR.");
		alreadyprocessed[pr_i_coarse2] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing coarse co-registration with correlation.";
		  getanswer();
		  }

		// ______ If requested: use estimated offsets from orbits as initial ______
		if (input_i_coarsecorr.initoffsetL == NaN && input_i_coarsecorr.initoffsetL == NaN) // flag for using answer of coarse1
		  {
		  input_i_coarsecorr.initoffsetL = slave.coarseoffsetL;
		  input_i_coarsecorr.initoffsetP = slave.coarseoffsetP;
		  }

		// ______Select method______
		switch (input_i_coarsecorr.method)
		  {
		  case cc_magfft:
			coarsecorrelfft(input_i_coarsecorr, master, slave); // infoimage
			break;
		  case cc_magspace:
			coarsecorrel(input_i_coarsecorr, master, slave); // infoimage
			break;
		  default:
			PRINT_ERROR("impossible error for checked input.")
			throw(input_error); // exit
		  }

		// ______ Update log files ______
		updatefile("scratchlogcoarse2", input_general.logfile);
		updatefile("scratchrescoarse2", input_general.i_resfile);
		PROGRESS.print("Finished COARSE_CORR.");
		DEBUG.print("Time spent for coarse coregistration with correlation:");
		printcpu();
		}

	// ______ Read estimates to interferogram info struct ______
	  if (alreadyprocessed[pr_i_coarse2])
		{
		INFO.print("");
		INFO.print("ifg: latest known processing stage: coarse (corr)");
		String c6initoffL = new String(new char[7]); // initial offset lines
		String c6initoffP = new String(new char[7]); // initial offset pixels
		readres(c6initoffL,sizeof(c6initoffL),input_general.i_resfile, "Coarse_correlation_translation_lines:", 0);
		readres(c6initoffP,sizeof(c6initoffP),input_general.i_resfile, "Coarse_correlation_translation_pixels:", 0);
		slave.coarseoffsetL = Integer.parseInt(c6initoffL); // used as initial value
		slave.coarseoffsetP = Integer.parseInt(c6initoffP); // used as initial value
		master.coarseoffsetL = -slave.coarseoffsetL; // (never used)
		master.coarseoffsetP = -slave.coarseoffsetP; // azifilt

		// ______ corners of current slave in master coordinate system ______
		final int32 sL0 = slave.currentwindow.linelo - slave.coarseoffsetL;
		final int32 sLN = slave.currentwindow.linehi - slave.coarseoffsetL;
		final int32 sP0 = slave.currentwindow.pixlo - slave.coarseoffsetP;
		final int32 sPN = slave.currentwindow.pixhi - slave.coarseoffsetP;

		// fill these as well to be sure ? (BK)
		interferogram.win.linelo = max(sL0,int32(master.currentwindow.linelo));
		interferogram.win.linehi = min(sLN,int32(master.currentwindow.linehi));
		interferogram.win.pixlo = max(sP0,int32(master.currentwindow.pixlo));
		interferogram.win.pixhi = min(sPN,int32(master.currentwindow.pixhi));
		}



	// ====== FILTER AZIMUTH: MASTER (and SLAVE directly after) ======
	// ______ Do this step here, use coarse corr for fDC polynomial ______
	  if (input_general.process[pr_m_filtazi])
		{
		PROGRESS.print("Start FILTAZI (master).");
		alreadyprocessed[pr_m_filtazi] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing filtering azimuth direction for master.";
		  getanswer();
		  }
		// ______ only master;  output file name of master ______
		input_ms_filtazi.foname = input_ms_filtazi.fomaster;
		RefObject<slcimage> TempRefObject15 = new RefObject<slcimage>(master);
		RefObject<slcimage> TempRefObject16 = new RefObject<slcimage>(slave);
		azimuthfilter(input_general, input_ms_filtazi, TempRefObject15, TempRefObject16);
		master = TempRefObject15.argvalue;
		slave = TempRefObject16.argvalue;
		updatefile("scratchresfiltazi", input_general.m_resfile); // update resfile
		updatefile("scratchlogfiltazi", input_general.logfile); // update logfile
		PROGRESS.print("Finished FILTAZI for master.");
		DEBUG.print("Time spent for azimuth filtering master:");
		printcpu();
		}
	  // ______ Update info of master, current filename, format, etc. ______
	  if (alreadyprocessed[pr_m_filtazi])
		{
		INFO.print("");
		INFO.print("master: latest known processing stage: azimuth filtered");
		String SECTIONID = new String(new char[ONE27]);
		SECTIONID = "*_Start_";
		SECTIONID += processcontrol[pr_m_filtazi];
		master.updateslcimage(input_general.m_resfile,SECTIONID);
		}


	// ====== FILTER AZIMUTH: SLAVE ======
	  if (input_general.process[pr_s_filtazi])
		{
		PROGRESS.print("Start FILTAZI (slave).");
		alreadyprocessed[pr_s_filtazi] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing filtering azimuth direction for slave.";
		  getanswer();
		  }
		// ______ only slave ... ______
		input_ms_filtazi.foname = input_ms_filtazi.foslave;
		RefObject<slcimage> TempRefObject17 = new RefObject<slcimage>(slave);
		RefObject<slcimage> TempRefObject18 = new RefObject<slcimage>(master);
		azimuthfilter(input_general, input_ms_filtazi, TempRefObject17, TempRefObject18);
		slave = TempRefObject17.argvalue;
		master = TempRefObject18.argvalue;
		updatefile("scratchresfiltazi", input_general.s_resfile); // update resfile
		updatefile("scratchlogfiltazi", input_general.logfile); // update logfile
		PROGRESS.print("Finished FILTAZI slave.");
		DEBUG.print("Time spent for azimuth filtering slave:");
		printcpu();
		}

	// ______ Update info of slave, filename current image, format, etc. ______
	  if (alreadyprocessed[pr_s_filtazi])
		{
		INFO.print("");
		INFO.print("slave: latest known processing stage: azimuth filtered");
		String SECTIONID = new String(new char[ONE27]);
		SECTIONID = "*_Start_";
		SECTIONID += processcontrol[pr_s_filtazi];
		slave.updateslcimage(input_general.s_resfile,SECTIONID);
		}



	// ====== FILTER RANGE: MASTER&SLAVE porbits ======
	// ______ Only one step for master and slave, after coarse coreg. ______
	// ______ Method based on orbits only. (For adaptive method see after resampling) ______
	// ______ Update of struct should only occur for method porbits, ______
	// ______ so first find that out (BK 26-mar-01) ______
	  if (input_general.process[pr_m_filtrange] && input_ms_filtrange.method==rf_porbits)
		{
		PROGRESS.print("Start FILTRANGE (porbits).");
		alreadyprocessed[pr_m_filtrange] =1; // update alreadypr.
		alreadyprocessed[pr_s_filtrange] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nRange filtering based on orbits master and slave.";
		  getanswer();
		  }

		// ______ Filter routine ______
		RefObject<orbit> TempRefObject19 = new RefObject<orbit>(masterorbit);
		RefObject<orbit> TempRefObject20 = new RefObject<orbit>(slaveorbit);
		rangefiltporbits(input_general, input_ms_filtrange, input_ellips, master, slave, TempRefObject19, TempRefObject20);
		masterorbit = TempRefObject19.argvalue;
		slaveorbit = TempRefObject20.argvalue;

		// ______ Update log files ______
		updatefile("scratchlogfiltrange", input_general.logfile);
		updatefile("scratchresMfiltrange", input_general.m_resfile);
		updatefile("scratchresSfiltrange", input_general.s_resfile);

		PROGRESS.print("Finished FILTRANGE.");
		DEBUG.print("Time spent for range filtering:");
		printcpu();
		}

	// ______ Update info of master, size interferogram ______
	// ______ Only update if method is porbits ______
	  if (alreadyprocessed[pr_m_filtrange])
		{
		String c10rfmethod = new String(new char[11]); // method
		readres(c10rfmethod,sizeof(c10rfmethod),input_general.m_resfile, "Method_rangefilt:", 0);
		if (!strcmp(c10rfmethod,"porbits"))
		  {
		  INFO.print("");
		  INFO.print("master: latest known processing stage: filtrange (orbits)");
		  String SECTIONID = new String(new char[ONE27]);
		  SECTIONID = "*_Start_";
		  SECTIONID += processcontrol[pr_m_filtrange];
		  master.updateslcimage(input_general.m_resfile,SECTIONID);
		  // ______ cut out overlay master slave exact ______
		  interferogram.win = master.currentwindow;
		  }
		}

	// ______ Update info of slave, size interferogram ______
	  if (alreadyprocessed[pr_s_filtrange]) // same time as _m_
		{
		String c10rfmethod = new String(new char[11]); // method
		readres(c10rfmethod,sizeof(c10rfmethod),input_general.s_resfile, "Method_rangefilt:", 0);
		if (!strcmp(c10rfmethod,"porbits"))
		  {
		  INFO.print("");
		  INFO.print("slave: latest known processing stage: filtrange (orbits)");
		  String SECTIONID = new String(new char[ONE27]);
		  SECTIONID = "*_Start_";
		  SECTIONID += processcontrol[pr_s_filtrange];
		  slave.updateslcimage(input_general.s_resfile,SECTIONID);
		  }
		}




	// ====== FINE CO-REGISTRATION ======
	  if (input_general.process[pr_i_fine])
		{
		PROGRESS.print("Start FINE.");
		alreadyprocessed[pr_i_fine] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing fine co-registration.";
		  getanswer();
		  }
		// ______ If requested: Read result coarse from file ______
		if (input_i_fine.initoffsetL == NaN && input_i_fine.initoffsetL == NaN) // asnwer of coarse2
		  {
		  input_i_fine.initoffsetL = slave.coarseoffsetL;
		  input_i_fine.initoffsetP = slave.coarseoffsetP;
		  }
		finecoreg(input_i_fine, master, slave); // selects other methods within

	//    // ______ Select method ______
	//    if (input_i_fine.method == fc_magfft  ||
	//        input_i_fine.method == fc_magspace  )
	//      finecoreg(input_i_fine, master, slave);    // selects other methods within
	//    else if (input_i_fine.method == fc_oversample) // still 2b implemented in place
	//      //finecoreg2(input_i_fine, master, slave)
	//      ;
	//    else
	//      {
	//      PRINT_ERROR("Method FC not implemented.")
	//      throw(input_error);// exit
	//      }

		// ______ Update log files ______
		updatefile("scratchlogfine", input_general.logfile);
		updatefile("scratchresfine", input_general.i_resfile);

		PROGRESS.print("Finished FINE.");
		DEBUG.print("Time spent for fine coregistration:");
		printcpu();

		// ______ Plot if requested, use gmt script: plotoffsets ______
		INFO.reset(); // make sure nothing in buffer
		INFO << "plotoffsets " << input_general.i_resfile << " " << master.currentwindow.linelo << " " << master.currentwindow.linehi << " " << master.currentwindow.pixlo << " " << master.currentwindow.pixhi << " " << input_i_fine.plotthreshold << " ";
		// ______ optional 7th argument for background ______
		if (input_i_fine.plotmagbg)
			INFO << master.file;
		INFO << "&" << "\n" << ends;
		String cmd = new String(new char[512]); // command string
		cmd = INFO.get_str();
		INFO.print("Following command plots the estimated offset vectors:");
		INFO.print(cmd);
		// ______ Actual make the call if requested ______
		if (input_i_fine.plotoffsets)
			system(cmd);

		}

	  // _______ Fill oversampling factor for coregpm routine ______
	  //uint osfactor = 32;// oversamplingsfactor
	  if (alreadyprocessed[pr_i_fine])
		{
		INFO.print("");
		INFO.print("ifg: latest known processing stage: fine");
		//char c4osfactor[4];
		//readres(c4osfactor,sizeof(c4osfactor),input_general.i_resfile,
		//        "Oversampling", 1);
		//osfactor = uint(atoi(c4osfactor));
		}


	// ____ start added by FvL ____

	// ====== RELATIVE TIMING ERROR ======
	  if (input_general.process[pr_i_timing])
		{
		PROGRESS.print("Start TIMING.");
		alreadyprocessed[pr_i_timing] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nComputation relative timing error.\n";
		  getanswer();
		  }

		// _____ Compute master-slave timing errors _____
		RefObject<Integer> TempRefObject21 = new RefObject<Integer>(slave.coarseorbitoffsetL);
		RefObject<Integer> TempRefObject22 = new RefObject<Integer>(slave.coarseorbitoffsetP);
		ms_timing_error(master, input_general.i_resfile, input_i_timing, TempRefObject21, TempRefObject22);
		slave.coarseorbitoffsetL = TempRefObject21.argvalue;
		slave.coarseorbitoffsetP = TempRefObject22.argvalue;

		// ______ Update log files ______
		updatefile("scratchlogtiming", input_general.logfile);
		updatefile("scratchrestiming", input_general.i_resfile);

		PROGRESS.print("Finished TIMING.");
		DEBUG.print("Time spent for computation relative timing error:");
		printcpu();
		}

	  // ______ Read estimates to slave info struct ______
	  if (alreadyprocessed[pr_i_timing])
		{
		String c6aztimingerrorlines = new String(new char[7]); // estimate azimuth timing error (lines)
		String c6rtimingerrorpixels = new String(new char[7]); // estimate range timing error (pixels)
		String c16aztimingerrorsec = new String(new char[17]); // estimate azimuth timing error (sec)
		String c16rtimingerrorsec = new String(new char[17]); // estimate range timing error (sec)
		readres(c6aztimingerrorlines,sizeof(c6aztimingerrorlines), input_general.i_resfile,"Estimated_azimuth_timing_error_lines", 1);
		readres(c6rtimingerrorpixels,sizeof(c6rtimingerrorpixels), input_general.i_resfile,"Estimated_range_timing_error_pixels", 1);
		readres(c16aztimingerrorsec,sizeof(c16aztimingerrorsec), input_general.i_resfile,"Estimated_azimuth_timing_error_sec", 1);
		readres(c16rtimingerrorsec,sizeof(c16rtimingerrorsec), input_general.i_resfile,"Estimated_range_timing_error_sec", 1);

		slave.add_az_t_error(Double.parseDouble(c16aztimingerrorsec));
		slave.add_rg_t_error(Double.parseDouble(c16rtimingerrorsec));
		slave.az_timing_error = Integer.parseInt(c6aztimingerrorlines);
		slave.r_timing_error = Integer.parseInt(c6rtimingerrorpixels);
		INFO << Double.parseDouble(c16aztimingerrorsec);
		INFO.print();
		INFO << Double.parseDouble(c16rtimingerrorsec);
		INFO.print();
		INFO << Integer.parseInt(c6aztimingerrorlines);
		INFO.print();
		INFO << Integer.parseInt(c6rtimingerrorpixels);
		INFO.print();

		}


	// ====== DEM ASSISTED COREGISTRATION ======
	  if (input_general.process[pr_i_demassist])
		{
		PROGRESS.print("Start DEMASSIST.");
		alreadyprocessed[pr_i_demassist] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nCoregistration using DEM.\n";
		  getanswer();
		  }
		demassist(input_general, input_ellips, input_i_demassist, master, slave, masterorbit, slaveorbit);

		// ______ Update log files ______
		updatefile("scratchlogdemassist", input_general.logfile);
		updatefile("scratchresdemassist", input_general.i_resfile);

		PROGRESS.print("Finished DEMASSIST.");
		DEBUG.print("Time spent for DEM assisted coregistration:");
		printcpu();
		}

	  // ______ Read estimates to slave info struct ______
	  if (alreadyprocessed[pr_i_demassist])
		{
		String c16slavel00 = new String(new char[17]); // delta line slave00
		String c16slavep00 = new String(new char[17]); // delta pixel slave00
		String c16slavel0N = new String(new char[17]); // delta line slave0N
		String c16slavep0N = new String(new char[17]); // delta pixel slave0N
		String c16slavelN0 = new String(new char[17]); // delta line slaveN0
		String c16slavepN0 = new String(new char[17]); // delta pixel slaveN0
		String c16slavelNN = new String(new char[17]); // delta line slaveNN
		String c16slavepNN = new String(new char[17]); // delta pixel slaveNN
		readres(c16slavel00,sizeof(c16slavel00), input_general.i_resfile,"Deltaline_slave00_dem:", 0);
		readres(c16slavep00,sizeof(c16slavep00), input_general.i_resfile,"Deltapixel_slave00_dem:", 0);
		readres(c16slavel0N,sizeof(c16slavel0N), input_general.i_resfile,"Deltaline_slave0N_dem:", 0);
		readres(c16slavep0N,sizeof(c16slavep0N), input_general.i_resfile,"Deltapixel_slave0N_dem:", 0);
		readres(c16slavelN0,sizeof(c16slavelN0), input_general.i_resfile,"Deltaline_slaveN0_dem:", 0);
		readres(c16slavepN0,sizeof(c16slavepN0), input_general.i_resfile,"Deltapixel_slaveN0_dem:", 0);
		readres(c16slavelNN,sizeof(c16slavelNN), input_general.i_resfile,"Deltaline_slaveNN_dem:", 0);
		readres(c16slavepNN,sizeof(c16slavepNN), input_general.i_resfile,"Deltapixel_slaveNN_dem:", 0);

		slave.add_offsetl00(Double.parseDouble(c16slavel00));
		slave.add_offsetp00(Double.parseDouble(c16slavep00));
		slave.add_offsetl0N(Double.parseDouble(c16slavel0N));
		slave.add_offsetp0N(Double.parseDouble(c16slavep0N));
		slave.add_offsetlN0(Double.parseDouble(c16slavelN0));
		slave.add_offsetpN0(Double.parseDouble(c16slavepN0));
		slave.add_offsetlNN(Double.parseDouble(c16slavelNN));
		slave.add_offsetpNN(Double.parseDouble(c16slavepNN));
		INFO << Double.parseDouble(c16slavel00);
		INFO.print();
		INFO << Double.parseDouble(c16slavep00);
		INFO.print();
		INFO << Double.parseDouble(c16slavel0N);
		INFO.print();
		INFO << Double.parseDouble(c16slavep0N);
		INFO.print();
		INFO << Double.parseDouble(c16slavelN0);
		INFO.print();
		INFO << Double.parseDouble(c16slavepN0);
		INFO.print();
		INFO << Double.parseDouble(c16slavelNN);
		INFO.print();
		INFO << Double.parseDouble(c16slavepNN);
		INFO.print();

		}

	// ____ end added by FvL ____


	// ====== COMPUTATION CO-REGISTRATION PARAMETERS ======
	  if (input_general.process[pr_i_coregpm])
		{
		PROGRESS.print("Start COREGPM.");
		alreadyprocessed[pr_i_coregpm] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing computation of co-registration parameters.";
		  getanswer();
		  }
		//coregpm(master, input_general.i_resfile, input_i_coregpm, osfactor);
		coregpm(master, slave, input_general.i_resfile, input_i_coregpm, alreadyprocessed[pr_i_demassist]);
		// ______ Update log files ______
		updatefile("scratchlogcpm", input_general.logfile);
		updatefile("scratchrescpm", input_general.i_resfile);
		PROGRESS.print("Finished COREGPM.");
		DEBUG.print("Time spent for computation of coregistration parameters:");
		printcpu();

		// ______ Plot if requested, use gmt script: plotcpm ______
		INFO.reset(); // make sure nothing in buffer
		INFO << "plotcpm CPM_Data " << master.currentwindow.linelo << " " << master.currentwindow.linehi << " " << master.currentwindow.pixlo << " " << master.currentwindow.pixhi << " ";
		// ______ optional 6th argument for background ______
		if (input_i_coregpm.plotmagbg)
			INFO << master.file;
		// not ok? BK 27-Jun-00: omem << "&" << ends;
		INFO << "&" << "\n" << ends;
		String cmd = new String(new char[512]); // command string
		cmd = INFO.get_str();
		INFO.print("Next command will plot solution:");
		INFO.print(cmd);
		// ______ Actual make the call if requested ______
		if (input_i_coregpm.plot)
			system(cmd);
		}

	  // ______ Fill matrix coefficients, interferogram info struct ______
	  // ______ Approximate overlap in master system in interf. window ______
	  if (alreadyprocessed[pr_i_coregpm])
		{
		INFO.print("");
		INFO.print("ifg: latest known processing stage: coregpm");
		String c10offL = new String(new char[11]);
		readres(c10offL,sizeof(c10offL),input_general.i_resfile,"Degree_cpm:");
		int32 degreecpm = Integer.parseInt(c10offL);
		coeff_cpmL = readcoeff(input_general.i_resfile, "Estimated_coefficientsL:", Ncoeffs(degreecpm));
		coeff_cpmP = readcoeff(input_general.i_resfile, "Estimated_coefficientsP:", Ncoeffs(degreecpm));
		// bk 1 sep 2000, req. for resample...
		//interferogram.win = getoverlap(master,slave,coeff_cpmL,coeff_cpmP);


		// ______ start added by FvL _____
		if (input_general.process[pr_s_resample]) // [MA] quick soln. if no resampling no need to read these parameters? To increase compatiblity with older result files
		  {
		  String c16slavel00 = new String(new char[17]); // delta line slave00
		  String c16slavep00 = new String(new char[17]); // delta pixel slave00
		  String c16slavel0N = new String(new char[17]); // delta line slave0N
		  String c16slavep0N = new String(new char[17]); // delta pixel slave0N
		  String c16slavelN0 = new String(new char[17]); // delta line slaveN0
		  String c16slavepN0 = new String(new char[17]); // delta pixel slaveN0
		  String c16slavelNN = new String(new char[17]); // delta line slaveNN
		  String c16slavepNN = new String(new char[17]); // delta pixel slaveNN
		  readres(c16slavel00,sizeof(c16slavel00), input_general.i_resfile,"Deltaline_slave00_poly:", 0);
		  readres(c16slavep00,sizeof(c16slavep00), input_general.i_resfile,"Deltapixel_slave00_poly:", 0);
		  readres(c16slavel0N,sizeof(c16slavel0N), input_general.i_resfile,"Deltaline_slave0N_poly:", 0);
		  readres(c16slavep0N,sizeof(c16slavep0N), input_general.i_resfile,"Deltapixel_slave0N_poly:", 0);
		  readres(c16slavelN0,sizeof(c16slavelN0), input_general.i_resfile,"Deltaline_slaveN0_poly:", 0);
		  readres(c16slavepN0,sizeof(c16slavepN0), input_general.i_resfile,"Deltapixel_slaveN0_poly:", 0);
		  readres(c16slavelNN,sizeof(c16slavelNN), input_general.i_resfile,"Deltaline_slaveNN_poly:", 0);
		  readres(c16slavepNN,sizeof(c16slavepNN), input_general.i_resfile,"Deltapixel_slaveNN_poly:", 0);

		  slave.add_offsetl00(Double.parseDouble(c16slavel00));
		  slave.add_offsetp00(Double.parseDouble(c16slavep00));
		  slave.add_offsetl0N(Double.parseDouble(c16slavel0N));
		  slave.add_offsetp0N(Double.parseDouble(c16slavep0N));
		  slave.add_offsetlN0(Double.parseDouble(c16slavelN0));
		  slave.add_offsetpN0(Double.parseDouble(c16slavepN0));
		  slave.add_offsetlNN(Double.parseDouble(c16slavelNN));
		  slave.add_offsetpNN(Double.parseDouble(c16slavepNN));
		  INFO << Double.parseDouble(c16slavel00);
		  INFO.print();
		  INFO << Double.parseDouble(c16slavep00);
		  INFO.print();
		  INFO << Double.parseDouble(c16slavel0N);
		  INFO.print();
		  INFO << Double.parseDouble(c16slavep0N);
		  INFO.print();
		  INFO << Double.parseDouble(c16slavelN0);
		  INFO.print();
		  INFO << Double.parseDouble(c16slavepN0);
		  INFO.print();
		  INFO << Double.parseDouble(c16slavelNN);
		  INFO.print();
		  INFO << Double.parseDouble(c16slavepNN);
		  } // [MA] end of quick soln.
		// ______ end added by FvL ______
		}


	  // ====== RESAMPLING OF SLAVE IMAGE ======
	  if (input_general.process[pr_s_resample])
		{
		PROGRESS.print("Start RESAMPLE.");
		alreadyprocessed[pr_s_resample] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing resampling of slave image. (might take some time.)";
		  getanswer();
		  }

		 if (input_s_resample.dbow_geo.pixhi != 0)
			{
			  PROGRESS.print("Computing RS_DBOW based on center lat/lon/height/width.");
			  real8 centerline;
			  real8 centerpix;
			  real8 centerphi = deg2rad(input_s_resample.dbow_geo.linelo/1e6-360.0);
			  real8 centerlambda = deg2rad(input_s_resample.dbow_geo.linehi/1e6-360.0);

			  real8 centerheight = input_general.terrain_height;
			  DEBUG.print("Using height of HEIGHT card");

			  // compute phi/lambda/h-->x,y,z
			  cn centerpos = input_ellips.ell2xyz(centerphi,centerlambda,centerheight);

			  DEBUG << "Converted RS_DBOW_GEO phi,lambda,hei --> x,y,z:" << centerphi << ", " << centerlambda << ", " << centerheight << " --> " << centerpos.x << " " << centerpos.y << " " << centerpos.z;
			  DEBUG.print();

	//          
	//            DEBUG << "Center master SLC: x,y,z: "
	//            << master.approxcentreoriginal.x << " " 
	//            << master.approxcentreoriginal.y << " "
	//            << master.approxcentreoriginal.z;
	//            DEBUG.print();
	//          

			  if (Math.abs(master.approxcentreoriginal.x - centerpos.x) > 60000.0) // 50km
				WARNING.print("RS_DBOW_GEO: coordinates seem to be outside SLC area? (X)");
			  if (Math.abs(master.approxcentreoriginal.y - centerpos.y) > 60000.0) // 50km
				WARNING.print("RS_DBOW_GEO: coordinates seem to be outside SLC area? (Y)");
			  if (Math.abs(master.approxcentreoriginal.z - centerpos.z) > 60000.0) // 50km
				WARNING.print("RS_DBOW_GEO: coordinates seem to be outside SLC area? (Z)");

			  RefObject<real8> TempRefObject23 = new RefObject<real8>(centerline);
			  RefObject<real8> TempRefObject24 = new RefObject<real8>(centerpix);
			  RefObject<orbit> TempRefObject25 = new RefObject<orbit>(masterorbit);
			  xyz2lp(TempRefObject23, TempRefObject24, master, TempRefObject25, centerpos, 10, 1e-3); // MAXITERATIONS, CONVERGENCE_TIME
			  centerline = TempRefObject23.argvalue;
			  centerpix = TempRefObject24.argvalue;
			  masterorbit = TempRefObject25.argvalue;

			  DEBUG << "RS_DBOW_GEO: center line/pix: " << centerline << " " << centerpix;
			  DEBUG.print();

			  int32 l0 = int32(centerline+0.5) - input_s_resample.dbow_geo.pixlo/2;
			  int32 lN = l0 + input_s_resample.dbow_geo.pixlo - 1;
			  int32 p0 = int32(centerpix+0.5) - input_s_resample.dbow_geo.pixhi/2;
			  int32 pN = p0 + input_s_resample.dbow_geo.pixhi - 1;

			  if (l0 < int32(master.currentwindow.linelo))
				  l0 = master.currentwindow.linelo;
			  if (lN > int32(master.currentwindow.linehi))
				  lN = master.currentwindow.linehi;
			  if (p0 < int32(master.currentwindow.pixlo))
				  p0 = master.currentwindow.pixlo;
			  if (pN > int32(master.currentwindow.pixhi))
				  pN = master.currentwindow.pixhi;
			  INFO << "RS_DBOW from GEO: " << l0 << " " << lN << " " << p0 << " " << pN;
			  INFO.print();

			  // ___ Simply fill dbow with it and it will be done correctly! ___
			  input_s_resample.dbow.linelo = uint(l0);
			  input_s_resample.dbow.linehi = uint(lN);
			  input_s_resample.dbow.pixlo = uint(p0);
			  input_s_resample.dbow.pixhi = uint(pN);
			}

		// ______ (interf.win contains approx. offset) ______
		resample(input_general, input_s_resample, master, slave, coeff_cpmL, coeff_cpmP, alreadyprocessed[pr_i_demassist]);
		// ______ Update log files ______
		updatefile("scratchlogresample", input_general.logfile);
		updatefile("scratchresresample", input_general.s_resfile);
		PROGRESS.print("Finished RESAMPLE.");
		DEBUG.print("Time spent for resampling of slave image:");
		printcpu();
		}

	  // ______ Update info of slave, size interferogram ______
	  if (alreadyprocessed[pr_s_resample])
		{
		INFO.print("");
		INFO.print("slave: latest known processing stage: resample");
		String SECTIONID = new String(new char[ONE27]);
		SECTIONID = "*_Start_";
		SECTIONID += processcontrol[pr_s_resample];
		slave.updateslcimage(input_general.s_resfile,SECTIONID);
		// (this is in master coord. system, as it should be) (BK 26-1-00)
		// in some routines this might be used, should be changed
		// wins are in same system
		interferogram.win = getoverlap(master.currentwindow, slave.currentwindow);
		}
	  // ___ generate a preview of amplitude if requested ______
	  if (input_general.process[pr_s_resample])
		{
		PROGRESS.print("calling preview for resampled slave");
		if (master.sensor!=SLC_RSAT && master.sensor!=SLC_ALOS && master.sensor!=SLC_TSX)
		  preview(input_general.preview, slave.currentwindow.pixels(), slave.formatflag, slave.file, "slave_rs_mag.ras", "-e 0.5 -s 1.0 -q mag -o sunraster -b -c gray -M 2/10");
		else if (master.sensor == SLC_ALOS || master.sensor == SLC_TSX)
		  preview(input_general.preview, slave.currentwindow.pixels(), slave.formatflag, slave.file, "slave_rs_mag.ras", "-e 0.5 -s 1.0 -q mag -o sunraster -b -c gray -M 10/10"); // [MA] higher mlook to reduce filesize
		else //multilooking depends on beam
		  preview(input_general.preview, slave.currentwindow.pixels(), slave.formatflag, slave.file, "slave_rs_mag.ras", "-e 0.5 -s 1.0 -q mag -o sunraster -b -c gray -M 5/6");
		}



	// ====== FILTER RANGE: MASTER&SLAVE adaptive ======
	// ______ Only one step for master and slave, after resampling. ______
	// ______ Method adaptive only, for method based on orbits see after coarse. ______
	  if (input_general.process[pr_m_filtrange] && input_ms_filtrange.method==rf_adaptive)
		{
		PROGRESS.print("Start FILTRANGE (adaptive).");
		alreadyprocessed[pr_m_filtrange] =1; // update alreadypr.
		alreadyprocessed[pr_s_filtrange] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nAdaptive range filtering master and slave.";
		  getanswer();
		  }

		// ______ Filter routine ______
		rangefilter(input_general, master, slave, interferogram, input_ms_filtrange);

		// ______ Update log files ______
		updatefile("scratchlogfiltrange", input_general.logfile);
		updatefile("scratchresMfiltrange", input_general.m_resfile);
		updatefile("scratchresSfiltrange", input_general.s_resfile);

		PROGRESS.print("Finished FILTRANGE.");
		DEBUG.print("Time spent for range filtering:");
		printcpu();
		}

	// ______ Update info of slave, size interferogram ______
	  if (alreadyprocessed[pr_m_filtrange])
		{
		String c10rfmethod = new String(new char[11]); // method
		readres(c10rfmethod,sizeof(c10rfmethod),input_general.m_resfile, "Method_rangefilt:", 0);
		if (!strcmp(c10rfmethod,"adaptive"))
		  {
		  INFO.print("");
		  INFO.print("master: latest known processing stage: adaptive range filtered");
		  String SECTIONID = new String(new char[ONE27]);
		  SECTIONID = "*_Start_";
		  SECTIONID += processcontrol[pr_m_filtrange];
		  master.updateslcimage(input_general.m_resfile,SECTIONID);
		  // ______ cut out overlay master slave exact ______
		  interferogram.win = master.currentwindow;
		  }
		}

	// ______ Update info of slave, size interferogram ______
	  if (alreadyprocessed[pr_s_filtrange]) // same time as _m_
		{
		String c10rfmethod = new String(new char[11]); // method
		readres(c10rfmethod,sizeof(c10rfmethod),input_general.s_resfile, "Method_rangefilt:", 0);
		if (!strcmp(c10rfmethod,"adaptive"))
		  {
		  INFO.print("");
		  INFO.print("slave: latest known processing stage: adaptive range filtered");
		  String SECTIONID = new String(new char[ONE27]);
		  SECTIONID = "*_Start_";
		  SECTIONID += processcontrol[pr_s_filtrange];
		  slave.updateslcimage(input_general.s_resfile,SECTIONID);
		  // ______ cut out overlay master slave exact, coord in master system ______
		  if (slave.currentwindow.linelo != master.currentwindow.linelo || slave.currentwindow.linehi != master.currentwindow.linehi || slave.currentwindow.pixlo != master.currentwindow.pixlo || slave.currentwindow.pixhi != master.currentwindow.pixhi)
			{
			WARNING.print("master/slave should overlap exactly after adaptive range filtering?");
			master.showdata();
			slave.showdata();
			}
		  }
		}






	// ______ Fill matrix coefficients flatearth ______
	// ______ because might be requested to b subtracted ______
	  if (alreadyprocessed[pr_i_comprefpha])
		{
		INFO.print("");
		INFO.print("ifg: latest known processing stage: comprefpha");
		String c10offL = new String(new char[11]);
		readres(c10offL,sizeof(c10offL),input_general.i_resfile,"Degree_flat:");
		int32 degreeflat = Integer.parseInt(c10offL);
		coeff_flat = readcoeff(input_general.i_resfile, "Estimated_coefficients_flatearth:", Ncoeffs(degreeflat));
		}


	// ====== COMPUTATION OF (COMPLEX) INTERFEROGRAM ======
	  if (input_general.process[pr_i_interfero])
		{
		PROGRESS.print("Start INTERFERO.");
		alreadyprocessed[pr_i_interfero] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing computation of interferogram.";
		  getanswer();
		  }
		// ______ Select method ______
		// ______ slave.currentwin is in master system (after resample) and is ______
		// ______ smaller than master.currentwin ______
		if (input_i_interfero.method==int_oldmethod)
		  compinterfero(master, slave, input_general, input_i_interfero, coeff_flat);
		else if (input_i_interfero.method==int_oversample)
		  {
		  PRINT_ERROR("NOT IMPLEMENTED IN THIS VERSION.")
		  throw(input_error); // exit
		  }
		else
		  {
		  PRINT_ERROR("unknown method interfero.")
		  throw(input_error); // exit
		  }

		// ______ Update log files ______
		updatefile("scratchloginterfero", input_general.logfile);
		updatefile("scratchresinterfero", input_general.i_resfile);

		PROGRESS.print("Finished INTERFERO.");
		DEBUG.print("Time spent for computation of interferogram:");
		printcpu();
		}

	// ______ Fill input struct for filenames ______
	  if (alreadyprocessed[pr_i_interfero]) // update productinfo
		{
		INFO.print("");
		INFO.print("ifg: latest known processing stage: interfero");
		String SECTIONID = new String(new char[ONE27]);
		SECTIONID = "*_Start_";
		SECTIONID += processcontrol[pr_i_interfero];
		interferogram.fillproductinfo(input_general.i_resfile,SECTIONID);
		}

	// ______ Generate preview if requested ______
	  if (input_general.process[pr_i_interfero])
		{
		// ___ floor is not required, we do int division (BK, april 2003) ___
		// int32(floor(interferogram.win.pixels()/interferogram.multilookP)),
		PROGRESS.print("calling preview for complex interferogram");
		preview(input_general.preview, int32((interferogram.win.pixels()/interferogram.multilookP)), interferogram.formatflag, interferogram.file, "interferogram_mag.ras", "-e 0.3 -s 1.0 -q mag -o sunraster -b -c gray -M 2/2");
		preview(input_general.preview, int32((interferogram.win.pixels()/interferogram.multilookP)), interferogram.formatflag, interferogram.file, "interferogram_pha.ras", "-q phase -o sunraster -b -c jet -M 2/2");
		preview(input_general.preview, int32((interferogram.win.pixels()/interferogram.multilookP)), interferogram.formatflag, interferogram.file, "interferogram_mix.ras", "-e 0.3 -s 1.2 -q mixed -o sunraster -b -c hsv -M 2/2");
		}



	// ====== COMPUTATION OF REFERENCE PHASE (FLATEARTH) ======
	  if (input_general.process[pr_i_comprefpha])
		{
		PROGRESS.print("Start COMPREFPHA.");
		alreadyprocessed[pr_i_comprefpha] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing computation of reference interferogram.";
		  getanswer();
		  }

		if (input_i_comprefpha.method == fe_porbits)
		  flatearth(input_i_comprefpha, input_ellips, master, slave, interferogram, masterorbit, slaveorbit);

		else if (input_i_comprefpha.method == fe_method2)
		  {
		  PRINT_ERROR("NOT IMPLEMENTED IN THIS VERSION.")
		  throw(input_error); // exit
		  }

		// ______ Update log files ______
		updatefile("scratchlogflat", input_general.logfile);
		updatefile("scratchresflat", input_general.i_resfile);

		PROGRESS.print("Finished COMPREFPHA.");
		DEBUG.print("Time spent for computation of reference interferogram 'flat earth':");
		printcpu();
		}

	// ______Fill matrix coefficients flatearth______
	  if (alreadyprocessed[pr_i_comprefpha])
		{
		INFO.print("");
		INFO.print("ifg: latest known processing stage: comprefpha");
		String c10offL = new String(new char[11]);
		readres(c10offL,sizeof(c10offL),input_general.i_resfile,"Degree_flat:");
		int32 degreeflat = Integer.parseInt(c10offL);

		coeff_flat = readcoeff(input_general.i_resfile, "Estimated_coefficients_flatearth:", Ncoeffs(degreeflat));

		//____ added by FvL ____________
		readres(c10offL,sizeof(c10offL),input_general.i_resfile,"Degree_h2ph:");
		int32 degreeh2ph = Integer.parseInt(c10offL);

		coeff_h2ph = readcoeff(input_general.i_resfile, "Estimated_coefficients_h2ph:", Ncoeffs(degreeh2ph));
		//_____ end added by FvL ________ 

		}



	// ====== SUBTRACT REFERENCE PHASE (FLATEARTH) ======
	  if (input_general.process[pr_i_subtrrefpha])
		{
		PROGRESS.print("Start SUBTRREFPHA.");
		alreadyprocessed[pr_i_subtrrefpha] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing subtraction of reference phase from interferogram.\n";
		  getanswer();
		  }

		// ______ Select appropriate module ______
		switch (input_i_subtrrefpha.method)
		  {
		  case srp_polynomial:
			subtrrefpha(master,interferogram, input_general,input_i_subtrrefpha, coeff_flat,coeff_h2ph); // coeff_h2ph added by FvL
			break;
		  case srp_exact:
			subtrrefpha(input_ellips, master, slave, interferogram, input_general, input_i_subtrrefpha, masterorbit, slaveorbit);
			break;
		  default:
			PRINT_ERROR("PANIC: not possible, bert.")
			throw(unhandled_case_error); // exit
		  }

		// ______ Update log files ______
		updatefile("scratchlogsubtrrefpha", input_general.logfile);
		updatefile("scratchressubtrrefpha", input_general.i_resfile);

		PROGRESS.print("Finished SUBTRREFPHA.");
		DEBUG.print("Time spent for subtraction of reference phase:");
		printcpu();
		}

	// ______ Fill input struct for filenames ______
	  if (alreadyprocessed[pr_i_subtrrefpha]) // update productinfo
		{
		INFO.print("");
		INFO.print("ifg: latest known processing stage: subtrrefpha");
		String SECTIONID = new String(new char[ONE27]);
		SECTIONID = "*_Start_";
		SECTIONID += processcontrol[pr_i_subtrrefpha];
		interferogram.fillproductinfo(input_general.i_resfile,SECTIONID);
		}

	// ______ Generate preview if requested ______
	  if (input_general.process[pr_i_subtrrefpha])
		{
		PROGRESS.print("calling preview for refpha corrected image");
		// ____ if refphase dumped, tmp change filename ___
		//if (!strcmp(interferogram.file,"NO_OUTPUT_ONLY_DUMPING_REF_PHA")
		if (input_i_subtrrefpha.dumponlyrefpha==true)
		  interferogram.file = input_i_subtrrefpha.forefpha;
		preview(input_general.preview, int32((interferogram.win.pixels()/interferogram.multilookP)), interferogram.formatflag, interferogram.file, "interferogram_srp_mag.ras", "-e 0.3 -s 1.0 -q mag -o sunraster -b -c gray -M 2/2");
		preview(input_general.preview, int32((interferogram.win.pixels()/interferogram.multilookP)), interferogram.formatflag, interferogram.file, "interferogram_srp_pha.ras", "-q phase -o sunraster -b -c jet -M 2/2");
		preview(input_general.preview, int32((interferogram.win.pixels()/interferogram.multilookP)), interferogram.formatflag, interferogram.file, "interferogram_srp_mix.ras", "-e 0.3 -s 1.2 -q mixed -o sunraster -b -c cool -M 2/2");
		if (input_i_subtrrefpha.dumponlyrefpha==true)
		  {
			  INFO.print("exiting, only dumped refpha");
			  exit(1);
		  }
		}



	// ====== COMPUTE REFERENCE PHASE (DEM) ======
	  if (input_general.process[pr_i_comprefdem])
		{
		PROGRESS.print("Start COMPREFDEM.");
		alreadyprocessed[pr_i_comprefdem] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing computation of reference phase from DEM.\n";
		  getanswer();
		  }
		// if (input_i_comprefdem.method == gtopo30)
		  radarcodedem(input_general, input_ellips, input_i_comprefdem, master, slave, interferogram, masterorbit, slaveorbit);

		// ______ Update log files ______
		updatefile("scratchlogcomprefdem", input_general.logfile);
		updatefile("scratchrescomprefdem", input_general.i_resfile);

		//WARNING("Not tested properly yet.");
		PROGRESS.print("Finished COMPREFDEM.");
		DEBUG.print("Time spent for computation of reference DEM:");
		printcpu();
		}

	// ______ Fill input struct for ______
	  if (alreadyprocessed[pr_i_comprefdem]) // update slcimage
		{
		INFO.print("");
		INFO.print("ifg: latest known processing stage: comprefdem");
		String SECTIONID = new String(new char[ONE27]);
		SECTIONID = "*_Start_";
		SECTIONID += processcontrol[pr_i_comprefdem];
		radarcodedrefdem.fillproductinfo(input_general.i_resfile,SECTIONID);
		}
	  // ___ generate preview if requested and just processed ___
	  if (input_general.process[pr_i_comprefdem])
		{
		// ______ See if we can trick cpxfiddle in believing this is cpx ______
		// ______ Should also work for real4 ______ // BK 16-Jun-2003 ______
		PROGRESS.print("calling preview for comprefdem (phase)");
		preview(input_general.preview, int32((radarcodedrefdem.win.pixels()/radarcodedrefdem.multilookP)), radarcodedrefdem.formatflag, radarcodedrefdem.file, "comprefdem_phase.ras", "-q normal -o sunraster -b -c jet -M 2/2");
		if (specified(input_i_comprefdem.forefdemhei))
		  {
		  PROGRESS.print("calling preview for comprefdem (height)");
		  preview(input_general.preview, int32((radarcodedrefdem.win.pixels()/radarcodedrefdem.multilookP)), FORMATR4, input_i_comprefdem.forefdemhei, "comprefdem_height.ras", "-q normal -o sunraster -b -c jet -M 2/2");
		  }
		}


	// ====== SUBTRACT REFERENCE PHASE (DEM) ======
	  if (input_general.process[pr_i_subtrrefdem])
		{
		PROGRESS.print("Start SUBTRREFDEM.");
		alreadyprocessed[pr_i_subtrrefdem] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing subtraction of reference DEM from interferogram.\n";
		  getanswer();
		  }
		subtrrefdem(interferogram,radarcodedrefdem,input_general,input_i_subtrrefdem);

		// ______ Update log files ______
		updatefile("scratchlogsubtrrefdem", input_general.logfile);
		updatefile("scratchressubtrrefdem", input_general.i_resfile);

		PROGRESS.print("Finished SUBTRREFDEM.");
		DEBUG.print("Time spent for subtraction of reference DEM:");
		printcpu();
		}

	// ______ Fill input struct for filenames ______
	  if (alreadyprocessed[pr_i_subtrrefdem]) // update slcimage
		{
		INFO.print("");
		INFO.print("ifg: latest known processing stage: subtrrefdem");
		//WARNING("check this ..., nothing changed except filenaem?");
		String SECTIONID = new String(new char[ONE27]);
		SECTIONID = "*_Start_";
		SECTIONID += processcontrol[pr_i_subtrrefdem];
		interferogram.fillproductinfo(input_general.i_resfile,SECTIONID);
		}

	// ______ Generate preview if requested ______
	  if (input_general.process[pr_i_subtrrefdem])
		{
		PROGRESS.print("calling preview for refdem corrected interferogram");
		preview(input_general.preview, int32((interferogram.win.pixels()/interferogram.multilookP)), interferogram.formatflag, interferogram.file, "interferogram_srd_mag.ras", "-e 0.3 -s 1.0 -q mag -o sunraster -b -c gray -M 2/2");
		preview(input_general.preview, int32((interferogram.win.pixels()/interferogram.multilookP)), interferogram.formatflag, interferogram.file, "interferogram_srd_pha.ras", "-q phase -o sunraster -b -c jet -M 2/2");
		preview(input_general.preview, int32((interferogram.win.pixels()/interferogram.multilookP)), interferogram.formatflag, interferogram.file, "interferogram_srd_mix.ras", "-e 0.3 -s 1.2 -q mixed -o sunraster -b -c hot -M 2/2");
		}


	// [don] the computation of the coherence image has been moved 
	//       after the computation of the reference phase (DEM) in 
	//       order to include the radarcodedrefdem_phase_correction 

	// ====== COMPUTATION OF COHERENCE IMAGE ======
	// (use also other external computed phase in future!) 
	  if (input_general.process[pr_i_coherence])
		{
		PROGRESS.print("Start COHERENCE.");
		alreadyprocessed[pr_i_coherence] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing computation of coherence.";
		  getanswer();
		  }

		// ______ Select method ______
		switch (input_i_coherence.method)
		  {
		  case coh_oldmethod:
			  // use flatearth, compute interferogram again
			  // compcoherence(master,slave,interferogram,
			  compcoherence(master, slave, input_general, input_i_coherence, coeff_flat);
			break;
		  case coh_newmethod:
							 // start_added_by_don
							 if (alreadyprocessed[pr_i_comprefdem])
				// use flatearth and radarcodedrefdem, compute interferogram again
				compcoherence(master, slave, radarcodedrefdem, input_general, input_i_coherence, coeff_flat);
							 else
				{
				  ERROR << "Requested step: " << pr_i_coherence << " (" << processcontrol[pr_i_coherence] << ") seems impossible, because step " << pr_i_comprefdem << " (" << processcontrol[pr_i_comprefdem] << ") is not in resultfile.";
				  PRINT_ERROR(ERROR.get_str())
				  throw(input_error); // exit
				}
							 // end_added_by_don
			break;
		  default:
			PRINT_ERROR("unknown method coherence.")
			throw(unhandled_case_error); // exit
		  }

		// ______ Update log files ______
		updatefile("scratchlogcoherence", input_general.logfile);
		updatefile("scratchrescoherence", input_general.i_resfile);

		PROGRESS.print("Finished COHERENCE.");
		DEBUG.print("Time spent for computation of coherence:");
		printcpu();
		}

	// ______ Fill input struct for filenames ______
	  if (alreadyprocessed[pr_i_coherence]) // update slcimage
		{
		INFO.print("");
		INFO.print("ifg: latest known processing stage: coherence");
		String SECTIONID = new String(new char[ONE27]);
		SECTIONID = "*_Start_";
		SECTIONID += processcontrol[pr_i_coherence];
		coherence.fillproductinfo(input_general.i_resfile,SECTIONID);
		}
	  // ___ generate preview if requested and just processed ___
	  if (input_general.process[pr_i_coherence])
		{
		// ______ See if we can trick cpxfiddle in believing this is cpx ______
		// ______ Should also work for real4 ______ // BK 16-Jun-2003 ______
		PROGRESS.print("calling preview for coherence");
		preview(input_general.preview, int32((coherence.win.pixels()/coherence.multilookP)), coherence.formatflag, coherence.file, "coherence.ras", "-q normal -o sunraster -b -c gray -M 2/2 -r 0.0/1.0 "); // coherence range 0 to 1.
		}


	// ====== PHASE FILTERING ======
	  if (input_general.process[pr_i_filtphase])
		{
		PROGRESS.print("Start FILTPHASE.");
		alreadyprocessed[pr_i_filtphase] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing phase filtering.";
		  getanswer();
		  }

		switch (input_i_filtphase.method)
		  {
		  case fp_goldstein:
			phasefilter(input_general, interferogram, input_i_filtphase);
			break;
		  case fp_spatialconv:
			spatialphasefilt(input_general, interferogram, input_i_filtphase);
			break;
		  case fp_spectral:
			phasefilterspectral(input_general, interferogram, input_i_filtphase);
			break;
		  default:
			PRINT_ERROR("PANIC: not possible.")
			throw(unhandled_case_error); // exit
		  }

	// ______ Update log files ______
		updatefile("scratchlogfiltphase", input_general.logfile);
		updatefile("scratchresfiltphase", input_general.i_resfile);

		PROGRESS.print("Finished FILTPHASE.");
		DEBUG.print("Time spent for phase filtering:");
		printcpu();
		}

	// ______ Fill info struct for filenames ______
	  if (alreadyprocessed[pr_i_filtphase]) // update productinfo
		{
		INFO.print("");
		INFO.print("ifg: latest known processing stage: filtphase");
		String SECTIONID = new String(new char[ONE27]);
		SECTIONID = "*_Start_";
		SECTIONID += processcontrol[pr_i_filtphase];
		interferogram.fillproductinfo(input_general.i_resfile,SECTIONID);
		}

	// ______ Generate preview if requested ______
	  if (input_general.process[pr_i_filtphase])
		{
		PROGRESS.print("calling preview for phase filtered corrected interferogram");
		preview(input_general.preview, int32((interferogram.win.pixels()/interferogram.multilookP)), interferogram.formatflag, interferogram.file, "interferogram_filt_mag.ras", "-e 0.3 -s 1.0 -q mag -o sunraster -b -c gray -M 2/2");
		preview(input_general.preview, int32((interferogram.win.pixels()/interferogram.multilookP)), interferogram.formatflag, interferogram.file, "interferogram_filt_pha.ras", "-q phase -o sunraster -b -c jet -M 2/2");
		preview(input_general.preview, int32((interferogram.win.pixels()/interferogram.multilookP)), interferogram.formatflag, interferogram.file, "interferogram_filt_mix.ras", "-e 0.3 -s 1.2 -q mixed -o sunraster -b -c hsv -M 2/2");
		}



	// ====== DIFFERENTIAL 3 PASS, DEFO PROCESSING ======
	  if (input_general.process[pr_i_dinsar])
		{
		PROGRESS.print("Start DINSAR.");
		alreadyprocessed[pr_i_dinsar] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing differential interferometry.";
		  getanswer();
		  }
		dinsar(input_general,input_i_dinsar,input_ellips, master,masterorbit, slave,slaveorbit,interferogram); // (defopair)

		// ______ Update log files ______
		updatefile("scratchlogdinsar", input_general.logfile);
		updatefile("scratchresdinsar", input_general.i_resfile);

		PROGRESS.print("Finished DINSAR.");
		DEBUG.print("Time spent for 3 pass differential interferometry:");
		printcpu();
		}

	  // ______ Fill info struct for filenames ______
	  // ______ use interferogram because you want to unwrap it? ______
	  if (alreadyprocessed[pr_i_dinsar]) // update productinfo
		{
		INFO.print("");
		INFO.print("ifg: latest known processing stage: dinsar");
		String SECTIONID = new String(new char[ONE27]);
		SECTIONID = "*_Start_";
		SECTIONID += processcontrol[pr_i_dinsar];
		interferogram.fillproductinfo(input_general.i_resfile,SECTIONID);
		}

	// ______ Generate preview if requested ______
	  if (input_general.process[pr_i_dinsar])
		{
		PROGRESS.print("calling preview for differential interferogram");
		preview(input_general.preview, int32((interferogram.win.pixels()/interferogram.multilookP)), interferogram.formatflag, interferogram.file, "interferogram_dinsar_mag.ras", "-e 0.3 -s 1.0 -q mag -o sunraster -b -c gray -M 2/2");
		preview(input_general.preview, int32((interferogram.win.pixels()/interferogram.multilookP)), interferogram.formatflag, interferogram.file, "interferogram_dinsar_pha.ras", "-q phase -o sunraster -b -c jet -M 2/2");
		preview(input_general.preview, int32((interferogram.win.pixels()/interferogram.multilookP)), interferogram.formatflag, interferogram.file, "interferogram_dinsar_mix.ras", "-e 0.3 -s 1.2 -q mixed -o sunraster -b -c jet -M 2/2");
		}




	// ====== UNWRAPPING INTERFEROGRAM ======
	  unwrappedinterf = interferogram;
	  if (input_general.process[pr_i_unwrap])
		{
		PROGRESS.print("Start UNWRAP.");
		alreadyprocessed[pr_i_unwrap] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing unwrapping of interferogram.";
		  getanswer();
		  }

		switch (input_i_unwrap.method)
		  {
		  // ___ Branch cut ___
		  case uw_method1:
			unwraptreeframon(input_general, input_i_unwrap, interferogram);
			break;

		  // ___ Statistic minimal cost flow ___
		  case uw_method2:
			snaphu_unwrap(input_general, input_i_unwrap, interferogram, master, slave, masterorbit, slaveorbit, input_ellips);
			if (specified(input_i_unwrap.snaphu_log))
			  updatefile(input_i_unwrap.snaphu_log, input_general.logfile);
			break;

		  // ___ Minimal cost flow algorithm ___
		  case uw_method3:
			PRINT_ERROR("NOT IMPLEMENTED IN THIS VERSION.")
			throw(unhandled_case_error); // exit
			break;

		  default:
			PRINT_ERROR("PANIC: NOT POSSIBLE.")
			throw(unhandled_case_error); // exit
		  }

	//    // ______ Apply bias to unwrapped interferogram, based on tiepoint ______
	//    // ______ This simply loads the created unwrapped ifg; ______
	//    // ______ computes the tiepoint, and saves it to the same file. ______
	//    if (abs(input_i_unwrap.tiepoint.x) < 0.001 &&
	//        abs(input_i_unwrap.tiepoint.y) < 0.001 &&
	//        abs(input_i_unwrap.tiepoint.z) < 0.001)
	//      {
	//      PROGRESS.print("Computating integration constant based on tiepoint");
	//      add_bias_unwrap(
	//              input_general, 
	//              input_i_unwrap,
	//              interferogram,
	//              master, slave,
	//              masterorbit, slaveorbit,
	//              input_ellips);
	//++++ also add in unwrap.h
	//      }

		// ______ Update log files ______
		updatefile("scratchlogunwrap", input_general.logfile);
		updatefile("scratchresunwrap", input_general.i_resfile);

		PROGRESS.print("Finished UNWRAP.");
		DEBUG.print("Time spent for unwrapping:");
		printcpu();
		}

	  // ______ Fill info struct for filenames ______
	  // ______ 2b sure, initialize unwrappedinterf with interferogram ______
	  unwrappedinterf = interferogram;
	  if (alreadyprocessed[pr_i_unwrap]) // update productinfo
		{
		INFO.print("");
		INFO.print("ifg: latest known processing stage: unwrap");
		String SECTIONID = new String(new char[ONE27]);
		SECTIONID = "*_Start_";
		SECTIONID += processcontrol[pr_i_unwrap];
		unwrappedinterf.fillproductinfo(input_general.i_resfile,SECTIONID);
		}
	  // ___ generate preview if requested and just now computed____
	  if (input_general.process[pr_i_unwrap])
		{
		// ______ Generate preview if requested and oformat was hgt?______
		// ______ See if we can trick cpxfiddle in believing this is cpx ______
		// ______ Should also work for real4 ______
		// BK 13-Apr-2003
		PROGRESS.print("calling preview for unwrapped interferogram");
		preview(input_general.preview, int32((unwrappedinterf.win.pixels()/unwrappedinterf.multilookP)), unwrappedinterf.formatflag, unwrappedinterf.file, "unwrapped_interferogram.ras", "-q normal -o sunraster -b -c jet -M 2/2"); // -V
		}


	// ====== SLANT TO HEIGHT CONVERSION ======
	  if (input_general.process[pr_i_slant2h])
		{
		PROGRESS.print("Start SLANT2H.");
		alreadyprocessed[pr_i_slant2h] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing slant to height conversion.";
		  getanswer();
		  }
		if (!strcmp(input_i_slant2h.fohei,unwrappedinterf.file))
		  {
		  PRINT_ERROR("slant2h: same filename input/output")
		  throw(input_error); // exit
		  }

		switch (input_i_slant2h.method)
		  {
		  case s2h_schwabisch:
			RefObject<orbit> TempRefObject26 = new RefObject<orbit>(masterorbit);
			RefObject<orbit> TempRefObject27 = new RefObject<orbit>(slaveorbit);
			slant2hschwabisch(input_general, input_i_slant2h, input_ellips, master, slave, unwrappedinterf, TempRefObject26, TempRefObject27);
			masterorbit = TempRefObject26.argvalue;
			slaveorbit = TempRefObject27.argvalue;
			break;

		  case s2h_rodriguez:
			RefObject<orbit> TempRefObject28 = new RefObject<orbit>(masterorbit);
			RefObject<orbit> TempRefObject29 = new RefObject<orbit>(slaveorbit);
			slant2hrodriguez(input_general, input_i_slant2h, input_ellips, master, slave, unwrappedinterf, coeff_flat, TempRefObject28, TempRefObject29, baseline);
			masterorbit = TempRefObject28.argvalue;
			slaveorbit = TempRefObject29.argvalue;
			break;

		  case s2h_ambiguity:
			RefObject<orbit> TempRefObject30 = new RefObject<orbit>(masterorbit);
			RefObject<orbit> TempRefObject31 = new RefObject<orbit>(slaveorbit);
			slant2hambiguity(input_general, input_i_slant2h, input_ellips, master, slave, unwrappedinterf, TempRefObject30, TempRefObject31, baseline);
			masterorbit = TempRefObject30.argvalue;
			slaveorbit = TempRefObject31.argvalue;
			// ______ This method includes geocoding, but probably contains trend ______
			alreadyprocessed[pr_i_geocoding] =1; // update alreadypr.
			break;

		  default:
			PRINT_ERROR("slant2h: unknown method.")
			throw(unhandled_case_error); // exit
		  }

		// ______ Update log files ______
		updatefile("scratchlogslant2h", input_general.logfile);
		updatefile("scratchresslant2h", input_general.i_resfile);

		PROGRESS.print("Finished SLANT2H.");
		DEBUG.print("Time spent for slant-to-height:");
		printcpu();
		}

	  // ______ Fill info struct for filenames ______
	  if (alreadyprocessed[pr_i_slant2h]) // update productinfo
		{
		INFO.print("");
		INFO.print("ifg: latest known processing stage: slant2h");
		String SECTIONID = new String(new char[ONE27]);
		SECTIONID = "*_Start_";
		SECTIONID += processcontrol[pr_i_slant2h];
		heightinradarsystem.fillproductinfo(input_general.i_resfile,SECTIONID);
		}
	  // ___ generate preview if requested and just processed ___
	  if (input_general.process[pr_i_slant2h])
		{
		// ______ Generate preview if requested and oformat was hgt?______
		// ______ See if we can trick cpxfiddle in believing this is cpx ______
		// ______ Should also work for real4 ______
		// ______ note that there is not point in preview of these heights? _____
		// BK 13-Apr-2003
		PROGRESS.print("calling preview for slant to height");
		preview(input_general.preview, int32((heightinradarsystem.win.pixels()/heightinradarsystem.multilookP)), heightinradarsystem.formatflag, heightinradarsystem.file, "height_in_radarcoordinates.ras", "-q normal -o sunraster -b -c jet -M 2/2"); // -V
		}



	// ====== GEOCODING ======
	  if (input_general.process[pr_i_geocoding])
		{
		PROGRESS.print("Start GEOCODE.");
		alreadyprocessed[pr_i_geocoding] =1; // update alreadypr.
		if (input_general.interactive)
		  {
		  cerr << "\nProcessing geocoding.";
		  getanswer();
		  }
		RefObject<orbit> TempRefObject32 = new RefObject<orbit>(masterorbit);
		geocode(input_general, input_i_geocode, input_ellips, master, heightinradarsystem, TempRefObject32);
		masterorbit = TempRefObject32.argvalue;

		// ______ Update log files ______
		updatefile("scratchloggeocode", input_general.logfile);
		updatefile("scratchresgeocode", input_general.i_resfile);
		PROGRESS.print("Finished GEOCODE.");
		DEBUG.print("Time spent for geocoding:");
		printcpu();
		}
	  // ______ No need fillproductinfo after last one ______
	//  // ______ but do it for consistentness/preview card one day _____
	//  // ______ Fill info struct for filenames ______
	  if (alreadyprocessed[pr_i_geocoding]) // update productinfo
		{
		INFO.print("");
		INFO.print("ifg: latest known processing stage: geocoding");
		}



	  // ====== TIDY UP ======
	  // ______ Update process control flags in resultfiles ______
	  if (processmaster)
		updateprocesscontrol(input_general.m_resfile, MASTERID);
	  if (processlave)
		updateprocesscontrol(input_general.s_resfile, SLAVEID);
	  if (processinterf)
		updateprocesscontrol(input_general.i_resfile, INTERFID);

	  DEBUG.print("Time spent for tidying up:");
	  printcpu();

	  System.out.print("\n\nProcessing results are in parameter file");
	  if ((processmaster && processlave) || (processmaster && processinterf) || (processlave && processinterf))
		{
		System.out.print("s:\n");
		}
	  else
		{
		System.out.print(":\n");
		}
	  if (processmaster)
		{
		System.out.print("   ");
		System.out.print(input_general.m_resfile);
		System.out.print("\n");
		}
	  if (processlave)
		{
		System.out.print("   ");
		System.out.print(input_general.s_resfile);
		System.out.print("\n");
		}
	  if (processinterf)
		{
		System.out.print("   ");
		System.out.print(input_general.i_resfile);
		System.out.print("\n");
		}
	  quote();
	  WARNING.summary();

	  } // --- end of try block; now catch thrown exceptions -------------
	  catch(USAGE_ERROR usage_error) // catch errors of EXCEPTION class
		{
		exit(-1); // normal exit;
		}
	  catch(EXCEPTION error) // catch errors of EXCEPTION class
		{
		cerr << "Caught error of EXCEPTION class!" << "\n";
		cerr << "It is: " << (String)error << "\n";
		printcpu();
		WARNING.summary();
		exit(1);
		}
	  catch(String error_string) // catch handled errors
		{
		cerr << "Caught error!" << "\n";
		cerr << "It is: " << error_string << "\n";
		printcpu();
		WARNING.summary();
		exit(2);
		}
	  catch (Throwable e) // catches other errors
		{
		cerr << "Caught an unhandled error!" << "\n";
		printcpu();
		WARNING.summary();
		exit(3);
		}

	  System.out.print("\n\nNormal termination.\nThank you for using Doris.\n\n");
	  return int32(0);
	  } // END main


	//***************************************************************
	// @brief  return extra libraries used with Doris                 *
	// @param  none                                                   *
	// * usage                                                        *
	// * MA 19-May-2009                                               *
	// * TODO merge with ioroutines usedinfo?                         *
	// ***************************************************************
	public static String checklibs()
	  {
	  TRACE_FUNCTION("checklib");
	  String libs = "none!";
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//	#if __USE_FFTW_LIBRARY__
//C++ TO JAVA CONVERTER NOTE: 'extern' variable declarations are not required in Java:
	//		  extern const char fftwf_version[];
			  //libs="fftw";
			  libs =fftwf_version;
	//	#endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//	#if __USE_VECLIB_LIBRARY__
			  if (libs.length() == 5) // len("none!") == 5
				{
				libs ="veclib";
				}
			  else
				{
				libs += ", veclib";
				}
	//	#endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//	#if __USE_LAPACK_LIBRARY__
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//	  #if !WIN32 && ! __CYGWIN__ && (!__MINGW32__)
		  //#warning "win32 platform"
			int32 major;
			int32 minor;
			int32 patch;
			ilaver(major, minor, patch); // get lapack version, imported from matrix.cc

			String libver = new String(new char[15]);
			//char* libver = "lapack-xx.xx.xx"; // didn't worked why?
			String.format(libver,"lapack-%d.%d.%d", major, minor, patch); // lapack-3.1.1
	//	  #else
			 String libver = "lapack";

	//	  #endif
			  if (libs.length() == 5)
				{
				libs =String(libver);
				}
			  else
				{
				libs += ", " + String(libver);
				}
	//	#endif
		return libs;
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
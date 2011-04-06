public class GlobalMembersCoregistration
{
//C++ TO JAVA CONVERTER TODO TASK: The implementation of the following method could not be found:
	//	void cgemv(char NamelessParameter1, int NamelessParameter2, int NamelessParameter3, complex<float> NamelessParameter4, complex<float> NamelessParameter5, int NamelessParameter6, complex<float> NamelessParameter7, int NamelessParameter8, complex<float> NamelessParameter9, complex<float> NamelessParameter10, int NamelessParameter11, int NamelessParameter12);
	}
//C++ TO JAVA CONVERTER TODO TASK: The implementation of the following method could not be found:
	//	STUPID_cr4 cdotu(int NamelessParameter1, complex<float> NamelessParameter2, int NamelessParameter3, complex<float> NamelessParameter4, int NamelessParameter5);
	}
	//#endif

//***************************************************************
// *    coarseporbit                                              *
// *                                                              *
// * computes translation of slave w.r.t. master                  *
// * slave(some point) = master(same point) + trans(l,p) =>       *
// *  trans = slavecoordinates - mastercoordinates                *
// * uses orbits to find coordinates of center of master          *
// *  then solves for lin,pix of these cn. for slave.             *
// *                                                              *
// * input:                                                       *
// *  -                                                           *
// * output:                                                      *
// *  -                                                           *
// *                                                              *
// *    Bert Kampes, 12-Dec-1998                                  *
// ***************************************************************



	// ====== Prototypes ======
	// ______ Compute coarse coregistration ______
	public static void coarseporbit(input_ell ell, slcimage master, slcimage slave, RefObject<orbit> masterorbit, RefObject<orbit> slaveorbit, BASELINE baseline)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"coarseporbit (BK 12-Dec-1998)"<<ends;
		  TRACE.print();
	  }
	  final short MAXITER = 10; // maximum number of iterations
	  final double CRITERPOS = 1e-6; // 1micrometer
	  final double CRITERTIM = 1e-10; // seconds (~10-6 m)

	  // ______Initial values______    :master.approxcentreoriginal.x .y .z
	  // ______Window______            :master.currentwindow.linelo/hi , pixlo/hi
	  // ______Time______              :master.t_azi0/N , t_range0/N

	  // ______Get (approx) center pixel of current window master______
	  final int cen_lin = (master.currentwindow.linelo + master.currentwindow.linehi)/2;
	  final int cen_pix = (master.currentwindow.pixlo  + master.currentwindow.pixhi) /2;
	  final double HEI = 0.0;

	  // ______ Compute x,y,z (fill P) ______
	  // ______ P.x/y/z contains (converged) solution ______
	  cn P;
	  final int lp2xyziter = lp2xyz(cen_lin,cen_pix,ell,master,masterorbit.argvalue,P,MAXITER,CRITERPOS);

	  // ______Compute line,pixel for slave of this xyz______
	  double lin;
	  double pix;
	  final int xyz2lpiter = xyz2lp(lin,pix,slave,slaveorbit.argvalue,P,MAXITER,CRITERTIM);

	  // ______ Some extra parameters (not used, just info) ______ // BK 19-Oct-2000
	  final int Bt = Btemp(master.utc1,slave.utc1);
	  // ______ Modeled quantities ______
	  final double Bperp = baseline.get_bperp(cen_lin,cen_pix,HEI);
	  final double Bpar = baseline.get_bpar(cen_lin,cen_pix,HEI);
	  final double theta = rad2deg(baseline.get_theta(cen_lin,cen_pix,HEI));
	  final double inc_angle = rad2deg(baseline.get_theta_inc(cen_lin,cen_pix,HEI));
	  // ______ Derived quantities ______
	  final double B = baseline.get_b(cen_lin,cen_pix,HEI);
	  final double alpha = rad2deg(baseline.get_alpha(cen_lin,cen_pix,HEI));
	  final double Bh = baseline.get_bhor(cen_lin,cen_pix,HEI);
	  final double Bv = baseline.get_bvert(cen_lin,cen_pix,HEI);

	  final double Hamb = baseline.get_hamb(cen_lin,cen_pix,HEI);
	  final double orb_conv = rad2deg(baseline.get_orb_conv(cen_lin,cen_pix,HEI));

	  // ______ offset = P_slave - P_master = lin - cen_lin ______
	  INFO << "Estimated translation (l,p): " << Math.floor(lin-cen_lin +.5) << ", " << Math.floor(pix-cen_pix +.5);
	  INFO.print();

	  // ______ Write to tmp files ______
	  ofstream scratchlogfile = new ofstream("scratchlogcoarse", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchlogfile,"coarseporbit: scratchlogcoarse",__FILE__,__LINE__);
	  scratchlogfile << "\n\n*******************************************************************" << "\n* COARSE_COREGISTRATION Orbits" << "\n*******************************************************************" << "\n(Approximate) center master (line,pixel,hei): " << cen_lin << ", " << cen_pix << ", " << HEI << "\nEllipsoid WGS84 coordinates of this pixel (x,y,z): (" << P.x << ", " << P.y << ", " << P.z << ")" << "\n(line,pixel) of these coordinates in slave: " << lin << ", " << pix << "\nEstimated translation slave w.r.t. master (l,p):" << rint(lin-cen_lin) << ", " << rint(pix-cen_pix) << "\nMaximum number of iterations: " << MAXITER << "\nCriterium for position (m): " << CRITERPOS << "\nCriterium for azimuth time (s): " << CRITERTIM << " (=~ " << CRITERTIM *7.e3 << "m)" << "\nNumber of iterations conversion line,pixel to xyz: " << lp2xyziter << "\nNumber of iterations conversion xyz to line,pixel: " << xyz2lpiter << "\n*******************************************************************\n";
	  scratchlogfile.close();

	  // ______ give some extra info in resfile: Bperp, Bpar, Bh, Bv, Btemp ______
	  ofstream scratchresfile = new ofstream("scratchrescoarse", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile,"coarseporbit: scratchrescoarse",__FILE__,__LINE__);
	  scratchresfile.setf(ios.right, ios.adjustfield);
		// ______ this is read/used: ______
	  scratchresfile << "\n\n*******************************************************************" << "\n*_Start_" << processcontrol[(int)AnonymousEnum.pr_i_coarse] << "\n*******************************************************************" << "\nSome info for pixel: " << cen_lin << ", " << cen_pix << " (not used):" << "\n  Btemp:     [days]:  " << setw(10) << setiosflags(ios.right) << Bt << "      \t// Temporal baseline" << "\n  Bperp      [m]:     " << setw(10) << setiosflags(ios.right) << onedecimal(Bperp) << "      \t// Perpendicular baseline" << "\n  Bpar       [m]:     " << setw(10) << setiosflags(ios.right) << onedecimal(Bpar) << "      \t// Parallel baseline" << "\n  Bh         [m]:     " << setw(10) << setiosflags(ios.right) << onedecimal(Bh) << "      \t// Horizontal baseline" << "\n  Bv         [m]:     " << setw(10) << setiosflags(ios.right) << onedecimal(Bv) << "      \t// Vertical baseline" << "\n  B          [m]:     " << setw(10) << setiosflags(ios.right) << onedecimal(B) << "      \t// Baseline (distance between sensors)" << "\n  alpha      [deg]:   " << setw(10) << setiosflags(ios.right) << onedecimal(alpha) << "      \t// Baseline orientation" << "\n  theta      [deg]:   " << setw(10) << setiosflags(ios.right) << onedecimal(theta) << "      \t// look angle" << "\n  inc_angle  [deg]:   " << setw(10) << setiosflags(ios.right) << onedecimal(inc_angle) << "      \t// incidence angle" << "\n  orbitconv  [deg]:   " << setw(10) << setiosflags(ios.right) << orb_conv << "      \t// angle between orbits" << "\n  Height_amb [m]:     " << setw(10) << setiosflags(ios.right) << onedecimal(Hamb) << "      \t// height = h_amb*phase/2pi (approximately)" << "\n  Control point master (line,pixel,hei) = (" << cen_lin << ", " << cen_pix << ", " << HEI << ")" << "\n  Control point slave  (line,pixel,hei) = (" << lin << ", " << pix << ", " << HEI << ")" << "\nEstimated translation slave w.r.t. master (slave-master):" << "\n  Positive offsetL: slave image is to the bottom" << "\n  Positive offsetP: slave image is to the right" << "\nCoarse_orbits_translation_lines:  \t" << rint(lin-cen_lin) << "\nCoarse_orbits_translation_pixels: \t" << rint(pix-cen_pix) << "\n*******************************************************************" << "\n* End_" << processcontrol[(int)AnonymousEnum.pr_i_coarse] << "_NORMAL" << "\n*******************************************************************\n";

	  // ______Tidy up______
	  scratchresfile.close();
	  PROGRESS.print("Coarse precise orbits coregistration finished.");
	  } // END coarseporbit

//***************************************************************
// *    coarsecorrel                                              *
// *                                                              *
// * computes translation of slave w.r.t. master                  *
// * slave(some point) = master(same point) + trans(l,p) =>       *
// *  trans = slavecoordinates - mastercoordinates                *
// * uses correlation between magnitude of slave/master image     *
// *                                                              *
// * requires things on disk, input                               *
// * input:                                                       *
// *  -                                                           *
// * output:                                                      *
// *  -                                                           *
// *                                                              *
// *    Bert Kampes, 12-Dec-1998                                  *
// ***************************************************************


	// ______ Coarse coregistration ______
	public static void coarsecorrel(input_coarsecorr coarsecorrinput, slcimage minfo, slcimage sinfo)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"coarsecorrel (BK 12-Dec-1998)"<<ends;
		  TRACE.print();
	  }

	  String dummyline = new String(new char[ONE27]); // for errormessages
	  //const uint Mfilelines   = minfo.currentwindow.lines();
	  //const uint Sfilelines   = sinfo.currentwindow.lines();
	  final int Nwin = coarsecorrinput.Nwin; // number of windows
	  int NwinNANrm = coarsecorrinput.Nwin; ///MA number of windows w/o -999
	  final int initoffsetL = coarsecorrinput.initoffsetL; // initila offset
	  final int initoffsetP = coarsecorrinput.initoffsetP; // initila offset
	  int MasksizeL = coarsecorrinput.MasksizeL; // size of correlation window
	  int MasksizeP = coarsecorrinput.MasksizeP; // size of correlation window
	  final int AccL = coarsecorrinput.AccL; // accuracy of initial offset
	  final int AccP = coarsecorrinput.AccP; // accuracy of initial offset
	  boolean pointsrandom = true;
	  if (specified(coarsecorrinput.ifpositions)) // filename specified
		pointsrandom = false; // only use those points


	//  INFO("Masksize ...

	// ______Only odd Masksize possible_____
	  boolean forceoddl = false;
	  boolean forceoddp = false;
	  if (!isodd(MasksizeL))
		{
		forceoddl = true;
		MasksizeL+=1; // force oddness
		}
	  if (!isodd(MasksizeP))
		{
		forceoddp = true;
		MasksizeP+=1; // force oddness
		}

	  // ______Corners of slave in master system______
	  // ______offset = A(slave system) - A(master system)______
	  final int sl0 = sinfo.currentwindow.linelo - initoffsetL;
	  final int slN = sinfo.currentwindow.linehi - initoffsetL;
	  final int sp0 = sinfo.currentwindow.pixlo - initoffsetP;
	  final int spN = sinfo.currentwindow.pixhi - initoffsetP;

	  // ______Corners of useful overlap master,slave in master system______
	  final int BORDER = 20; // slightly smaller
	  final int l0 = max((int)minfo.currentwindow.linelo,sl0) + 0.5 *MasksizeL + AccL + BORDER;
	  final int lN = min((int)minfo.currentwindow.linehi,slN) - 0.5 *MasksizeL - AccL - BORDER;
	  final int p0 = max((int)minfo.currentwindow.pixlo,sp0) + 0.5 *MasksizeP + AccP + BORDER;
	  final int pN = min((int)minfo.currentwindow.pixhi,spN) - 0.5 *MasksizeP - AccP - BORDER;
	  final window overlap = new window(l0,lN,p0,pN);

	  // ______Distribute Nwin points over window______
	  // ______Centers(i,0): line, (i,1): pixel, (i,2) flagfromdisk______
	  matrix<Integer> Centers;
	  if (pointsrandom) // no filename specified
		{
		Centers = distributepoints((float)Nwin, overlap);
		}

	  else // read in points (center of windows) from file
		{
		Centers.resize(Nwin,3);
		ifstream ifpos;
		openfstream(ifpos,coarsecorrinput.ifpositions);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		bk_assert(ifpos,coarsecorrinput.ifpositions,__FILE__,__LINE__);
		int ll;
		int pp;
		for (int i =0; i<Nwin; ++i)
		  {
		  ifpos >> ll >> pp;
		  Centers(i,0) = (int)ll; // correct for lower left corner
		  Centers(i,1) = (int)pp; // correct for lower left corner
		  Centers(i,2) = (int)1; // flag from file
		  ifpos.getline(dummyline,ONE27,'\n'); // goto next line.
		  }
		ifpos.close();

		// ______ Check last point ivm. EOL after last position in file ______
		if (Centers(Nwin-1,0) == Centers(Nwin-2,0) && Centers(Nwin-1,1) == Centers(Nwin-2,1))
		  {
		  Centers(Nwin-1,0) = (int)(.5*(lN + l0) + 27); // random
		  Centers(Nwin-1,1) = (int)(.5*(pN + p0) + 37); // random
		  WARNING << "CC: there should be no EOL after last point in file: " << coarsecorrinput.ifpositions;
		  WARNING.print();
		  }

		// ______ Check if points are in overlap ______
		// ______ no check for uniqueness of points ______
		boolean troubleoverlap = false;
		for (int i =0; i<Nwin; ++i)
		  {
		  if (Centers(i,0) < l0)
			{
			troubleoverlap =true;
			WARNING << "COARSE_CORR: point from file: " << i+1 << " " << Centers(i,0) << " " << Centers(i,1) << " outside overlap master, slave. New position: ";
			Centers(i,0) = l0 + l0-Centers(i,0);
			WARNING << Centers(i,0) << " " << Centers(i,1);
			WARNING.print();
			}
		  if (Centers(i,0) > lN)
			{
			troubleoverlap =true;
			WARNING << "COARSE_CORR: point from file: " << i+1 << " " << Centers(i,0) << " " << Centers(i,1) << " outside overlap master, slave. New position: ";
			Centers(i,0) = lN + lN-Centers(i,0);
			WARNING << Centers(i,0) << " " << Centers(i,1);
			WARNING.print();
			}
		  if (Centers(i,1) < p0)
			{
			troubleoverlap =true;
			WARNING << "COARSE_CORR: point from file: " << i+1 << " " << Centers(i,0) << " " << Centers(i,1) << " outside overlap master, slave. New position: ";
			Centers(i,1) = p0 + p0-Centers(i,1);
			WARNING << Centers(i,0) << " " << Centers(i,1);
			WARNING.print();
			}
		  if (Centers(i,1) > pN)
			{
			troubleoverlap =true;
			WARNING << "COARSE_CORR: point from file: " << i+1 << " " << Centers(i,0) << " " << Centers(i,1) << " outside overlap master, slave. New position: ";
			Centers(i,1) = pN + pN-Centers(i,1);
			WARNING << Centers(i,0) << " " << Centers(i,1);
			WARNING.print();
			}
		  }
		if (troubleoverlap) // give some additional info
		  {
		  WARNING << "FINE: there were points from file outside overlap (l0,lN,p0,pN): " << l0 << " " << lN << " " << p0 << " " << pN << ends;
		  WARNING.print();
		  }
		}


	  // ______Compute correlation of these points______
	  matrix<complex<Float>> Mcmpl;
	  matrix<complex<Float>> Scmpl;
	  matrix<Float> Master; // amplitude master
	  matrix<Float> Mask; // amplitude slave
	  matrix<Float> Correl; // matrix with correlations
	  matrix<Float> Result = new matrix(Nwin,3); // R(i,0)=correlation; (i,1)=delta l; (i,2)=delta p;

	  // ______ Progress messages ______
	  int percent = 0;
	  int tenpercent = rint(Nwin/10.0); // round
	  if (tenpercent ==0) // avoid error: x%0
		  tenpercent = 1000;
	  for (int i =0;i<Nwin;i++)
		{
		if (i%tenpercent ==0)
		  {
		  PROGRESS << "COARSE_CORR: " << setw(3) << percent << "%" << ends;
		  PROGRESS.print();
		  percent += 10;
		  }

		// ______Center of window in master system______
		int cenMwinL = Centers(i,0);
		int cenMwinP = Centers(i,1);

		window master; // size=masksize+2*acc.
		master.linelo = cenMwinL - (MasksizeL-1)/2 -AccL; // ML is forced odd
		master.linehi = master.linelo + MasksizeL +2 *AccL - 1;
		master.pixlo = cenMwinP - (MasksizeP-1)/2 - AccP; // MP is forced odd
		master.pixhi = master.pixlo + MasksizeP +2 *AccP - 1;

		// ______Same points in slave system (disk)______
		window slavemask; // size=masksize
		int cenSwinL = cenMwinL + initoffsetL; // adjust initoffset
		int cenSwinP = cenMwinP + initoffsetP; // adjust initoffset
		slavemask.linelo = cenSwinL - (MasksizeL-1)/2; // ML is forced odd
		slavemask.linehi = slavemask.linelo + MasksizeL - 1;
		slavemask.pixlo = cenSwinP - (MasksizeP-1)/2; // MP is forced odd
		slavemask.pixhi = slavemask.pixlo + MasksizeP - 1;

		// ______Read windows from files, compute magnitude______
		Mcmpl = minfo.readdata(master);
		Scmpl = sinfo.readdata(slavemask);
		Master = magnitude(Mcmpl);
		Mask = magnitude(Scmpl);

		// ______Compute correlation matrix and find maximum______
		Correl = correlate(Master,Mask);
		int L;
		int P;
	//    MA: if maximum correlation is 0, which is due to NaNs, assign -999
	//    so in getoffset they are disregarded as in magfft. See getoffset.
	//    real4 corr = max(Correl, L, P);             // returns also L,P
		float corr = (max(Correl, L, P) == 0) ? -999 : max(Correl, L, P); // returns also L,P

		int relcenML = master.linehi - cenMwinL; // system of matrix
		int relcenMP = master.pixhi - cenMwinP; // system of matrix
		int reloffsetL = relcenML - L;
		int reloffsetP = relcenMP - P;
		int offsetL = reloffsetL + initoffsetL; // estimated offset lines
		int offsetP = reloffsetP + initoffsetP; // estimated offset pixels

		Result(i,0) = corr;
		Result(i,1) = offsetL;
		Result(i,2) = offsetP;
		}

	  // ______Get correct offsetL, offsetP______
	  int offsetLines = -999;
	  int offsetPixels = -999;
	  RefObject<Integer> TempRefObject = new RefObject<Integer>(offsetLines);
	  RefObject<Integer> TempRefObject2 = new RefObject<Integer>(offsetPixels);
	  getoffset(Result, TempRefObject, TempRefObject2);
	  offsetLines = TempRefObject.argvalue;
	  offsetPixels = TempRefObject2.argvalue;


	  // ______Write to files______
	  ofstream scratchlogfile = new ofstream("scratchlogcoarse2", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchlogfile,"coarsecorrel: scratchlogcoarse2",__FILE__,__LINE__);
	  scratchlogfile << "\n\n*******************************************************************" << "\n* COARSE_COREGISTRATION: Correlation" << "\n*******************************************************************" << "\nNumber of correlation windows: \t" << Nwin << "\nCorrelation window size (l,p): \t" << MasksizeL << ", " << MasksizeP;
		if (forceoddl)
			scratchlogfile << "(l forced odd) ";
		if (forceoddp)
			scratchlogfile << "(p forced odd)";
	  scratchlogfile << "\nSearchwindow size (l,p): \t\t" << MasksizeL + 2 *AccL << ", " << MasksizeP + 2 *AccP << "\nNumber \tposl \tposp \toffsetl offsetp \tcorrelation\n";
	  for (int k =0; k<Nwin; k++)
		{
		// MA remove NaN valued coh windows from  Nwin, to be used in resfile
		if (Result(k,0) == -999)
			NwinNANrm = NwinNANrm - 1;
		scratchlogfile << k << "\t" << Centers(k,0) << "\t" << Centers(k,1) << "\t" << Result(k,1) << "\t" << Result(k,2) << "\t" << Result(k,0) << "\n";
		 }
	  scratchlogfile << "Estimated total offset (l,p): \t" << offsetLines << ", " << offsetPixels << "\n*******************************************************************\n";
	  scratchlogfile.close();

	  ofstream scratchresfile = new ofstream("scratchrescoarse2", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile,"coarsecorrel: scratchrescoarse2",__FILE__,__LINE__);
	  scratchresfile << "\n\n*******************************************************************" << "\n*_Start_" << processcontrol[(int)AnonymousEnum.pr_i_coarse2] << "\n*******************************************************************" << "\nEstimated translation slave w.r.t. master:" << "\nCoarse_correlation_translation_lines: \t" << offsetLines << "\nCoarse_correlation_translation_pixels: \t" << offsetPixels << "\nNumber of correlation windows: \t\t" << NwinNANrm << " of " << Nwin;
	  scratchresfile << "\n\n#     center(l,p)   coherence   offsetL   offsetP\n";
		for (int k =0; k<Nwin; k++)
		 {
		  //MA remove/skip -999 values before writing resfile. For magspace.
		  // All the values are kept in  doris.log 
		  if (Result(k,0) == -999)
			  continue;
		  scratchresfile << k << " \t" << Centers(k,0) << " \t" << Centers(k,1) << " \t" << Result(k,0) << " \t" << Result(k,1) << " \t" << Result(k,2) << "\n";
		 }
					 //<< "\n* End_coarse_correlation:_NORMAL"
	  scratchresfile << "\n*******************************************************************" << "\n* End_" << processcontrol[(int)AnonymousEnum.pr_i_coarse2] << "_NORMAL" << "\n*******************************************************************\n";
	  scratchresfile.close();

	// ______Tidy up______
	  INFO << "Individually estimated translations (#, l, p, corr, offl, offp): ";
	  INFO.print();
	  for (int k =0; k<Nwin; k++)
		{
		INFO << k << " \t" << Centers(k,0) << " \t" << Centers(k,1) << " \t" << Result(k,0) << " \t" << Result(k,1) << " \t" << Result(k,2);
		INFO.print();
		}
	  INFO << "Estimated translation (l,p): " << offsetLines << ", " << offsetPixels << ends;
	  INFO.print();
	  PROGRESS.print("Coarse coregistration based on correlation finished.");
	  } // END coarsecorrel

//***************************************************************
// *    coarsecorrelfft                                           *
// *                                                              *
// * computes translation of slave w.r.t. master                  *
// * slave(some point) = master(same point) + trans(l,p) =>       *
// *  trans = slavecoordinates - mastercoordinates                *
// * uses correlation between magnitude of slave/master image     *
// * uses fft to compute coherence, no subtraction of mean        *
// *                                                              *
// * requires thingsa on disk, input                              *
// * input:                                                       *
// *  -                                                           *
// * output:                                                      *
// *  -                                                           *
// *                                                              *
// *    Bert Kampes, 12-Dec-1998                                  *
// ***************************************************************


	// ______ Corr by fft ______
	public static void coarsecorrelfft(input_coarsecorr coarsecorrinput, slcimage minfo, slcimage sinfo)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"coarsecorrelfft (BK 12-Dec-1998)"<<ends;
		  TRACE.print();
	  }
	  if (coarsecorrinput.method != cc_magfft)
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "unknown method, This routine is only for cc_magfft method.";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(argument_error);
		}

	  String dummyline = new String(new char[ONE27]); // for errormessages
	  //const uint Mfilelines   = minfo.currentwindow.lines();
	  //const uint Sfilelines   = sinfo.currentwindow.lines();
	  final int Nwin = coarsecorrinput.Nwin; // number of windows
	  int NwinNANrm = coarsecorrinput.Nwin; ///MA number of windows w/o -999
	  final int initoffsetL = coarsecorrinput.initoffsetL; // initial offset
	  final int initoffsetP = coarsecorrinput.initoffsetP; // initial offset
	  final int MasksizeL = coarsecorrinput.MasksizeL; // size of correlation window
	  final int MasksizeP = coarsecorrinput.MasksizeP; // size of correlation window

	  boolean pointsrandom = true;
	  if (specified(coarsecorrinput.ifpositions)) // filename specified
		pointsrandom = false; // only use these points

	  // ______Only pow2 Masksize possible_____
	  if (!ispower2(MasksizeL))
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "coarse correl fft: MasksizeL should be 2^n";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}
	  if (!ispower2(MasksizeP))
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "coarse correl fft: MasksizeP should be 2^n";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}

	  // ______Corners of slave in master system______
	  // ______offset = [A](slave system) - [A](master system)______
	  final int sl0 = sinfo.currentwindow.linelo - initoffsetL;
	  final int slN = sinfo.currentwindow.linehi - initoffsetL;
	  final int sp0 = sinfo.currentwindow.pixlo - initoffsetP;
	  final int spN = sinfo.currentwindow.pixhi - initoffsetP;

	  // ______Corners of useful overlap master,slave in master system______
	  final int BORDER = 20; // slightly smaller
	  final int l0 = max((int)minfo.currentwindow.linelo,sl0) + BORDER;
	  final int lN = min((int)minfo.currentwindow.linehi,slN) - MasksizeL - BORDER;
	  final int p0 = max((int)minfo.currentwindow.pixlo,sp0) + BORDER;
	  final int pN = min((int)minfo.currentwindow.pixhi,spN) - MasksizeP - BORDER;
	  final window overlap = new window(l0,lN,p0,pN);

	  // ______Distribute Nwin points over window______
	  // ______Minlminp(i,0): line, (i,1): pixel, (i,2) flagfromdisk______
	  matrix<Integer> Minlminp;
	  if (pointsrandom) // no filename specified
		{
		Minlminp = distributepoints((float)Nwin, overlap);
		}
	  else // read in points (center of windows) from file
		{
		Minlminp.resize(Nwin,3);
		ifstream ifpos;
		openfstream(ifpos,coarsecorrinput.ifpositions);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		bk_assert(ifpos,coarsecorrinput.ifpositions,__FILE__,__LINE__);
		int ll;
		int pp;
		for (int i =0; i<Nwin; ++i)
		  {
		  ifpos >> ll >> pp;
		  Minlminp(i,0) = (int)(ll-0.5 *MasksizeL); // correct for lower left corner
		  Minlminp(i,1) = (int)(pp-0.5 *MasksizeP); // correct for lower left corner
		  Minlminp(i,2) = (int)1; // flag from file
		  ifpos.getline(dummyline,ONE27,'\n'); // goto next line.
		  }
		ifpos.close();

		// ______ Check last point ivm. EOL after last position in file ______
		if (Minlminp(Nwin-1,0) == Minlminp(Nwin-2,0) && Minlminp(Nwin-1,1) == Minlminp(Nwin-2,1))
		  {
		  Minlminp(Nwin-1,0) = (int)(.5*(lN + l0) + 27); // random
		  Minlminp(Nwin-1,1) = (int)(.5*(pN + p0) + 37); // random
		  }

		// ______ Check if points are in overlap ______
		// ______ no check for uniqueness of points ______
		boolean troubleoverlap = false;
		for (int i =0; i<Nwin; ++i)
		  {
		  if (Minlminp(i,0) < l0)
			{
			troubleoverlap =true;
			WARNING << "COARSECORR: point from file: " << i+1 << " " << Minlminp(i,0) +.5 *MasksizeL << " " << Minlminp(i,1) +.5 *MasksizeP << " outside overlap master, slave. New position: ";
			Minlminp(i,0) = l0 + l0-Minlminp(i,0);
			WARNING << Minlminp(i,0) << " " << Minlminp(i,1);
			WARNING.print();
			}
		  if (Minlminp(i,0) > lN)
			{
			troubleoverlap =true;
			WARNING << "COARSECORR: point from file: " << i+1 << " " << Minlminp(i,0) +.5 *MasksizeL << " " << Minlminp(i,1) +.5 *MasksizeP << " outside overlap master, slave. New position: ";
			Minlminp(i,0) = lN + lN-Minlminp(i,0);
			WARNING << Minlminp(i,0) << " " << Minlminp(i,1);
			WARNING.print();
			}
		  if (Minlminp(i,1) < p0)
			{
			troubleoverlap =true;
			WARNING << "COARSECORR: point from file: " << i+1 << " " << Minlminp(i,0) +.5 *MasksizeL << " " << Minlminp(i,1) +.5 *MasksizeP << " outside overlap master, slave. New position: ";
			Minlminp(i,1) = p0 + p0-Minlminp(i,1);
			WARNING << Minlminp(i,0) << " " << Minlminp(i,1);
			WARNING.print();
			}
		  if (Minlminp(i,1) > pN)
			{
			troubleoverlap =true;
			WARNING << "COARSECORR: point from file: " << i+1 << " " << Minlminp(i,0) + 0.5 *MasksizeL << " " << Minlminp(i,1) + 0.5 *MasksizeP << " outside overlap master, slave. New position: ";
			Minlminp(i,1) = pN + pN-Minlminp(i,1);
			WARNING << Minlminp(i,0) << " " << Minlminp(i,1);
			WARNING.print();
			}
		  }
		if (troubleoverlap) // give some additional info
		  {
		  WARNING << "COARSECORR: point in input file outside overlap (l0,lN,p0,pN): " << l0 << " " << lN << " " << p0 << " " << pN;
		  WARNING.print();
		  }
		}

	  // ______Compute coherence of these points______
	  matrix<complex<Float>> Master;
	  matrix<complex<Float>> Mask;
	  matrix<Float> Result = new matrix(Nwin,3); // R(i,0):delta l;
													//  R(i,1):delta p; R(i,2):correl
	  // ______ Progress messages ______
	  int percent = 0;
	  int tenpercent = rint(Nwin/10.0); // round
	  if (tenpercent ==0) // avoid error: x%0
		  tenpercent = 1000;
	  for (int i =0; i<Nwin; ++i)
		{
		if (i%tenpercent ==0)
		  {
		  PROGRESS << "COARSE_CORR: " << setw(3) << percent << "%";
		  PROGRESS.print();
		  percent += 10;
		  }

		// ______Minlminp (lower left corners) of window in master system______
		final int minMwinL = Minlminp(i,0);
		final int minMwinP = Minlminp(i,1);
		DEBUG.print(" ");
		DEBUG << "Window: " << i << " [" << minMwinL << ", " << minMwinP << "]";
		DEBUG.print();
		window master = new window(minMwinL, minMwinL+MasksizeL-1, minMwinP, minMwinP+MasksizeP-1); // size=masksize
		// ______Same points in slave system (disk)______
		window mask = new window(minMwinL+initoffsetL, minMwinL+initoffsetL+MasksizeL-1, minMwinP+initoffsetP, minMwinP+initoffsetP+MasksizeP-1);

		// ______Read windows from files______
		Master = minfo.readdata(master);
		Mask = sinfo.readdata(mask);

		// ______ Coherence/max correlation ______
		float offsetL;
		float offsetP;
		//const real4 coheren = corrfft(absMaster,absMask,offsetL,offsetP);
		//const real4 coheren = coherencefft(Master, Mask, 
		//  1, MasksizeL/2, MasksizeP/2, //do not ovs, search full matrix for max
		//  offsetL,offsetP);// returned
		RefObject<Float> TempRefObject = new RefObject<Float>(offsetL);
		RefObject<Float> TempRefObject2 = new RefObject<Float>(offsetP);
		final float coheren = crosscorrelate(Master, Mask, 1, MasksizeL/2, MasksizeP/2, TempRefObject, TempRefObject2); // returned
		offsetL = TempRefObject.argvalue;
		offsetP = TempRefObject2.argvalue;
		DEBUG << "Offset between chips (l,p)    = " << offsetL << ", " << offsetP;
		DEBUG.print();

		// ______ Store result of this patch ______
		Result(i,0) = coheren;
		Result(i,1) = initoffsetL + offsetL; // total estimated offset
		Result(i,2) = initoffsetP + offsetP; // total estimated offset
		DEBUG << "Offset between images on disk = " << Result(i,1) << ", " << Result(i,2) << " (corr=" << coheren << ")";
		DEBUG.print();
		} // for nwin

	  // ______ Position approx. with respect to center of window ______
	  // ______ correct position array for center instead of lower left ______
	  for (int i =0; i<Nwin; i++)
		{
		Minlminp(i,0) += (int)(0.5 *MasksizeL);
		Minlminp(i,1) += (int)(0.5 *MasksizeP);
		}

	  // ______ Get good general estimate for offsetL, offsetP ______
	  int offsetLines = -999;
	  int offsetPixels = -999;
	  RefObject<Integer> TempRefObject3 = new RefObject<Integer>(offsetLines);
	  RefObject<Integer> TempRefObject4 = new RefObject<Integer>(offsetPixels);
	  getoffset(Result, TempRefObject3, TempRefObject4);
	  offsetLines = TempRefObject3.argvalue;
	  offsetPixels = TempRefObject4.argvalue;

	  // ______ Write to files ______
	  ofstream scratchlogfile = new ofstream("scratchlogcoarse2", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchlogfile,"coarsecorrelfft: scratchlogcoarse2",__FILE__,__LINE__);
	  scratchlogfile << "\n\n*******************************************************************" << "\n* COARSE_COREGISTRATION: Correlation" << "\n*******************************************************************" << "\nNumber of correlation windows: \t" << Nwin << "\nwindow size (l,p):             \t" << MasksizeL << ", " << MasksizeP << "\n\nNumber \tposL \tposP \toffsetL offsetP\tcorrelation\n";
	  for (int k =0; k<Nwin; k++)
		{
		// MA remove NaN valued coh windows from  Nwin, to be used in resfile
		if (Result(k,0) == -999)
			NwinNANrm = NwinNANrm - 1;
		scratchlogfile << k << "\t" << Minlminp(k,0) << "\t" << Minlminp(k,1) << "\t" << Result(k,1) << "\t" << Result(k,2) << "\t" << Result(k,0) << "\n";
		 }
	  scratchlogfile << "Estimated total offset (l,p): \t" << offsetLines << ", " << offsetPixels << "\nCoherence -999 values are disregarded in the analysis." << "\n*******************************************************************\n";
	  scratchlogfile.close();

	  ofstream scratchresfile = new ofstream("scratchrescoarse2", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile,"coarsecorrelfft: scratchrescoarse2",__FILE__,__LINE__);
	  scratchresfile << "\n\n*******************************************************************" << "\n*_Start_" << processcontrol[(int)AnonymousEnum.pr_i_coarse2] << "\n*******************************************************************" << "\nEstimated translation slave w.r.t. master:" << "\nCoarse_correlation_translation_lines: \t" << offsetLines << "\nCoarse_correlation_translation_pixels: \t" << offsetPixels << "\nNumber of correlation windows: \t\t" << NwinNANrm << " of " << Nwin;
	  scratchresfile << "\n\n#     center(l,p)   coherence   offsetL   offsetP\n";
		for (int k =0; k<Nwin; k++)
		 {
		  //MA remove/skip NaN -999 values before writing resfile. For magfft.
		  // All the values are kept in  doris.log 
		  if (Result(k,0) == -999)
			  continue;
		  scratchresfile << k << " \t" << Minlminp(k,0) << " \t" << Minlminp(k,1) << " \t" << Result(k,0) << " \t" << Result(k,1) << " \t" << Result(k,2) << "\n";
		 }
					 //<< "\n* End_coarse_correlation:_NORMAL"
	  scratchresfile << "\n*******************************************************************" << "\n* End_" << processcontrol[(int)AnonymousEnum.pr_i_coarse2] << "_NORMAL" << "\n*******************************************************************\n";

	// ______Tidy up______
	  scratchresfile.close();
	  INFO << "Individual estimated translations (#, l, p, corr, offl, offp):";
	  INFO.print();
	  for (int k =0; k<Nwin; k++)
		{
		INFO << k << " \t" << Minlminp(k,0) << " \t" << Minlminp(k,1) << " \t" << Result(k,0) << " \t" << Result(k,1) << " \t" << Result(k,2) << " \t";
		INFO.print();
		}

	  INFO << "Estimated overall translation (l,p): " << offsetLines << ", " << offsetPixels;
	  INFO.print();
	  INFO << "Coherence -999 values are disregarded in the analysis."; //MA see getoffset
	  INFO.print();
	  PROGRESS.print("Coarse coregistration based on correlation finished.");
	  } // END coarsecorrelfft

//***************************************************************
// *    mtiming_correl (coarse (ok) + fine ?)                     *
// *                                                              *
// * computes translation of  master w.r.t. DEM (sim. amplitude)  *
// * master(some point) = simamp (same point) + trans(l,p) =>     *
// *  trans = mastercoordinates - simamp.coordinates              *
// * uses correlation between magnitude of master and simamp      *
// * image to estimate overall shift.                             *
// *                                                              *
// * requires things on disk, input                               *
// * input:                                                       *
// *  - input settings,                                           *
// *  - master info                                               *
// *  - simulated amplitude info                                  *
// * output:                                                      *
// *  - coarse offsets between dem and the master                 *
// *                                                              *
// *    Bert Kampes, 12-Dec-1998  (coarsecorr)                    *
// *    Batuhan Osmanoglu, 30-JUL-2007 (demcorr for phase)        *
// *    Mahmut Arikan, 12-Nov-2008                                *
// ***************************************************************

	// ______ Sim. Amplitude coregistration (magspace) ______
			// const input_coarsecorr  &input,
	public static void mtiming_correl(input_mtiming mtiminginput, slcimage minfo, productinfo sinfo)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"mtiming_correl (MA,BO 12-Nov-2008)"<<ends;
		  TRACE.print();
	  }

	  final String STEP ="MTIMING: "; // step name
	  String dummyline = new String(new char[ONE27]); // for errormessages
	  //const uint Mfilelines   = minfo.currentwindow.lines();
	  //const uint Sfilelines   = sinfo.currentwindow.lines();
	  final int Nwin = mtiminginput.Nwin; // number of windows
	  int NwinNANrm = mtiminginput.Nwin; ///MA number of windows w/o -999
	  final int initoffsetL = mtiminginput.initoffsetL; // initial offset not nec for simamp
	  final int initoffsetP = mtiminginput.initoffsetP; // initial offset
	  int MasksizeL = mtiminginput.MasksizeL; // size of correlation window
	  int MasksizeP = mtiminginput.MasksizeP; // size of correlation window
	  final int AccL = mtiminginput.AccL; // accuracy of initial offset
	  final int AccP = mtiminginput.AccP; // accuracy of initial offset
	  boolean pointsrandom = true;
	  if (specified(mtiminginput.ifpositions)) // filename specified
		pointsrandom = false; // only use these points


	//  INFO("Masksize ...

	// ______Only odd Masksize possible_____
	  boolean forceoddl = false;
	  boolean forceoddp = false;
	  if (!isodd(MasksizeL))
		{
		forceoddl = true;
		MasksizeL+=1; // force oddness
		}
	  if (!isodd(MasksizeP))
		{
		forceoddp = true;
		MasksizeP+=1; // force oddness
		}

	  // ______Corners of simamp(dem) in master system______
	  // ______offset = A(master system) - A(slave system)______
	  final int sl0 = sinfo.win.linelo - initoffsetL; // [MA] sim. ampl. image extend should be
	  final int slN = sinfo.win.linehi - initoffsetL; // the same as master crop extend. Kept for the convience.
	  final int sp0 = sinfo.win.pixlo - initoffsetP;
	  final int spN = sinfo.win.pixhi - initoffsetP;
	  DEBUG << "slave l0: " << sl0 << " slN " << slN << " sp0 " << sp0 << " spN " << spN;
	  DEBUG.print();

	  // ______Corners of useful overlap master,slave in master system______
	  final int BORDER = 20; // slightly smaller
	  final int l0 = max((int)minfo.currentwindow.linelo,sl0) + 0.5 *MasksizeL + AccL + BORDER;
	  final int lN = min((int)minfo.currentwindow.linehi,slN) - 0.5 *MasksizeL - AccL - BORDER;
	  final int p0 = max((int)minfo.currentwindow.pixlo,sp0) + 0.5 *MasksizeP + AccP + BORDER;
	  final int pN = min((int)minfo.currentwindow.pixhi,spN) - 0.5 *MasksizeP - AccP - BORDER;
	//
	//  // ______Check masksize against height and width of the crop______
	//  if( int32(MasksizeL) > int32(lN-l0) || int32(MasksizeP) > int32(pN-p0) )
	//    {
	//     ERROR << "MTE: Impossible to continue! Masksize larger than the overlapping crop width or height. Please check.";
	//     ERROR.print();     
	//     ERROR << "MTE: MasksizeL [" << MasksizeL << "] > crop height [" << int32(lN-l0) << "] ?";
	//     ERROR.print();     
	//     ERROR << "MTE: MasksizeP [" << MasksizeP << "] >  crop width [" << int32(pN-p0) << "] ?";
	//     ERROR.print();     
	//    throw(input_error) ;
	//    }
	//
	  DEBUG << "mastercurrentwinl0: " << minfo.currentwindow.linelo << " lN " << minfo.currentwindow.linehi << " p0 " << minfo.currentwindow.pixlo << " pN " << minfo.currentwindow.pixhi;
	  DEBUG.print();
	  DEBUG << "         master l0: " << l0 << " lN " << lN << " p0 " << p0 << " pN " << pN;
	  DEBUG.print();
	  final window overlap = new window(l0,lN,p0,pN);

	  DEBUG << "overlap l0: " << l0 << " lN " << lN << " p0 " << p0 << " pN " << pN;
	  DEBUG.print();

	  // ______Distribute Nwin points over window______
	  // ______Centers(i,0): line, (i,1): pixel, (i,2) flagfromdisk______
	  matrix<Integer> Centers;
	  if (pointsrandom) // no filename specified
		{
		Centers = distributepoints((float)Nwin, overlap);
		}

	  else // read in points (center of windows) from file
		{
		Centers.resize(Nwin,3);
		ifstream ifpos;
		openfstream(ifpos,mtiminginput.ifpositions);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		bk_assert(ifpos,mtiminginput.ifpositions,__FILE__,__LINE__);
		int ll;
		int pp;
		for (int i =0; i<Nwin; ++i)
		  {
		  ifpos >> ll >> pp;
		  Centers(i,0) = (int)ll; // correct for lower left corner
		  Centers(i,1) = (int)pp; // correct for lower left corner
		  Centers(i,2) = (int)1; // flag from file
		  ifpos.getline(dummyline,ONE27,'\n'); // goto next line.
		  }
		ifpos.close();

		// ______ Check last point ivm. EOL after last position in file ______
		if (Centers(Nwin-1,0) == Centers(Nwin-2,0) && Centers(Nwin-1,1) == Centers(Nwin-2,1))
		  {
		  Centers(Nwin-1,0) = (int)(.5*(lN + l0) + 27); // random
		  Centers(Nwin-1,1) = (int)(.5*(pN + p0) + 37); // random
		  WARNING << "MTE: there should be no EOL after last point in file: " << mtiminginput.ifpositions;
		  WARNING.print();
		  }

		// ______ Check if points are in overlap ______
		// ______ no check for uniqueness of points ______
		boolean troubleoverlap = false;
		for (int i =0; i<Nwin; ++i)
		  {
		  if (Centers(i,0) < l0)
			{
			troubleoverlap =true;
			WARNING << STEP << "point from file: " << i+1 << " " << Centers(i,0) << " " << Centers(i,1) << " outside overlap master, slave. New position: ";
			Centers(i,0) = l0 + l0-Centers(i,0);
			WARNING << Centers(i,0) << " " << Centers(i,1);
			WARNING.print();
			}
		  if (Centers(i,0) > lN)
			{
			troubleoverlap =true;
			WARNING << STEP << "point from file: " << i+1 << " " << Centers(i,0) << " " << Centers(i,1) << " outside overlap master, slave. New position: ";
			Centers(i,0) = lN + lN-Centers(i,0);
			WARNING << Centers(i,0) << " " << Centers(i,1);
			WARNING.print();
			}
		  if (Centers(i,1) < p0)
			{
			troubleoverlap =true;
			WARNING << STEP << "point from file: " << i+1 << " " << Centers(i,0) << " " << Centers(i,1) << " outside overlap master, slave. New position: ";
			Centers(i,1) = p0 + p0-Centers(i,1);
			WARNING << Centers(i,0) << " " << Centers(i,1);
			WARNING.print();
			}
		  if (Centers(i,1) > pN)
			{
			troubleoverlap =true;
			WARNING << STEP << "point from file: " << i+1 << " " << Centers(i,0) << " " << Centers(i,1) << " outside overlap master, slave. New position: ";
			Centers(i,1) = pN + pN-Centers(i,1);
			WARNING << Centers(i,0) << " " << Centers(i,1);
			WARNING.print();
			}
		  }
		if (troubleoverlap) // give some additional info
		  {
		  WARNING << STEP << "there were points in input file which lie outside overlap (l0,lN,p0,pN): " << l0 << " " << lN << " " << p0 << " " << pN << ends;
		  WARNING.print();
		  }
		}

	  // ______Compute correlation of these points______
	  matrix<complex<Float>> Mcmpl; // Master complex image
	  matrix<Float> Sampl; // Simulated amplitude
	  matrix<Float> mMag; // amplitude master
	  matrix<Float> Correl; // matrix with correlations
	  matrix<Float> Result = new matrix(Nwin,3); // R(i,0)=correlation; (i,1)=delta l; (i,2)=delta p;

	  // ______ Progress messages ______
	  int percent = 0;
	  int tenpercent = rint(Nwin/10.0); // round
	  if (tenpercent ==0) // avoid error: x%0
		  tenpercent = 1000;
	  for (int i =0; i<Nwin; ++i)
		{
		if (i%tenpercent ==0)
		  {
		  PROGRESS << STEP << setw(3) << percent << "%" << ends;
		  PROGRESS.print();
		  percent += 10;
		  }

		// ______Center of window in master system______
		 int cenMwinL = Centers(i,0);
		 int cenMwinP = Centers(i,1);

		DEBUG.print(" ");
		DEBUG << "Window(cen): " << i << " [" << cenMwinL << ", " << cenMwinP << "]";
		DEBUG.print();

		window mwin; // big patch: size=masksize+2*acc.
		mwin.linelo = cenMwinL - (MasksizeL-1)/2 -AccL; // ML is forced odd
		mwin.linehi = mwin.linelo + MasksizeL +2 *AccL - 1;
		mwin.pixlo = cenMwinP - (MasksizeP-1)/2 - AccP; // MP is forced odd
		mwin.pixhi = mwin.pixlo + MasksizeP +2 *AccP - 1;

	  // Products actually only hold data within the window. 
	  // Therefore we need to convert back to file's(x,y) before reading data.
	  // Batu 2007 08 01
	  // uint cenSwinL    = cenMwinL + initoffsetL - sinfo.win.linelo +1 ;          // adjust initoffset
	  // [MA] this is fixed in products::readr4
		// ______Same points in slave system (disk)______
		window swin; // small patch: size=masksize
		int cenSwinL = cenMwinL + initoffsetL; // adjust initoffset
		int cenSwinP = cenMwinP + initoffsetP; // adjust initoffset
		swin.linelo = cenSwinL - (MasksizeL-1)/2; // ML is forced odd
		swin.linehi = swin.linelo + MasksizeL - 1;
		swin.pixlo = cenSwinP - (MasksizeP-1)/2; // MP is forced odd
		swin.pixhi = swin.pixlo + MasksizeP - 1;
	//    DEBUG << "   cenSwinL " << cenSwinL << " cenSwinP " << cenSwinL; 
	//    DEBUG.print();
	//    DEBUG << "sl0 " << swin.linelo << " slN " << swin.linehi  << " sp0 " << swin.pixlo << " spN " << swin.pixhi;
	//    DEBUG.print();

		// ______Read windows from files, compute magnitude______
		// Sampl  = sinfo.readdatar4(master); // readfile(Sampl,master,numberoflatpixels?,winfromfile?,zerooffset)
		Mcmpl = minfo.readdata(swin); // small patch
		Sampl = sinfo.readdatar4(mwin); // big patch
		mMag = magnitude(Mcmpl);
		matrix<Float> sMask = mMag; // amplitude small patch from master that shifts over
		matrix<Float> mMask = Sampl; // amplitude big patch from simamp

		// ______Compute correlation matrix and find maximum______
		//#Correl = correlate(Master,Mask); 
		Correl = correlate(mMask,sMask); // correlate(simamp,masteramp)
		int L;
		int P;
	//    MA: if maximum correlation is 0, which is due to NaNs, assign -999
	//    so in getoffset they are disregarded.
	//    real4 corr = max(Correl, L, P);             // returns also L,P
		float corr = (max(Correl, L, P) == 0) ? -999 : max(Correl, L, P); // returns also L,P

		int relcenML = mwin.linehi - cenMwinL; // system of matrix
		int relcenMP = mwin.pixhi - cenMwinP; // system of matrix
		int reloffsetL = relcenML - L;
		int reloffsetP = relcenMP - P;
		DEBUG << "Offset between chips (l,p)    = " << reloffsetL << ", " << reloffsetP;
		DEBUG.print();

		// ______ Store result of this patch ______
		Result(i,0) = corr;
		Result(i,1) = initoffsetL + reloffsetL; // total estimated offset lines
		Result(i,2) = initoffsetP + reloffsetP; // total estimated offset pixels
		DEBUG << "Offset between images on disk = " << Result(i,1) << ", " << Result(i,2) << " (corr=" << corr << ")";
		DEBUG.print();
		} // for nwin

	  // ______ Get good general estimate for offsetL, offsetP ______
	  int offsetLines = -999; // NaN
	  int offsetPixels = -999;
	  //getoffset(Result,offsetLines,offsetPixels);   // getoffsets based on Mean
	  RefObject<Integer> TempRefObject = new RefObject<Integer>(offsetLines);
	  RefObject<Integer> TempRefObject2 = new RefObject<Integer>(offsetPixels);
	  getmodeoffset(Result, TempRefObject, TempRefObject2); // [MA] max occurence
	  offsetLines = TempRefObject.argvalue;
	  offsetPixels = TempRefObject2.argvalue;


	  // ______ Convert offsets to seconds and write  master time offset to res file ______
	  // using overall coarse offsets determing master timing error

	  // ______ Initialize Variables ______
	  double masterAztime = -999;
	  double masterRatime = -999;

	  // ______ Compute Time ______
	  // minus sign is due to the offsets being reference to DEM (offset = master-dem)
	  offsets2timing(minfo, -offsetLines, -offsetPixels, masterAztime, masterRatime); // using overall offsets to
																					  // determine master timing error 

	  INFO << "Estimated master azimuth timing error [sec]: " << masterAztime << " sec.";
	  INFO.print();
	  INFO << "Estimated master range timing error   [sec]: " << masterRatime << " sec.";
	  INFO.print();

	  // azimuth and range time are later updated at proccess.cc
	  // ______ End of conversion offset to timing errors______

	  // ______ Write to files ______
	  ofstream scratchlogfile = new ofstream("scratchlogmtiming", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchlogfile,"mtiming_correl: scratchlogmtiming",__FILE__,__LINE__);
	  scratchlogfile << "\n\n*******************************************************************" << "\n* MTIMING_COREGISTRATION: Correlation" << "\n*******************************************************************" << "\nCorrelation method: \t\t" << " magspace" << "\nNumber of correlation windows: \t" << Nwin << "\nCorrelation window size (l,p): \t" << MasksizeL << ", " << MasksizeP;
		if (forceoddl)
			scratchlogfile << " (l forced odd)";
		if (forceoddp)
			scratchlogfile << " (p forced odd)";
	  scratchlogfile << "\nSearchwindow size (l,p): \t\t" << MasksizeL + 2 *AccL << ", " << MasksizeP + 2 *AccP << "\nNumber \tposl \tposp \toffsetl offsetp \tcorrelation\n";
	  for (int k =0; k<Nwin; k++)
		{
		// MA remove NaN valued coh windows from  Nwin, to be used in resfile
		if (Result(k,0) == -999)
			NwinNANrm = NwinNANrm - 1;
		scratchlogfile << k << "\t" << Centers(k,0) << "\t" << Centers(k,1) << "\t" << Result(k,1) << "\t" << Result(k,2) << "\t" << Result(k,0) << "\n";
		 }
	  scratchlogfile << "Estimated total offset (l,p): \t" << offsetLines << ", " << offsetPixels << "\nCoherence NaN values are disregarded in the analysis." << "\n*******************************************************************\n";
	  scratchlogfile.close();

	  ofstream scratchresfile = new ofstream("scratchresmtiming", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile,"mtiming_correl: scratchresmtiming",__FILE__,__LINE__);
	  scratchresfile << "\n\n*******************************************************************" << "\n*_Start_" << processcontrol[(int)AnonymousEnum.pr_m_mtiming] << " " << "" << "\n*******************************************************************" << "\nCorrelation method \t\t\t: \t" << "magspace " << "(" << MasksizeL + 2 *AccL << "," << MasksizeP + 2 *AccP << ")" << "\nNumber of correlation windows used \t: \t" << NwinNANrm << " of " << Nwin << "\nEstimated translation master w.r.t. synthetic amplitude (master-dem):" << "\n  Positive offsetL: master image is to the bottom" << "\n  Positive offsetP: master image is to the right" << "\nCoarse_correlation_translation_lines    : \t" << offsetLines << "\nCoarse_correlation_translation_pixels   : \t" << offsetPixels << "\nMaster_azimuth_timing_error             : \t" << masterAztime << " sec." << "\nMaster_range_timing_error               : \t" << masterRatime << " sec."; // in seconds
	  scratchresfile << "\n*******************************************************************" << "\n* End_" << processcontrol[(int)AnonymousEnum.pr_m_mtiming] << "_NORMAL" << "\n*******************************************************************\n";
	  scratchresfile.close();

	// ______Tidy up______
	  INFO << "Individual estimated translations (#, l, p, corr, offl, offp):";
	  INFO.print();
	  for (int k =0; k<Nwin; k++)
		{
		INFO << k << " \t" << Centers(k,0) << " \t" << Centers(k,1) << " \t" << Result(k,0) << " \t" << Result(k,1) << " \t" << Result(k,2) << " \t";
		INFO.print();
		}

	  PROGRESS << "Estimated overall translation (l,p): " << offsetLines << ", " << offsetPixels << " (used)" << ends;
	  PROGRESS.print();
	  INFO << "Coherence NaN values are disregarded in the analysis."; //MA see getoffset
	  INFO.print();
	  PROGRESS.print("MASTER TIMING Error estimation finished.");
	  } // END mtiming_correl

//***************************************************************
// *    mtiming_correlfft                                         *
// *                                                              *
// * computes translation of  master w.r.t. DEM (sim. amplitude)  *
// * master(some point) = simamp (same point) + trans(l,p) =>     *
// *  trans = mastercoordinates - simamp.coordinates              *
// * uses correlation between magnitude of master and simamp      *
// * image to estimate overall shift.                             *
// * uses fft to compute coherence                                *
// *                                                              *
// * requires things on disk, input                               *
// * input:                                                       *
// *  - input settings,                                           *
// *  - master info                                               *
// *  - simulated amplitude info                                  *
// * output:                                                      *
// *  - coarse offsets between dem and the master                 *
// *                                                              *
// *    Bert Kampes, 12-Dec-1998 (coarsecorrelfft)                *
// *    Mahmut Arikan, 04-Dec-2008 
// ***************************************************************

	// ______ Sim. Amplitude coregistration (magfft) ______
	public static void mtiming_correlfft(input_mtiming mtiminginput, slcimage minfo, productinfo sinfo)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"mtiming_correlfft (MA 04-Dec-2008)"<<ends;
		  TRACE.print();
	  }
	  if (mtiminginput.method != cc_magfft)
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "unknown method, This routine is only for cc_magfft method.";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(argument_error);
		}

	  final String STEP ="MTIMING: "; // step name
	  String dummyline = new String(new char[ONE27]); // for errormessages
	  //const uint Mfilelines   = minfo.currentwindow.lines();
	  //const uint Sfilelines   = sinfo.currentwindow.lines();
	  final int Nwin = mtiminginput.Nwin; // number of windows
	  int NwinNANrm = mtiminginput.Nwin; ///MA number of windows w/o -999
	  final int initoffsetL = mtiminginput.initoffsetL; // initial offset
	  final int initoffsetP = mtiminginput.initoffsetP; // initial offset
	  final int MasksizeL = mtiminginput.MasksizeL; // size of correlation window
	  final int MasksizeP = mtiminginput.MasksizeP; // size of correlation window

	  boolean pointsrandom = true;
	  if (specified(mtiminginput.ifpositions)) // filename specified
		pointsrandom = false; // only use these points

	  // ______Only pow2 Masksize possible_____
	  if (!ispower2(MasksizeL))
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "mtiming correl fft: MasksizeL should be 2^n";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}
	  if (!ispower2(MasksizeP))
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "mtiming correl fft: MasksizeP should be 2^n";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}

	  // ______Corners of simamp(dem) in master system______
	  // ______offset = [A](slave system) - [A](master system)______
	  final int sl0 = sinfo.win.linelo - initoffsetL;
	  final int slN = sinfo.win.linehi - initoffsetL;
	  final int sp0 = sinfo.win.pixlo - initoffsetP;
	  final int spN = sinfo.win.pixhi - initoffsetP;
	  DEBUG << "slave l0: " << sl0 << " slN " << slN << " sp0 " << sp0 << " spN " << spN;
	  DEBUG.print();

	  // ______Corners of useful overlap master,slave in master system______
	  final int BORDER = 20; // slightly smaller
	  final int l0 = max((int)minfo.currentwindow.linelo,sl0) + BORDER;
	  final int lN = min((int)minfo.currentwindow.linehi,slN) - MasksizeL - BORDER;
	  final int p0 = max((int)minfo.currentwindow.pixlo,sp0) + BORDER;
	  final int pN = min((int)minfo.currentwindow.pixhi,spN) - MasksizeP - BORDER;
	  final window overlap = new window(l0,lN,p0,pN);

	  DEBUG << "overlap l0: " << l0 << " lN " << lN << " p0 " << p0 << " pN " << pN;
	  DEBUG.print();

	  // ______Distribute Nwin points over window______
	  // ______Minlminp(i,0): line, (i,1): pixel, (i,2) flagfromdisk______
	  matrix<Integer> Minlminp;
	  if (pointsrandom) // no filename specified
		{
		Minlminp = distributepoints((float)Nwin, overlap);
		}
	  else // read in points (center of windows) from file
		{
		Minlminp.resize(Nwin,3);
		ifstream ifpos;
		openfstream(ifpos,mtiminginput.ifpositions);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		bk_assert(ifpos,mtiminginput.ifpositions,__FILE__,__LINE__);
		int ll;
		int pp;
		for (int i =0; i<Nwin; ++i)
		  {
		  ifpos >> ll >> pp;
		  Minlminp(i,0) = (int)(ll-0.5 *MasksizeL); // correct for lower left corner
		  Minlminp(i,1) = (int)(pp-0.5 *MasksizeP); // correct for lower left corner
		  Minlminp(i,2) = (int)1; // flag from file
		  ifpos.getline(dummyline,ONE27,'\n'); // goto next line.
		  }
		ifpos.close();

		// ______ Check last point ivm. EOL after last position in file ______
		if (Minlminp(Nwin-1,0) == Minlminp(Nwin-2,0) && Minlminp(Nwin-1,1) == Minlminp(Nwin-2,1))
		  {
		  Minlminp(Nwin-1,0) = (int)(.5*(lN + l0) + 27); // random
		  Minlminp(Nwin-1,1) = (int)(.5*(pN + p0) + 37); // random
		  WARNING << "MTE: there should be no EOL after last point in file: " << mtiminginput.ifpositions;
		  WARNING.print();
		  }

		// ______ Check if points are in overlap ______
		// ______ no check for uniqueness of points ______
		boolean troubleoverlap = false;
		for (int i =0; i<Nwin; ++i)
		  {
		  if (Minlminp(i,0) < l0)
			{
			troubleoverlap =true;
			WARNING << STEP << "point from file: " << i+1 << " " << Minlminp(i,0) +.5 *MasksizeL << " " << Minlminp(i,1) +.5 *MasksizeP << " outside overlap master, slave. New position: ";
			Minlminp(i,0) = l0 + l0-Minlminp(i,0);
			WARNING << Minlminp(i,0) << " " << Minlminp(i,1);
			WARNING.print();
			}
		  if (Minlminp(i,0) > lN)
			{
			troubleoverlap =true;
			WARNING << STEP << "point from file: " << i+1 << " " << Minlminp(i,0) +.5 *MasksizeL << " " << Minlminp(i,1) +.5 *MasksizeP << " outside overlap master, slave. New position: ";
			Minlminp(i,0) = lN + lN-Minlminp(i,0);
			WARNING << Minlminp(i,0) << " " << Minlminp(i,1);
			WARNING.print();
			}
		  if (Minlminp(i,1) < p0)
			{
			troubleoverlap =true;
			WARNING << STEP << "point from file: " << i+1 << " " << Minlminp(i,0) +.5 *MasksizeL << " " << Minlminp(i,1) +.5 *MasksizeP << " outside overlap master, slave. New position: ";
			Minlminp(i,1) = p0 + p0-Minlminp(i,1);
			WARNING << Minlminp(i,0) << " " << Minlminp(i,1);
			WARNING.print();
			}
		  if (Minlminp(i,1) > pN)
			{
			troubleoverlap =true;
			WARNING << STEP << "point from file: " << i+1 << " " << Minlminp(i,0) + 0.5 *MasksizeL << " " << Minlminp(i,1) + 0.5 *MasksizeP << " outside overlap master, slave. New position: ";
			Minlminp(i,1) = pN + pN-Minlminp(i,1);
			WARNING << Minlminp(i,0) << " " << Minlminp(i,1);
			WARNING.print();
			}
		  }
		if (troubleoverlap) // give some additional info
		  {
		  WARNING << STEP << "there were points in input file which lie outside overlap (l0,lN,p0,pN): " << l0 << " " << lN << " " << p0 << " " << pN << ends;
		  WARNING.print();
		  }
		}

	  // ______Compute coherence of these points______
	  matrix<complex<Float>> Mcmpl; // Master complex image
	  matrix<Float> Sampl; // Simulated amplitude
	  matrix<complex<Float>> Scmpl; // real4 simamp --> creal4 simamp
	  matrix<Float> Result = new matrix(Nwin,3); // R(i,0):delta l;
									  //  R(i,1):delta p; R(i,2):correl

	  // ______ Progress messages ______
	  int percent = 0;
	  int tenpercent = rint(Nwin/10.0); // round
	  if (tenpercent ==0) // avoid error: x%0
		  tenpercent = 1000;
	  for (int i =0; i<Nwin; ++i)
		{
		if (i%tenpercent ==0)
		  {
		  PROGRESS << STEP << setw(3) << percent << "%" << ends;
		  PROGRESS.print();
		  percent += 10;
		  }

		// ______Minlminp (lower left corners) of window in master system______
		final int minMwinL = Minlminp(i,0);
		final int minMwinP = Minlminp(i,1);
		DEBUG.print(" ");
		DEBUG << "Window(ll): " << i << " [" << minMwinL << ", " << minMwinP << "]";
		DEBUG.print();
		window mwin = new window(minMwinL, minMwinL+MasksizeL-1, minMwinP, minMwinP+MasksizeP-1); // size=mask window size
		// ______Same points in slave system (disk)______
		window swin = new window(minMwinL+initoffsetL, minMwinL+initoffsetL+MasksizeL-1, minMwinP+initoffsetP, minMwinP+initoffsetP+MasksizeP-1);

		// ______Read windows from files______
		Mcmpl = minfo.readdata(swin); // master read patch
		Sampl = sinfo.readdatar4(mwin); // simamp (DEM) read patch
		Scmpl = mat2cr4(Sampl);
		Sampl.resize(1,1); // dealloc...
		matrix<complex<Float>> sMask = Mcmpl; // complex patch from the master that shifts over
		matrix<complex<Float>> mMask = Scmpl; // complex patch from the simamp
												 // patch sizes are equal but
												 // shifted patch can have initial
												 // offset

		// ______ Coherence/max correlation ______
		float offsetL;
		float offsetP;
		//const real4 coheren = corrfft(absMaster,absMask,offsetL,offsetP);
		//const real4 coheren = coherencefft(Master, Mask, 
		//  1, MasksizeL/2, MasksizeP/2, //do not ovs, search full matrix for max
		//  offsetL,offsetP);// returned
		RefObject<Float> TempRefObject = new RefObject<Float>(offsetL);
		RefObject<Float> TempRefObject2 = new RefObject<Float>(offsetP);
		final float coheren = crosscorrelate(mMask, sMask, 1, MasksizeL/2, MasksizeP/2, TempRefObject, TempRefObject2); // returned
		offsetL = TempRefObject.argvalue;
		offsetP = TempRefObject2.argvalue;
		DEBUG << "Offset between chips (l,p)    = " << offsetL << ", " << offsetP;
		DEBUG.print();

		// ______ Store result of this patch ______
		Result(i,0) = coheren;
		Result(i,1) = initoffsetL + offsetL; // total estimated offset
		Result(i,2) = initoffsetP + offsetP; // total estimated offset
		DEBUG << "Offset between images on disk = " << Result(i,1) << ", " << Result(i,2) << " (corr=" << coheren << ")";
		DEBUG.print();
		} // for nwin

	  // ______ Position approx. with respect to center of window ______
	  // ______ correct position array for center instead of lower left ______
	  for (int i =0; i<Nwin; i++)
		{
		Minlminp(i,0) += (int)(0.5 *MasksizeL);
		Minlminp(i,1) += (int)(0.5 *MasksizeP);
		}

	  // ______ Get good general estimate for offsetL, offsetP ______
	  int offsetLines = -999; // NaN
	  int offsetPixels = -999;
	  //getoffset(Result,offsetLines,offsetPixels);   // getoffsets based on Mean
	  RefObject<Integer> TempRefObject3 = new RefObject<Integer>(offsetLines);
	  RefObject<Integer> TempRefObject4 = new RefObject<Integer>(offsetPixels);
	  getmodeoffset(Result, TempRefObject3, TempRefObject4); // [MA] max occurence
	  offsetLines = TempRefObject3.argvalue;
	  offsetPixels = TempRefObject4.argvalue;


	  // ______ Convert offsets to seconds and write master time offset to res file ______
	  // using overall coarse offsets determing master timing error

	  // check if some timing card are already defined: do this in processor.cc: see timingerror_flag

	  // ______ Initialize Variables ______
	  double masterAztime = -999;
	  double masterRatime = -999;

	  // ______ Compute Time ______
	  // minus sign is due to the offsets being reference to DEM (offset = master-dem)
	  offsets2timing(minfo, -offsetLines, -offsetPixels, masterAztime, masterRatime); // using overall offsets to
																					  // determine master timing error 

	  INFO << "Estimated master azimuth timing error [sec]: " << masterAztime << " sec.";
	  INFO.print();
	  INFO << "Estimated master range timing error   [sec]: " << masterRatime << " sec.";
	  INFO.print();

	  // azimuth and range time are later updated at proccess.cc
	  // ______ End of conversion offset to timing errors______

	  // ______ Write to files ______
	  ofstream scratchlogfile = new ofstream("scratchlogmtiming", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchlogfile,"mtiming_correlfft: scratchlogmtiming",__FILE__,__LINE__);
	  scratchlogfile << "\n\n*******************************************************************" << "\n* MTIMING_COREGISTRATION: Correlation" << "\n*******************************************************************" << "\nCorrelation method: \t\t" << " magfft" << "\nNumber of correlation windows: \t" << Nwin << "\nCorrelation window size (l,p):             \t" << MasksizeL << ", " << MasksizeP << "\n\nNumber \tposL \tposP \toffsetL offsetP\tcorrelation\n";
	  for (int k =0; k<Nwin; k++)
		{
		// MA remove NaN valued coh windows from  Nwin, to be used in resfile
		if (Result(k,0) == -999)
			NwinNANrm = NwinNANrm - 1;
		scratchlogfile << k << "\t" << Minlminp(k,0) << "\t" << Minlminp(k,1) << "\t" << Result(k,1) << "\t" << Result(k,2) << "\t" << Result(k,0) << "\n";
		 }
	  scratchlogfile << "Estimated total offset (l,p): \t" << offsetLines << ", " << offsetPixels << "\nCoherence NaN values are disregarded in the analysis." << "\n*******************************************************************\n";
	  scratchlogfile.close();

	  ofstream scratchresfile = new ofstream("scratchresmtiming", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile,"mtiming_correlfft: scratchresmtiming",__FILE__,__LINE__);
	  scratchresfile << "\n\n*******************************************************************" << "\n*_Start_" << processcontrol[(int)AnonymousEnum.pr_m_mtiming] << "\n*******************************************************************" << "\nCorrelation method \t\t\t: \t" << "magfft " << "(" << MasksizeL << "," << MasksizeP << ")" << "\nNumber of correlation windows used \t: \t" << NwinNANrm << " of " << Nwin << "\nEstimated translation master w.r.t. synthetic amplitude (master-dem):" << "\n  Positive offsetL: master image is to the bottom" << "\n  Positive offsetP: master image is to the right" << "\nCoarse_correlation_translation_lines    : \t" << offsetLines << "\nCoarse_correlation_translation_pixels   : \t" << offsetPixels << "\nMaster_azimuth_timing_error             : \t" << masterAztime << " sec." << "\nMaster_range_timing_error               : \t" << masterRatime << " sec."; // in seconds
	//  scratchresfile << "\n\n#     center(l,p)   coherence   offsetL   offsetP\n";
	//    for (uint k=0; k<Nwin; k++)
	//     { 
	//      //MA remove/skip NaN: -999 values before writing resfile. For magfft.
	//      // All the values are kept in  doris.log 
	//      if (  Result(k,0) == -999 ) continue;
	//      scratchresfile << k  << " \t" << Minlminp(k,0) << " \t" << Minlminp(k,1) << " \t"
	//           << Result(k,0)  << " \t" << Result(k,1)  << " \t" << Result(k,2)  << "\n";
	//     }
	  scratchresfile << "\n*******************************************************************" << "\n* End_" << processcontrol[(int)AnonymousEnum.pr_m_mtiming] << "_NORMAL" << "\n*******************************************************************\n";
	  scratchresfile.close();

	// ______Tidy up______
	  INFO << "Individual estimated translations (#, l, p, corr, offl, offp):";
	  INFO.print();
	  for (int k =0; k<Nwin; k++)
		{
		INFO << k << " \t" << Minlminp(k,0) << " \t" << Minlminp(k,1) << " \t" << Result(k,0) << " \t" << Result(k,1) << " \t" << Result(k,2) << " \t";
		INFO.print();
		}

	  PROGRESS << "Estimated overall translation (l,p): " << offsetLines << ", " << offsetPixels << " (used)" << ends;
	  PROGRESS.print();
	  INFO << "Coherence NaN values are disregarded in the analysis."; //MA see getoffset
	  INFO.print();
	  //PROGRESS.print("SIMAMP coregistration based on correlation finished.");
	  PROGRESS.print("MASTER TIMING Error estimation finished.");
	  } // END mtiming_correlfft

//***************************************************************
// * corrfft                                                      *
// *                                                              *
// * coherence in spectral domain by fft's                        *
// *  uses extension with zeros                                   *
// *  pixel level                                                 *
// *                                                              *
// * input:                                                       *
// *  - Master                                                    *
// *  - Mask (size Master)                                        *
// * output:                                                      *
// *  - coherence value                                           *
// *  - updated offsetL, P                                        *
// *    positive offsetL: Mask is shifted up                      *
// *    positive offsetP: Mask is shifted left                    *
// *                                                              *
// *    Bert Kampes, 16-Feb-1999                                  *
// * note: this routine can be speeded up by removing matrices    *
// * powerma* and by using *= instead of dotmult                  *
// * for now this is not done because it requires only little time*
// * note also that for coarse coregistration division by powers  *
// * is not really required, cross products are good enough.      *
// *    Bert Kampes, 18-Oct-1999                                  *
// ***************************************************************
//
//real4 corrfft(
//         const matrix<real4> &Master,                   // magnitude image
//         const matrix<real4> &Mask,                     // magnitude image
//         real4 &offsetL,                                // updated
//         real4 &offsetP)                                // updated 
//  {
//  TRACE_FUNCTION("corrfft (BK 18-Oct-1999)");
//  // ______ Internal variables ______
//  const int32 L     = Master.lines();
//  const int32 P     = Master.pixels();
//  const int32 twoL  = 2*L;
//  const int32 twoP  = 2*P;
//  const int32 halfL = L/2;
//  const int32 halfP = P/2;
//
//  // ______ Check input ______
//  if (L != Mask.lines() || P != Mask.pixels())
//    {
//    PRINT_ERROR("Mask, Master not same size.")
//    throw(input_error);
//    }
//  if (!(ispower2(L) || ispower2(P)))
//    {
//    PRINT_ERROR("Mask, Master size not power of 2.")
//    throw(input_error);
//    }
//
//  // ======Compute powers for submatrices======
//  register int32 i;
//  register int32 j;
//  const complr4 ONE(1.0);
//  matrix<complr4> Master2(twoL,twoP);// init 0
//  matrix<complr4> Mask2(twoL,twoP);  // init 0
//  matrix<complr4> blok2(twoL,twoP);  // init 0
//
//  // ====== Powers, misuse Master2, blok2 ======
//  // ______ First powers to use one matrix less; 3 is minimum ______
//  const real4 meanMaster = mean(Master);
//  const real4 meanMask   = mean(Mask);
//  for (i=0; i<L; i++)
//    {
//    for (j=0; j<P; j++)
//      {
//      Master2(i,j)   = complr4(sqr(Master(i,j)-meanMaster));// intensity image
//      Mask2(i+L,j+P) = complr4(sqr(Mask(i,j)-meanMask));// only real part
//      blok2(i+halfL,j+halfP) = ONE;                     // only real part
//      }
//    }
//  fft2d(Master2);
//  fft2d(Mask2);
//  fft2d(blok2);
//
//  // ______ new way: test this ______
//  blok2.conj();                // conjugated in blok2
//  Master2 *= blok2;
//  Master2.conj();              // M2=original(b2)*conj(m2)
//  ifft2d(Master2);             // norms master in Master2
//  Mask2   *= blok2;            // powers in spectral domain
//  ifft2d(Mask2);               // norms slave in Mask2
//
//  // ______ Use real(block2) to store sqrt(norms1*norms2) ______
//  for (i=0; i<=L; i++)         // all shifts
//    for (j=0; j<=P; j++)       // all shifts
//      blok2(i,j) = complr4(sqrt(real(Master2(i,j))*real(Mask2(i,j))));
//
//  // ====== Cross products covariance master/slave ======
//  Master2.clean();             // init 0
//  Mask2.clean();               // init 0
//  for (i=0;i<L;i++)
//    {
//    for (j=0;j<P;j++)
//      {
//      Master2(i,j)           = complr4(Master(i,j)-meanMaster); // only real
//      Mask2(i+halfL,j+halfP) = complr4(Mask(i,j)-meanMask);     // part mag. image
//      }
//    }
//
//  // ======FFT's of master/mask======
//  // padded with N zeros to prevent periodical convolution
//  fft2d(Master2);
//  fft2d(Mask2);
//
//  // ______ Store in Mask2 crossproducts in spectral/space domain ______
//  Master2.conj();
//  Mask2 *= Master2;            // corr. by zero padding
//  ifft2d(Mask2);               // space domain (real only)
//
//  // ====== Correlation in space domain for all shifts [-N/2,N/2] ======
//  //real4 coher;
//  real4 maxcoher = -999.0;
//  for (i=0; i<=L; i++)         // all shifts
//    {
//    for (j=0; j<=P; j++)       // all shifts
//      {
//      const real4 coher = real(Mask2(i,j)) / real(blok2(i,j));
//      if (coher > maxcoher)
//        {
//        maxcoher = coher;
//        offsetL  = -halfL + i; // update by reference
//        offsetP  = -halfP + j; // update by reference
//        }
//      }
//    }
//  return maxcoher;
//  } // END corrfft
//


//***************************************************************
// *    distributepoints                                          *
// *                                                              *
// * Returns matrix with distributed points in of input.          *
// * First window at (win.linelo,win.pixlo),                      *
// *  divided over wl lines,                                      *
// *  with dp distance in pixel direction.                        *
// *                                                              *
// * input:                                                       *
// *  - number of windows                                         *
// *  - window which should be divided                            *
// * output:                                                      *
// *  - matrix <uint> (NW,3) =(l,p, flagfromdisk==0)              *
// *                                                              *
// *    Bert Kampes, 21-Jan-1999                                  *
// ***************************************************************

	// ______ Corr by fft ______
	// superseded by coherencefft with factor=1
	//real4 corrfft( 
	//      const matrix<real4>     &magnitudeMaster,
	//      const matrix<real4>     &magnitudeMask,
	//      real4                   &offsetL,
	//      real4                   &offsetP);


	// ______ Distribute nW windows over win ______
	public static matrix<Integer> distributepoints(float nW, window win)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"distributepoints (BK 21-Jan-1999)"<<ends;
		  TRACE.print();
	  }
	  float lines = win.linehi - win.linelo + 1;
	  float pixels = win.pixhi - win.pixlo + 1;

	  int numw = nW;
	  matrix<Integer> Result = new matrix(numw,(int)3);
	  // ______ Distribution for dl=dp ______
	  float wp = Math.sqrt(nW/(lines/pixels)); // wl: #windows in line direction
	  float wl = nW / wp; // wp: #windows in pixel direction
	  if (wl < wp) // switch wl,wp : later back
		wl = wp;
	  int wlint = rint(wl); // round largest
	  float deltal = (lines-1) / ((float)(wlint-1));
	  int totp = pixels *wlint;
	  float deltap = ((float)(totp-1)) / ((float)(nW-1));
	  float p = -deltap;
	  float l = 0.;
	  int lcnt = 0;
//C++ TO JAVA CONVERTER NOTE: 'register' variable declarations are not supported in Java:
//ORIGINAL LINE: register int i;
	  int i;
	  for (i =0; i<nW; i++)
		{
		p += deltap;
		while (rint(p)>=pixels) // round
		  {
		  p -= pixels;
		  lcnt++;
		  }
		l = lcnt * deltal;
		Result(i,0) = (int)(rint(l));
		Result(i,1) = (int)(rint(p));
		}

	  // ______ Correct distribution to window ______
	  for (i =0; i<nW; i++)
		{
		Result(i,0) += win.linelo;
		Result(i,1) += win.pixlo;
		}

	  return Result;
	  } // END distributepoints

//***************************************************************
// *    getoffset                                                 *
// *                                                              *
// * Returns offset in line and pixel direction                   *
// *  based on matrix with estimated offests                      *
// *  by correlation                                              *
// * Checks on consistency, THRESHOLD 0.4 for correlation         *
// *                                                              *
// * input:                                                       *
// *  - matrix<real4> with corr,offsets(l,p)                      *
// * output:                                                      *
// *  - (offL,offP)                                               *
// *                                                              *
// * See also Documentation page 6.                               *
// *                                                              *
// *    Bert Kampes, 21-Jan-1999                                  *
// ***************************************************************


	// ______ Estimate offset based on consistency ______
	public static void getoffset(matrix<Float> Result, RefObject<Integer> offsetLines, RefObject<Integer> offsetPixels)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"getoffset (BK 21-Jan-1999)"<<ends;
		  TRACE.print();
	  }
	  if (Result.pixels() != 3)
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "code 901: input not 3 width";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}
	  // --- First sort estimated offsets on coherence ascending! ---
	  DEBUG.print("sorting on coherence.");
	  //DEBUG.print("unsorted input matrix:");
	  //Result.showdata();
	  matrix<Float> sortResult = -Result;
	  mysort2(sortResult); // sort matrix on first column (coh)
	  sortResult = -sortResult;
	  //DEBUG.print("sorted matrix:");
	  //sortResult.showdata();

	  // --- Set offset to highest coherence estimate ---
	  offsetLines.argvalue = (int)(rint(sortResult(0,1))); //rounds negative too
	  offsetPixels.argvalue = (int)(rint(sortResult(0,2))); //rounds negative too
	  final int nW = sortResult.lines();
	  int nWNANrm = sortResult.lines(); //MA added for removal of -999 values
	  if (nW ==1)
		  return;

	  // --- Threshold on coherence ---
	  float var_coh = 0.0;
	  float mean_coh = 0.0;
	  for (int i =0; i<nW; i++)
	  { //MA fix to ignore -999 values from statistics
		if (sortResult(i,0) == -999)
		{
		nWNANrm = nWNANrm - 1;
		continue;
		}
		 mean_coh+=sortResult(i,0);
	  }
	  //mean_coh /= real4(nW);
	  mean_coh /= (float)nWNANrm;
	  for (int i =0; i<nW; i++)
	  { //MA fix to ignore -999 values from statistics
	   if (sortResult(i,0) == -999)
		   continue;
		var_coh +=sqr(sortResult(i,0)-mean_coh);
	  }
	  //var_coh /= real4(nW-1);
	  var_coh /= (float)(nWNANrm-1);
	  INFO << "Mean coherence at estimated positions: " << mean_coh;
	  INFO.print();
	  final float std_coh = Math.sqrt(var_coh);
	  INFO << "Standard deviation coherence:          " << std_coh;
	  INFO.print();
	  final float thresh_coh = mean_coh;
	  INFO << "Using as threshold:                    " << thresh_coh;
	  INFO.print();
	  int cnt = 1; // estimates above threshold
	  mean_coh = sortResult(0,0); // mean above threshold
	  INFO.print("Using following data to determine coarse image offset:");
	  INFO.print("coherence    offset_L    offset_P");
	  INFO.print("------------------------------------------------------");
	  INFO << sortResult(0,0) << "      " << sortResult(0,1) << "        " << sortResult(0,2);
	  INFO.print();
	  for (int i =1; i<nW; i++)
		{
		if (sortResult(i,0)>=thresh_coh)
		  {
		  cnt++;
		  mean_coh += sortResult(i,0);
		  offsetLines.argvalue += (int)(rint(sortResult(i,1))); // round
		  offsetPixels.argvalue += (int)(rint(sortResult(i,2))); // round
		  INFO << sortResult(i,0) << "      " << sortResult(i,1) << "        " << sortResult(i,2);
		  INFO.print();
		  }
		}

	  // ___ Report stats ___
	  if (cnt > 1)
		{
		mean_coh /= (float)cnt;
		final float meanL = offsetLines.argvalue; // float mean
		final float meanP = offsetPixels.argvalue; // float mean
		offsetLines.argvalue = (int)(rint((double)offsetLines.argvalue/(double)cnt)); // round
		offsetPixels.argvalue = (int)(rint((double)offsetPixels.argvalue/(double)cnt)); // round
		float var_L = 0.0;
		float var_P = 0.0;
		for (int i =0; i<cnt; i++)
			var_L+=sqr(sortResult(i,1)-meanL);
		for (int i =0; i<cnt; i++)
			var_P+=sqr(sortResult(i,2)-meanP);
		var_L /= (float)(cnt-1);
		var_P /= (float)(cnt-1);
		INFO << "Standard deviation offset L = " << Math.sqrt(var_L);
		INFO.print();
		INFO << "Standard deviation offset P = " << Math.sqrt(var_P);
		INFO.print();
		if (Math.sqrt(var_L)>6.0 || Math.sqrt(var_P)>6.0)
		  WARNING.print("Check estimated offset coarse corr: it seems unreliable.");
		}

	  // ___ Warn if appropriate ___
	  if (mean_coh < 0.2)
		{
		WARNING.print("getoffset: mean coherence of estimates used < 0.2");
		WARNING.print("(please check bottom of LOGFILE to see if offset is OK)");
		}
	  if (nW < 6)
		{
		WARNING.print("getoffset: number of windows to estimate offset < 6");
		WARNING.print("(please check bottom of LOGFILE to see if offset is OK)");
		}

	//
	//  int32 cnt;
	//  int32 valueL;
	//  int32 valueP;
	//  real4 correl;
	//  int32 highestcnt    = 0;
	//  real4 highestcorrel = 0.0;
	//  for (i=0; i<nW; i++)
	//    {
	//    valueL = int32(Result(i,0)+0.5);
	//    valueP = int32(Result(i,1)+0.5);
	//    correl = Result(i,2);
	//    if (correl > highestcorrel) 
	//      highestcorrel = correl;
	//    cnt = 0;
	//    for (j=0; j<nW; j++)
	//      {
	//      if (abs(Result(j,0) - valueL) < 2  &&  
	//          abs(Result(j,1) - valueP) < 2)
	//        cnt++;
	//      }
	//    if (cnt > highestcnt)
	//      {
	//      highestcnt   = cnt;
	//      offsetLines  = valueL;                    // Return offsetLines
	//      offsetPixels = valueP;                    // Return offsetPixels
	//      }
	//    }
	//
	//  // ______ Check result ______
	//  real4 THRESHOLD = 0.3;
	//  if (nW < 6)
	//    {
	//    WARNING.print("getoffset: number of windows to estimate offset < 6");
	//    WARNING.print("(please check bottom of LOGFILE)");
	//    }
	//  if (highestcnt < 0.2*nW)
	//    {
	//    WARNING.print("getoffset: estimated offset not consistent with other estimates.");
	//    WARNING.print("(check bottom of LOGFILE)");
	//    }
	//  if (highestcorrel < THRESHOLD)
	//    {
	//    WARNING << "getoffset: estimated translation has correlation of: "
	//         << highestcorrel;
	//    WARNING.print();
	//    WARNING.print("(please check bottom of LOGFILE)");
	//    }
	//
	  } // END getoffset

//***************************************************************
// *    getmodeoffset                                             *
// *                                                              *
// * Returns offset in line and pixel direction                   *
// *  based on matrix with estimated offests                      *
// *  by correlation                                              *
// * Checks on consistency, THRESHOLD 0.4 for correlation         *
// *                                                              *
// * input:                                                       *
// *  - matrix<real4> with corr,offsets(l,p)                      *
// * output:                                                      *
// *  - (offL,offP)                                               *
// *                                                              *
// * See also Documentation page 6.                               *
// *                                                              *
// *    Bert Kampes,   21-Jan-1999 (getoffset)                    *
// *    Mahmut Arikan, 09-Dec-2008                                *
// ***************************************************************

	// ______ Estimate offset based on consistency ______
	public static void getmodeoffset(matrix<Float> Result, RefObject<Integer> offsetLines, RefObject<Integer> offsetPixels)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"getmodeoffset (MA 09-Dec-2008)"<<ends;
		  TRACE.print();
	  }
	  if (Result.pixels() != 3)
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "code 901: input not 3 width";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}

	  // --- First sort estimated offsets on coherence ascending! ---
	  DEBUG.print("sorting on coherence.");
	  //DEBUG.print("unsorted input matrix:");
	  //Result.showdata();
	//
	//  Result(0,0)=0.1   ; Result(0,1)=3 ; Result(0,2)=2  ;   
	//  Result(1,0)=0.2   ; Result(1,1)=1 ; Result(1,2)=4  ;   
	//  Result(2,0)=0.3   ; Result(2,1)=4 ; Result(2,2)=1  ;   
	//  Result(3,0)=0.2   ; Result(3,1)=1 ; Result(3,2)=3  ;   
	//  Result(4,0)=0.2   ; Result(4,1)=3 ; Result(4,2)=1  ;   
	//  Result(5,0)=0.3   ; Result(5,1)=1 ; Result(5,2)=-1  ;   
	//  Result(6,0)=0.2   ; Result(6,1)=4 ; Result(6,2)=0  ;   
	//  Result(7,0)=0.1   ; Result(7,1)=1 ; Result(7,2)=-2  ;   
	//
	//Result.showdata();
	//cerr << endl;
	//
	  matrix<Float> sortResult = -Result;
	  mysort2(sortResult); // sort matrix on first column (coh)
											// sorts ascending
	//  sortResult.showdata(); cout << endl;
	//  mysort2selcol(sortResult, 1);
	  sortResult = -sortResult; // max coh at top
	  //DEBUG.print("sorted matrix:");
	  //sortResult.showdata();

	  // --- Set offset to highest coherence estimate ---
	  offsetLines.argvalue = (int)(rint(sortResult(0,1))); // rounds negative too, was -999
	  offsetPixels.argvalue = (int)(rint(sortResult(0,2))); // rounds negative too
													// [ why set to highes coherence mean 
													// loop index could start from i==0.]

	  // ______ Remove window offests with -999 (NaN) coherence values _____ 
	  // added by [MA]
	  final int nW = sortResult.lines(); // Number of windows
	  int nWNANrm = nW; // Number of windows without NAN values
	  if (nW ==1)
		  return;

	  // --- Threshold on coherence ---
	  float var_coh = 0.0;
	  float mean_coh = 0.0;
	  for (int i =0; i<nW; i++) // [MA] fix to ignore -999 values from statistics
		{
		 if (sortResult(i,0) == -999) // if NaN
		   {
			nWNANrm -= 1; // determine number of windows without NaN
			continue;
		   }
		 mean_coh+=sortResult(i,0);
		}
	  //mean_coh /= real4(nW);
	  mean_coh /= (float)nWNANrm; // mean coherence

	  for (int i =0; i<nW; i++) // [MA fix to ignore -999 values from statistics
		{
		 if (sortResult(i,0) == -999)
			 continue;
		 var_coh +=sqr(sortResult(i,0)-mean_coh);
		}
	  //var_coh /= real4(nW-1);
	  var_coh /= (float)(nWNANrm-1); // mean variance

	  INFO << "Mean coherence at estimated positions: " << mean_coh;
	  INFO.print();
	  final float std_coh = Math.sqrt(var_coh);
	  INFO << "Standard deviation coherence:          " << std_coh;
	  INFO.print();

	  // ______ Statistics about threshold ______ 
	  final float thresh_coh = mean_coh;
	  INFO << "Using as threshold:                    " << thresh_coh;
	  INFO.print();

	  DEBUG.print("Using following data to determine coarse image offset:");
	  DEBUG.print("coherence    offset_L    offset_P");
	  DEBUG.print("------------------------------------------------------");
	  DEBUG << sortResult(0,0) << "      " << sortResult(0,1) << "        " << sortResult(0,2); // print the line w/ max. coh.
	  DEBUG.print();

	  int cnt = 1; // estimates above threshold
	  mean_coh = sortResult(0,0); // new mean above threshold
	  for (register int i =1; i<nW; i++)
		{
		if (sortResult(i,0)>=thresh_coh)
		  {
		  cnt++;
		  mean_coh += sortResult(i,0);
		  offsetLines.argvalue += (int)(rint(sortResult(i,1))); // round
		  offsetPixels.argvalue += (int)(rint(sortResult(i,2))); // round
		  DEBUG << sortResult(i,0) << "      " << sortResult(i,1) << "        " << sortResult(i,2);
		  DEBUG.print();

		  } // values above threshold
		} // end loop and print

	  // ___ Report stats ___
	  if (cnt > 1)
		{
		mean_coh /= (float)cnt;
		final float meanL = offsetLines.argvalue; // float mean
		final float meanP = offsetPixels.argvalue; // float mean
		offsetLines.argvalue = (int)(rint((double)offsetLines.argvalue/(double)cnt)); // round
		offsetPixels.argvalue = (int)(rint((double)offsetPixels.argvalue/(double)cnt)); // round
		float var_L = 0.0;
		float var_P = 0.0;
		for (register int i =0; i<cnt; i++)
			var_L+=sqr(sortResult(i,1)-meanL);
		for (register int i =0; i<cnt; i++)
			var_P+=sqr(sortResult(i,2)-meanP);
		var_L /= (float)(cnt-1);
		var_P /= (float)(cnt-1);
		INFO << "Standard deviation offset L = " << Math.sqrt(var_L);
		INFO.print();
		INFO << "Standard deviation offset P = " << Math.sqrt(var_P);
		INFO.print();
		if (Math.sqrt(var_L)>6.0 || Math.sqrt(var_P)>6.0)
		  WARNING.print("Check estimated offset coarse corr: it seems unreliable.");
		}

	  INFO << "Estimated overall mean translation (l,p): " << offsetLines.argvalue << ", " << offsetPixels.argvalue << " (not used)" << ends;
	  INFO.print();

	// ofstream scratchlogfile("scratchlogmtiming", ios::out | ios::trunc);
	// bk_assert(scratchlogfile,"mtiming_correl: scratchlogmtiming",__FILE__,__LINE__);
	// // INFO.rdbuf(scratchlogfile.rdbuf());
	// INFO.print();


	   // _____ Mode of offsets _____  [MA]
	  PROGRESS.print("getmodeoffset: Start mode analysis ");
	  PROGRESS << "Using as threshold:  " << thresh_coh << " and checking for mode value";
	  PROGRESS.print();
	  INFO.print("Using following data to determine coarse image offset:");
	  INFO.print("avg. coh    offset_L    offset_P  occurence  index");
	  INFO.print("------------------------------------------------------");

	  mysort231(sortResult); // re-sort on 2nd, 3rd than 1st column
	  // sortResult.showdata();       
	  int mode_val = 0; // mode count, mode index
	  int mode_idx = -1;
	  int evenmode_val = 0; // check for equal values of mode
	  int nEven = 0;
	  int L =NaN; // Line, Pixel, frequency
	  int P =NaN;
	  int offset_freq =0;
	  float offset_mcoh =0.0; // avg. coherence for each set of offsets
	  for (register int i =0; i<nW; i++) // Major reason of this main loop is individual stdout request.
		{
		if (sortResult(i,0)>=thresh_coh)
		  {
		  // _____ frequency of offsets _____  [MA]
		  if (L != (int)(rint(sortResult(i,1))) || P != (int)(rint(sortResult(i,2)))) // the same offset multiple times
			{
			  L =(int)(rint(sortResult(i,1))); // get initial values
			  P =(int)(rint(sortResult(i,2)));
			}
		  else
			 {
			  continue; // L, P equal to previous values then skip counting
						   // since matrix is sorted on L,P
			 }
		  offset_freq =0; // reset
		  offset_mcoh =0;
		  for (register int j =0; j<nW; j++) // scan data for occurences of an offset
			{ // for all offsets
			 if (L == (int)(rint(sortResult(j,1))) && P == (int)(rint(sortResult(j,2))))
			   {
				 offset_freq++;
				 offset_mcoh += sortResult(j,0); // for decission on even mode values
			   } // at different L,P pair.
			} // end scan data

		  if (offset_freq > mode_val)
			{
			 mode_val =offset_freq;
			 mode_idx =i; // index of mode value
													  // in magfft if you correlate the same
													  // slc patches. index get a value other than
													  // 1. that's okay when all offset are zero.
			}
		  else if (mode_val == offset_freq)
			{
			 if (evenmode_val != offset_freq) // initialize with one
				 nEven =1;
			 evenmode_val =offset_freq;
			 nEven++;
			}

		  offset_mcoh /= (float)offset_freq;

		  // _____ for each offset pair above threshold list frequency _____
		  INFO << offset_mcoh << "\t " << L << "\t   " << P << "\t\t" << offset_freq << "\t " << mode_idx;
		  INFO.print();

		  } // above threshold
		} // end mode

		// _____ Even occurence check _____
		if (mode_val == evenmode_val) // there are even values of mode.
		  {
			WARNING << "There are " << nEven << " offset pairs which has equal mode values are equal.";
			WARNING.print();
			WARNING << "Check offset results and logs, and increase the number and/or the size of the correlation windows.";
			WARNING.print();
		  }


	  offsetLines.argvalue = (int)(rint(sortResult(mode_idx,1))); // update mode offsets
	  offsetPixels.argvalue = (int)(rint(sortResult(mode_idx,2)));
	  PROGRESS.print("getmodeoffset: End of mode analysis ");

	  // ___ Warn if appropriate ___
	  if (mean_coh < 0.2)
		{
		WARNING.print("getmodeoffset: mean coherence of estimates used < 0.2");
		WARNING.print("(please check bottom of LOGFILE to see if offset is OK)");
		}
	  if (mode_val == 1)
		{
		WARNING.print("getmodeoffset: all the offset occurence == 1. There is no mode value. ");
		WARNING.print("(please check bottom of LOGFILE to see if offset is OK or change window size.)");
		}
	  if (nW < 6)
		{
		WARNING.print("getmodeoffset: number of windows to estimate offset < 6");
		WARNING.print("(please check bottom of LOGFILE to see if offset is OK)");
		}


	  } // END getmodeoffset

//***************************************************************
// *    finecoreg                                                 *
// *                                                              *
// * computes translation of slave w.r.t. master                  *
// * slave(some point) = master(same point) + trans(l,p) =>       *
// *  trans = slavecoordinates - mastercoordinates                *
// * in NWIN windows.                                             *
// * Then solves polynomial for best transformation to master     *
// * with coregpm routine/step                                    *
// *                                                              *
// * input:                                                       *
// *  -                                                           *
// * output:                                                      *
// *  -                                                           *
// *    Bert Kampes, 12-Dec-1998                                  *
// * distribute points can be with input file as well, besides    *
// * letting Doris randomly distribute npoints.                   *
// *    BK 29-Oct-99                                              *
// ***************************************************************


	// ______ Fine coregistration ______
	public static void finecoreg(input_fine fineinput, slcimage minfo, slcimage sinfo)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"finecoreg (BK 29-Oct-99)"<<ends;
		  TRACE.print();
	  }
	  String dummyline = new String(new char[ONE27]);
	  //const uint Mfilelines   = minfo.currentwindow.lines();
	  //const uint Sfilelines   = sinfo.currentwindow.lines();
	  final int Nwin = fineinput.Nwin; // n windows, from file or random
	  int NwinNANrm = fineinput.Nwin; // [MA] number of windows w/o NaN
	  final int initoffsetL = fineinput.initoffsetL; // initial offset
	  final int initoffsetP = fineinput.initoffsetP; // initial offset
	  int MasksizeL = fineinput.MasksizeL; // size of correlation window
	  int MasksizeP = fineinput.MasksizeP; // size of correlation window
	  int AccL = fineinput.AccL; // size of small chip
	  int AccP = fineinput.AccP; // size of small chip
	  final int OVS = fineinput.osfactor; // factor
	  boolean pointsrandom = true;
	  if (specified(fineinput.ifpositions)) // filename specified
		pointsrandom = false; // only use these points

	  // ______Correct sizes if in space domain______
	  if (fineinput.method == fc_magspace || fineinput.method == fc_cmplxspace)
		{
		INFO.print("Adapting size of window for space method");
		MasksizeL += 2 *fineinput.AccL;
		MasksizeP += 2 *fineinput.AccP;
		}

	  // ______Corners of slave in master system______
	  // ______offset = [A](slave system) - [A](master system)______
	  final int sl0 = sinfo.currentwindow.linelo - initoffsetL;
	  final int slN = sinfo.currentwindow.linehi - initoffsetL;
	  final int sp0 = sinfo.currentwindow.pixlo - initoffsetP;
	  final int spN = sinfo.currentwindow.pixhi - initoffsetP;

	  // ______Corners of useful overlap master,slave in master system______
	  final int BORDER = 20; // make slightly smaller
	  final int l0 = max((int)minfo.currentwindow.linelo,sl0) + BORDER;
	  final int lN = min((int)minfo.currentwindow.linehi,slN) - MasksizeL - BORDER;
	  final int p0 = max((int)minfo.currentwindow.pixlo,sp0) + BORDER;
	  final int pN = min((int)minfo.currentwindow.pixhi,spN) - MasksizeP - BORDER;
	  final window overlap = new window(l0,lN,p0,pN);

	  // ______ Distribute Nwin points over window, or read from file ______
	  // ______ Minlminp(i,0): line, (i,1): pixel, (i,2) flagfromdisk ______
	  matrix<Integer> Minlminp;
	  if (pointsrandom) // no filename specified
		{
		Minlminp = distributepoints((float)Nwin, overlap);
		}

	  else // read in points (center of windows) from file
		{
		Minlminp.resize(Nwin,3);
		ifstream ifpos = new ifstream(fineinput.ifpositions, ios.in);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		bk_assert(ifpos,fineinput.ifpositions,__FILE__,__LINE__);
		int ll;
		int pp;
		for (int i =0; i<Nwin; ++i)
		  {
		  ifpos >> ll >> pp;
		  Minlminp(i,0) = (int)(ll - 0.5 *MasksizeL); // correct for lower left corner
		  Minlminp(i,1) = (int)(pp - 0.5 *MasksizeP); // correct for lower left corner
		  Minlminp(i,2) = (int)1; // flag from file
		  ifpos.getline(dummyline,ONE27,'\n'); // goto next line.
		  }
		ifpos.close();
		// ______ Check last point for possible EOL after last position in file ______
		if (Minlminp(Nwin-1,0) == Minlminp(Nwin-2,0) && Minlminp(Nwin-1,1) == Minlminp(Nwin-2,1))
		  {
		  Minlminp(Nwin-1,0) = (int)(0.5*(lN + l0) + 27); // random
		  Minlminp(Nwin-1,1) = (int)(0.5*(pN + p0) + 37); // random
		  }
		// ______ Check if points are in overlap ______
		// ______ no check for uniqueness of points ______
		boolean troubleoverlap = false;
		for (int i =0; i<Nwin; ++i)
		  {
		  if (Minlminp(i,0) < l0)
			{
			troubleoverlap =true;
			WARNING << "FINE: point from file: " << i+1 << " " << Minlminp(i,0) + 0.5 *MasksizeL << " " << Minlminp(i,1) + 0.5 *MasksizeP << " outside overlap master, slave. New position: ";
			Minlminp(i,0) = l0 + l0-Minlminp(i,0);
			WARNING << Minlminp(i,0) << " " << Minlminp(i,1);
			WARNING.print();
			}
		  if (Minlminp(i,0) > lN)
			{
			troubleoverlap =true;
			WARNING << "FINE: point from file: " << i+1 << " " << Minlminp(i,0) + 0.5 *MasksizeL << " " << Minlminp(i,1) + 0.5 *MasksizeP << " outside overlap master, slave. New position: ";
			Minlminp(i,0) = lN + lN-Minlminp(i,0);
			WARNING << Minlminp(i,0) << " " << Minlminp(i,1);
			WARNING.print();
			}
		  if (Minlminp(i,1) < p0)
			{
			troubleoverlap =true;
			WARNING << "FINE: point from file: " << i+1 << " " << Minlminp(i,0) + 0.5 *MasksizeL << " " << Minlminp(i,1) + 0.5 *MasksizeP << " outside overlap master, slave. New position: ";
			Minlminp(i,1) = p0 + p0-Minlminp(i,1);
			WARNING << Minlminp(i,0) << " " << Minlminp(i,1);
			WARNING.print();
			}
		  if (Minlminp(i,1) > pN)
			{
			troubleoverlap =true;
			WARNING << "FINE: point from file: " << i+1 << " " << Minlminp(i,0) + 0.5 *MasksizeL << " " << Minlminp(i,1) + 0.5 *MasksizeP << " outside overlap master, slave. New position: ";
			Minlminp(i,1) = pN + pN-Minlminp(i,1);
			WARNING << Minlminp(i,0) << " " << Minlminp(i,1);
			WARNING.print();
			}
		  }
		if (troubleoverlap) // give some additional info
		  {
		  WARNING << "FINE: there were points from file outside overlap (l0,lN,p0,pN): " << l0 << " " << lN << " " << p0 << " " << pN;
		  WARNING.print();
		  }
		}

	  // ______Compute coherence of these points______
	  matrix<complex<Float>> Master;
	  matrix<complex<Float>> Mask;
	  matrix<Float> Result = new matrix(Nwin,3); // R(i,0):delta l;
											// R(i,1):delta p; R(i,2):correl

	  // ______ Progress message ______
	  int tenpercent = rint(Nwin/10.0);
	  if (tenpercent ==0)
		  tenpercent = 1000;
	  int percent = 0;

	  // ====== Compute for all locations ======
	  for (int i =0;i<Nwin;i++)
		{
		// ______ Give progress message ______
		if (i%tenpercent ==0)
		  {
		  PROGRESS << "FINE: " << setw(3) << percent << "%";
		  PROGRESS.print();
		  percent += 10;
		  }

		// ______Minlminp (lower left corners) of window in master system______
		final int minMwinL = Minlminp(i,0);
		final int minMwinP = Minlminp(i,1);
		DEBUG.print(" ");
		DEBUG << "Window: " << i << " [" << minMwinL << ", " << minMwinP << "]";
		DEBUG.print();
		window master = new window(minMwinL, minMwinL+MasksizeL-1, minMwinP, minMwinP+MasksizeP-1); // size=masksize
		// ______Same points in slave system (disk)______
		window mask = new window(minMwinL+initoffsetL, minMwinL+initoffsetL+MasksizeL-1, minMwinP+initoffsetP, minMwinP+initoffsetP+MasksizeP-1); // size=masksize
		// ______Read windows from files______
		Master = minfo.readdata(master);
		Mask = sinfo.readdata(mask);

		// ______Coherence______
		// ______update offsetL/P______
		float offsetL;
		float offsetP;
		float coheren;
		switch (fineinput.method)
		  {
		  //case fc_cmplxfft:
		  //WARNING("THIS METHOD IS NOT OK YET, I RECOMMEND MAGNITUDE.");
		  //coheren = coherencefft(fineinput, Master, Mask, offsetL, offsetP);
		  //break;
		  //case fc_cmplxspace:
		  //WARNING("THIS METHOD IS NOT OK YET, I RECOMMEND MAGNITUDE.");
		  //coheren = coherencespace(fineinput, Master, Mask, offsetL, offsetP);
		  //break;
		  case fc_magfft: // fast: oversample coherence
			{
			//coheren     = coherencefft(Master, Mask, OVS, AccL, AccP,
			//                           offsetL, offsetP);// returned
			if (AccL > MasksizeL/2) // [MA] fix for Acc being half of Masksize at max
			  {
			   AccL = MasksizeL/2;
			   WARNING << "FINE: AccL for magfft can be half of the window size at max, changing to " << AccL;
			   WARNING.print();
			  }
			else if (AccP > MasksizeP/2)
			  {
			   AccP = MasksizeP/2;
			   WARNING << "FINE: AccP for magfft can be half of the window size at max, changing to " << AccP;
			   WARNING.print();
			  }

			RefObject<Float> TempRefObject = new RefObject<Float>(offsetL);
			RefObject<Float> TempRefObject2 = new RefObject<Float>(offsetP);
			coheren = crosscorrelate(Master, Mask, OVS, AccL, AccP, TempRefObject, TempRefObject2); // returned
			offsetL = TempRefObject.argvalue;
			offsetP = TempRefObject2.argvalue;
			break;
			}
		  // ====== New method (BK 13 Aug 2005) ======
		  // ====== This should work for ERS/N1; different PRFs ======
		  case fc_oversample: // slow (better): oversample complex data first
			{

			if (AccL > MasksizeL/2) // [MA] fix for Acc being half of Masksize at max
			  {
			   AccL = MasksizeL/2;
			   WARNING << "FINE: AccL for magfft can be half of the window size at max, changing to " << AccL;
			   WARNING.print();
			  }
			else if (AccP > MasksizeP/2)
			  {
			   AccP = MasksizeP/2;
			   WARNING << "FINE: AccP for magfft can be half of the window size at max, changing to " << AccP;
			   WARNING.print();
			  }

			// ______ Oversample complex chips by factor two ______
			// ______ neg.shift input shifts to -> 0
			DEBUG.print("Centering azimuth spectrum patches around 0");
			final float m_pixlo = master.pixlo; // neg.shift -> 0
			final float s_pixlo = mask.pixlo; // neg.shift -> 0
			shiftazispectrum(Master,minfo,-m_pixlo); // shift from fDC to zero
			shiftazispectrum(Mask, sinfo,-s_pixlo); // shift from fDC to zero
			DEBUG.print("Oversampling patches with factor two using zero padding");
			final matrix<complex<Float>> m_ovs_chip = oversample(Master,2,2);
			final matrix<complex<Float>> s_ovs_chip = oversample(Mask, 2,2);
			// ______ Peak in cross-corr of magnitude of ovs data ______
			DEBUG.print("Cross-correlating magnitude of ovs patches");
			DEBUG.print("(no need to shift spectrum back)"); // (else account for ovs..)
			//coheren = coherencefft(m_ovs_chip, s_ovs_chip, 
			//                       OVS/2, 2*AccL, 2*AccP, 
			//                       offsetL,offsetP);
			RefObject<Float> TempRefObject3 = new RefObject<Float>(offsetL);
			RefObject<Float> TempRefObject4 = new RefObject<Float>(offsetP);
			coheren = crosscorrelate(m_ovs_chip, s_ovs_chip, OVS/2, 2 *AccL, 2 *AccP, TempRefObject3, TempRefObject4);
			offsetL = TempRefObject3.argvalue;
			offsetP = TempRefObject4.argvalue;
			offsetL /= 2.0; // orig data oversampled by factor 2
			offsetP /= 2.0; // orig data oversampled by factor 2
			break;
			}
		  case fc_magspace:
			RefObject<Float> TempRefObject5 = new RefObject<Float>(offsetL);
			RefObject<Float> TempRefObject6 = new RefObject<Float>(offsetP);
			coheren = coherencespace(fineinput, Master, Mask, TempRefObject5, TempRefObject6);
			offsetL = TempRefObject5.argvalue;
			offsetP = TempRefObject6.argvalue;
			break;
		  default:
			{
				ERROR.terminate();
				String cp_s = new String(new char[256]);
				cp_s = "unknown method for fine coregistration.";
				ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
				ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
				ERROR.print();
			}
			throw(unhandled_case_error);
		  } // switch method
		  Result(i,0) = initoffsetL + offsetL;
		  Result(i,1) = initoffsetP + offsetP;
		  Result(i,2) = coheren;
		  INFO << "Fine offset between small patches:   " << Result(i,0) << ", " << Result(i,1) << " (coh="<<coheren<<")";
		  INFO.print();
		} // for nwin


	  // ______ Position approx. with respect to center of window ______
	  // ______ correct position array for center instead of lower left ______
	  for (int i =0; i<Nwin; i++)
		{
		Minlminp(i,0) += (int)(0.5 *MasksizeL);
		Minlminp(i,1) += (int)(0.5 *MasksizeP);
		}


	  // ______Write to files______
	  ofstream scratchlogfile = new ofstream("scratchlogfine", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchlogfile,"finecoreg: scratchlogfine",__FILE__,__LINE__);
	  scratchlogfile << "\n\n*******************************************************************" << "\n* FINE_COREGISTRATION" << "\n*******************************************************************" << "\nNumber of correlation windows: \t" << Nwin << "\nwindow size (l,p):             \t" << MasksizeL << ", " << MasksizeP << "\nInitial offsets:               \t" << initoffsetL << ", " << initoffsetP << "\nOversampling factor:           \t" << OVS << "\n\nNumber \tposl \tposp \toffsetl offsetp\tcorrelation\n";
	  for (int i =0;i<Nwin;i++)
		{ // MA remove NaN valued coh windows from Nwin, to be used in resfile
		if (isnan(Result(i,2)))
			NwinNANrm = NwinNANrm - 1;
		scratchlogfile << setiosflags(ios.fixed) << setiosflags(ios.showpoint) << setiosflags(ios.right) << setw(8) << setprecision(0) << i << " " << Minlminp(i,0) << " " << Minlminp(i,1) << " " << setprecision(3) << Result(i,0) << " " << Result(i,1) << " " << setprecision(2) << Result(i,2) << "\n";
		 }
	  scratchlogfile << "\n*******************************************************************\n";
	  scratchlogfile.close();

	  ofstream scratchresfile = new ofstream("scratchresfine", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile,"finecoreg: scratchresfine",__FILE__,__LINE__);

		//Changed by MA <<  Nwin
	  scratchresfile << "\n\n*******************************************************************" << "\n*_Start_" << processcontrol[(int)AnonymousEnum.pr_i_fine] << "\n*******************************************************************" << "\nInitial offsets (l,p):             \t" << initoffsetL << ", " << initoffsetP << "\nWindow_size_L_for_correlation:     \t" << MasksizeL << "\nWindow_size_P_for_correlation:     \t" << MasksizeP << "\nMax. offset that can be estimated: \t" << MasksizeL/2 << "\nPeak search ovs window (l,p):      \t" << 2 *AccL << " , " << 2 *AccP << "\nOversampling factor:               \t" << OVS << "\nNumber_of_correlation_windows:     \t" << NwinNANrm << "\nNumber \tposL \tposP \toffsetL offsetP\tcorrelation\n";
	  scratchresfile.close();

	  FILE resfile;
	  resfile =fopen("scratchresfine","a");
	  for (int i =0; i<Nwin; i++)
	   { //MA remove/skip NaN values before writing resfile.
	   if (isnan(Result(i,2)))
		   continue;
		fprintf(resfile,"%4.0f %5.0f %5.0f %# 9.2f %# 9.2f %# 6.2f\n", (float)i, (float)(Minlminp(i,0)), (float)(Minlminp(i,1)), Result(i,0), Result(i,1), Result(i,2));
	  }

	  fprintf(resfile, "\n*******************************************************************");
	  fprintf(resfile,"%s%s%s", "\n* End_", processcontrol[(int)AnonymousEnum.pr_i_fine], "_NORMAL");
	  fprintf(resfile, "\n*******************************************************************\n");

	  // ______Tidy up______
	  fclose(resfile);
	  PROGRESS.print("Fine coregistration finished.");
	  } // END finecoreg

//***************************************************************
// * coherencefft                                                 *
// *                                                              *
// * coherence in spectral domain by fft's based on magnitude     *
// *  uses extension with zeros.  returns relative shift between  *
// *  two patched and the estimated correlation.                  *
// *                                                              *
// * input:                                                       *
// *  - Master                                                    *
// *  - Mask (size Master)                                        *
// * output:                                                      *
// *  - coherence value [-1 1]                                    *
// *  - updated offsetL, P                                        *
// *    positive offsetL: Mask is shifted up                      *
// *    positive offsetP: Mask is shifted left                    *
// *                                                              *
// *    Bert Kampes, 03-Feb-1999                                  *
// * bugfix? streamlined, based on magnitude forced               *
// *    Bert Kampes, 16-Nov-1999                                  *
// * 1) should find max at pixel level, then oversample sub-pixel *
// * but it seems to be implemented strangely                     *
// * 2) oversampling should be performed on complex images with   *
// * factor 2, so to avoid aliasing of spectrum (shift azi).      *
// *    Bert Kampes, 12-Aug-2005                                  *
// ***************************************************************
		//const input_fine &fineinput,


	// ______ Correlation with FFT ______
			//const input_fine      &fineinput, 
	public static float coherencefft(matrix<complex<Float>> Master, matrix<complex<Float>> Mask, int ovsfactor, int AccL, int AccP, RefObject<Float> offsetL, RefObject<Float> offsetP)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"coherencefft (BK 16-Nov-1999)"<<ends;
		  TRACE.print();
	  }
	  // ______ Internal variables ______
	  final int L = Master.lines();
	  final int P = Master.pixels();
	  final int twoL = 2 *L;
	  final int twoP = 2 *P;
	  final int halfL = L/2;
	  final int halfP = P/2;

	  // ______ Check input ______
	  if (!ispower2(ovsfactor))
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "coherencefft factor not power of 2";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}
	  if (Master.lines() != Mask.lines() || Master.pixels() != Mask.pixels())
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "Mask, Master not same size.";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}
	  if (!(ispower2(L) || ispower2(P)))
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "Mask, Master size not power of 2.";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}

	  // ______ Zero mean magnitude images ______
	  DEBUG.print("Using de-meaned magnitude patches for incoherent cross-correlation");
	  matrix<Float> magMaster = magnitude(Master);
	  matrix<Float> magMask = magnitude(Mask);
	  magMaster -= mean(magMaster);
	  magMask -= mean(magMask);

	  // ====== FFT's of master/mask ======
	  // ______ Pad with N zeros to prevent periodical convolution ______
	  matrix<complex<Float>> Master2 = new matrix(twoL,twoP); // initial 0
	  matrix<complex<Float>> Mask2 = new matrix(twoL,twoP); // initial 0
	  window windef = new window(0,0,0,0); // defaults to total matrix
	  window win1 = new window(0, L-1, 0, P-1);
	  window win2 = new window(halfL, halfL+L-1, halfP, halfP+P-1);
	  Master2.setdata(win1,mat2cr4(magMaster),windef); // zero-mean magnitude
	  Mask2.setdata(win2,mat2cr4(magMask),windef); // zero-mean magnitude

	  // ______ Crossproducts in spectral/space domain ______
	  // ______ Use Mask2 to store cross products temporarly ______
	  fft2d(Master2);
	  fft2d(Mask2);
	  Master2.conj();
	  Mask2 *= Master2; // corr = conj(M).*S
	  ifft2d(Mask2); // cross prod. in space

	  // ______ keep cross-products for shifts [-AccL,+AccL) ______
	  window wintmp = new window(halfL-AccL, halfL+AccL-1, halfP-AccP, halfP+AccP-1);
	  matrix<complex<Float>> TMP = new matrix(wintmp,Mask2);
	  matrix<Float> Covar = real(TMP); // imag==0

	  // ====== Compute norms, zero padded matrices ======
	  final complex<Float> ONE = new complex(1.0);
	  matrix<complex<Float>> blok = new matrix(L,P);
	  blok.setdata(ONE); // only real part

	  Master2.clean(); // reset to zeros
	  Mask2.clean(); // reset to zeros
	  Master2.setdata(win1,mat2cr4(sqr(magMaster)),windef); // use Master2 for intensity
	  Mask2.setdata(win2,blok,windef); // use Mask2 for padded Block

	  fft2d(Master2); // (intensity of master)
	  fft2d(Mask2); // (block)
	  Mask2.conj(); // conj(block)
	  Master2 *= Mask2;
	  Master2.conj(); // Master2 == conj(Master)*block
	  ifft2d(Master2);
	  // ______ Master2 now contains norms of master image in space domain ______
	  // ______ Resize to shifts [-AccL,+AccL) ______
	  TMP.setdata(Master2,wintmp); // fill TMP
	  matrix<Float> pmaster = real(TMP); // norms in pmaster

	  // ====== Now compute norms for slave image ======
	  Master2.clean(); // reset to zeros
	  window win5 = new window(L,twoL-1,P,twoP-1);
	  Master2.setdata(win5,mat2cr4(sqr(magMask)),windef);
	  fft2d(Master2); // (intensity of slave)
	  Master2 *= Mask2; // Master2 == conj(block)*Slave
	  ifft2d(Master2);
	  // ______ Master2 now contains norms of slave image in space domain ______
	  // ______ Resize to shifts [-AccL,+AccL) ______
	  TMP.setdata(Master2,wintmp); // fill TMP
	  final matrix<Float> pmask = real(TMP); // norms in pmask
	  pmaster *= pmask;
	  Covar /= Math.sqrt(pmaster);

	  // ====== Estimate shift by oversampling estimated correlation ======
	  int offL;
	  int offP;
	  final float maxcorr = (ovsfactor ==1) ? max(Covar,offL,offP) : max(oversample(Covar,ovsfactor,ovsfactor),offL,offP);
	  offsetL.argvalue = -(float)AccL + (float)offL/(float)ovsfactor; // update by reference
	  offsetP.argvalue = -(float)AccP + (float)offP/(float)ovsfactor; // update by reference
	  return maxcorr;
	  } // END coherencefft

//***************************************************************
// * crosscorrelate                                               *
// *                                                              *
// * cross correlation of zero-meaned magnitude of two patches    *
// *  uses ffts, some tricks for speed-up.                        *
// *  optionally improves peak position to sub-pixel.             *
// * This is an improvement upon coherencefft: faster and local peak *
// * Better to put this in matrixspecs                            *
// *                                                              *
// * input:                                                       *
// *  - Master                                                    *
// *  - Mask (same size as Master)                                *
// * output:                                                      *
// *  - peak correlation value [-1 1]                             *
// *  - updated offsetL, P                                        *
// *    positive offsetL: Mask is shifted up                      *
// *    positive offsetP: Mask is shifted left                    *
// *                                                              *
// * Bert Kampes, 12-Aug-2005                                     *
// ***************************************************************
//C++ TO JAVA CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in Java):
private matrix<complex<Float>> BLOCK;


	// ______ Correlation with FFT ______
	public static float crosscorrelate(matrix<complex<Float>> Master, matrix<complex<Float>> Mask, int ovsfactor, int AccL, int AccP, RefObject<Float> offsetL, RefObject<Float> offsetP)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"crosscorrelate (BK 12-Aug-2005)"<<ends;
		  TRACE.print();
	  }
	  // ______ Internal variables ______
	  final int L = Master.lines();
	  final int P = Master.pixels();
	  final int twoL = 2 *L;
	  final int twoP = 2 *P;
	  final int halfL = L/2;
	  final int halfP = P/2;

	  // ______ Check input ______
	  if (Master.lines() != Mask.lines() || Master.pixels() != Mask.pixels())
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "Mask, Master not same size.";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}
	  if (!(ispower2(L) || ispower2(P)))
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "Mask, Master size not power of 2.";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}
	  if (!ispower2(ovsfactor))
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "coherencefft factor not power of 2";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}

	  // ______ Zero mean magnitude images ______
	  DEBUG.print("Using de-meaned magnitude patches for incoherent cross-correlation");
	  matrix<Float> magMaster = magnitude(Master);
	  matrix<Float> magMask = magnitude(Mask);
	  magMaster -= mean(magMaster);
	  magMask -= mean(magMask);

	  // ====== (1) Compute cross-products of Master/Mask ======
	  // ______ Pad with N zeros to prevent periodical convolution ______
	  matrix<complex<Float>> Master2 = new matrix(twoL,twoP); // initial 0
	  matrix<complex<Float>> Mask2 = new matrix(twoL,twoP); // initial 0
	  window windef = new window(0,0,0,0); // defaults to total matrix
	  window win1 = new window(0, L-1, 0, P-1);
	  window win2 = new window(halfL, halfL+L-1, halfP, halfP+P-1);
	  Master2.setdata(win1,mat2cr4(magMaster),windef); // zero-mean magnitude
	  Mask2.setdata(win2,mat2cr4(magMask),windef); // zero-mean magnitude
	  // ______ Crossproducts in spectral/space domain ______
	  // ______ Use Mask2 to store cross products temporarly ______
	  fft2d(Master2);
	  fft2d(Mask2);
	  Master2.conj();
	  Mask2 *= Master2; // corr = conj(M).*S
	  ifft2d(Mask2); // real(Mask2): cross prod. in space

	  // ====== (2) compute norms for all shifts ======
	  // ______ use tricks to do this efficient ______
	  // ______ real(Mask2) contains cross-products ______
	  // ______ Mask2(0,0):Mask2(N,N) for shifts = -N/2:N/2 ______
	  // ______ rest of this matrix should not be used ______
	  // ______ Use Master2 to store intensity here in re,im ______
	  Master2.clean(); // reset to zeros
//C++ TO JAVA CONVERTER NOTE: 'register' variable declarations are not supported in Java:
//ORIGINAL LINE: register int l,p;
	  int l;
	  int p;
	  // --- flipud(fliplr(master^2) in real ---
	  // --- mask^2 in imag part; this saves a fft ---
	  // --- automatically the real/imag parts contain the norms ---
	  for (L =L; l<twoL; ++l)
		for (P =P; p<twoP; ++p)
		  Master2(l,p) = complex<Float>(sqr(magMaster(twoL-1-l,twoP-1-p)), sqr(magMask(l-L,p-P)));
	  // --- use a static block for fast computation ---
	//C++ TO JAVA CONVERTER NOTE: This static local variable declaration (not allowed in Java) has been moved just prior to the method:
	//  static matrix<complex<float>> BLOCK; // initial 0
	  if ((int)(BLOCK.lines())!=twoL || (int)(BLOCK.pixels())!=twoP)
		{
		DEBUG << "crosscorrelate:changing static block to size [" << twoL << ", " << twoP << "]";
		DEBUG.print();
		BLOCK.resize(twoL,twoP);
		for (L =halfL; l<halfL+L; ++l)
		  for (P =halfP; p<halfP+P; ++p)
			BLOCK(l,p) = complex<Float>(1.0);
		fft2d(BLOCK);
		BLOCK.conj(); // static variable: keep this for re-use
		}
	  // _____ Compute the cross-products, i.e., the norms for each shift ---
	  // ______ Master2(0,0):Master2(N,N) for shifts = -N/2:N/2 ______
	  fft2d(Master2);
	  Master2 *= BLOCK;
	  ifft2d(Master2); // real(Master2): powers of Master; imag(Master2): Mask


	  // ====== (3) find maximum correlation at pixel level ======
	  matrix<Float> Covar = new matrix(L+1,P+1); // correlation for each shift
	  float maxcorr = -999.0;
	  int maxcorrL = 0; // local index in Covar of maxcorr
	  int maxcorrP = 0; // local index in Covar of maxcorr
	  for (L =0; l<=L; ++l) // all shifts
		{
		for (P =0; p<=P; ++p) // all shifts
		  {
		  Covar(l,p) = real(Mask2(l,p)) / Math.sqrt(real(Master2(l,p))*imag(Master2(l,p)));
		  if (Covar(l,p) > maxcorr)
			{
			maxcorr = Covar(l,p);
			maxcorrL = l; // local index in Covar of maxcorr
			maxcorrP = p; // local index in Covar of maxcorr
			}
		  }
		}
	  offsetL.argvalue = -halfL + maxcorrL; // update by reference
	  offsetP.argvalue = -halfP + maxcorrP; // update by reference
	  DEBUG << "Pixel level offset:     " << offsetL.argvalue << ", " << offsetP.argvalue << " (corr=" << maxcorr << ")";
	  DEBUG.print();

	  // ====== (4) oversample to find peak sub-pixel ======
	  // ====== Estimate shift by oversampling estimated correlation ======
	  if (ovsfactor>1)
		{
		// --- (4a) get little chip around max. corr, if possible ---
		// --- make sure that we can copy the data ---
		if (maxcorrL<AccL)
		  {
		  DEBUG << "Careful, decrease AccL or increase winsizeL";
		  DEBUG.print();
		  maxcorrL = AccL;
		  }
		if (maxcorrP<AccP)
		  {
		  DEBUG << "Careful, decrease AccP or increase winsizeP";
		  DEBUG.print();
		  maxcorrP = AccP;
		  }
		if (maxcorrL>(L-AccL))
		  {
		  DEBUG << "Careful, decrease AccL or increase winsizeL";
		  DEBUG.print();
		  maxcorrL = L-AccL;
		  }
		if (maxcorrP>(P-AccP))
		  {
		  DEBUG << "Careful, decrease AccP or increase winsizeP";
		  DEBUG.print();
		  maxcorrP = P-AccP;
		  }
		// --- Now get the chip around max corr ---
		//matrix<real4> chip(2*AccL,2*AccP);// locally oversample corr
		//for (l=maxcorrL-AccL; l<maxcorrL+AccL; ++l)
		//  for (p=maxcorrP-AccP; p<maxcorrP+AccP; ++p)
		//    chip(l-(maxcorrL-AccL),p-(maxcorrP-AccP)) = Covar(l,p);
		window win3 = new window(maxcorrL-AccL,maxcorrL+AccL-1, maxcorrP-AccP,maxcorrP+AccP-1);
		final matrix<Float> chip = new matrix(win3,Covar); // construct as part
		// --- (4b) oversample chip to obtain sub-pixel max ---
		int offL;
		int offP;
		maxcorr = max(oversample(chip, ovsfactor, ovsfactor), offL,offP);
		offsetL.argvalue = -halfL + maxcorrL - AccL + (float)offL/(float)ovsfactor;
		offsetP.argvalue = -halfP + maxcorrP - AccP + (float)offP/(float)ovsfactor;
		DEBUG << "Sub-pixel level offset: " << offsetL.argvalue << ", " << offsetP.argvalue << " (corr=" << maxcorr << ")";
		DEBUG.print();
		}
	  return maxcorr;
	  } // END crosscorrelate

//***************************************************************
// * coherencespace                                               *
// *                                                              *
// * coherence in space domain based on magnitude                 *
// *  uses extension with zeros                                   *
// *                                                              *
// * input:                                                       *
// *  - Master                                                    *
// *  - Mask (size Master)                                        *
// * output:                                                      *
// *  - coherence value                                           *
// *  - updated offsetL, P                                        *
// *                                                              *
// *    Bert Kampes, 03-Feb-1999                                  *
// ***************************************************************


	// ______ Correlation in space domain ______
	public static float coherencespace(input_fine fineinput, matrix<complex<Float>> Master, matrix<complex<Float>> Mask, RefObject<Float> offsetL, RefObject<Float> offsetP)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"coherencespace (BK 03-Feb-1999)"<<ends;
		  TRACE.print();
	  }

	  // ______ Internal variables ______
	  final int L = Master.lines();
	  final int P = Master.pixels();
	  final int AccL = fineinput.AccL;
	  final int AccP = fineinput.AccP;
	  final int factor = fineinput.osfactor;

	  // ______Select parts of Master/slave______
	  final int MasksizeL = L - 2 *AccL;
	  final int MasksizeP = P - 2 *AccP;

	  // ______ Check input ______
	  if (!ispower2(AccL) || !ispower2(AccP))
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "AccL should be power of 2 for oversampling.";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}
	  if (MasksizeL < 4 || MasksizeP < 4)
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "Correlationwindow size too small (<4; size= FC_winsize-2*FC_Acc).";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}

	  // ______Shift center of Slave over Master______
	  window winmask = new window(AccL, AccL+MasksizeL-1, AccP, AccP+MasksizeP-1);

	  matrix<Float> coher = new matrix(2 *AccL,2 *AccP); // store result
													// 1st element: shift==AccL
	  window windef = new window(0, 0, 0, 0); // defaults to total

	  switch (fineinput.method)
		{
		case fc_cmplxspace:
		  {
		  {
			  ERROR.terminate();
			  String cp_s = new String(new char[256]);
			  cp_s = "not implemented in v1.0";
			  ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			  ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			  ERROR.print();
		  }
		  throw(unhandled_case_error);
		  break;
		  }

		case fc_magspace:
		  {
		  matrix<Float> magMask = magnitude(Mask); // magnitude
		  magMask -= mean(magMask); // subtract mean
		  matrix<Float> Mask2 = new matrix(winmask,magMask); // construct as part
		  float normmask = norm2(Mask2);
		  matrix<Float> Master2 = new matrix(MasksizeL, MasksizeP);
		  matrix<Float> magMaster = magnitude(Master);
		  magMaster -= mean(magMaster);
		  window winmaster;
		  for (register int i =0;i<2 *AccL;i++)
			{
			winmaster.linelo = i;
			winmaster.linehi = i+MasksizeL-1;
			for (register int j =0;j<2 *AccP;j++)
			  {
			  winmaster.pixlo = j;
			  winmaster.pixhi = j+MasksizeP-1;
			  Master2.setdata(windef,magMaster,winmaster);
			  // ______Coherence for this position______
			  float cohs1s2 = 0.;
			  float cohs1s1 = 0.;
			  for (register int k =0;k<MasksizeL;k++)
				{
				for (register int l =0;l<MasksizeP;l++)
				  {
				  cohs1s2 += (Master2(k,l) * Mask2(k,l));
				  cohs1s1 += sqr(Master2(k,l));
				  }
				}
			  coher(i,j) = cohs1s2 / Math.sqrt(cohs1s1 * normmask); // [-1 1]
			  }
			}
		  break;
		  }

		default:
		  {
			  ERROR.terminate();
			  String cp_s = new String(new char[256]);
			  cp_s = "unknown method";
			  ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			  ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			  ERROR.print();
		  }
		  throw(unhandled_case_error);
		} // switch method


	  // ______ Correlation in space domain ______
	  int offL;
	  int offP;
	  final matrix<Float> coher8 = oversample(coher,factor,factor);
	  final float maxcor = max(coher8,offL,offP);
	  offsetL.argvalue = AccL - offL/(float)factor; // update by reference
	  offsetP.argvalue = AccP - offP/(float)factor; // update by reference
	  return maxcor;
	  } // END coherencespace

//***************************************************************
// * coregpm                                                      *
// *                                                              *
// * Compute coregistration parameters (least squares)            *
// *                                                              *
// * input:                                                       *
// *  - Position of windows                                       *
// *  - Computed offsets                                          *
// *                                                              *
// * output:                                                      *
// *  - coregistration parameters to file                         *
// *    (wrt. normalized master grid)                             *
// *                                                              *
// *    Bert Kampes, 22-Feb-1999                                  *
// *    Bert Kampes, 26-Oct-1999 normalized coordinates           *
// * changed matrices from real4 to real8,                        *
// #%// BK 22-Mar-2001                                            *
// ***************************************************************


	// ______ Compute coregistration parameters ______
	//      const window            &originalmaster,
	public static void coregpm(slcimage master, slcimage slave, String i_resfile, input_coregpm coregpminput, short demassist)
			//const uint          oversamplingsfactorfine)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"coregpm (BK 26-Oct-1999)"<<ends;
		  TRACE.print();
	  }
	  // ______Names of variables in this routine______
	  // unknowns:          x
	  // solution unknowns: placed in rhsL/P
	  // observations:      y
	  // covariance obs.:   Qy_1
	  // designmatrix:      A
	  // normalmatrix:      N = At. Qy-1. A
	  // covariance unkn.:  N_1 (Qx_hat, inverse normalmatrix)
	  // estimates:         *_hat

	  final float THRESHOLD = coregpminput.threshold; // threshold ...
	  final int DEGREE = coregpminput.degree; // degree of polynomial
	  final int MAX_ITERATIONS = coregpminput.maxiter; // max. of pnts to remove
	  final float CRIT_VALUE = coregpminput.k_alpha; // crit. value outlier removal
	  final int Nunk = Ncoeffs(DEGREE); // Number of unknowns/direction

	  // ______ Normalize data for polynomial ______
	  final double minL = master.originalwindow.linelo;
	  final double maxL = master.originalwindow.linehi;
	  final double minP = master.originalwindow.pixlo;
	  final double maxP = master.originalwindow.pixhi;

	  // ______ A priori sigma of  offset ______
	  // ______ Read this factor from the result file 
	  // ______ "Oversampling factor: 32"
	  // ______ "Window_size_L_for_correlation: 4"
	  // ______ "Window_size_P_for_correlation: 121"
	  DEBUG.print("Reading oversampling factor from result file");
	  int osfactor = 32; // oversamplingsfactor
	  int corrwinL = 64; // window size to compute FINE correlation
	  int corrwinP = 64; // window size to compute FINE correlation
	  String c4osfactor = new String(new char[4]);
	  String c10corrwinL = new String(new char[10]);
	  String c10corrwinP = new String(new char[10]);
	  boolean found = readres(c4osfactor,sizeof(c4osfactor),i_resfile, "Oversampling", 1);
	  if (found)
		  osfactor = (int)(Integer.parseInt(c4osfactor));
	  found = readres(c10corrwinL,sizeof(c10corrwinL),i_resfile, "Window_size_L_for_correlation:", 0);
	  if (found)
		  corrwinL = (int)(Integer.parseInt(c10corrwinL));
	  found = readres(c10corrwinP,sizeof(c10corrwinP),i_resfile, "Window_size_P_for_correlation:", 0);
	  if (found)
		  corrwinP = (int)(Integer.parseInt(c10corrwinP));
	  corrwinL = max(10,corrwinL-8); // if fft method peak is not at center
	  corrwinP = max(10,corrwinP-8); // +then effective number of samples is smaller
	  // _____ oversampling factor is bin in which maximum can be found _____
	  // _____ ovsf=16-->apriorisigma=0.03
	  final float ACCURACY = 0.5 * (1.0/((float)osfactor));

	  // but we need coreg accuracy of 0.1 pixel about.  therefore use a priori
	  // based on experience here, and different for azimuth and range
	  // this also helps our automated outlier detection and testing hopefully.
	  // BK 15-Apr-2003
	  // if the image is oversampled, then still use orig spacing
	  float SIGMAL =-999.9; // sigma in orig pixels
	  float SIGMAP =-999.9; // seems range direction is better???
	  if (coregpminput.weightflag!=3)
		{
		SIGMAL = 0.15/master.ovs_az; // sigma in orig pixels
		SIGMAP = 0.10/master.ovs_rg; // seems range direction is better???
		DEBUG.print("Using a smaller sigma in range, because it seems that can be estimated better");
		INFO << "a priori std.dev offset vectors line direction [samples]:  " << SIGMAL;
		INFO.print();
		INFO << "a priori std.dev offset vectors pixel direction [samples]: " << SIGMAP;
		INFO.print();
		}

	  // ______ Find #points > threshold ______
	  matrix<Float> Data = getofffile(i_resfile, THRESHOLD);
	  // ______ Data contains the following: ______
	  // Data(i,0) = winnumber; Data(i,1) = posL; Data(i,2) = posP; 
	  // Data(i,3) = offL;      Data(i,4) = offP; Data(i,5) = corr;

	  // ______ start added by FvL ______
	  ifstream DeltaLfile;
	  ifstream DeltaPfile;
	  streampos pos;

	  if (demassist != 0)
		{
		  openfstream(DeltaLfile,"delta_line.raw");
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  bk_assert(DeltaLfile,"delta_line.raw",__FILE__,__LINE__);
		  openfstream(DeltaPfile,"delta_pixel.raw");
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  bk_assert(DeltaPfile,"delta_pixel.raw",__FILE__,__LINE__);

		  int posL;
		  int posP;
		  float offL;
		  float offP;
		  double deltaL;
		  double deltaP;
		  final int sizer8 = sizeof(double);
		  float ms_az_timing_error_L = slave.az_timing_error; // ms = masterslave: relative timing error
		  float ms_r_timing_error_P = slave.r_timing_error;

		  for (register int ii =0; ii<Data.lines(); ii++)
			{
			  posL = (int)(Data(ii,1));
			  posP = (int)(Data(ii,2));
			  offL = Data(ii,3);
			  offP = Data(ii,4);
			  pos = (streampos)((posL-master.currentwindow.linelo)* master.currentwindow.pixels() + posP - master.currentwindow.pixlo);
			  pos = (streampos)(pos * sizer8);

			  DeltaLfile.seekg(pos,ios.beg);
			  DeltaPfile.seekg(pos,ios.beg);

			  DeltaLfile.read((char)&deltaL, sizer8);
			  DeltaPfile.read((char)&deltaP, sizer8);

			  Data(ii,3) = offL-(float)deltaL-ms_az_timing_error_L;
			  Data(ii,4) = offP-(float)deltaP-ms_r_timing_error_P;
			}
		}

	  // ______ end added by FvL ______

	  int ITERATION = 0;
	  int DONE = 0;
	  // sqr: level significance: alpha=0.001; power of test: gamma=0.80
	  //real4 CRIT_VALUE = sqrt(3.29);
	  INFO << "Critical value for outlier test: " << CRIT_VALUE;
	  INFO.print();
	  int winL = 0; // window number to be removed
	  int winP = 0; // window number of largest w -test in range
	  matrix<Double> eL_hat;
	  matrix<Double> eP_hat;
	  matrix<Double> wtestL;
	  matrix<Double> wtestP;
	  matrix<Double> rhsL;
	  matrix<Double> rhsP;
	  matrix<Double> Qx_hat;
	  double maxdev = 0.0;
	  double overallmodeltestL = 0.0;
	  double overallmodeltestP = 0.0;
	  double maxwL;
	  double maxwP;
//C++ TO JAVA CONVERTER NOTE: 'register' variable declarations are not supported in Java:
//ORIGINAL LINE: register int i,j,k,index;
	  int i;
	  int j;
	  int k;
	  int index;
	  while (DONE != 1)
		{
		PROGRESS << "Start iteration " << ITERATION;
		PROGRESS.print();
		// ______ Remove identified outlier from previous estimation ______
		if (ITERATION != 0)
		  {
		  matrix<Float> tmp_DATA = Data; //(remove_observation_i,*);
		  Data.resize(Data.lines()-1, Data.pixels());
		  j = 0; // counter over reduced obs.vector
		  for (i =0; i<tmp_DATA.lines(); i++) // counter over original window numbers
			{
			if (i != winL) // do not copy the one to be removed.
			  {
			  Data.setrow(j,tmp_DATA.getrow(i)); // copy back without removed obs.
			  j++; // fill next row of Data
			  }
			else
			  {
			  INFO << "Removing observation " << i << " from observation vector.";
			  INFO.print();
			  }
			}
		  }

		// ______Check redundancy______
		int Nobs = Data.lines(); // Number of points > threshold
		if (Nobs < Nunk)
		  {
		  {
			  ERROR.terminate();
			  String cp_s = new String(new char[256]);
			  cp_s = "coregpm: Number of windows > threshold is smaller than parameters solved for.";
			  ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			  ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			  ERROR.print();
		  }
		  throw(input_error);
		  }

		// ______Set up system of equations______
		// ______Order unknowns: A00 A10 A01 A20 A11 A02 A30 A21 A12 A03 for degree=3______
		matrix<Double> yL = new matrix(Nobs,1); // observation
		matrix<Double> yP = new matrix(Nobs,1); // observation
		matrix<Double> A = new matrix(Nobs,Nunk); // designmatrix
		matrix<Double> Qy_1 = new matrix(Nobs,1); // a priori covariance matrix (diag)

		// ______ Normalize data for polynomial ______
		INFO << "coregpm: polynomial normalized by factors: " << minL << " " << maxL << " " << minP << " " << maxP << " to [-2,2]";
		INFO.print();

		// ______Fill matrices______
		DEBUG.print("Setting up design matrix for LS adjustment");
		for (i =0; i<Nobs; i++)
		  {
		  double posL = normalize((double)(Data(i,1)),minL,maxL);
		  double posP = normalize((double)(Data(i,2)),minP,maxP);
		  yL(i,0) = (double)(Data(i,3));
		  yP(i,0) = (double)(Data(i,4));
		  DEBUG << "coregpm: (" << posL << ", "<< posP << "): yL=" << yL(i,0) << " yP=" << yP(i,0);
		  DEBUG.print();
		  // ______Set up designmatrix______
		  index = 0;
		  for (j =0; j<=DEGREE; j++)
			{
			for (k =0; k<=j; k++)
			  {
			  A(i,index) = Math.pow(posL,(double)(j-k)) * Math.pow(posP,(double)k);
			  index++;
			  }
			}
		  }


		// ______Weight matrix data______
		DEBUG.print("Setting up (inverse of) covariance matrix for LS adjustment");
		switch(coregpminput.weightflag)
		  {
		  case 0:
			for (i =0; i<Nobs; i++)
			  Qy_1(i,0) = (double)1.0;
			break;
		  case 1:
			DEBUG.print("Using sqrt(coherence) as weights.");
			for (i =0; i<Nobs; i++)
			  Qy_1(i,0) = (double)(Data(i,5)); // more weight to higher correlation
			// ______ Normalize weights to avoid influence on estimated var.factor ______
			INFO.print("Normalizing covariance matrix for LS estimation.");
			Qy_1 = Qy_1 / mean(Qy_1); // normalize weights (for tests!)
			break;
		  case 2:
			DEBUG.print("Using coherence as weights.");
			for (i =0; i<Nobs; i++)
			  Qy_1(i,0) = (double)(Data(i,5))*(double)(Data(i,5)); // more weight to higher correlation
			// ______ Normalize weights to avoid influence on estimated var.factor ______
			INFO.print("Normalizing covariance matrix for LS estimation.");
			Qy_1 = Qy_1 / mean(Qy_1); // normalize weights (for tests!)
			break;
		  // --- Bamler paper igarss 2000 and 2004; Bert Kampes, 16-Aug-2005 ---
		  case 3:
			// for coherent cross-correlation the precision of the shift is
			// sigma_cc = sqrt(3/(2N))*sqrt(1-coh^2)/(pi*coh) in units of pixels
			// for incoherent cross-correlation as we do, sigma seems approx. [BK]
			// sigma_ic = sqrt(2/coh)*sigma_cc
			// actually with osf^1.5 (but we will ignore that here)
			// it seems for large N this is to optimistic, maybe because of a bias
			// in the coherence estimator, or some other reason;  in any case,
			// the result is a large number of warnings.
			DEBUG.print("Using expression Bamler04 as weights.");
			for (i =0; i<Nobs; i++)
			  {
			  // N_corr: number of samples for cross-corr; approx. FC_WINSIZE
			  // number of effictive samples depends on data ovs factor 
			  // Bamler 2000: also on oversampling ratio of data, but ignored here.
			  final float N_corr = corrwinL *corrwinP;
			  final float coh = Data(i,5); // estimated correlation; assume unbiased?
			  final float sigma_cc = Math.sqrt(3.0/(2.0 *N_corr))*Math.sqrt(1.0-sqr(coh))/(PI *coh);
			  final float sigma_ic = Math.sqrt(2.0/coh)*sigma_cc;
			  DEBUG << "Window " << i << ": estimated coherence   = " << coh;
			  DEBUG.print();
			  DEBUG << "Window " << i << ": sigma(estimated shift) for coherent cross-correlation = " << sigma_cc << " [pixel]";
			  DEBUG.print();
			  DEBUG << "Window " << i << ": sigma(estimated shift) = " << sigma_ic << " [pixel]";
			  DEBUG.print();
			  Qy_1(i,0) = 1.0/sqr(sigma_ic); // Qy_1=diag(inverse(Qy));
			  SIGMAL = 1.0; // remove this factor effectively
			  SIGMAP = 1.0; // remove this factor effectively
			  }
			break;
		  default:
			{
				ERROR.terminate();
				String cp_s = new String(new char[256]);
				cp_s = "Panic, not possible with checked input.";
				ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
				ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
				ERROR.print();
			}
			throw(unhandled_case_error);
		  }


		// ______Compute Normalmatrix, rghthandside______
		matrix<Double> N = matTxmat(A,diagxmat(Qy_1,A));
		//matrix<real8> rhsL = matTxmat(A,diagxmat(Qy_1,yL));
		//matrix<real8> rhsP = matTxmat(A,diagxmat(Qy_1,yP));
		//matrix<real8> Qx_hat = N;
		rhsL = matTxmat(A,diagxmat(Qy_1,yL));
		rhsP = matTxmat(A,diagxmat(Qy_1,yP));
		Qx_hat = N;
		// ______Compute solution______
		choles(Qx_hat); // Cholesky factorisation normalmatrix
		solvechol(Qx_hat,rhsL); // Solution unknowns in rhs
		solvechol(Qx_hat,rhsP); // Solution unknowns in rhs
		invertchol(Qx_hat); // Covariance matrix of unknowns
		// ______Test inverse______
		for (i =0; i<Qx_hat.lines(); i++)
		  for (j =0; j<i; j++)
			Qx_hat(j,i) = Qx_hat(i,j); // repair Qx
		maxdev = max(Math.abs(N *Qx_hat-eye((double)(Qx_hat.lines()))));
		INFO << "coregpm: max(abs(N*inv(N)-I)) = " << maxdev;
		INFO.print();
		// ___ use trace buffer to store string, remember to rewind it ___
		if (maxdev > .01)
		  {
		  ERROR << "coregpm: maximum deviation N*inv(N) from unity = " << maxdev << ". This is larger than 0.01";
		  ERROR.print(ERROR.get_str());
		  throw(some_error);
		  }
		else if (maxdev > .001)
		  {
		  WARNING << "coregpm: maximum deviation N*inv(N) from unity = " << maxdev << ". This is between 0.01 and 0.001";
		  WARNING.print();
		  }


		// ______Some other stuff, scale is ok______
		matrix<Double> Qy_hat = A * (matxmatT(Qx_hat,A));
		matrix<Double> yL_hat = A * rhsL;
		matrix<Double> yP_hat = A * rhsP;
		//matrix<real8> eL_hat      = yL - yL_hat;
		//matrix<real8> eP_hat      = yP - yP_hat;
		eL_hat = yL - yL_hat;
		eP_hat = yP - yP_hat;
		//  matrix<real4> Qe_hat    = Qy - Qy_hat;
		matrix<Double> Qe_hat = -Qy_hat;
		for (i =0; i<Nobs; i++)
		  Qe_hat(i,i) += (1. / Qy_1(i,0));

		// ______Overall model test (variance factor)______
		overallmodeltestL = 0.;
		overallmodeltestP = 0.;
		for (i =0; i<Nobs; i++)
		  {
		  overallmodeltestL += sqr(eL_hat(i,0))*Qy_1(i,0);
		  overallmodeltestP += sqr(eP_hat(i,0))*Qy_1(i,0);
		  }
		overallmodeltestL = (overallmodeltestL/sqr(SIGMAL)) /(Nobs-Nunk); // this is sigma hat!
		overallmodeltestP = (overallmodeltestP/sqr(SIGMAP)) /(Nobs-Nunk); // not OMT!
		INFO << "coregpm: overallmodeltest Lines = " << overallmodeltestL;
		INFO.print();
		INFO << "coregpm: overallmodeltest Pixels = " << overallmodeltestP;
		INFO.print();

		// ______Datasnooping, assume Qy diag______
		wtestL.resize(Nobs,1);
		wtestP.resize(Nobs,1);
		for (i =0; i<Nobs; i++)
		  {
		  wtestL(i,0) = eL_hat(i,0) / (Math.sqrt(Qe_hat(i,i))*SIGMAL); // computed excl.var.factor
		  wtestP(i,0) = eP_hat(i,0) / (Math.sqrt(Qe_hat(i,i))*SIGMAP);
		  }

		int dumm = 0;
		maxwL = max(Math.abs(wtestL),winL,dumm); // returns winL
		maxwP = max(Math.abs(wtestP),winP,dumm); // returns winP
		INFO << "maximum wtest statistic azimuth = " << maxwL << " for window number: " << Data(winL,0);
		INFO.print();
		INFO << "maximum wtest statistic range   = " << maxwP << " for window number: " << Data(winP,0);
		INFO.print();
		// --- use summed wtest for outlier detection ---
		// #%// BK 21-Oct-2003
		matrix<Double> wtestsum = sqr(wtestL)+sqr(wtestP); // (Nobs,1)
		double maxwsum = max(wtestsum,winL,dumm); // idx to remove
		INFO << "Detected outlier:  summed sqr.wtest = " << maxwsum << "; observation: " << winL << "; window number: " << Data(winL,0);
		INFO.print();


		// ______ Test if we are done yet ______
		if (Nobs <= Nunk)
		  {
		  WARNING.print("NO redundancy!  Exiting iterations.");
		  DONE = 1; // cannot remove more than this
		  }
		// seems something fishy here..., b-method of testing delft
		//    if (max(overallmodeltestL,overallmodeltestP) < 1.0)
		//      {
		//      INFO.print("OMTs accepted, not iterating anymore (final solution reached).");
		//      DONE = 1;// ok (?).
		//      }
		if (max(maxwL,maxwP) <= CRIT_VALUE) // all tests accepted?
		  {
		  INFO.print("All outlier tests accepted! (final solution computed)");
		  DONE = 1; // yeah!
		  }
		if (ITERATION >= MAX_ITERATIONS)
		  {
		  INFO.print("max. number of iterations reached (exiting loop).");
		  DONE = 1; // we reached max. (or no max_iter specified)
		  }

		// ______ Only warn if last iteration has been done ______
		if (DONE == 1)
		  {
		  // ___ use trace buffer to store string, remember to rewind it ___
		  if (overallmodeltestL > 10)
			{
			WARNING << "coregpm: overallmodeltest Lines = " << overallmodeltestL << ends;
			WARNING.print();
			WARNING << " is larger than 10. (Suggest model or a priori sigma not correct.)";
			WARNING.print();
			}
		  // ___ use trace buffer to store string, remember to rewind it ___
		  if (overallmodeltestP > 10)
			{
			WARNING << "coregpm: overallmodeltest Pixels = " << overallmodeltestP;
			WARNING.print();
			WARNING << " is larger than 10.\n(suggests a priori sigma not correct.)";
			WARNING.print();
			}
		  // if a priori sigma is correct, max wtest should be something like 1.96
		  if (max(maxwL,maxwP)>200.0)
			{
			WARNING << "Recommendation: remove window number: " << Data(winL,0) << " and re-run step COREGPM.  max. wtest is: " << max(maxwL,maxwP) << ".";
			WARNING.print();
			}

		  // this test seems to generate too many warnings...
		  // ______Test of Jaron Samson's thesis: not ok/ depends on SIGMA ...______
		  // //real4 expected = 1. / 400.;   // WRONG, 1./sqr(28) ???
		  // real8 expected = 1.0 / sqr(28); //784.;
		  // expected /= sqr(SIGMAL);        // correct for variance factor
		  // for (i=0; i<Nobs; i++)
		  //   if (Qy_hat(i,i) > expected)
		  //     {
		  //     WARNING << "coregpm: Qy_hat too large for window: "
		  //          << Data(i,0);
		  //     WARNING.print();
		  //     }

		  } // Only warn when done iterating.
		ITERATION++; // update counter here!
		} // iterations remove outliers


	  // ____ start added by FvL _________
	  // Determine inverse transformation
	  // (slave corners only, needed for overlap)

	  // ______ Normalize data for polynomial ______
	  final double sminL = slave.originalwindow.linelo;
	  final double smaxL = slave.originalwindow.linehi;
	  final double sminP = slave.originalwindow.pixlo;
	  final double smaxP = slave.originalwindow.pixhi;

	  // ______Check redundancy______
	  int Nobs = Data.lines(); // Number of points > threshold

	  // ______Set up system of equations for slave______
	  // ______Order unknowns: A00 A10 A01 A20 A11 A02 A30 A21 A12 A03 for degree=3______
	  matrix<Double> srhsL;
	  matrix<Double> srhsP;
	  matrix<Double> yL = new matrix(Nobs,1); // observation
	  matrix<Double> yP = new matrix(Nobs,1); // observation
	  matrix<Double> A = new matrix(Nobs,Nunk); // designmatrix
	  matrix<Double> Qy_1 = new matrix(Nobs,1); // a priori covariance matrix (diag)

	  // ______ Normalize data for polynomial ______
	  INFO << "coregpm: slave polynomial normalized by factors: " << sminL << " " << smaxL << " " << sminP << " " << smaxP << " to [-2,2]";
	  INFO.print();

	  // ______Fill matrices______
	  DEBUG.print("Setting up design matrix for LS adjustment");
	  for (i =0; i<Nobs; i++)
		{
		  double posL = normalize((double)(Data(i,1)+Data(i,3)),sminL,smaxL);
		  double posP = normalize((double)(Data(i,2)+Data(i,4)),sminP,smaxP);
		  yL(i,0) = (double)(-Data(i,3));
		  yP(i,0) = (double)(-Data(i,4));
		  DEBUG << "coregpm: (" << posL << ", "<< posP << "): yL=" << yL(i,0) << " yP=" << yP(i,0);
		  DEBUG.print();
		  // ______Set up designmatrix______
		  index = 0;
		  for (j =0; j<=DEGREE; j++)
			{
			  for (k =0; k<=j; k++)
				{
				  A(i,index) = Math.pow(posL,(double)(j-k)) * Math.pow(posP,(double)k);
				  index++;
				}
			}
		}

		// ______Weight matrix data______
		DEBUG.print("Setting up (inverse of) covariance matrix for LS adjustment");
		switch(coregpminput.weightflag)
		  {
		  case 0:
			for (i =0; i<Nobs; i++)
			  Qy_1(i,0) = (double)1.0;
			break;
		  case 1:
			DEBUG.print("Using sqrt(coherence) as weights.");
			for (i =0; i<Nobs; i++)
			  Qy_1(i,0) = (double)(Data(i,5)); // more weight to higher correlation
			// ______ Normalize weights to avoid influence on estimated var.factor ______
			INFO.print("Normalizing covariance matrix for LS estimation.");
			Qy_1 = Qy_1 / mean(Qy_1); // normalize weights (for tests!)
			break;
		  case 2:
			DEBUG.print("Using coherence as weights.");
			for (i =0; i<Nobs; i++)
			  Qy_1(i,0) = (double)(Data(i,5))*(double)(Data(i,5)); // more weight to higher correlation
			// ______ Normalize weights to avoid influence on estimated var.factor ______
			INFO.print("Normalizing covariance matrix for LS estimation.");
			Qy_1 = Qy_1 / mean(Qy_1); // normalize weights (for tests!)
			break;
		  // --- Bamler paper igarss 2000 and 2004; Bert Kampes, 16-Aug-2005 ---
		  case 3:
			// for coherent cross-correlation the precision of the shift is
			// sigma_cc = sqrt(3/(2N))*sqrt(1-coh^2)/(pi*coh) in units of pixels
			// for incoherent cross-correlation as we do, sigma seems approx. [BK]
			// sigma_ic = sqrt(2/coh)*sigma_cc
			// actually with osf^1.5 (but we will ignore that here)
			// it seems for large N this is to optimistic, maybe because of a bias
			// in the coherence estimator, or some other reason;  in any case,
			// the result is a large number of warnings.
			DEBUG.print("Using expression Bamler04 as weights.");
			for (i =0; i<Nobs; i++)
			  {
			  // N_corr: number of samples for cross-corr; approx. FC_WINSIZE
			  // number of effictive samples depends on data ovs factor 
			  // Bamler 2000: also on oversampling ratio of data, but ignored here.
			  final float N_corr = corrwinL *corrwinP;
			  final float coh = Data(i,5); // estimated correlation; assume unbiased?
			  final float sigma_cc = Math.sqrt(3.0/(2.0 *N_corr))*Math.sqrt(1.0-sqr(coh))/(PI *coh);
			  final float sigma_ic = Math.sqrt(2.0/coh)*sigma_cc;
			  DEBUG << "Window " << i << ": estimated coherence   = " << coh;
			  DEBUG.print();
			  DEBUG << "Window " << i << ": sigma(estimated shift) for coherent cross-correlation = " << sigma_cc << " [pixel]";
			  DEBUG.print();
			  DEBUG << "Window " << i << ": sigma(estimated shift) = " << sigma_ic << " [pixel]";
			  DEBUG.print();
			  Qy_1(i,0) = 1.0/sqr(sigma_ic); // Qy_1=diag(inverse(Qy));
			  SIGMAL = 1.0; // remove this factor effectively
			  SIGMAP = 1.0; // remove this factor effectively
			  }
			break;
		  default:
			{
				ERROR.terminate();
				String cp_s = new String(new char[256]);
				cp_s = "Panic, not possible with checked input.";
				ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
				ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
				ERROR.print();
			}
			throw(unhandled_case_error);
		  }

	  // ______Compute Normalmatrix, rghthandside______
	  matrix<Double> N = matTxmat(A,diagxmat(Qy_1,A)); //use same Qy_1
	  srhsL = matTxmat(A,diagxmat(Qy_1,yL));
	  srhsP = matTxmat(A,diagxmat(Qy_1,yP));
	  Qx_hat = N;

	  // ______Compute solution______
	  choles(Qx_hat); // Cholesky factorisation normalmatrix
	  solvechol(Qx_hat,srhsL); // Solution unknowns in rhs
	  solvechol(Qx_hat,srhsP); // Solution unknowns in rhs
	  invertchol(Qx_hat); // Covariance matrix of unknowns

	  double slave_l0 = slave.currentwindow.linelo;
	  double slave_lN = slave.currentwindow.linehi;
	  double slave_p0 = slave.currentwindow.pixlo;
	  double slave_pN = slave.currentwindow.pixhi;

	  double deltaline_slave00;
	  double deltapixel_slave00;
	  double deltaline_slave0N;
	  double deltapixel_slave0N;
	  double deltaline_slaveN0;
	  double deltapixel_slaveN0;
	  double deltaline_slaveNN;
	  double deltapixel_slaveNN;

	  deltaline_slave00 = polyval(normalize(slave_l0,sminL,smaxL), normalize(slave_p0,sminP,smaxP), srhsL,DEGREE);
	  deltapixel_slave00 = polyval(normalize(slave_l0,sminL,smaxL), normalize(slave_p0,sminP,smaxP), srhsP,DEGREE);
	  deltaline_slave0N = polyval(normalize(slave_l0,sminL,smaxL), normalize(slave_pN,sminP,smaxP), srhsL,DEGREE);
	  deltapixel_slave0N = polyval(normalize(slave_l0,sminL,smaxL), normalize(slave_pN,sminP,smaxP), srhsP,DEGREE);
	  deltaline_slaveN0 = polyval(normalize(slave_lN,sminL,smaxL), normalize(slave_p0,sminP,smaxP), srhsL,DEGREE);
	  deltapixel_slaveN0 = polyval(normalize(slave_lN,sminL,smaxL), normalize(slave_p0,sminP,smaxP), srhsP,DEGREE);
	  deltaline_slaveNN = polyval(normalize(slave_lN,sminL,smaxL), normalize(slave_pN,sminP,smaxP), srhsL,DEGREE);
	  deltapixel_slaveNN = polyval(normalize(slave_lN,sminL,smaxL), normalize(slave_pN,sminP,smaxP), srhsP,DEGREE);

	  // ____ end added by FvL _________

	  // ______ Create dump file for making plots ______
	  ofstream cpmdata = new ofstream("CPM_Data", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(cpmdata,"coregpm: CPM_DATA",__FILE__,__LINE__);
	  cpmdata << "File: CPM_Data" << "\nThis file contains information on the least squares" << "\n estimation of the coregistration parameters." << "\nThis info is used in the plotcmp script." << "\nThere are 10 columns with:" << "\nWindow number, position L, position P, " << "\n offsetL (observation), offsetP (observation), correlation," << "\n estimated errorL, errorP, w-test statistics for L, P." << "\nwin   posL  posP      offL      offP  corr      eL     eP  wtstL  wtstP" << "\n------------------------------------------------------------\n";
	  cpmdata.close();

	  // ______ Only way to format in c++ since stupid iomanip dont work? ______
	  FILE cpm;
	  cpm =fopen("CPM_Data","a");
	  //for (i=0; i<Nobs; i++)
	  for (i =0; i<Data.lines(); i++)
		fprintf(cpm, "%4.0f %5.0f %5.0f %# 9.2f %# 9.2f %# 6.2f %6.2f %6.2f %6.2f %6.2f\n", Data(i,0), Data(i,1), Data(i,2), Data(i,3), Data(i,4), Data(i,5), eL_hat(i,0), eP_hat(i,0), Math.abs(wtestL(i,0)), Math.abs(wtestP(i,0)));
	  fclose(cpm);


	  // ====== Write results to scratch files ======
	  ofstream scratchlogfile = new ofstream("scratchlogcpm", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchlogfile,"coregpm: scratchlogcpm",__FILE__,__LINE__);
	  scratchlogfile << "\n\n*******************************************************************" << "\n* COMP_COREGPM:" << "\n*******************************************************************" << "\nA polynomial model is weighted least squares estimated" << "\nfor azimuth and range through the FINE offset vectors." << "\nThe number of coefficient are the unknowns, the number of" << "\nobservations are the offset vectors above the THRESHOLD" << "\nspecified in the input file.  To estimate the unknowns, at" << "\nleast the number of observations must equal the number of unknowns." << "\nIf there are more observations, we can statistically test" << "\nwhether the observations fit the model, and whether there are" << "\noutliers in the observations, which we like to remove." << "\nThe overall model test does the first.  Wtest the second." << "\nWe advice to remove some bad estimated offsets by hand based" << "\nthe largest w-test, and to iterate running this step until" << "\nno outlier is identified anymore.  A great tool is plotting" << "\nthe observations and errors, which can be done with the utility" << "\nscripts provided by Doris (calls to GMT)." << "\nAlso see any book on LS methods." << "\n\nDegree of model:\t\t\t\t" << DEGREE << "\nThreshold on data (correlation):\t\t\t" << THRESHOLD << "\nOversmaplings factor used in fine:           \t" << osfactor << "\nThis means maximum can be found within [samples]: \t" << ACCURACY << "\nA priori sigma azimuth (based on experience): \t" << SIGMAL << "\nA priori sigma range (based on experience): \t" << SIGMAP << "\nNumber of observations: \t\t\t" << Data.lines() << "\nNumber of rejected observations: \t\t\t" << ITERATION << "\nNumber of unknowns: \t\t\t\t" << Nunk << "\nOverall model test in Azimuth direction: \t" << overallmodeltestL << "\nOverall model test in Range direction: \t\t" << overallmodeltestP << "\nLargest w test statistic in Azimuth direction: \t" << maxwL << "\n  for window number: \t\t\t\t" << Data(winL,0) << "\nLargest w test statistic in Range direction: \t" << maxwP << "\n  for window number: \t\t\t\t" << Data(winP,0) << "\nMaximum deviation from unity Normalmatrix*Covar(unknowns): \t" << maxdev << "\nEstimated parameters in Azimuth direction" << "\nx_hat \tstd" << "\n(a00 | a10 a01 | a20 a11 a02 | a30 a21 a12 a03 | ...)\n";
	  for (i =0; i<Nunk; i++)
		scratchlogfile << setiosflags(ios.fixed) << setiosflags(ios.showpoint) << setiosflags(ios.right) << setw(8) << setprecision(4) << rhsL(i,0) << " \t" << Math.sqrt(Qx_hat(i,i)) << "\n";
	  scratchlogfile << "\nEstimated parameters in Range direction" << "\n(b00 | b10 b01 | b20 b11 b02 | b30 b21 b12 b03 | ...)\n";
	  for (i =0; i<Nunk; i++)
		scratchlogfile << setiosflags(ios.fixed) << setiosflags(ios.showpoint) << setiosflags(ios.right) << setw(8) << setprecision(4) << rhsP(i,0) << " \t" << Qx_hat(i,i) << "\n";

	  scratchlogfile << "\nCovariance matrix estimated parameters:" << "\n---------------------------------------\n";
	  for (i =0; i<Nunk; i++)
		{
		for (j =0; j<Nunk; j++)
		  {
		  scratchlogfile << setiosflags(ios.fixed) << setiosflags(ios.showpoint) << setiosflags(ios.right) << setw(8) << setprecision(4) << Qx_hat(i,j) << " ";
		  }
		scratchlogfile << "\n";
		}
	  scratchlogfile << "\n" << "\nDeltaline_slave00_poly:                    \t" << deltaline_slave00 << "\nDeltapixel_slave00_poly:                   \t" << deltapixel_slave00 << "\nDeltaline_slave0N_poly:                    \t" << deltaline_slave0N << "\nDeltapixel_slave0N_poly:                   \t" << deltapixel_slave0N << "\nDeltaline_slaveN0_poly:                    \t" << deltaline_slaveN0 << "\nDeltapixel_slaveN0_poly:                   \t" << deltapixel_slaveN0 << "\nDeltaline_slaveNN_poly:                    \t" << deltaline_slaveNN << "\nDeltapixel_slaveNN_poly:                   \t" << deltapixel_slaveNN;

	  scratchlogfile << "\n*******************************************************************\n";
	  scratchlogfile.close();


	  ofstream scratchresfile = new ofstream("scratchrescpm", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile,"coregpm: scratchrescpm",__FILE__,__LINE__);

	  scratchresfile.setf(ios.scientific, ios.floatfield);
	  scratchresfile.setf(ios.right, ios.adjustfield);
	  scratchresfile.precision(8);
	  scratchresfile.width(18);

	  scratchresfile << "\n\n*******************************************************************" << "\n*_Start_" << processcontrol[(int)AnonymousEnum.pr_i_coregpm] << "\n*******************************************************************" << "\nDegree_cpm:\t" << DEGREE << "\nEstimated_coefficientsL:\n";
	  int coeffL = 0;
	  int coeffP = 0;
	  for (i =0; i<Nunk; i++)
		{
		if (rhsL(i,0) < 0.)
		  scratchresfile << rhsL(i,0);
		else
		  scratchresfile << " " << rhsL(i,0);

		// ______ Add coefficient number behind value ______
		scratchresfile << " \t" << coeffL << " " << coeffP << "\n";
		coeffL--;
		coeffP++;
		if (coeffL == -1)
		  {
		  coeffL = coeffP;
		  coeffP = 0;
		  }
		}

	  coeffL = 0;
	  coeffP = 0;
	  scratchresfile << "\nEstimated_coefficientsP:\n";
	  for (i =0; i<Nunk; i++)
		{
		if (rhsP(i,0) < 0.)
		  scratchresfile << rhsP(i,0);
		else
		  scratchresfile << " " << rhsP(i,0);

		// ______ Add coefficient number behind value ______
		scratchresfile << " \t" << coeffL << " " << coeffP << "\n";
		coeffL--;
		coeffP++;
		if (coeffL == -1)
		  {
		  coeffL = coeffP;
		  coeffP = 0;
		  }
		}
	   scratchresfile << "\nDeltaline_slave00_poly:                    \t" << deltaline_slave00 << "\nDeltapixel_slave00_poly:                   \t" << deltapixel_slave00 << "\nDeltaline_slave0N_poly:                    \t" << deltaline_slave0N << "\nDeltapixel_slave0N_poly:                   \t" << deltapixel_slave0N << "\nDeltaline_slaveN0_poly:                    \t" << deltaline_slaveN0 << "\nDeltapixel_slaveN0_poly:                   \t" << deltapixel_slaveN0 << "\nDeltaline_slaveNN_poly:                    \t" << deltaline_slaveNN << "\nDeltapixel_slaveNN_poly:                   \t" << deltapixel_slaveNN;
					 //<< "\n* End_coregpm:_NORMAL"
	 scratchresfile << "\n*******************************************************************" << "\n* End_" << processcontrol[(int)AnonymousEnum.pr_i_coregpm] << "_NORMAL" << "\n*******************************************************************\n";
	  scratchresfile.close();


	// ====== Compute offsets for corners ======
	// BK 18-May-2000
	  // read rhsL from file due top format... double
	  matrix<Double> Lcoeff = readcoeff("scratchrescpm", "Estimated_coefficientsL:",Ncoeffs(DEGREE));
	  // read rhsP from file due top format... double
	  matrix<Double> Pcoeff = readcoeff("scratchrescpm", "Estimated_coefficientsP:",Ncoeffs(DEGREE));
	  matrix<Float> x_axis = new matrix(2,1);
	  matrix<Float> y_axis = new matrix(2,1);
	  x_axis(0,0) = minL;
	  x_axis(1,0) = maxL;
	  y_axis(0,0) = minP;
	  y_axis(1,0) = maxP;
	  normalize(x_axis,minL,maxL);
	  normalize(y_axis,minP,maxP);
	  matrix<Float> offsetcornersL = polyval<Float>(x_axis,y_axis,Lcoeff); // MA
	  matrix<Float> offsetcornersP = polyval<Float>(x_axis,y_axis,Pcoeff);
	  INFO.print(" ");
	  INFO.print("Modeled transformation in azimuth:");
	  INFO.print("-------------------------------------------------");
	  INFO << "  First line:    " << offsetcornersL(0,0) << " ... " << offsetcornersL(0,1);
	  INFO.print();
	  INFO.print("                    :           :");
	  INFO << "  Last line:     " << offsetcornersL(1,0) << " ... " << offsetcornersL(1,1);
	  INFO.print();
	  INFO.print("\n");
	  INFO.print("Modeled transformation in range:");
	  INFO.print("-------------------------------------------------");
	  INFO << "  First line:    " << offsetcornersP(0,0) << " ... " << offsetcornersP(0,1);
	  INFO.print();
	  INFO.print("                    :           :");
	  INFO << "  Last line:     " << offsetcornersP(1,0) << " ... " << offsetcornersP(1,1);
	  INFO.print();
	  INFO.print(" ");


	  // ====== Dump evaluated polynomial if requested ======
	  // BK 17-May-2000
	  if (coregpminput.dumpmodel)
		{
		DEBUG.print("Do evaluation of coreg model with stepsize 100 pixels or so...");
		DEBUG.print("And account for currentwindow, not orig window...");
		PROGRESS.print("Started dumping evaluated model azimuth.");
		TRACE.print(); // empty buffer to be sure
		TRACE << "offsetazi_" << master.originalwindow.lines() << "_" << master.originalwindow.pixels() << ".r4";
		String fileazi = new String(new char[ONE27]);
		fileazi = TRACE.get_str();
		TRACE.print(); // empty buffer to be sure

		ofstream dumpfile;
		openfstream(dumpfile,fileazi,true);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		bk_assert(dumpfile,fileazi,__FILE__,__LINE__);

		// polyval both standing x,y... (?)
		// matrix<real4> p_axis(1,master.originalwindow.pixels());
		matrix<Float> l_axis = new matrix(1,1); // ...
		matrix<Float> p_axis = new matrix(master.originalwindow.pixels(),1);
		for (i =0; i<p_axis.pixels(); ++i)
		  p_axis(i,0) = i+master.originalwindow.pixlo; // multilook==1 ?
		normalize(p_axis,minP,maxP);

		// azimuth
		for (i =master.originalwindow.linelo; i<=master.originalwindow.linehi; ++i) // all lines
		  {
		  l_axis(0,0) = normalize((double)i,minL,maxL);
		  matrix<Float> MODEL = polyval<Float>(l_axis,p_axis,Lcoeff);
		  dumpfile << MODEL;
		  }
		dumpfile.close();
		INFO << "Dumped model azimuth offset to file: " << fileazi << " format: real4; number of lines: " << master.originalwindow.lines() << " number of pixels: " << master.originalwindow.pixels();
		INFO.print();

		// ______ same for range ______
		PROGRESS.print("Started dumping evaluated model range.");
		TRACE.print(); // empty buffer to be sure
		TRACE << "offsetrange_" << master.originalwindow.lines() << "_" << master.originalwindow.pixels() << ".r4";
		String filerange = new String(new char[ONE27]);
		filerange = TRACE.get_str();
		TRACE.print(); // empty buffer to be sure
		ofstream dumpfile2;
		openfstream(dumpfile2,filerange,true);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		bk_assert(dumpfile2,filerange,__FILE__,__LINE__);

		for (i =master.originalwindow.linelo; i<=master.originalwindow.linehi; ++i) // all lines
		  {
		  l_axis(0,0) = normalize((double)i,minL,maxL);
		  matrix<Float> MODEL = polyval<Float>(l_axis,p_axis,Pcoeff);
		  dumpfile2 << MODEL;
		  }
		INFO << "Dumped model range offset to file: " << filerange << " format: real4; number of lines: " << master.originalwindow.lines() << " number of pixels: " << master.originalwindow.pixels();
		INFO.print();
		dumpfile2.close();
		}

	  // ====== Tidy up ======
	  PROGRESS.print("finished computation of coregistration parameters.");
	  } // END coregpm

//***************************************************************
// *    getofffile                                                *
// *                                                              *
// * Returns matrix (real4) with data of fine coreg from file     *
// *   mat(i,0)=window number                                     *
// *   mat(i,1)=position: line coordinate                         *
// *   mat(i,2)=position: pixels coordinate                       *
// *   mat(i,3)=offset: line direction                            *
// *   mat(i,4)=offset: pixle direction                           *
// *   mat(i,5)=correlation:                                      *
// * searches for "Number_of_correlation_windows:"                *
// *                                                              *
// *    Bert Kampes, 24-Feb-1999                                  *
// ***************************************************************
			//const uint             oversamplingsfactor);


	// ______ Read observations from file ______
	public static matrix<Float> getofffile(String file, float threshold)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"getofffile (BK 24-Feb-1999)"<<ends;
		  TRACE.print();
	  }
	  String dummyline = new String(new char[ONE27]);
	  String word = new String(new char[EIGHTY]);
	  boolean foundsection = false;

	  ifstream infile;
	  openfstream(infile,file);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(infile,file,__FILE__,__LINE__);


	  // ======Search file for data section======
	  while (infile != null)
		{
		infile >> word;
		if (strcmp("Number_of_correlation_windows:",word)) // no pattern match.
		  {
		  infile.getline(dummyline,ONE27,'\n'); // goto next line.
		  }
		else // in data section
		  {
		  foundsection =true;
		  int N; // number of points
		  infile >> N;
		  infile.getline(dummyline,ONE27,'\n'); // next line
		  infile.getline(dummyline,ONE27,'\n'); // skip line with info
		  int pos = infile.tellg(); // position of start data
		  int Nobs = 0; // number points > threshold
		  float winnumber; // on file
		  float posL;
		  float posP;
		  float offL;
		  float offP;
		  float corr;
//C++ TO JAVA CONVERTER NOTE: 'register' variable declarations are not supported in Java:
//ORIGINAL LINE: register int i;
		  int i;
		  for (i =0;i<N;i++)
			{
			infile >> winnumber >> posL >> posP >> offL >> offP >> corr;
			infile.getline(dummyline,ONE27,'\n'); // goto next data record
			if (corr > threshold)
			  Nobs++;
			}

		  if (Nobs == 0)
			{
			{
				ERROR.terminate();
				String cp_s = new String(new char[256]);
				cp_s = "code ???: No data found > threshold.";
				ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
				ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
				ERROR.print();
			}
			throw(some_error);
			}

		  matrix<Float> Data = new matrix(Nobs,6);
		  infile.seekg(pos); // return to start data
		  int cnti = -1;
		  for (i =0;i<N;i++)
			{
			infile >> winnumber >> posL >> posP >> offL >> offP >> corr;
			infile.getline(dummyline,ONE27,'\n'); // goto next data record
			if (corr > threshold)
			  {
			  cnti++;
			  Data(cnti,0) = winnumber;
			  Data(cnti,1) = posL;
			  Data(cnti,2) = posP;
			  Data(cnti,3) = offL;
			  Data(cnti,4) = offP;
			  Data(cnti,5) = corr;
			  }
			}

		  infile.close();
		  return Data;
		  } // else
		} // file

	  // ______Tidy up______
	  if (!foundsection)
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "code 401: getofffile: couldn't find data section in file.";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(some_error);
		}
	  infile.close();

	  // --- return a dummy here since some compiler like that ---
	  return matrix<Float>(999,999); // BK 07-Apr-2003
	  } // END getofffile

//***************************************************************
// *    resample                                                  *
// *                                                              *
// * Resample slave to master grid based on coregistration        *
// *  parameters.                                                 *
// * if dbow==0 then default to overlap, else dbow,               *
// *  write 0's where it does not overlap                         *
// * (later at interf.comp. if master<slave then doris exits! bug)*
// *                                                              *
// * input:                                                       *
// *  - inputoptions                                              *
// * output:                                                      *
// *  - void (file)                                               *
// *                                                              *
// *    Bert Kampes, 16-Mar-1999                                  *
// * added DBOW master add zeros.                                 *
// #%// BK 21-Aug-2000BOW master add zeros.                       *
// * shift data to center of spectrum before resampling, shift    *
// * back afterwards. (see e.g. thesis Geudtner)                  *
// #%// BK 09-Nov-2000                                            *
// * Seems to be a bug in shifting the data spectrum if more      *
// * buffers are used, working on it.                             *
// * (Increase FORSURE variable if crash)                         *
// #%// BK 19-Nov-2000                                            *
// ***************************************************************


	// ______ Resample slave ______
	public static void resample(input_gen generalinput, input_resample resampleinput, slcimage master, slcimage slave, matrix<Double> cpmL, matrix<Double> cpmP, short demassist)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"resample (BK 16-Mar-1999; BK 09-Nov-2000)"<<ends;
		  TRACE.print();
	  }
	  if (resampleinput.shiftazi==true)
		DEBUG.print("shifting kernelL to data fDC BK 26-Oct-2002");
	  // ___ Handle input ___
	  final int BUFFERMEMSIZE = generalinput.memory; // Bytes 500MB --> 500 000 000 bytes
	  final int Npoints = resampleinput.method%100; // #pnts interpolator
	  if (isodd(Npoints))
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "resample only even point interpolators.";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}
	  final int Npointsd2 = Npoints/2;
	  final int Npointsd2m1 = Npointsd2-1;
	  //const uint  Sfilelines   = slave.currentwindow.lines();
	  final int sizeofci16 = sizeof(complex<Short>);
	  final int sizeofcr4 = sizeof(complex<Float>);

	  // ______ Normalize data for polynomial ______
	  final double minL = master.originalwindow.linelo;
	  final double maxL = master.originalwindow.linehi;
	  final double minP = master.originalwindow.pixlo;
	  final double maxP = master.originalwindow.pixhi;
	  INFO << "resample: polynomial normalized by factors: " << minL << " " << maxL << " " << minP << " " << maxP << " to [-2,2]";
	  INFO.print();

	  // ______ For KNAB/Raised Cosine kernel if requested ______
	  // ______ Because kernel is same in az. and rg. min. must be used.
	  final float CHI_az = slave.prf/slave.abw; // oversampling factor az
	  final float CHI_rg = (slave.rsr2x/2.0)/slave.rbw; // oversampling factor rg
	  final float CHI = min(CHI_az,CHI_rg); // min. oversampling factor of data
	  INFO << "Oversampling ratio azimuth (PRF/ABW):    " << CHI_az;
	  INFO.print();
	  INFO << "Oversampling ratio azimuth (RSR/RBW):    " << CHI_rg;
	  INFO.print();
	  INFO << "KNAB/RC kernel uses: oversampling ratio: " << CHI;
	  INFO.print();
	  if (CHI < 1.1)
		{
		WARNING << "Oversampling ratio: " << CHI << " not optimal for KNAB/RC";
		WARNING.print();
		}


	  // ====== Create lookup table ======
	  // ______ e.g. four point interpolator
	  // ______ interpolating point: p=6.4925
	  // ______ required points: 5, 6, 7, 8
	  // ______ kernel number from lookup table: floor(.4925*INTERVAL+.5)
	  // ______  table[0]= 0 1 0 0 ;table[INTERVAL]= 0 0 1 0
	  // ______ intervals in lookup table: dx
	  // ______ for high doppler 100 is OK (fdc=3prf; 6pi --> 10deg error?)
	  final int INTERVAL = 127; // precision: 1./INTERVAL [pixel]
	  final int Ninterval = INTERVAL + 1; // size of lookup table
	  final double dx = 1.0/INTERVAL; // interval look up table
	  INFO << "resample: lookup table size: " << Ninterval;
	  INFO.print();

//C++ TO JAVA CONVERTER NOTE: 'register' variable declarations are not supported in Java:
//ORIGINAL LINE: register int i;
	  int i;
	  matrix<Float> x_axis = new matrix(Npoints,1);
	  for (i =0; i<Npoints; ++i)
		x_axis(i,0) = 1.0 - Npointsd2 + i; // start at [-1 0 1 2]

	  // ______ Lookup table complex because of multiplication with complex ______
	  // ______ Loopkup table for azimuth and range and ______
	  // ______ shift spectrum of azi kernel with doppler centroid ______
	  // ______ kernel in azimuth should be sampled higer ______
	  // ______ and may be different from range due to different ______
	  // ______ oversampling ratio and spectral shift (const) ______
	  matrix<complex<Float>>[] pntKernelAz = new matrix[Ninterval]; // kernel in azimuth
	  matrix<complex<Float>>[] pntKernelRg = new matrix[Ninterval]; // kernel in range
	  // ______ same axis required for shift azimuth spectrum as used ______
	  // ______ for kernel to avoid phase shift (Raffaele Nutricato) ______
	  matrix<Float>[] pntAxis = new matrix[Ninterval];

	  for (i =0; i<Ninterval; ++i)
		{
		pntKernelAz[i] = new matrix<complex<Float>> (Npoints,1);
		pntKernelRg[i] = new matrix<complex<Float>> (Npoints,1);
		pntAxis[i] = new matrix<Float> (Npoints,1); // only used for azishift
		switch(resampleinput.method)
		  {
		  // --- Extremely simple kernels (not good, but fast) ---
		  case rs_rect:
			(pntKernelAz[i]) = mat2cr4rect(x_axis);
			(pntKernelRg[i]) = mat2cr4rect(x_axis);
			break;
		  case rs_tri:
			(pntKernelAz[i]) = mat2cr4tri(x_axis);
			(pntKernelRg[i]) = mat2cr4tri(x_axis);
			break;
		  // --- Truncated sinc ---
		  case rs_ts6p:
			(pntKernelAz[i]) = mat2cr4ts6(x_axis);
			(pntKernelRg[i]) = mat2cr4ts6(x_axis);
			break;
		  case rs_ts8p:
			(pntKernelAz[i]) = mat2cr4ts8(x_axis);
			(pntKernelRg[i]) = mat2cr4ts8(x_axis);
			break;
		  case rs_ts16p:
			(pntKernelAz[i]) = mat2cr4ts16(x_axis);
			(pntKernelRg[i]) = mat2cr4ts16(x_axis);
			break;
		  // --- Cubic Convolution kernel: theoretical better than truncated sinc. ---
		  case rs_cc4p:
			(pntKernelAz[i]) = mat2cr4cc4(x_axis);
			(pntKernelRg[i]) = mat2cr4cc4(x_axis);
			break;
		  case rs_cc6p:
			(pntKernelAz[i]) = mat2cr4cc6(x_axis);
			(pntKernelRg[i]) = mat2cr4cc6(x_axis);
			break;
		  // --- KNAB kernel: theoretical better than cubic conv. ---
		  case rs_knab4p:
			(pntKernelAz[i]) = mat2cr4knab(x_axis, CHI_az, 4);
			(pntKernelRg[i]) = mat2cr4knab(x_axis, CHI_rg, 4);
			break;
		  case rs_knab6p:
			(pntKernelAz[i]) = mat2cr4knab(x_axis, CHI_az, 6);
			(pntKernelRg[i]) = mat2cr4knab(x_axis, CHI_rg, 6);
			break;
		  case rs_knab8p:
			(pntKernelAz[i]) = mat2cr4knab(x_axis, CHI_az, 8);
			(pntKernelRg[i]) = mat2cr4knab(x_axis, CHI_rg, 8);
			break;
		  case rs_knab10p:
			(pntKernelAz[i]) = mat2cr4knab(x_axis, CHI_az, 10);
			(pntKernelRg[i]) = mat2cr4knab(x_axis, CHI_rg, 10);
			break;
		  case rs_knab16p:
			(pntKernelAz[i]) = mat2cr4knab(x_axis, CHI_az, 16);
			(pntKernelRg[i]) = mat2cr4knab(x_axis, CHI_rg, 16);
			break;
		  // --- Raised cosine: theoretical best ---
		  case rs_rc6p:
			(pntKernelAz[i]) = mat2cr4rc_kernel(x_axis, CHI_az, 6);
			(pntKernelRg[i]) = mat2cr4rc_kernel(x_axis, CHI_rg, 6);
			break;
		  case rs_rc12p:
			(pntKernelAz[i]) = mat2cr4rc_kernel(x_axis, CHI_az, 12);
			(pntKernelRg[i]) = mat2cr4rc_kernel(x_axis, CHI_rg, 12);
			break;
		  default:
			{
				ERROR.terminate();
				String cp_s = new String(new char[256]);
				cp_s = "impossible.";
				ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
				ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
				ERROR.print();
			}
			throw(unhandled_case_error);
		  } //kernel selector
		(pntAxis[i]) = x_axis; // to shift kernelL use: k*=exp(-i*2pi*axis*fdc/prf)
		x_axis -= dx; // Note: 'wrong' way (mirrored)
		}
	  // ====== Usage: pntKernelAz[0]->showdata(); or (*pntKernelAz[0][0]).showdata(); ======
	  // ______ Log kernels to check sum, etc. ______
	  DEBUG.print("Overview of LUT for interpolation kernel follows:");
	  DEBUG.print("-------------------------------------------------");

	  for (i =0; i<Ninterval; ++i)
		{
		for (int x =0; x<Npoints; ++x)
//C++ TO JAVA CONVERTER TODO TASK: There are no simple equivalents to function pointers in Java:
	//	 DEBUG << ((pntAxis[i])(x,0)) << "      ";
		DEBUG.print();
		float sum_az = 0.0;
		float sum_rg = 0.0;
		for (int x =0; x<Npoints; ++x)
		  {
//C++ TO JAVA CONVERTER TODO TASK: There are no simple equivalents to function pointers in Java:
	//	 DEBUG << real((pntKernelAz[i])(x,0)) << " "; // complex kernel
//C++ TO JAVA CONVERTER TODO TASK: There are no simple equivalents to function pointers in Java:
	//	 sum_az += real((pntKernelAz[i])(x,0));
//C++ TO JAVA CONVERTER TODO TASK: There are no simple equivalents to function pointers in Java:
	//	 sum_rg += real((pntKernelRg[i])(x,0));
		  }
		DEBUG << "(sum=" << sum_az << ")";
		DEBUG.print();
		DEBUG.print("Normalizing kernel by dividing LUT elements by sum:");
		(pntKernelAz[i]) /= sum_az;
		(pntKernelRg[i]) /= sum_rg;
		// ______ Only show azimuth kernel ______
		for (int x =0; x<Npoints; ++x)
//C++ TO JAVA CONVERTER TODO TASK: There are no simple equivalents to function pointers in Java:
	//	 DEBUG << real((pntKernelAz[i])(x,0)) << " "; // complex kernel; normalized
		DEBUG.print();
		}
	  PROGRESS.print("Resample: normalized lookup table created (kernel and axis).");

	  // ______Save some time by computing degree here______
	  final int degree_cpmL = degree(cpmL.size());
	  final int degree_cpmP = degree(cpmP.size());

	  // ______ Initialization (needed for DEM assist) [FvL] _____
	  double deltaL_dem;
	  double deltaP_dem;
	  float deltaL_poly;
	  float deltaP_poly;
	  float interpL;
	  float interpP;
	  float ms_az_timing_error_L = slave.az_timing_error;
	  float ms_r_timing_error_P = slave.r_timing_error;
	  final int sizer8 = sizeof(double);
	  ifstream DeltaLfile;
	  ifstream DeltaPfile;

	  if (demassist != 0)
		{
		  openfstream(DeltaLfile,"delta_line.raw");
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  bk_assert(DeltaLfile,"delta_line.raw",__FILE__,__LINE__);
		  openfstream(DeltaPfile,"delta_pixel.raw");
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  bk_assert(DeltaPfile,"delta_pixel.raw",__FILE__,__LINE__);
		}

	  streampos pos;

	  // ______Corners of overlap in master system______
	  // changed by FvL

	  window overlap;
	  if (demassist != 0)
		overlap = getoverlap(master,slave,(double)Npointsd2,(double)ms_az_timing_error_L,(double)ms_r_timing_error_P);
	  else
		overlap = getoverlap(master,slave,(double)Npointsd2,(double)0,(double)0);


	  // ====== Adjust overlap possibly for RS_DBOW card ======
	  int write0lines1 = 0; // DBOW card, 0's at start
	  int write0linesN = 0;
	  int write0pixels1 = 0;
	  int write0pixelsN = 0;
	  if (!(resampleinput.dbow.linelo == 0 && resampleinput.dbow.linehi == 0 && resampleinput.dbow.pixlo == 0 && resampleinput.dbow.pixhi == 0))
		{
		// ______ Check if overlap is large enough to contain DBOW ______
		if (resampleinput.dbow.linelo > overlap.linehi)
		  {
		  {
			  ERROR.terminate();
			  String cp_s = new String(new char[256]);
			  cp_s = "RS_DBOW: specified min. line larger than max. line of overlap.";
			  ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			  ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			  ERROR.print();
		  }
		  throw(input_error);
		  }
		if (resampleinput.dbow.linehi < overlap.linelo)
		  {
		  {
			  ERROR.terminate();
			  String cp_s = new String(new char[256]);
			  cp_s = "RS_DBOW: specified max. line smaller than min. line of overlap.";
			  ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			  ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			  ERROR.print();
		  }
		  throw(input_error);
		  }
		if (resampleinput.dbow.pixlo > overlap.pixhi)
		  {
		  {
			  ERROR.terminate();
			  String cp_s = new String(new char[256]);
			  cp_s = "RS_DBOW: specified min. pixel larger than max. pixel of overlap.";
			  ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			  ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			  ERROR.print();
		  }
		  throw(input_error);
		  }
		if (resampleinput.dbow.pixhi < overlap.pixlo)
		  {
		  {
			  ERROR.terminate();
			  String cp_s = new String(new char[256]);
			  cp_s = "RS_DBOW: specified max. pixel smaller than min. pixel of overlap.";
			  ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			  ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			  ERROR.print();
		  }
		  throw(input_error);
		  }

		write0lines1 = overlap.linelo - resampleinput.dbow.linelo;

		if (write0lines1 < 0) // smaller window selected
			write0lines1 = 0;
		write0linesN = -overlap.linehi + resampleinput.dbow.linehi;

		if (write0linesN < 0) // smaller window selected
			write0linesN = 0;
		write0pixels1 = overlap.pixlo - resampleinput.dbow.pixlo;

		if (write0pixels1 < 0) // smaller window selected
			write0pixels1 = 0;
		write0pixelsN = -overlap.pixhi + resampleinput.dbow.pixhi;

		if (write0pixelsN < 0) // smaller window selected
			write0pixelsN = 0;

		if (resampleinput.dbow.linelo < overlap.linelo)
		  {
		  WARNING << "RS_DBOW: min. line < overlap (writing: " << write0lines1 << " lines with zeros before first resampled line).";
		  WARNING.print();
		  }
		else
		  overlap.linelo = resampleinput.dbow.linelo; // correct it
		if (resampleinput.dbow.linehi > overlap.linehi)
		  {
		  WARNING << "RS_DBOW: max. line > overlap (writing: " << write0linesN << " lines with zeros after last resampled line).";
		  WARNING.print();
		  }
		else
		  overlap.linehi = resampleinput.dbow.linehi; // correct it

		if (resampleinput.dbow.pixlo < overlap.pixlo)
		  {
		  WARNING << "RS_DBOW: min. pixel < overlap (writing: " << write0pixels1 << " columns with zeros before first resampled column).";
		  WARNING.print();
		  }
		else
		  overlap.pixlo = resampleinput.dbow.pixlo; // correct it

		if (resampleinput.dbow.pixhi > overlap.pixhi)
		  {
		  WARNING << "RS_DBOW: max. pixel > overlap (writing: " << write0pixelsN << " columns with zeros after last resampled column).";
		  WARNING.print();
		  }
		else
		  overlap.pixhi = resampleinput.dbow.pixhi; // correct it

		} // adjust overlap


	  // ______ Buffersize output matrix ______
	  final int Npointsxsize = Npoints *sizeofcr4; // size for memcpy (fill PART)
	  final int npixels = slave.currentwindow.pixels();
	  final double bytesperline = sizeofcr4 * npixels;
	  // ___ COMMENTED OUT, OLD WAY SHIFT DATA, now shiftkernel ___
	  final double bigmatrices = 2.5; // BUFFER, RESULT & PART buffers
	  final int nlines = (BUFFERMEMSIZE/bigmatrices)/bytesperline; // buffer nlines

	  // ______ Declare/allocate matrices ______
	  matrix<complex<Float>> BUFFER; // load after output is written
	  matrix<complex<Float>> RESULT = new matrix(nlines,overlap.pixhi-overlap.pixlo+1);
	  matrix<complex<Float>> PART = new matrix(Npoints,Npoints);

	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if __USE_VECLIB_LIBRARY__
	  matrix<complex<Float>> TMPRES = new matrix(Npoints,1);
	  int Np = Npoints; // must be non-constant
	  int ONEint = 1; // must have pointer to 1
	  complex<Float> c4alpha = new complex(1.,0.0);
	  complex<Float> c4beta = new complex(0.0,0.0);
	  STUPID_cr4 ANS = new STUPID_cr4(); // VECLIB struct return type ofstream ofile;
	//#else
	  ofstream ofile;
	//#endif


	  // ====== Open output file ======
	  openfstream(ofile,resampleinput.fileout,generalinput.overwrit);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(ofile,resampleinput.fileout,__FILE__,__LINE__);

	  // ________ First write zero lines if appropriate (DBOW) ______
	  switch (resampleinput.oformatflag)
		{
		case FORMATCR4:
		  {
		  final complex<Float> zerocr4 = new complex(0,0);
		  for (int thisline =0; thisline<write0lines1; ++thisline)
			for (int thispixel =0; thispixel<(int)(RESULT.pixels())+write0pixels1+write0pixelsN; ++thispixel)
			  ofile.write((char) zerocr4, sizeofcr4);
		  break;
		  }
		case FORMATCI2:
		  {
		  final complex<Short> zeroci16 = new complex(0,0);
		  for (int thisline =0; thisline<write0lines1; ++thisline)
			for (int thispixel =0; thispixel<(int)(RESULT.pixels())+write0pixels1+write0pixelsN; ++thispixel)
			  ofile.write((char) zeroci16, sizeofci16);
		  break;
		  }
		default:
		  {
			  ERROR.terminate();
			  String cp_s = new String(new char[256]);
			  cp_s = "impossible format";
			  ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			  ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			  ERROR.print();
		  }
		  throw(unhandled_case_error);
		}

	  // ______ Info ______
	  INFO << "Overlap window: " << overlap.linelo << ":" << overlap.linehi << ", " << overlap.pixlo << ":" << overlap.pixhi;
	  INFO.print();


	  // ______ Progress messages ______
	  int percent = 0;
	  int tenpercent = rint(overlap.lines()/10.0); // round
	  if (tenpercent ==0) // avoid error: x%0
		  tenpercent = 1000;

	  // ====== Resample all lines that are requested ======
	  boolean newbufferrequired = true; // read initial slave buffer
//C++ TO JAVA CONVERTER NOTE: 'register' variable declarations are not supported in Java:
//ORIGINAL LINE: register int linecnt = -1;
	  int linecnt = -1; // indicate output buffer full
	  int firstline = 0; // slave system
	  int lastline = 0;

//C++ TO JAVA CONVERTER NOTE: 'register' variable declarations are not supported in Java:
//ORIGINAL LINE: register int line;
	  int line; // loop counter master system
//C++ TO JAVA CONVERTER NOTE: 'register' variable declarations are not supported in Java:
//ORIGINAL LINE: register int pixel;
	  int pixel; // loop counter master system
	  for (line =overlap.linelo; line<=overlap.linehi; line++)
		{
		// ______ Progress messages ______
		if (((line-overlap.linelo)%tenpercent)==0)
		  {
		  PROGRESS << "RESAMPLE: " << setw(3) << percent << "%";
		  PROGRESS.print();
		  percent += 10;
		  }

		// ====== Write RESULT to disk if it is full (write last bit at end) ======
		if (linecnt ==(int)(RESULT.lines())-1) // ==nlines
		  {
		  newbufferrequired = true; // do load slave from file
		  DEBUG << "Writing slave: [" << line-RESULT.lines() << ":" << line-1 << ", " << overlap.pixlo << ":" << overlap.pixhi << "] (master coord. system)";
		  DEBUG.print();
		  linecnt = 0;
		  // ______ Actually write ______
		  switch (resampleinput.oformatflag)
			{
			case FORMATCR4:
			  {
			  // old, now first write zeropixels...: ofile << RESULT;
			  final complex<Float> zerocr4 = new complex(0.0, 0.0);
			  for (int thisline =0; thisline<(int)(RESULT.lines()); ++thisline)
				{
				// ______ Write zero pixels at start ______
				for (int thispixel =0; thispixel<write0pixels1; ++thispixel)
				  {
				  ofile.write((char) zerocr4, sizeofcr4);
				  }
				// ______ WRITE the interpolated data per row ______
				ofile.write((char) RESULT[thisline][0], RESULT.pixels()*sizeof(RESULT(0,0)));
				// ______ Write zero pixels at end ______
				for (int thispixel =0; thispixel<write0pixelsN; ++thispixel)
				  {
				  ofile.write((char) zerocr4, sizeofcr4);
				  }
				}
			  break;
			  }
			case FORMATCI2:
			  {
			  final complex<Short> zeroci16 = new complex(0,0);
			  complex<Short> castedresult;
			  for (int thisline =0; thisline<(int)(RESULT.lines()); ++thisline)
				{
				// ______ Write zero pixels at start ______
				for (int thispixel =0; thispixel<write0pixels1; ++thispixel)
				  {
				  ofile.write((char) zeroci16, sizeofci16);
				  }
				// ______ Write the interpolated data per row ______
				for (int thispixel =0; thispixel<(int)(RESULT.pixels()); ++thispixel)
				  {
				  // no default conversion, this seems slow, test this (BK)
				  castedresult = cr4toci2(RESULT(thisline,thispixel));
				  ofile.write((char) castedresult, sizeofci16);
				  }
				// ______ Write zero pixels at end ______
				for (int thispixel =0; thispixel<write0pixelsN; ++thispixel)
				  {
				  ofile.write((char) zeroci16, sizeofci16);
				  }
				}
			  break;
			  }
			default:
			  {
				  ERROR.terminate();
				  String cp_s = new String(new char[256]);
				  cp_s = "impossible format";
				  ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
				  ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
				  ERROR.print();
			  }
			  throw(unhandled_case_error);
			}
		  }

		else // output buffer not full yet
		  {
		  linecnt++;
		  }

		// ====== Read slave buffer if justwritten || firstblock ======
		if (newbufferrequired ==true)
		  {
		  newbufferrequired = false; // only load after output
											// written
		  if (demassist != 0)
			{
			  pos = (streampos)((line-master.currentwindow.linelo)*master.currentwindow.pixels() + overlap.pixlo - master.currentwindow.pixlo);
			  pos = (streampos)(pos * sizer8);
			  DeltaLfile.seekg(pos,ios.beg); // [MA] better to check for failbit
			  DeltaLfile.read((char)&deltaL_dem, sizer8);

			  deltaL_poly = polyval(normalize((float)line,minL,maxL), normalize((float)overlap.pixlo,minP,maxP), cpmL,degree_cpmL);

			  float firstline_pixlo = line + deltaL_dem + deltaL_poly + ms_az_timing_error_L;

			  pos = (streampos)((line-master.currentwindow.linelo)*master.currentwindow.pixels() + overlap.pixhi - master.currentwindow.pixlo);
			  pos = (streampos)(pos * sizer8);
			  DeltaLfile.seekg(pos,ios.beg);
			  DeltaLfile.read((char)&deltaL_dem, sizer8);

			  deltaL_poly = polyval(normalize((float)line,minL,maxL), normalize((float)overlap.pixhi,minP,maxP), cpmL,degree_cpmL);

			  float firstline_pixhi = line + deltaL_dem + deltaL_poly + ms_az_timing_error_L;


			  int line2 = line + nlines - 1;

			  // LAST BUFFER FIX
			  // [DON] Davide Nitti,  the overrun of last line due to buffer nlines.
			  // start added by don
			  if (line2 > (int)master.currentwindow.linehi)
			  {
				 DEBUG << "Variable line2: [ACTUAL Value: " << line2 << " - NEW Value: " << master.currentwindow.linehi << "]";
				 DEBUG.print();
				 line2 = master.currentwindow.linehi;
			  }
			  // end added by don

			  pos = (streampos)((line2-master.currentwindow.linelo)*master.currentwindow.pixels() + overlap.pixlo - master.currentwindow.pixlo);
			  pos = (streampos)(pos * sizer8);
			  DeltaLfile.seekg(pos,ios.beg);
			  DeltaLfile.read((char)&deltaL_dem, sizeof(deltaL_dem)); // [MA] sizer8 --> sizeof(deltaL_dem)

			  if (DeltaLfile.fail()) // [MA] put it to a proper class
			  {
			  WARNING << "Failed to read position: " << pos; // coherence will be lost in lastbuffer
			  WARNING.print();
			  // exit(1) 
			  }

			  deltaL_poly = polyval(normalize((float)line2,minL,maxL), normalize((float)overlap.pixlo,minP,maxP), cpmL,degree_cpmL);

			  float lastline_pixlo = (float)(line2 + deltaL_dem + deltaL_poly + ms_az_timing_error_L);

			  pos = (streampos)((line2-master.currentwindow.linelo)*master.currentwindow.pixels() + overlap.pixhi - master.currentwindow.pixlo);
			  pos = (streampos)(pos * sizer8);
			  DeltaLfile.seekg(pos,ios.beg);
			  DeltaLfile.read((char)&deltaL_dem, sizer8);

			  deltaL_poly = polyval(normalize((float)line2,minL,maxL), normalize((float)overlap.pixhi,minP,maxP), cpmL,degree_cpmL);

			  float lastline_pixhi = (float)(line2 + deltaL_dem + deltaL_poly + ms_az_timing_error_L);

			  firstline = (int)(Math.ceil(min(firstline_pixlo,firstline_pixhi)))-Npoints;
			  lastline = (int)(Math.ceil(min(lastline_pixlo,lastline_pixhi)))+Npoints;
			  }
		  else
			{
			  firstline = (int)(Math.ceil(min(line + polyval(normalize((float)line,minL,maxL), normalize((float)overlap.pixlo,minP,maxP), cpmL,degree_cpmL), line + polyval(normalize((float)line,minL,maxL), normalize((float)overlap.pixhi,minP,maxP), cpmL,degree_cpmL)))) - Npoints;
			  int line2 = line + nlines - 1;
			  lastline = (int)(Math.ceil(min(line2 + polyval(normalize((float)line2,minL,maxL), normalize((float)overlap.pixlo,minP,maxP), cpmL,degree_cpmL), line2 + polyval(normalize((float)line2,minL,maxL), normalize((float)overlap.pixhi,minP,maxP), cpmL,degree_cpmL)))) + Npoints;
			}

		  //const int32 FORSURE = 5;                // buffer larger 2*FORSURE start/end
		  final int FORSURE = 25; // buffer larger 2*FORSURE start/end
		  firstline -= FORSURE;
		  lastline += FORSURE;

		  // ______ Don't compare apples with pears, uint<->int! ______
		  if (firstline < (int)slave.currentwindow.linelo)
			firstline = slave.currentwindow.linelo;

		  if (lastline > (int)slave.currentwindow.linehi)
			lastline = slave.currentwindow.linehi;
		  // ______ Fill slave BUFFER from disk ______
		  window winslavefile = new window(firstline, lastline, slave.currentwindow.pixlo, slave.currentwindow.pixhi);
		  DEBUG << "Reading slave: [" << winslavefile.linelo << ":" << winslavefile.linehi << ", " << winslavefile.pixlo << ":" << winslavefile.pixhi << "]";
		  DEBUG.print();
		  BUFFER = slave.readdata(winslavefile);
		  } // ___end: Read new slave buffer to resample outputbuffer


		// ====== Actual resample all pixels this output line ======
		for (pixel =overlap.pixlo; pixel<=(int)overlap.pixhi; pixel++)
		  {
			if (demassist != 0)
			  {

				//pos = overlap.pixels() * ( line - overlap.linelo ) + pixel - overlap.pixlo;
				pos = (streampos)((line-master.currentwindow.linelo)*master.currentwindow.pixels() + pixel - master.currentwindow.pixlo);
				pos = (streampos)(pos * sizer8);

				DeltaLfile.seekg(pos,ios.beg);
				DeltaPfile.seekg(pos,ios.beg);

				DeltaLfile.read((char)&deltaL_dem, sizer8);
				DeltaPfile.read((char)&deltaP_dem, sizer8);

				deltaL_poly = polyval(normalize((float)line,minL,maxL), normalize((float)pixel,minP,maxP), cpmL,degree_cpmL);
				deltaP_poly = polyval(normalize((float)line,minL,maxL), normalize((float)pixel,minP,maxP), cpmP,degree_cpmP);

				interpL = (float)(line + deltaL_dem + deltaL_poly + ms_az_timing_error_L);
				interpP = (float)(pixel + deltaP_dem + deltaP_poly + ms_r_timing_error_P);

			  }
			else
			  {

				// ______ Evaluate coregistration polynomial ______
				// bk 25-10-99 why don't i do this per buffer, that's faster. (but more mem)
				//interpL = line  + polyval(line,pixel,cpmL,degree_cpmL); // e.g. 255.35432
				//interpP = pixel + polyval(line,pixel,cpmP,degree_cpmP); // e.g. 2.5232
				// ______ BK USE normalized coordinates, do this smarter .... !!!!
				interpL = line + polyval(normalize((float)line,minL,maxL), normalize((float)pixel,minP,maxP), cpmL,degree_cpmL); // e.g. 255.35432
				interpP = pixel + polyval(normalize((float)line,minL,maxL), normalize((float)pixel,minP,maxP), cpmP,degree_cpmP); // e.g. 2.5232
			  }


		  // ______ Get correct lines for interpolation ______
		  final int fl_interpL = interpL;
		  final int fl_interpP = interpP;
		  final int firstL = fl_interpL - Npointsd2m1; // e.g. 254 (5 6 7)
		  final int firstP = fl_interpP - Npointsd2m1; // e.g. 1 (2 3 4)
		  final float interpLdec = interpL - fl_interpL; // e.g. .35432
		  final float interpPdec = interpP - fl_interpP; // e.g. .5232

		  // ______ Copy kernels here, change kernelL if required _ // BK 26-Oct-2002
		  // ______ Faster to have two kernel lookup tables ! _____
		  // ______ I have that now, but still make copy (slow) ______
		  final int kernelnoL = interpLdec *INTERVAL+0.5; // lookup table index
		  final int kernelnoP = interpPdec *INTERVAL+0.5; // lookup table index
		  matrix<complex<Float>> kernelL = (pntKernelAz[kernelnoL]); // local copy to change
		  final matrix<complex<Float>> kernelP = (pntKernelRg[kernelnoP]); // local copy

	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if __DEBUG
				   //maybe modify... [FvL]
		  // ______This shouldn't be possible...______
		  final int Npointsm1 = Npoints-1;
		  if (firstL < slave.currentwindow.linelo)
			{
			WARNING.print("firstL smaller than on disk (required for interpolation). continuing");
			RESULT(linecnt,pixel-overlap.pixlo) = complex<Float>(0.,0.);
			continue; // with next pixel
			}
		  if (firstL+Npointsm1 > slave.currentwindow.linehi)
			{
			WARNING << "lastL larger than on disk (required for interpolation). continuing" << "lineL: " << firstL+Npointsm1 << " > " << slave.currentwindow.linehi;
			WARNING.print();
			RESULT(linecnt,pixel-overlap.pixlo) = complex<Float>(0.,0.);
			continue; // with next pixel
			}
		  if (firstP < slave.currentwindow.pixlo)
			{
			WARNING.print("firstP smaller than on disk (required for interpolation). continuing");
			RESULT(linecnt,pixel-overlap.pixlo) = complex<Float>(0.,0.);
			continue; // with next pixel
			}
		  if (firstP+Npointsm1 > slave.currentwindow.pixhi)
			{
			WARNING.print("lastP larger than on disk (required for interpolation). continuing");
			RESULT(linecnt,pixel-overlap.pixlo) = complex<Float>(0.,0.);
			continue; // with next pixel
			}
	//#endif

		  // ______ Shift azimuth kernel with fDC before interpolation ______
		  if (resampleinput.shiftazi==true)
			{
			// ___ Doppler centroid is function of range only ____
			final float tmp = 2.0 *PI *slave.pix2fdc(interpP)/slave.prf;
			// ___ to shift spectrum of convolution kernel to fDC of data, multiply
			// ___ in the space domain with a phase trend of -2pi*t*fdc/prf
			// ___ (to shift back (no need) you would use +fdc), see manual;
			for (i =0; i<Npoints; ++i)
			  {
			  // ___ Modify kernel, shift spectrum to fDC ___
//C++ TO JAVA CONVERTER TODO TASK: There are no simple equivalents to function pointers in Java:
	//		 final float t = ((pntAxis[kernelnoL])(i,0))*tmp;
			  //kernelL(i,0)  *= complr4(cos(t),-sin(t));// note '-' (see manual)
			  kernelL(i,0) *= complex<Float>(Math.cos(t),-Math.sin(t)); // note '-' (see manual)
			  }
			}

			// ______ For speed: define setdata internally (memcpy) ______
			for (i =0; i<Npoints; i++)
//C++ TO JAVA CONVERTER TODO TASK: The memory management function 'memcpy' has no equivalent in Java:
			  memcpy(PART[i], BUFFER[i+firstL-firstline]+ firstP-slave.currentwindow.pixlo,Npointsxsize);

			// ====== Some speed considerations ======
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//		#if __USE_VECLIB_LIBRARY__
			// ______Compute PART * kernelP______
			cgemv("T", Np, Np, c4alpha, PART[0], Np, kernelP[0], ONEint, c4beta, TMPRES[0], ONEint, 1);
			// ______Compute Result * kernelL; put in matrix RESULT______
			ANS = cdotu(Np, TMPRES[0], ONEint, kernelL[0], ONEint);
			RESULT(linecnt,pixel-overlap.pixlo) = complex<Float>(ANS.re,ANS.im);
	//		#else
			// ______ NO VECLIB: slower, but works ______
			RESULT(linecnt,pixel-overlap.pixlo) = ((matTxmat(PART *kernelP, kernelL))(0,0));
	//		#endif
		  } // for all pixels in overlap
		} // for all lines in overlap


	  // ______ Write last lines of Result to disk (filled upto linecnt) ______
	  DEBUG << "Writing slave: [" << overlap.linehi-linecnt << ":" << overlap.linehi << ", " << overlap.pixlo << ":" << overlap.pixhi << "] (master coord. system)";
	  DEBUG.print();

	  // ______ Actually write ______
	  switch (resampleinput.oformatflag)
		{
		case FORMATCR4:
		  {
		  final complex<Float> zerocr4 = new complex(0,0);
		  for (int thisline =0; thisline<=linecnt; thisline++)
			{
			// ______ Write zero pixels at start ______
			for (int thispixel =0; thispixel<write0pixels1; ++thispixel)
			  {
			  ofile.write((char) zerocr4, sizeofcr4);
			  }
			// ______ WRITE the interpolated data per row ______
			ofile.write((char) RESULT[thisline][0], RESULT.pixels()*sizeofcr4);
			// ______ Write zero pixels at end ______
			for (int thispixel =0; thispixel<write0pixelsN; ++thispixel)
			  {
			  ofile.write((char) zerocr4, sizeofcr4);
			  }
			}
		  break;
		  }
		case FORMATCI2:
		  {
		  final complex<Short> zeroci16 = new complex(0,0);
		  complex<Short> castedresult;
		  for (int thisline =0; thisline<=linecnt; thisline++)
			{
			// ______ Write zero pixels at start ______
			for (int thispixel =0; thispixel<write0pixels1; ++thispixel)
			  {
			  ofile.write((char) zeroci16, sizeofci16);
			  }
			// ______ Write the interpolated data per row ______
			for (int thispixel =0; thispixel<(int)(RESULT.pixels()); thispixel++)
			  {
			  castedresult = cr4toci2(RESULT(thisline,thispixel));
			  ofile.write((char) castedresult, sizeofci16);
			  }
			// ______ Write zero pixels at end ______
			for (int thispixel =0; thispixel<write0pixelsN; ++thispixel)
			  {
			  ofile.write((char) zeroci16, sizeofci16);
			  }
			}
		  break;
		  }
		default:
		  {
			  ERROR.terminate();
			  String cp_s = new String(new char[256]);
			  cp_s = "impossible format";
			  ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			  ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			  ERROR.print();
		  }
		  throw(unhandled_case_error);
		}


	  // ====== Write last zero lines if appropriate (DBOW card) ======
	  switch (resampleinput.oformatflag)
		{
		case FORMATCR4:
		  {
		  complex<Float> zerocr4 = new complex(0,0);
		  for (int thisline =0; thisline<write0linesN; ++thisline)
			for (int thispixel =0; thispixel<(int)(RESULT.pixels())+write0pixels1+write0pixelsN; ++thispixel)
			  ofile.write((char) zerocr4, sizeofcr4);
		  break;
		  }
		case FORMATCI2:
		  {
		  complex<Short> zeroci16 = new complex(0,0);
		  for (int thisline =0; thisline<write0linesN; ++thisline)
			for (int thispixel =0; thispixel<(int)(RESULT.pixels())+write0pixels1+write0pixelsN; ++thispixel)
			  ofile.write((char) zeroci16, sizeofci16);
		  break;
		  }
		default:
		  {
			  ERROR.terminate();
			  String cp_s = new String(new char[256]);
			  cp_s = "impossible format";
			  ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			  ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			  ERROR.print();
		  }
		  throw(unhandled_case_error);
		}
	  ofile.close();



	  // ====== Write results to slave resfile ======
	  String rsmethod = new String(new char[EIGHTY]);
	  switch(resampleinput.method)
		{
		case rs_rect:
		  rsmethod = "nearest neighbour";
		  break;
		case rs_tri:
		  rsmethod = "piecewise linear";
		  break;
		case rs_cc4p:
		  rsmethod = "4 point cubic convolution";
		  break;
		case rs_cc6p:
		  rsmethod = "6 point cubic convolution";
		  break;
		case rs_ts6p:
		  rsmethod = "6 point truncated sinc";
		  break;
		case rs_ts8p:
		  rsmethod = "8 point truncated sinc";
		  break;
		case rs_ts16p:
		  rsmethod = "16 point truncated sinc";
		  break;
		case rs_knab4p:
		  rsmethod = "4 point knab kernel";
		  break;
		case rs_knab6p:
		  rsmethod = "6 point knab kernel";
		  break;
		case rs_knab8p:
		  rsmethod = "8 point knab kernel";
		  break;
		case rs_knab10p:
		  rsmethod = "10 point knab kernel";
		  break;
		case rs_knab16p:
		  rsmethod = "16 point knab kernel";
		  break;
		case rs_rc6p:
		  rsmethod = "6 point raised cosine kernel";
		  break;
		case rs_rc12p:
		  rsmethod = "12 point raised cosine kernel";
		  break;
		default:
		  {
			  ERROR.terminate();
			  String cp_s = new String(new char[256]);
			  cp_s = "impossible.";
			  ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			  ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			  ERROR.print();
		  }
		  throw(unhandled_case_error);
		}

	  String rsoformat = new String(new char[EIGHTY]);
	  switch(resampleinput.oformatflag)
		{
		case FORMATCR4:
		  rsoformat = "complex_real4";
		  break;
		case FORMATCI2:
		  rsoformat = "complex_short";
		  break;
		default:
		  {
			  ERROR.terminate();
			  String cp_s = new String(new char[256]);
			  cp_s = "impossible.";
			  ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			  ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			  ERROR.print();
		  }
		  throw(unhandled_case_error);
		}


	  // --- Write result file ---
	  ofstream scratchlogfile = new ofstream("scratchlogresample", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchlogfile,"resample: scratchlogresample",__FILE__,__LINE__);
	  scratchlogfile << "\n\n*******************************************************************" << "\n* RESAMPLE:" << "\n*******************************************************************" << "\nData_output_file: \t\t\t" << resampleinput.fileout << "\nData_output_format: \t\t\t" << rsoformat << "\nInterpolation kernel: \t\t\t" << rsmethod << "\nResampled slave size in master system: \t" << overlap.linelo - write0lines1 << ", " << overlap.linehi + write0linesN << ", " << overlap.pixlo - write0pixels1 << ", " << overlap.pixhi + write0pixelsN << "\n*******************************************************************\n";
	  scratchlogfile.close();

	  ofstream scratchresfile = new ofstream("scratchresresample", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile,"resample: scratchresresample",__FILE__,__LINE__);
	  scratchresfile << "\n\n*******************************************************************" << "\n*_Start_" << processcontrol[(int)AnonymousEnum.pr_s_resample] << "\n*******************************************************************" << "\nShifted azimuth spectrum:             \t\t" << resampleinput.shiftazi << "\nData_output_file:                     \t\t" << resampleinput.fileout << "\nData_output_format:                   \t\t" << rsoformat << "\nInterpolation kernel:                 \t\t" << rsmethod << "\nFirst_line (w.r.t. original_master):  \t\t" << overlap.linelo - write0lines1 << "\nLast_line (w.r.t. original_master):   \t\t" << overlap.linehi + write0linesN << "\nFirst_pixel (w.r.t. original_master): \t\t" << overlap.pixlo - write0pixels1 << "\nLast_pixel (w.r.t. original_master):  \t\t" << overlap.pixhi + write0pixelsN << "\n*******************************************************************" << "\n* End_" << processcontrol[(int)AnonymousEnum.pr_s_resample] << "_NORMAL" << "\n*******************************************************************\n";
	  scratchresfile.close();


	  // ______Tidy up______
	  DEBUG.print("deleting new matrix, memory errors could be caused by this");
	  for (i =0;i<Ninterval;i++) // like this ???
		{
		pntKernelAz[i] = null;
		pntKernelRg[i] = null;
		pntAxis[i] = null;
		//    delete [] pntKernelAz[i];
		}
	  DEBUG.print("Exiting resample.");
	  } // END resample

//***************************************************************
// * ms_timing_error                                              *
// *                                                              *
// * relative timing error between master and slave               *
// *                                                              *
// * input:                                                       *
// *  - master                                                    *
// *  - interferogram result file                                 *
// *  - timing input                                              *
// *  - coarse_orbit_offsetL                                      *
// *  - coarse_orbit_offsetP                                      *
// * output:                                                      *
// *  - ms_az_timing_error_L                                      *
// *  - ms_r_timing_error_P                                       *
// *  - ms_az_timing_error                                        *
// *  - ms_r_timing_error                                         *
// *                                                              *
// *    Freek van Leijen, 06-SEP-2007                             *
// ***************************************************************

	// ______ Compute master-slave timing error ______
	public static void ms_timing_error(slcimage master, String i_resfile, input_reltiming timinginput, RefObject<Integer> coarse_orbit_offsetL, RefObject<Integer> coarse_orbit_offsetP)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"ms_timing_error (FvL 6-SEP-2007)"<<ends;
		  TRACE.print();
	  }

		INFO << coarse_orbit_offsetL.argvalue;
		INFO.print();
		INFO << coarse_orbit_offsetP.argvalue;
		INFO.print();

	  final float THRESHOLD = timinginput.threshold; // threshold ...
	  final int MAX_ITERATIONS = timinginput.maxiter; // max. of pnts to remove
	  final float CRIT_VALUE = timinginput.k_alpha; // crit. value outlier removal
	  final int DEGREE = 0; // degree of polynomial
	  final int Nunk = Ncoeffs(DEGREE); // Number of unknowns/direction

	  // ______ Normalize data for polynomial ______
	  final double minL = master.originalwindow.linelo;
	  final double maxL = master.originalwindow.linehi;
	  final double minP = master.originalwindow.pixlo;
	  final double maxP = master.originalwindow.pixhi;

	  // ______ A priori sigma of  offset ______
	  // ______ Read this factor from the result file 
	  // ______ "Oversampling factor: 32"
	  // ______ "Window_size_L_for_correlation: 4"
	  // ______ "Window_size_P_for_correlation: 121"
	  DEBUG.print("Reading oversampling factor from result file");
	  int osfactor = 32; // oversamplingsfactor
	  int corrwinL = 64; // window size to compute FINE correlation
	  int corrwinP = 64; // window size to compute FINE correlation
	  String c4osfactor = new String(new char[4]);
	  String c10corrwinL = new String(new char[10]);
	  String c10corrwinP = new String(new char[10]);
	  boolean found = readres(c4osfactor,sizeof(c4osfactor),i_resfile, "Oversampling", 1);
	  if (found)
		  osfactor = (int)(Integer.parseInt(c4osfactor));
	  found = readres(c10corrwinL,sizeof(c10corrwinL),i_resfile, "Window_size_L_for_correlation:", 0);
	  if (found)
		  corrwinL = (int)(Integer.parseInt(c10corrwinL));
	  found = readres(c10corrwinP,sizeof(c10corrwinP),i_resfile, "Window_size_P_for_correlation:", 0);
	  if (found)
		  corrwinP = (int)(Integer.parseInt(c10corrwinP));
	  corrwinL = max(10,corrwinL-8); // if fft method peak is not at center
	  corrwinP = max(10,corrwinP-8); // +then effective number of samples is smaller
	  // _____ oversampling factor is bin in which maximum can be found _____
	  // _____ ovsf=16-->apriorisigma=0.03
	  final float ACCURACY = 0.5 * (1.0/((float)osfactor));

	  // but we need coreg accuracy of 0.1 pixel about.  therefore use a priori
	  // based on experience here, and different for azimuth and range
	  // this also helps our automated outlier detection and testing hopefully.
	  // BK 15-Apr-2003
	  // if the image is oversampled, then still use orig spacing
	  float SIGMAL = 0.15/master.ovs_az; // sigma in orig pixels
	  float SIGMAP = 0.10/master.ovs_rg; // seems range direction is better???
	  INFO.print("Using a smaller sigma in range, because it seems that can be estimated better");
	  INFO << "a priori std.dev offset vectors line direction [samples]:  " << SIGMAL;
	  INFO.print();
	  INFO << "a priori std.dev offset vectors pixel direction [samples]: " << SIGMAP;
	  INFO.print();

	  // ______ Find #points > threshold ______
	  matrix<Float> Data = getofffile(i_resfile, THRESHOLD);
	  // ______ Data contains the following: ______
	  // Data(i,0) = winnumber; Data(i,1) = posL; Data(i,2) = posP; 
	  // Data(i,3) = offL;      Data(i,4) = offP; Data(i,5) = corr;


	  int ITERATION = 0;
	  int DONE = 0;
	  // sqr: level significance: alpha=0.001; power of test: gamma=0.80
	  //real4 CRIT_VALUE = sqrt(3.29);
	  INFO << "Critical value for outlier test: " << CRIT_VALUE;
	  INFO.print();
	  int winL = 0; // window number to be removed
	  int winP = 0; // window number of largest w -test in range
	  matrix<Double> eL_hat;
	  matrix<Double> eP_hat;
	  matrix<Double> wtestL;
	  matrix<Double> wtestP;
	  matrix<Double> rhsL;
	  matrix<Double> rhsP;
	  matrix<Double> Qx_hat;
	  double maxdev = 0.0;
	  double overallmodeltestL = 0.0;
	  double overallmodeltestP = 0.0;
	  double maxwL;
	  double maxwP;
//C++ TO JAVA CONVERTER NOTE: 'register' variable declarations are not supported in Java:
//ORIGINAL LINE: register int i,j,k,index;
	  int i;
	  int j;
	  int k;
	  int index;
	  while (DONE != 1)
		{
		DEBUG << "Start iteration " << ITERATION;
		DEBUG.print();
		// ______ Remove identified outlier from previous estimation ______
		if (ITERATION != 0)
		  {
		  matrix<Float> tmp_DATA = Data; //(remove_observation_i,*);
		  Data.resize(Data.lines()-1, Data.pixels());
		  j = 0; // counter over reduced obs.vector
		  for (i =0; i<tmp_DATA.lines(); i++) // counter over original window numbers
			{
			if (i != winL) // do not copy the one to be removed.
			  {
			  Data.setrow(j,tmp_DATA.getrow(i)); // copy back without removed obs.
			  j++; // fill next row of Data
			  }
			else
			  {
			  DEBUG << "Removing observation " << i << " from observation vector.";
			  DEBUG.print();
			  }
			}
		  }

		// ______Check redundancy______
		int Nobs = Data.lines(); // Number of points > threshold
		if (Nobs < Nunk)
		  {
		  {
			  ERROR.terminate();
			  String cp_s = new String(new char[256]);
			  cp_s = "ms_timing_error: Number of windows > threshold is smaller than parameters solved for.";
			  ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			  ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			  ERROR.print();
		  }
		  throw(input_error);
		  }

		// ______Set up system of equations______
		// ______Order unknowns: A00 A10 A01 A20 A11 A02 A30 A21 A12 A03 for degree=3______
		matrix<Double> yL = new matrix(Nobs,1); // observation
		matrix<Double> yP = new matrix(Nobs,1); // observation
		matrix<Double> A = new matrix(Nobs,Nunk); // designmatrix
		matrix<Double> Qy_1 = new matrix(Nobs,1); // a priori covariance matrix (diag)

		// ______ Normalize data for polynomial ______
		DEBUG << "ms_timing_error: polynomial normalized by factors: " << minL << " " << maxL << " " << minP << " " << maxP << " to [-2,2]";
		DEBUG.print();

		// ______Fill matrices______
		DEBUG.print("Setting up design matrix for LS adjustment");
		for (i =0; i<Nobs; i++)
		  {
		  double posL = normalize((double)(Data(i,1)),minL,maxL);
		  double posP = normalize((double)(Data(i,2)),minP,maxP);
		  yL(i,0) = (double)(Data(i,3));
		  yP(i,0) = (double)(Data(i,4));
		  DEBUG << "ms_timing_error: (" << posL << ", "<< posP << "): yL=" << yL(i,0) << " yP=" << yP(i,0);
		  DEBUG.print();
		  // ______Set up designmatrix______
		  index = 0;
		  for (j =0; j<=DEGREE; j++)
			{
			for (k =0; k<=j; k++)
			  {
			  A(i,index) = Math.pow(posL,(double)(j-k)) * Math.pow(posP,(double)k);
			  index++;
			  }
			}
		  }


		// ______Weight matrix data______
		DEBUG.print("Setting up (inverse of) covariance matrix for LS adjustment");
		for (i =0; i<Nobs; i++)
		  Qy_1(i,0) = (double)1.0; //unweighted, could be changed later


		// ______Compute Normalmatrix, rghthandside______
		matrix<Double> N = matTxmat(A,diagxmat(Qy_1,A));
		rhsL = matTxmat(A,diagxmat(Qy_1,yL));
		rhsP = matTxmat(A,diagxmat(Qy_1,yP));
		Qx_hat = N;
		// ______Compute solution______
		choles(Qx_hat); // Cholesky factorisation normalmatrix
		solvechol(Qx_hat,rhsL); // Solution unknowns in rhs
		solvechol(Qx_hat,rhsP); // Solution unknowns in rhs
		invertchol(Qx_hat); // Covariance matrix of unknowns
		// ______Test inverse______
		for (i =0; i<Qx_hat.lines(); i++)
		  for (j =0; j<i; j++)
			Qx_hat(j,i) = Qx_hat(i,j); // repair Qx
		maxdev = max(Math.abs(N *Qx_hat-eye((double)(Qx_hat.lines()))));
		DEBUG << "ms_timing_error: max(abs(N*inv(N)-I)) = " << maxdev;
		DEBUG.print();
		// ___ use trace buffer to store string, remember to rewind it ___
		if (maxdev > .01)
		  {
		  ERROR << "ms_timing_error: maximum deviation N*inv(N) from unity = " << maxdev << ". This is larger than 0.01";
		  ERROR.print(ERROR.get_str());
		  throw(some_error);
		  }
		else if (maxdev > .001)
		  {
		  WARNING << "ms_timing_error: maximum deviation N*inv(N) from unity = " << maxdev << ". This is between 0.01 and 0.001";
		  WARNING.print();
		  }


		// ______Some other stuff, scale is ok______
		matrix<Double> Qy_hat = A * (matxmatT(Qx_hat,A));
		matrix<Double> yL_hat = A * rhsL;
		matrix<Double> yP_hat = A * rhsP;
		eL_hat = yL - yL_hat;
		eP_hat = yP - yP_hat;
		matrix<Double> Qe_hat = -Qy_hat;
		for (i =0; i<Nobs; i++)
		  Qe_hat(i,i) += (1. / Qy_1(i,0));

		// ______Overall model test (variance factor)______
		overallmodeltestL = 0.;
		overallmodeltestP = 0.;
		for (i =0; i<Nobs; i++)
		  {
		  overallmodeltestL += sqr(eL_hat(i,0))*Qy_1(i,0);
		  overallmodeltestP += sqr(eP_hat(i,0))*Qy_1(i,0);
		  }
		overallmodeltestL = (overallmodeltestL/sqr(SIGMAL)) /(Nobs-Nunk); // this is sigma hat!
		overallmodeltestP = (overallmodeltestP/sqr(SIGMAP)) /(Nobs-Nunk); // not OMT!
		DEBUG << "ms_timing_error: overallmodeltest Lines = " << overallmodeltestL;
		DEBUG.print();
		DEBUG << "ms_timing_error: overallmodeltest Pixels = " << overallmodeltestP;
		DEBUG.print();

		// ______Datasnooping, assume Qy diag______
		wtestL.resize(Nobs,1);
		wtestP.resize(Nobs,1);
		for (i =0; i<Nobs; i++)
		  {
		  wtestL(i,0) = eL_hat(i,0) / (Math.sqrt(Qe_hat(i,i))*SIGMAL); // computed excl.var.factor
		  wtestP(i,0) = eP_hat(i,0) / (Math.sqrt(Qe_hat(i,i))*SIGMAP);
		  }

		int dumm = 0;
		maxwL = max(Math.abs(wtestL),winL,dumm); // returns winL
		maxwP = max(Math.abs(wtestP),winP,dumm); // returns winP
		DEBUG << "maximum wtest statistic azimuth = " << maxwL << " for window number: " << Data(winL,0);
		DEBUG.print();
		DEBUG << "maximum wtest statistic range   = " << maxwP << " for window number: " << Data(winP,0);
		DEBUG.print();
		// --- use summed wtest for outlier detection ---
		// #%// BK 21-Oct-2003
		matrix<Double> wtestsum = sqr(wtestL)+sqr(wtestP); // (Nobs,1)
		double maxwsum = max(wtestsum,winL,dumm); // idx to remove
		DEBUG << "Detected outlier:  summed sqr.wtest = " << maxwsum << "; observation: " << winL << "; window number: " << Data(winL,0);
		DEBUG.print();


		// ______ Test if we are done yet ______
		if (Nobs <= Nunk)
		  {
		  WARNING.print("NO redundancy!  Exiting iterations.");
		  DONE = 1; // cannot remove more than this
		  }
		// seems something fishy here..., b-method of testing delft
		//    if (max(overallmodeltestL,overallmodeltestP) < 1.0)
		//      {
		//      INFO.print("OMTs accepted, not iterating anymore (final solution reached).");
		//      DONE = 1;// ok (?).
		//      }
		if (max(maxwL,maxwP) <= CRIT_VALUE) // all tests accepted?
		  {
		  INFO.print("All outlier tests accepted! (final solution computed)");
		  DONE = 1; // yeah!
		  }
		if (ITERATION >= MAX_ITERATIONS)
		  {
		  INFO.print("max. number of iterations reached (exiting loop).");
		  DONE = 1; // we reached max. (or no max_iter specified)
		  }

		// ______ Only warn if last iteration has been done ______
		if (DONE == 1)
		  {
		  // ___ use trace buffer to store string, remember to rewind it ___
		  if (overallmodeltestL > 10)
			{
			WARNING << "ms_timing_error: overallmodeltest Lines = " << overallmodeltestL << ends;
			WARNING.print();
			WARNING << " is larger than 10. (Suggest model or a priori sigma not correct.)";
			WARNING.print();
			}
		  // ___ use trace buffer to store string, remember to rewind it ___
		  if (overallmodeltestP > 10)
			{
			WARNING << "ms_timing_error: overallmodeltest Pixels = " << overallmodeltestP;
			WARNING.print();
			WARNING << " is larger than 10.\n(suggests a priori sigma not correct.)";
			WARNING.print();
			}

		  } // Only warn when done iterating.
		ITERATION++; // update counter here!
		} // iterations remove outliers


	  // Calculate master-slave timing errors
	  int ms_az_timing_error_L = coarse_orbit_offsetL.argvalue-(int)(rint(rhsL(0,0)));
	  int ms_r_timing_error_P = coarse_orbit_offsetP.argvalue-(int)(rint(rhsP(0,0)));

	  double ms_az_timing_error = (double)ms_az_timing_error_L/master.prf;
	  double ms_r_timing_error = (double)ms_r_timing_error_P/master.rsr2x;

	  INFO << "Orbit azimuth offset (master-slave): " << coarse_orbit_offsetL.argvalue << " lines.";
	  INFO.print();
	  INFO << "Orbit range offset (master-slave): " << coarse_orbit_offsetP.argvalue << " pixels.";
	  INFO.print();
	  INFO << "Estimated azimuth offset (master-slave): " << rhsL(0,0) << " lines.";
	  INFO.print();
	  INFO << "Estimated range offset (master-slave): " << rhsP(0,0) << " pixels.";
	  INFO.print();

	  INFO << "Estimated azimuth timing error (master-slave): " << ms_az_timing_error_L << " lines.";
	  INFO.print();
	  INFO << "Estimated range timing error (master-slave): " << ms_r_timing_error_P << " pixels.";
	  INFO.print();

	  INFO << "Estimated azimuth timing error (master-slave) [sec]: " << ms_az_timing_error << " sec.";
	  INFO.print();
	  INFO << "Estimated range timing error (master-slave) [sec]: " << ms_r_timing_error << " sec.";
	  INFO.print();


	  // ______ Write to tmp files ______
	  ofstream scratchlogfile = new ofstream("scratchlogtiming", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchlogfile,"timing: scratchlogtiming",__FILE__,__LINE__);
	  scratchlogfile << "\n\n*******************************************************************" << "\n* RELATIVE_TIMING_ERROR" << "\n*******************************************************************" << "\nOrbit_azimuth_offset (master-slave):    \t" << coarse_orbit_offsetL.argvalue << " lines." << "\nOrbit_range_offset (master-slave):      \t" << coarse_orbit_offsetP.argvalue << " pixels." << "\nEstimated_azimuth_offset (master-slave): " << rhsL(0,0) << " lines." << "\nEstimated_range_offset (master-slave): " << rhsP(0,0) << " pixels." << "\nEstimated_azimuth_timing_error_lines (master-slave): " << ms_az_timing_error_L << " lines." << "\nEstimated_range_timing_error_pixels (master-slave): " << ms_r_timing_error_P << " pixels." << "\nEstimated_azimuth_timing_error_sec (master-slave): " << ms_az_timing_error << " sec." << "\nEstimated_range_timing_error_sec (master-slave): " << ms_r_timing_error << " sec." << "\n*******************************************************************\n";
	  scratchlogfile.close();

	  ofstream scratchresfile = new ofstream("scratchrestiming", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile,"timing: scratchrestiming",__FILE__,__LINE__);
	  scratchresfile.setf(ios.right, ios.adjustfield);
	  scratchresfile << "\n\n*******************************************************************" << "\n*_Start_" << processcontrol[(int)AnonymousEnum.pr_i_timing] << "\n*******************************************************************" << "\nOrbit_azimuth_offset (master-slave):    \t" << coarse_orbit_offsetL.argvalue << " lines." << "\nOrbit_range_offset (master-slave):      \t" << coarse_orbit_offsetP.argvalue << " pixels." << "\nEstimated_azimuth_offset (master-slave): " << rhsL(0,0) << " lines." << "\nEstimated_range_offset (master-slave): " << rhsP(0,0) << " pixels." << "\nEstimated_azimuth_timing_error_lines (master-slave): " << ms_az_timing_error_L << " lines." << "\nEstimated_range_timing_error_pixels (master-slave): " << ms_r_timing_error_P << " pixels." << "\nEstimated_azimuth_timing_error_sec (master-slave): " << ms_az_timing_error << " sec." << "\nEstimated_range_timing_error_sec (master-slave): " << ms_r_timing_error << " sec." << "\n*******************************************************************" << "\n* End_" << processcontrol[(int)AnonymousEnum.pr_i_timing] << "_NORMAL" << "\n*******************************************************************\n";


	  // ====== Tidy up ======
	  scratchresfile.close();
	  PROGRESS.print("Finished computation of master-slave timing error.");

	  } // END rel_timing_error

//***************************************************************
// *    cc4                                                       *
// *                                                              *
// * cubic convolution 4 points                                   *
// *                                                              *
// * input:                                                       *
// *  - x-axis                                                    *
// * output:                                                      *
// *  - y=f(x); function evaluated at x                           *
// *                                                              *
// *    Bert Kampes, 16-Mar-1999                                  *
// ***************************************************************

	// ______Interpolation kernals______
	public static matrix<Float> cc4(matrix<Float> x)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"cc4 (BK 16-Mar-1999)"<<ends;
		  TRACE.print();
	  }
	  if (x.pixels() != 1)
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "cc4: standing vectors only.";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}

	  float alpha = -1.0;
	  matrix<Float> y = new matrix(x.lines(),1);
	  for (register int i =0;i<y.lines();i++)
		{
		float xx2 = sqr(x(i,0));
		float xx = Math.sqrt(xx2);
		if (xx < 1)
		  y(i,0) = (alpha+2)*xx2 *xx - (alpha+3)*xx2 + 1;
		else if (xx < 2)
		  y(i,0) = alpha *xx2 *xx - 5 *alpha *xx2 + 8 *alpha *xx - 4 *alpha;
		else
		  y(i,0) = 0.0;
		}
	  return y;
	  } // END cc4

//***************************************************************
// *    cc6                                                       *
// *                                                              *
// * cubic convolution 6 points                                   *
// *                                                              *
// * input:                                                       *
// *  - x-axis                                                    *
// * output:                                                      *
// *  - y=f(x); function evaluated at x                           *
// *                                                              *
// *    Bert Kampes, 16-Mar-1999                                  *
// * corrected (alfa+beta)->(alfa-beta) after correction in paper *
// * by Ramon Hanssen                                             *
// *    Bert Kampes, 16-Mar-1999                                  *
// ***************************************************************
	public static matrix<Float> cc6(matrix<Float> x)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"cc6 (BK 16-Mar-1999)"<<ends;
		  TRACE.print();
	  }
	  if (x.pixels() != 1)
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "cc6: standing vectors only.";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}

	  float alpha = -.5;
	  float beta = .5;
	  matrix<Float> y = new matrix(x.lines(),1);
	  for (register int i =0;i<y.lines();i++)
		{
		float xx2 = sqr(x(i,0));
		float xx = Math.sqrt(xx2);
		if (xx < 1)
		  y(i,0) = (alpha-beta+2)*xx2 *xx - (alpha-beta+3)*xx2 + 1;
		//y(i,0) = (alpha+beta+2)*xx2*xx - (alpha+beta+3)*xx2 + 1;??wrong in paper?
		else if (xx < 2)
		  y(i,0) = alpha *xx2 *xx - (5 *alpha-beta)*xx2 + (8 *alpha-3 *beta)*xx - (4 *alpha-2 *beta);
		else if (xx < 3)
		  y(i,0) = beta *xx2 *xx - 8 *beta *xx2 + 21 *beta *xx - 18 *beta;
		else
		  y(i,0) = 0.;
		}
	  return y;
	  } // END cc6

//***************************************************************
// *    ts6                                                       *
// *                                                              *
// * truncated sinc 6 points                                      *
// *                                                              *
// * input:                                                       *
// *  - x-axis                                                    *
// * output:                                                      *
// *  - y=f(x); function evaluated at x                           *
// *                                                              *
// *    Bert Kampes, 16-Mar-1999                                  *
// ***************************************************************
	public static matrix<Float> ts6(matrix<Float> x)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"ts6 (BK 16-Mar-1999)"<<ends;
		  TRACE.print();
	  }
	  if (x.pixels() != 1)
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "ts6: standing vectors only.";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}

	  matrix<Float> y = new matrix(x.lines(),1);
	  for (register int i =0;i<y.lines();i++)
		y(i,0) = sinc(x(i,0)) * rect(x(i,0)/6.0);
	  return y;
	  } // END ts6

//***************************************************************
// *    ts8                                                       *
// *                                                              *
// * truncated sinc 8 points                                      *
// *                                                              *
// * input:                                                       *
// *  - x-axis                                                    *
// * output:                                                      *
// *  - y=f(x); function evaluated at x                           *
// *                                                              *
// *    Bert Kampes, 16-Mar-1999                                  *
// ***************************************************************
	public static matrix<Float> ts8(matrix<Float> x)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"ts8 (BK 16-Mar-1999)"<<ends;
		  TRACE.print();
	  }
	  if (x.pixels() != 1)
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "ts8: standing vectors only.";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}

	  matrix<Float> y = new matrix(x.lines(),1);
	  for (register int i =0;i<y.lines();i++)
		y(i,0) = sinc(x(i,0)) * rect(x(i,0)/8.0);
	  return y;
	  } // END ts8

//***************************************************************
// *    ts16                                                      *
// *                                                              *
// * truncated sinc 16 points                                     *
// *                                                              *
// * input:                                                       *
// *  - x-axis                                                    *
// * output:                                                      *
// *  - y=f(x); function evaluated at x                           *
// *                                                              *
// *    Bert Kampes, 16-Mar-1999                                  *
// ***************************************************************
	public static matrix<Float> ts16(matrix<Float> x)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"ts16 (BK 16-Mar-1999)"<<ends;
		  TRACE.print();
	  }
	  if (x.pixels() != 1)
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "ts16: standing vectors only.";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}
	  matrix<Float> y = new matrix(x.lines(),1);
	  for (register int i =0;i<y.lines();i++)
		y(i,0) = sinc(x(i,0)) * rect(x(i,0)/16.0);
	  return y;
	  } // END ts16

//***************************************************************
// *    rect                                                      *
// *                                                              *
// * rect function for matrix (stepping function?)                *
// *                                                              *
// * input:                                                       *
// *  - x-axis                                                    *
// * output:                                                      *
// *  - y=f(x); function evaluated at x                           *
// *                                                              *
// *    Bert Kampes, 16-Mar-1999                                  *
// ***************************************************************
	public static matrix<Float> rect(matrix<Float> x)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"rect (BK 16-Mar-1999)"<<ends;
		  TRACE.print();
	  }
	  if (x.pixels() != 1)
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "rect: standing vectors only.";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}

	  matrix<Float> y = new matrix(x.lines(),1);
	  for (register int i =0;i<y.lines();i++)
		y(i,0) = rect(x(i,0));
	  return y;
	  } // END rect

//***************************************************************
// *    tri                                                       *
// *                                                              *
// * tri function for matrix (piecewize linear?, triangle)        *
// *                                                              *
// * input:                                                       *
// *  - x-axis                                                    *
// * output:                                                      *
// *  - y=f(x); function evaluated at x                           *
// *                                                              *
// *    Bert Kampes, 16-Mar-1999                                  *
// ***************************************************************
	public static matrix<Float> tri(matrix<Float> x)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"tri (BK 16-Mar-1999)"<<ends;
		  TRACE.print();
	  }
	  if (x.pixels() != 1)
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "tri: standing vectors only.";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}
	  matrix<Float> y = new matrix(x.lines(),1);
	  for (register int i =0;i<y.lines();i++)
		y(i,0) = tri(x(i,0));
	  return y;
	  } // END tri

//***************************************************************
// *    knab                                                      *
// *                                                              *
// * KNAB window of N points, oversampling factor CHI             *
// *                                                              *
// * defined by: Migliaccio IEEE letters vol41,no5, pp1105,1110, 2003 *
// * k = sinc(x).*(cosh((pi*v*L/2)*sqrt(1-(2.*x./L).^2))/cosh(pi*v*L/2));
// *                                                              *
// * input:                                                       *
// *  - x-axis                                                    *
// *  - oversampling factor of bandlimited sigal CHI              *
// *  - N points of kernel size                                   *
// * output:                                                      *
// *  - y=f(x); function evaluated at x                           *
// *                                                              *
// *    Bert Kampes, 22-DEC-2003                                  *
// ***************************************************************
	// ___ knab: oversampling factor of signal CHI, number of points N ___
	public static matrix<Float> knab(matrix<Float> x, float CHI, int N)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"knab (BK 22-Dec-2003)"<<ends;
		  TRACE.print();
	  }
	  if (x.pixels() != 1)
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "knab: standing vectors only.";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}
	  matrix<Float> y = new matrix(x.lines(),1);
	  float v = 1.0-1.0/CHI;
	  float vv = PI *v *(float)N/2.0;
	  for (register int i =0;i<y.lines();i++)
		y(i,0) = sinc(x(i,0))*Math.cosh(vv *Math.sqrt(1.0-sqr(2.0 *x(i,0)/(float)N)))/Math.cosh(vv);
	  return y;
	  } // END knab

//***************************************************************
// *    rc_kernel                                                 *
// *                                                              *
// * Raised Cosine window of N points, oversampling factor CHI    *
// *                                                              *
// * defined by: Cho, Kong and Kim, J.Elektromagn.Waves and appl  *
// *  vol19, no.1, pp, 129-135, 2005;                             *
// * claimed to be best, 0.9999 for 6 points kernel.              *
// * k(x) = sinc(x).*[cos(v*pi*x)/(1-4*v^2*x^2)]*rect(x/L)     *
// *  where v = 1-B/fs = 1-1/Chi (roll-off factor; ERS: 15.55/18.96)*
// *        L = 6 (window size)                                   *
// *                                                              *
// * input:                                                       *
// *  - x-axis                                                    *
// *  - oversampling factor of bandlimited sigal CHI              *
// *  - N points of kernel size                                   *
// * output:                                                      *
// *  - y=f(x); function evaluated at x                           *
// *                                                              *
// #%// Bert Kampes, 28-Jul-2005
// ***************************************************************
	// ___ raised cosine: oversampling factor of signal CHI, number of points N ___
	public static matrix<Float> rc_kernel(matrix<Float> x, float CHI, int N)
	  {
	  {
		  TRACE.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<"rc_kernel (BK 28-Jul-2005)"<<ends;
		  TRACE.print();
	  }
	  if (x.pixels() != 1)
		{
		{
			ERROR.terminate();
			String cp_s = new String(new char[256]);
			cp_s = "rc_kernel: standing vectors only.";
			ERROR.reset();
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends;
			ERROR.print();
		}
		throw(input_error);
		}
	  matrix<Float> y = new matrix(x.lines(),1);
	  float v = 1.0-1.0/CHI; // alpha in paper cho05
	  for (register int i =0;i<y.lines();i++)
		y(i,0) = sinc(x(i,0)) * rect(x(i,0)/(float)N)* Math.cos(v *PI *x(i,0)) / (1.0-sqr(2.0 *v *x(i,0)));
	  return y;
	  } // END rc_kernel

	//#endif // COREGISTRATION_H





	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if WIN32
	// Jia defined min max here, but I did this in constants.hh
	//#define max _MAX
	//#define min _MIN
	//#endif
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
// * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/coregistration.cc,v $     *
// * $Revision: 3.40 $                                            *
// * $Date: 2005/10/18 13:46:51 $                                 *
// * $Author: kampes $                                            *
// *                                                              *
// * -coarse based on orbits.                                     *
// * -coarse based on correlation with fft/in space domain.       *
// * -fine coregistration offset vector computation.              *
// * -computation coregistration parameters.                      *
// * -computation flat earth correction.                          *
// * -resampling of slave to master grid.                         *
// ***************************************************************


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
// * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/coregistration.hh,v $
// * $Revision: 3.13 $
// * $Date: 2005/08/24 10:03:18 $
// * $Author: kampes $
// *
// * Routines for coregistration.
// ***************************************************************


//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
//#if ! COREGISTRATION_H
//#define COREGISTRATION_H


// Jia defined this for compilation under windows
// Bert Kampes, 24-Aug-2005
//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
//#if _MSC_VER > 1000
//#endif // _MSC_VER > 1000




// ______ Prototypes ______
//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
//#if __USE_VECLIB_LIBRARY__
// ______This typedef is needed for output of veclib routine cdot?______
public class STUPID_cr4
{
	public float re;
	public float im;
}

//C++ TO JAVA CONVERTER TODO TASK: Extern blocks are not supported in Java.
extern "C"
{

//C++ TO JAVA CONVERTER TODO TASK: Extern blocks are not supported in Java.
extern "C"
{



final class DefineConstantsCoregistration
{
	public static final String SWNAME = "Doris (Delft o-o Radar Interferometric Software)";
	public static final String SWVERSION = "Version  4.02 (08-Jun-2009)";
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
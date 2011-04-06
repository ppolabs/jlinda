public class GlobalMembersReferencephase
{

	// Using A Two-Dimensional Quality Mesh Generator and Delaunay Triangulator
	// from  Jonathan Richard Shewchuk
	// Some definition for triangulate call
	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//#define VOID int
	//C++ TO JAVA CONVERTER NOTE: The following #define macro was replaced in-line:
	//#define REAL double
	//#define ANSI_DECLARATORS



	//***************************************************************
	// *    flatearth                                                 *
	// *                                                              *
	// * Compute polynomial model for 'flat earth' correction.        *
	// *  fie(l,p) = sumj=0:d sumk=0:d Ajk l^j p^k (NOT bert 8sept99) *
	// * precise orbits are used to compute delta range for Npoints   *
	// * after which the polynomial model is fitted (LS).             *
	// *                                                              *
	// * input:                                                       *
	// *  - inputoptions                                              *
	// *  - info structs                                              *
	// *  - platform data points                                      *
	// * output:                                                      *
	// *  - void (result to file "scratchresflat")                    *
	// *  - coefficients normalized wrt. original window of master    *
	// *                                                              *
	// *    Bert Kampes, 09-Mar-1999                                  *
	// *    Bert Kampes, 26-Oct-1999 normalization of coeff.,         *
	// *    dump to logfile: var(unknowns) == diag(inv(AtA))          *
	// ***************************************************************
	public static void flatearth(input_comprefpha comprefphainput, input_ell ellips, slcimage master, slcimage slave, productinfo interferogram, RefObject<orbit> masterorbit, RefObject<orbit> slaveorbit)
	  {
	  TRACE_FUNCTION("flatearth (BK 26-Oct-1999)")
	  final int32 MAXITER = 10;
	  final real8 CRITERPOS = 1e-6;
	  final real8 CRITERTIM = 1e-10;
	  String dummyline = new String(new char[ONE27]);

	  INFO << "FLATEARTH: MAXITER: " << MAXITER << "; " << "CRITERPOS: " << CRITERPOS << " m; " << "CRITERTIM: " << CRITERTIM << " s";
	  INFO.print();

	  // ______ Normalization factors for polynomial ______
	  final real8 minL = master.originalwindow.linelo;
	  final real8 maxL = master.originalwindow.linehi;
	  final real8 minP = master.originalwindow.pixlo;
	  final real8 maxP = master.originalwindow.pixhi;
	  INFO << "flatearth: polynomial normalized by factors: " << minL << " " << maxL << " " << minP << " " << maxP << " to [-2,2]";
	  INFO.print();

	  // ______Handling of input______
	  final real8 m_minpi4cdivlam = (-4.0 *PI *SOL)/master.wavelength;
	  final real8 s_minpi4cdivlam = (-4.0 *PI *SOL)/slave.wavelength;
	  DEBUG << "master wavelength = " << master.wavelength;
	  DEBUG.print();
	  DEBUG << "slave  wavelength = " << slave.wavelength;
	  DEBUG.print();
	  final int32 DEGREE = comprefphainput.degree;
	  final int32 Nunk = Ncoeffs(DEGREE); // Number of unknowns
	  boolean pointsrandom = true;
	  if (specified(comprefphainput.ifpositions))
		pointsrandom = false; // only use those points



	  // ______ Distribute points wel distributed over win ______
	  // ______ or read from ascii file ______
	  // ______(i,0): line, (i,1): pixel, (i,2) flagfromdisk______
	  matrix<uint> Position;
	  final uint Npoints = comprefphainput.Npoints;
//C++ TO JAVA CONVERTER NOTE: 'register' variable declarations are not supported in Java:
//ORIGINAL LINE: register int32 i,j,k,index;
	  int32 i;
	  int32 j;
	  int32 k;
	  int32 index;

	  if (pointsrandom) // no filename specified
		{
		Position = distributepoints(Npoints, interferogram.win);
		}
	  else // read from file
		{
		Position.resize(Npoints, 3);
		//ifstream ifpos(comprefphainput.ifpositions, ios::in);
		ifstream ifpos;
		RefObject<ifstream> TempRefObject = new RefObject<ifstream>(ifpos);
		openfstream(TempRefObject, comprefphainput.ifpositions);
		ifpos = TempRefObject.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		bk_assert(ifpos, comprefphainput.ifpositions, __FILE__, __LINE__);
		uint ll;
		uint pp;
		for (i =0; i<Npoints; ++i)
		  {
		  ifpos >> ll >> pp;
		  Position(i,0) = uint(ll);
		  Position(i,1) = uint(pp);
		  Position(i,2) = uint(1); // flag from file
		  ifpos.getline(dummyline,ONE27,'\n'); // goto next line.
		  }
		ifpos.close();

		// ______ Check last point ivm. EOL after last position in file ______
		if (Position(Npoints-1,0) == Position(Npoints-2,0) && Position(Npoints-1,1) == Position(Npoints-2,1))
		  {
		  Position(Npoints-1,0) = uint(.5*(minL + maxL) + 27); // random
		  Position(Npoints-1,1) = uint(.5*(minP + maxP) + 37); // random
		  WARNING << "refpha: there should be no EOL after last point in file: " << comprefphainput.ifpositions;
		  WARNING.print();
		  }

		// ______ Check if points are in overlap ______
		// ______ no check for uniqueness of points ______
		}

	  matrix<real8> y = new matrix(Npoints, 1); // observation
	  matrix<real8> y_h2ph = new matrix(Npoints, 1); // observation, h2ph factors, added by FvL
	  matrix<real8> A = new matrix(Npoints, Nunk); // designmatrix

	  // ______Check redundancy______
	  if (Npoints < Nunk)
		{
		PRINT_ERROR("flatearth: Number of points is smaller than parameters solved for.")
		throw(input_error);
		}



	  // ======Compute delta r for all points======
	  for (i =0; i<Npoints; ++i)
		{
		final real8 line = Position(i,0);
		final real8 pixel = Position(i,1);

		// ______ Compute azimuth/range time of this pixel______
		//const real8 m_trange = pix2tr(pixel,master.t_range1,master.rsr2x);
		final real8 m_trange = master.pix2tr(pixel);
		final real8 m_tazi = master.line2ta(line); // added by FvL

		// ______ Compute xyz of this point P from position in image ______
		cn P; // point, returned by lp2xyz
		RefObject<cn> TempRefObject2 = new RefObject<cn>(P);
		lp2xyz(line, pixel, ellips, master, masterorbit, TempRefObject2, MAXITER, CRITERPOS);
		P = TempRefObject2.argvalue;

		// ______ Compute xyz for slave satelite from P ______
		real8 s_tazi; // returned
		real8 s_trange; // returned
		RefObject<real8> TempRefObject3 = new RefObject<real8>(s_tazi);
		RefObject<real8> TempRefObject4 = new RefObject<real8>(s_trange);
		xyz2t(TempRefObject3, TempRefObject4, slave, slaveorbit, P, MAXITER, CRITERTIM);
		s_tazi = TempRefObject3.argvalue;
		s_trange = TempRefObject4.argvalue;


	// ______Compute delta range ~= phase______
	// ______ real8 dr = dist(m_possat,pospoint) - dist(s_possat,pospoint);
	// ______ real8 phase = -pi4*(dr/LAMBDA);
	// ______  dr    == M-S         want if no flatearth M-S - flatearth = M-S-(M-S)=0
	// ______  phase == -4pi*dr/lambda == 4pi*(S-M)/lambda
	// BK: 24-9: actually defined as: phi = +pi4/lambda * (r1-r2) ???
		// real8 phase = pi4*((dist(s_possat,pospoint)-dist(m_possat,pospoint))/LAMBDA);
		//y(i,0) = pi4*((dist(s_possat,pospoint)-dist(m_possat,pospoint))/LAMBDA);
		//y(i,0) = pi4divlam*(s_possat.dist(pospoint)-m_possat.dist(pospoint));
		//y(i,0) = minpi4cdivlam * (m_trange - s_trange);
		y(i,0) = m_minpi4cdivlam *m_trange - s_minpi4cdivlam *s_trange;
		DEBUG << "l=" << line << " p=" << pixel << " t1=" << m_trange << " t2=" << s_trange << " fe=" << y(i,0) << " [rad]";
		DEBUG.print();

		// ____________________________________________________________________________________
		// _____________ Vector with h2ph factors for random number of points by FvL __________
		//_____________________________________________________________________________________

		cn Psat_master = masterorbit.argvalue.getxyz(m_tazi);
		cn Psat_slave = slaveorbit.argvalue.getxyz(s_tazi);
		real8 B = Psat_master.dist(Psat_slave); // abs. value
		// const real8 Bpar = P.dist(M) - P.dist(S);    // sign ok
		real8 Bpar = SOL*(m_trange-s_trange); // sign ok

		// ______ if (MP>SP) then S is to the right of slant line, then B perp is positive.
		cn r1 = Psat_master.min(P);
		cn r2 = Psat_slave.min(P);
		// real8 theta = Psat_master.angle(r1);  // look angle
		real8 theta = P.angle(r1); // incidence angle
		real8 theta_slave = P.angle(r2); // incidence angle slave
		real8 Bperp = (theta > theta_slave) ? Math.sqrt(sqr(B)-sqr(Bpar)) : -Math.sqrt(sqr(B)-sqr(Bpar));

		y_h2ph(i,0) = Bperp/(m_trange *SOL *Math.sin(theta));

		// ____________________________________________________________________________________
		// _____________ End added part by FvL ________________________________________________
		//_____________________________________________________________________________________

		// ______Set up system of equations______
		// ______Order unknowns: A00 A10 A01 A20 A11 A02 A30 A21 A12 A03 for degree=3______
		// ______  normalize data [-2,2] ______
		real8 posL = normalize(line,minL,maxL);
		real8 posP = normalize(pixel,minP,maxP);
		index = 0;
		for (j =0; j<=DEGREE; j++)
		  {
		  for (k =0; k<=j; k++)
			{
			A(i,index) = Math.pow(posL,real8(j-k)) * Math.pow(posP,real8(k));
			index++;
			}
		  }
		}


	  // ======Compute polynomial for these phases (LS)======
	  // ______Compute Normalmatrix, rghthandside______
	  matrix<real8> N = matTxmat(A, A);
	  matrix<real8> rhs = matTxmat(A, y);
	  matrix<real8> rhs_h2ph = matTxmat(A, y_h2ph); // Added by FvL, same A matrix can be used


	  // ______Compute solution______
	  matrix<real8> Qx_hat = N;
	  RefObject<matrix<real4>> TempRefObject5 = new RefObject<matrix<real4>>(Qx_hat);
	  choles(TempRefObject5); // Cholesky factorisation normalmatrix
	  Qx_hat = TempRefObject5.argvalue;
	  RefObject<matrix<real4>> TempRefObject6 = new RefObject<matrix<real4>>(rhs);
	  solvechol(Qx_hat, TempRefObject6); // Estimate of unknowns in rhs
	  rhs = TempRefObject6.argvalue;
	  RefObject<matrix<real4>> TempRefObject7 = new RefObject<matrix<real4>>(rhs_h2ph);
	  solvechol(Qx_hat, TempRefObject7); // Estimate of unknowns in rhs_h2ph, added by FvL
	  rhs_h2ph = TempRefObject7.argvalue;
	  RefObject<matrix<real4>> TempRefObject8 = new RefObject<matrix<real4>>(Qx_hat);
	  invertchol(TempRefObject8); // Covariance matrix
	  Qx_hat = TempRefObject8.argvalue;


	  // ______Test inverse______
	  for (i =0; i<Qx_hat.lines(); i++)
		for (j =0; j<i; j++)
		  Qx_hat(j,i) = Qx_hat(i,j); // repair Qx_hat
	  final real8 maxdev = maxMath.abs(N *Qx_hat-eye(real8Qx_hat.lines()));
	  INFO << "flatearth: max(abs(N*inv(N)-I)) = " << maxdev;
	  INFO.print();
	  if (maxdev > .01)
		{
		ERROR << "Deviation too large. Decrease degree or number of points?";
		PRINT_ERROR(ERROR.get_str())
		throw(some_error);
		}
	  else if (maxdev > .001)
		{
		WARNING << "Deviation quite large. Decrease degree or number of points?";
		WARNING.print();
		}
	  else
		{
		INFO.print("Deviation is OK.");
		}

	  // ______Some other stuff, scale is ok______
	  //  matrix<real8> Qy_hat        = A * (matxmatT(Qx_hat,A));
	  matrix<real8> y_hat = A * rhs;
	  matrix<real8> y_hat_h2ph = A * rhs_h2ph; // added by FvL
	  matrix<real8> e_hat = y - y_hat;
	  matrix<real8> e_hat_h2ph = y_h2ph - y_hat_h2ph; // added by FvL


	  // ______Overall model test (variance factor)______
	  // ... ?





	  // ______ Wrap offset ______ 
	  // BK 30/9/99 do not do this, later absolute ref. phase is used.
	  // in s2h rodriguez for example.
	  // it does not change anything for compinterfero etc.
	  //  rhs(0,0) = remainder(rhs(0,0),2*PI);



	  // ______Write results to file______
	  ofstream scratchlogfile = new ofstream("scratchlogflat", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchlogfile, "flatearth: scratchlogflat", __FILE__, __LINE__);

					 //<< "\n*_Start_" << processcontrol[pr_i_comprefpha]
	  scratchlogfile << "\n\n*******************************************************************" << "\n* FLATEARTH: " << "\n*******************************************************************" << "\nDegree_flat:\t" << DEGREE << "\nEstimated coefficients:\n" << "\nx_hat \tstd:\n";
	  for (i =0; i<Nunk; i++)
		scratchlogfile << setiosflags(ios.fixed) << setiosflags(ios.showpoint) << setiosflags(ios.right) << setw(8) << setprecision(4) << rhs(i,0) << " \t" << Math.sqrt(Qx_hat(i,i)) << "\n";

	  // ___________________ added by FvL _________________________________________________________

	  scratchlogfile << "\n" << "\nDegree_h2ph:\t" << DEGREE << "\nEstimated coefficients:\n" << "\nx_hat \tstd:\n";
	  for (i =0; i<Nunk; i++)
		scratchlogfile << setiosflags(ios.fixed) << setiosflags(ios.showpoint) << setiosflags(ios.right) << setw(8) << setprecision(4) << rhs_h2ph(i,0) << " \t" << Math.sqrt(Qx_hat(i,i)) << "\n";
	  // ___________________ end added by FvL _________________________________________________________

	  scratchlogfile << "\nCovariance matrix estimated parameters:" << "\n---------------------------------------\n";
	  for (i =0; i<Nunk; i++)
		{
		for (j =0; j<Nunk; j++)
		  {
		  scratchlogfile << setiosflags(ios.fixed) << setiosflags(ios.showpoint) << setiosflags(ios.right) << setw(8) << setprecision(4) << Qx_hat(i,j) << " ";
		  }
		scratchlogfile << "\n";
		}

	  scratchlogfile << "\nMaximum deviation N*inv(N):" << setiosflags(ios.scientific) << maxdev << "\nSome more info for each observation:" << "\nline \tpixel \tobs \t\tobs_hat \t\t err_hat\n";
	  for (i =0; i<Npoints; i++)
		scratchlogfile << Position(i,0) << "\t" << Position(i,1) << "\t" << y(i,0) << "\t" << y_hat(i,0) << "\t" << setiosflags(ios.fixed) << setiosflags(ios.showpoint) << setiosflags(ios.right) << setw(8) << setprecision(4) << e_hat(i,0) << "\n";
	  scratchlogfile << "\nMaximum absolute error: \t\t\t" << maxMath.abs(e_hat) << "\n*******************************************************************\n";
	  //  scratchlogfile.close(); // commented out for coordinates dumped below

	  ofstream scratchresfile = new ofstream("scratchresflat", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile, "flatearth: scratchresflat", __FILE__, __LINE__);

	  scratchresfile.setf(ios.scientific, ios.floatfield);
	  scratchresfile.setf(ios.right, ios.adjustfield);
	  scratchresfile.precision(8);
	  scratchresfile.width(18);

					 //<< "\n*_Start_flat_earth"
	  scratchresfile << "\n\n*******************************************************************" << "\n*_Start_" << processcontrol[pr_i_comprefpha] << "\n*******************************************************************" << "\nDegree_flat:\t" << DEGREE << "\nEstimated_coefficients_flatearth:\n";
	  int32 coeffL = 0;
	  int32 coeffP = 0;
	  for (i =0; i<Nunk; i++)
		{
		if (rhs(i,0) < 0.)
		  scratchresfile << rhs(i,0);
		else
		  scratchresfile << " " << rhs(i,0);

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

	  //_________ added by FvL _______________________________________  
	  scratchresfile << "\n" << "\nDegree_h2ph:\t" << DEGREE << "\nEstimated_coefficients_h2ph:\n";
	  coeffL = 0;
	  coeffP = 0;
	  for (i =0; i<Nunk; i++)
		{
		if (rhs_h2ph(i,0) < 0.)
		  scratchresfile << rhs_h2ph(i,0);
		else
		  scratchresfile << " " << rhs_h2ph(i,0);

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
	  //_________ end added by FvL _______________________________________  

	  scratchresfile << "*******************************************************************" << "\n* End_" << processcontrol[pr_i_comprefpha] << "_NORMAL" << "\n*******************************************************************\n";


	  // ====== Compute coordinates of corners of interferogram here ======
	  // ______ (though better place this somewhere else ....
	  real8 phi; // rad, returned
	  real8 lambda;
	  real8 height;
	  RefObject<real8> TempRefObject9 = new RefObject<real8>(phi);
	  RefObject<real8> TempRefObject10 = new RefObject<real8>(lambda);
	  RefObject<real8> TempRefObject11 = new RefObject<real8>(height);
	  lp2ell(interferogram.win.linelo, interferogram.win.pixlo, ellips, master, masterorbit, TempRefObject9, TempRefObject10, TempRefObject11, MAXITER, CRITERPOS);
	  phi = TempRefObject9.argvalue;
	  lambda = TempRefObject10.argvalue;
	  height = TempRefObject11.argvalue;
	  INFO << "Coordinates of corner interferogram: " << interferogram.win.linelo << ", " << interferogram.win.pixlo << " = " << rad2deg(phi) << ", " << rad2deg(lambda);
	  INFO.print();

	  scratchlogfile << "\n\n********************************************" << "\n* [Lat_Long] coordinates of crop [START]   *" << "\n********************************************\n";

	  scratchlogfile << "\nCoords_of_ifg_corner [l,p] : [phi,lam]: " << interferogram.win.linelo << " , " << interferogram.win.pixlo << " = " << rad2deg(phi) << " , " << rad2deg(lambda);

	  RefObject<real8> TempRefObject12 = new RefObject<real8>(phi);
	  RefObject<real8> TempRefObject13 = new RefObject<real8>(lambda);
	  RefObject<real8> TempRefObject14 = new RefObject<real8>(height);
	  lp2ell(interferogram.win.linehi, interferogram.win.pixlo, ellips, master, masterorbit, TempRefObject12, TempRefObject13, TempRefObject14, MAXITER, CRITERPOS);
	  phi = TempRefObject12.argvalue;
	  lambda = TempRefObject13.argvalue;
	  height = TempRefObject14.argvalue;
	  INFO << "\nCoordinates of corner interferogram: " << interferogram.win.linehi << ", " << interferogram.win.pixlo << " = " << rad2deg(phi) << ", " << rad2deg(lambda);
	  INFO.print();

	  scratchlogfile << "\nCoords_of_ifg_corner [l,p] : [phi,lam]: " << interferogram.win.linehi << " , " << interferogram.win.pixlo << " = " << rad2deg(phi) << " , " << rad2deg(lambda);

	  RefObject<real8> TempRefObject15 = new RefObject<real8>(phi);
	  RefObject<real8> TempRefObject16 = new RefObject<real8>(lambda);
	  RefObject<real8> TempRefObject17 = new RefObject<real8>(height);
	  lp2ell(interferogram.win.linelo, interferogram.win.pixhi, ellips, master, masterorbit, TempRefObject15, TempRefObject16, TempRefObject17, MAXITER, CRITERPOS);
	  phi = TempRefObject15.argvalue;
	  lambda = TempRefObject16.argvalue;
	  height = TempRefObject17.argvalue;
	  INFO << "\nCoordinates of corner interferogram: " << interferogram.win.linelo << ", " << interferogram.win.pixhi << " = " << rad2deg(phi) << ", " << rad2deg(lambda);
	  INFO.print();

	  scratchlogfile << "\nCoordinates of ifgs  [low,hi] corner [l,p] : [phi,lam]: " << interferogram.win.linelo << " , " << interferogram.win.pixhi << " = " << rad2deg(phi) << " , " << rad2deg(lambda);

	  RefObject<real8> TempRefObject18 = new RefObject<real8>(phi);
	  RefObject<real8> TempRefObject19 = new RefObject<real8>(lambda);
	  RefObject<real8> TempRefObject20 = new RefObject<real8>(height);
	  lp2ell(interferogram.win.linehi, interferogram.win.pixhi, ellips, master, masterorbit, TempRefObject18, TempRefObject19, TempRefObject20, MAXITER, CRITERPOS);
	  phi = TempRefObject18.argvalue;
	  lambda = TempRefObject19.argvalue;
	  height = TempRefObject20.argvalue;
	  INFO << "\nCoordinates of corner interferogram: " << interferogram.win.linehi << ", " << interferogram.win.pixhi << " = " << rad2deg(phi) << ", " << rad2deg(lambda);
	  INFO.print();

	  scratchlogfile << "\nCoords_of_ifg_corner [l,p] : [phi,lam]: " << interferogram.win.linehi << " , " << interferogram.win.pixhi << " = " << rad2deg(phi) << " , " << rad2deg(lambda);

	  scratchlogfile << "\n\n******************************************" << "\n* [Lat_Long] coordinates of crop [STOP]  *" << "\n******************************************\n";
	  // ______Tidy up______
	  scratchresfile.close(); // close res file
	  scratchlogfile.close(); // close log.out file
	  } // END flatearth



	//***************************************************************
	// *    demassist                                                 * 
	// *                                                              *
	// * Coregistration based on DEM (SRTM)                          *
	// * DEM on equiangular grid (lat/lon) assumed                    *
	// * DEM seems stored from North to South                         *
	// *                                                              *
	// * Freek van Leijen, Liu Guang, Mahmut Arikan, 21-Sep-2007      *
	// ***************************************************************
	public static void demassist(input_gen generalinput, input_ell ellips, input_demassist demassistinput, slcimage master, slcimage slave, RefObject<orbit> masterorbit, RefObject<orbit> slaveorbit)
	  {
	  TRACE_FUNCTION("demassist (FvL 21-SEP-2007)")

	  final int32 MAXITER = 10;
	  final real8 CRITERPOS = 1e-6;
	  final real8 CRITERTIM = 1e-10;

	  final real8 lat0file = demassistinput.demlatleftupper; // first pix on disk w02090
	  final real8 lon0file = demassistinput.demlonleftupper; // first pix on disk
	  final real8 DEMdeltalat = demassistinput.demdeltalat; // in radians
	  final real8 DEMdeltalon = demassistinput.demdeltalon; // in radians
	  final int32 numberoflonpixels = demassistinput.demcols; // NCOLS on file
	  final int32 numberoflatpixels = demassistinput.demrows; // NROWS on file
	  final real8 NODATA = demassistinput.demnodata; // (BK 4 may 2001)
	  final boolean outputdemi = specified(demassistinput.fodemi); // if spec. then output
	  final boolean outputrefdemhei = specified(demassistinput.forefdemhei);

	  ////const real8 m_min4picdivlam = (-4.0*PI*SOL)/master.wavelength;
	  ////const real8 s_min4picdivlam = (-4.0*PI*SOL)/slave.wavelength;

	  final real8 latNfile = lat0file-DEMdeltalat*(numberoflatpixels-1); // upper=max. lat value
	  final real8 lonNfile = lon0file+DEMdeltalon*(numberoflonpixels-1); // left=min. lon value

	  // ______ Extra info ______
	  INFO << "DEM input: w/e/s/n:          \t" << rad2deg(lon0file) << "/" << rad2deg(lonNfile) << "/" << rad2deg(latNfile) << "/" << rad2deg(lat0file);
	  INFO.print();

	  // ______ Get corners of master (approx) to select DEM ______
	  // ______ in radians (if height were zero)______
	  real8 extralat = (1.5 *DEMdeltalat + deg2rad(4.0/25.0));
	  real8 extralong = (1.5 *DEMdeltalon + deg2rad(4.0/25.0));

	  real8 phimin;
	  real8 phimax;
	  real8 lambdamin;
	  real8 lambdamax;
	  int32 indexphi0DEM;
	  int32 indexphiNDEM;
	  int32 indexlambda0DEM;
	  int32 indexlambdaNDEM;
	  RefObject<real8> TempRefObject = new RefObject<real8>(phimin);
	  RefObject<real8> TempRefObject2 = new RefObject<real8>(phimax);
	  RefObject<real8> TempRefObject3 = new RefObject<real8>(lambdamin);
	  RefObject<real8> TempRefObject4 = new RefObject<real8>(lambdamax);
	  RefObject<int32> TempRefObject5 = new RefObject<int32>(indexphi0DEM);
	  RefObject<int32> TempRefObject6 = new RefObject<int32>(indexphiNDEM);
	  RefObject<int32> TempRefObject7 = new RefObject<int32>(indexlambda0DEM);
	  RefObject<int32> TempRefObject8 = new RefObject<int32>(indexlambdaNDEM);
	  getcorners(master.currentwindow.linelo, master.currentwindow.linehi, master.currentwindow.pixlo, master.currentwindow.pixhi, extralat, extralong, lat0file, lon0file, DEMdeltalat, DEMdeltalon, numberoflatpixels, numberoflonpixels, ellips, master, masterorbit, TempRefObject, TempRefObject2, TempRefObject3, TempRefObject4, TempRefObject5, TempRefObject6, TempRefObject7, TempRefObject8);
	  phimin = TempRefObject.argvalue;
	  phimax = TempRefObject2.argvalue;
	  lambdamin = TempRefObject3.argvalue;
	  lambdamax = TempRefObject4.argvalue;
	  indexphi0DEM = TempRefObject5.argvalue;
	  indexphiNDEM = TempRefObject6.argvalue;
	  indexlambda0DEM = TempRefObject7.argvalue;
	  indexlambdaNDEM = TempRefObject8.argvalue;

	  // ______ Extra info ______
	  INFO << "DEM input required: w/e/s/n: \t" << rad2deg(lambdamin) << "/" << rad2deg(lambdamax) << "/" << rad2deg(phimin) << "/" << rad2deg(phimax);
	  INFO.print();
	  INFO << "For window (l0,lN,p0,pN):    \t" << master.currentwindow.linelo << " " << master.currentwindow.linehi << " " << master.currentwindow.pixlo << " " << master.currentwindow.pixhi;
	  INFO.print();


	  // ______ Check corners of DEM ______
	  // check if DEM is appropriate for master crop
	  // DEM should at least partially cover master crop
	  // note: phi is [90:-90]
	  if (phimax <= latNfile) // DEM is more north than master
		{
		ERROR << "master crop outside DEM: most South latitude: " << rad2deg(latNfile) << " [deg]; master crop requires: " << rad2deg(phimax) << " [deg]";
		PRINT_ERROR(ERROR.get_str())
		//throw(some_error);
		}
	  // DEM is more south than master crop
	  if (phimin >= lat0file) // largest latitude at first line of file
		{
		ERROR << "master crop outside DEM: most North latitude: " << rad2deg(lat0file) << " [deg]; master crop requires: " << rad2deg(phimax) << " [deg]";
		PRINT_ERROR(ERROR.get_str())
		//throw(some_error);
		}
	  if (lambdamax <= lon0file)
		{
		ERROR << "master crop outside DEM: most West longitude: " << rad2deg(lon0file) << " [deg]; master crop window requires: " << rad2deg(lambdamax) << " [deg]";
		PRINT_ERROR(ERROR.get_str())
		//throw(some_error);
		}
	  if (lambdamin >= lonNfile)
		{
		ERROR << "master crop outside DEM: most East longitude: " << rad2deg(lonNfile) << " [deg]; master crop window requires: " << rad2deg(lambdamin) << " [deg]";
		PRINT_ERROR(ERROR.get_str())
		//throw(some_error);
		}


	  //===================================================================
	  //============ First loop: radarcode DEM ============================
	  //============ (DEM geometry)            ============================
	  //===================================================================

	  int32 numvalid = 0; // number of good values, not NODATA in buffer
	  int32 numNODATA = 0; // number of NODATA values in buffer
	  real8 meancroppedDEM = 0.0; // to detect byte order problems, formats
	  real8 min_input_dem = 100000.0; // stats
	  real8 max_input_dem = -100000.0; // stats

	  // ______ Compute buffer size radarcoding DEM______
	  final real8 BUFFERMEMSIZE = generalinput.memory; // Bytes
	  int32 NcolsDEM = indexlambdaNDEM-indexlambda0DEM+1;
	  int32 NrowsDEM = indexphiNDEM-indexphi0DEM+1;
	  final real8 Nrows_possible_DEM = BUFFERMEMSIZE / (5 *8 *NcolsDEM);
	  int32 bufferlines = Math.ceil(Nrows_possible_DEM); // [MA] checked ok. Sinces SLC is not multilooked, see comprefdem for solution
	  if (bufferlines>NrowsDEM)
		  bufferlines =NrowsDEM;
	  int32 numfullbuffers = NrowsDEM / bufferlines;
	  int32 restlines = NrowsDEM % bufferlines;
	  int32 extrabuffer = (restlines == 0) ? 0 : 1;

	  // ______ Extra info ______
	  INFO << "Radar coding of DEM in: " << numfullbuffers << " buffers of " << bufferlines << " lines and " << extrabuffer << " extra buffer of " << restlines << " lines.";
	  INFO.print();

	  INFO << "NcolsDEM: " << NcolsDEM;
	  INFO.print();
	  INFO << "NrowsDEM: " << NrowsDEM;
	  INFO.print();



	  // ______ Open (temporary) output files ______
	  // DEM heights 
	  ofstream demofile;
	  RefObject<ofstream> TempRefObject9 = new RefObject<ofstream>(demofile);
	  openfstream(TempRefObject9, demassistinput.fodem, generalinput.overwrit); // dem_crop radarcoded
	  demofile = TempRefObject9.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(demofile, demassistinput.fodem, __FILE__, __LINE__);

	  // master line coordinates of DEM
	  ofstream masterdemlineoutfile = new ofstream("master_demline.temp", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(masterdemlineoutfile, "master_demline.temp", __FILE__, __LINE__);

	  // master pixel coordinates of DEM
	  ofstream masterdempixeloutfile = new ofstream("master_dempixel.temp", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(masterdempixeloutfile, "master_dempixel.temp", __FILE__, __LINE__);

	  // delta line coordinates of DEM
	  ofstream deltademlineoutfile = new ofstream("delta_demline.temp", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(deltademlineoutfile, "delta_demline.temp", __FILE__, __LINE__);

	  // master pixel coordinates of DEM
	  ofstream deltadempixeloutfile = new ofstream("delta_dempixel.temp", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(deltadempixeloutfile, "delta_dempixel.temp", __FILE__, __LINE__);



	  // ______ DEM loop per buffer ______
//C++ TO JAVA CONVERTER NOTE: 'register' variable declarations are not supported in Java:
//ORIGINAL LINE: register int32 j,i;
	  int32 j; // DEM index grid counter, register j first to ensure allocation
	  int32 i;
	  for (register int32 buffer =0; buffer<numfullbuffers+extrabuffer; ++buffer)
		{

		 // Determine indices for buffer
		final int32 indexphi0BUFFER = indexphi0DEM+buffer *bufferlines;
		final int32 blines = (buffer == numfullbuffers) ? restlines : bufferlines;
		final int32 indexphiNBUFFER = indexphi0BUFFER+(blines-1);
		matrix<real4> DEM = new matrix(blines, NcolsDEM);

		// ______ Extra info ______
		PROGRESS << "Buffer# [l0:lN, p0:pN]: " << buffer+1 << " [" << indexphi0BUFFER << ": " << indexphiNBUFFER << ", " << indexlambda0DEM << ": " << indexlambdaNDEM << "]";
		PROGRESS.print();

		// ______ lat/lon for first pixel in matrix read from file ______
		// ______ upper is max. latitude, left is min. longitude ______
		final real8 upperleftphi = lat0file-indexphi0BUFFER *DEMdeltalat;
		final real8 upperleftlambda = lon0file+indexlambda0DEM *DEMdeltalon;

		window zerooffset = new window(0,0,0,0);
		window winfromfile = new window(indexphi0BUFFER,indexphiNBUFFER, indexlambda0DEM,indexlambdaNDEM);

		// ______ Read in grdfile of DEM in matrix R4 (raw data, no header) _______
		// ______ added formats (BK 4-May-2001) ______
		PROGRESS << "Reading crop of DEM for buffer: " << buffer+1;
		PROGRESS.print();
		DEBUG.print("Reading input DEM into real4 matrix (buffer).");
		switch (demassistinput.iformatflag)
		  {
		  // ______ Read as short BE, then convert to host order ______
		  case FORMATI2_BIGENDIAN:
			{
			matrix<int16> DEMi2 = new matrix(blines, NcolsDEM);
			readfile(DEMi2,demassistinput.firefdem,numberoflatpixels,winfromfile,zerooffset);
			for (int32 iii =0; iii<DEM.lines(); ++iii)
			  for (int32 jjj =0; jjj<DEM.pixels(); ++jjj)
				DEM(iii,jjj) = real4(ntohs(DEMi2(iii,jjj))); // cast to real4
			DEMi2.resize(1, 1); // dealloc...
			INFO.print("Read crop of input DEM: format: SHORT SIGNED INT.");
			break;
			}

		  case FORMATI2:
			{
			matrix<int16> DEMi2 = new matrix(blines, NcolsDEM);
			readfile(DEMi2,demassistinput.firefdem,numberoflatpixels,winfromfile,zerooffset);
			for (int32 iii =0; iii<DEM.lines(); ++iii)
			  for (int32 jjj =0; jjj<DEM.pixels(); ++jjj)
				DEM(iii,jjj) = DEMi2(iii,jjj); // cast to real4
			DEMi2.resize(1, 1); // dealloc...
			INFO.print("Read crop of input DEM: format: SHORT SIGNED INT.");
			break;
			}

		  case FORMATR4:
			readfile(DEM,demassistinput.firefdem,numberoflatpixels,winfromfile,zerooffset);
			INFO.print("Read crop of input DEM: format: REAL4.");
			break;
		  case FORMATR8:
			{
			matrix<real8> DEMr8 = new matrix(blines, NcolsDEM);
			readfile(DEMr8,demassistinput.firefdem,numberoflatpixels,winfromfile,zerooffset);
			for (int32 iii =0; iii<DEM.lines(); ++iii)
			  for (int32 jjj =0; jjj<DEM.pixels(); ++jjj)
				DEM(iii,jjj) = DEMr8(iii,jjj); // cast to real4
			DEMr8.resize(1, 1); // dealloc...
			INFO.print("Read crop of input DEM: format: REAL8.");
			break;
			}
		  default:
			PRINT_ERROR("totally impossible, checked input.")
			//throw(unhandled_case_error);
		  }


		// ----- Loop over DEM for stats ------------------------
		real8 min_dem_buffer = 100000.0;
		real8 max_dem_buffer = -100000.0;
		for (i =0; i<DEM.lines(); ++i)
		  {
		  // ----- Loop over oversampled matrix in x ------
		  for (j =0; j<DEM.pixels(); ++j)
			{
			if(DEM(i,j)!=NODATA)
			  {
			  numvalid++;
			  meancroppedDEM += DEM(i,j); // divide by numvalid later
			  if (DEM(i,j)<min_dem_buffer) //buffer
				  min_dem_buffer =DEM(i,j);
			  if (DEM(i,j)>max_dem_buffer) // stats
				  max_dem_buffer =DEM(i,j);
			  }
			else
			  {
			  numNODATA++;
			  }
			} //loop dem for stats
		  } //loop dem for stats
		min_input_dem = min(min_input_dem,min_dem_buffer); //global stats
		max_input_dem = max(max_input_dem,max_dem_buffer); //global stats


		// ====== Radarcoding DEM ==============================
		// ______ DEM contains values from leftupper with ______
		// ______ spacing (DEMdeltalat,DEMdeltalon) ______
		// ______ Transform DEM to l,p,refphase ______
		PROGRESS.print("Converting DEM to radar system for this buffer.");
		final int32 NpointsDEM = DEM.size();
		final int32 NpixelsDEM = DEM.pixels();
		// ______ Extra info ______
		INFO << "Number of points in DEM: " << NpointsDEM;
		INFO.print();

		matrix<real8> masterDEMline = new matrix(DEM.lines(), DEM.pixels());
		matrix<real8> masterDEMpixel = new matrix(DEM.lines(), DEM.pixels());
		matrix<real8> deltaDEMline = new matrix(DEM.lines(), DEM.pixels());
		matrix<real8> deltaDEMpixel = new matrix(DEM.lines(), DEM.pixels());

		// --- Loop DEM ---
		real8 phi;
		real8 lambda;
		real8 height;
		real8 m_l;
		real8 m_p;
		real8 s_l;
		real8 s_p;


		phi = upperleftphi;
		for (i =0; i<DEM.lines(); ++i)
		  {
		  if ((i%100)==0)
			{
			// ______ Extra info ______
			PROGRESS << "Radarcoding DEM line: " << i << " (" << Math.floor(.5+(100.*real8(i)/real8DEM.lines())) << "%)";
			PROGRESS.print();
			}

		  lambda = upperleftlambda;
		  for (j =0; j<DEM.pixels(); ++j)
			{
		height = DEM(i,j);
		RefObject<real8> TempRefObject10 = new RefObject<real8>(m_l);
		RefObject<real8> TempRefObject11 = new RefObject<real8>(m_p);
		ell2lp(TempRefObject10, TempRefObject11, ellips, master, masterorbit, phi, lambda, height, MAXITER, CRITERTIM);
		m_l = TempRefObject10.argvalue;
		m_p = TempRefObject11.argvalue;
		RefObject<real8> TempRefObject12 = new RefObject<real8>(s_l);
		RefObject<real8> TempRefObject13 = new RefObject<real8>(s_p);
		ell2lp(TempRefObject12, TempRefObject13, ellips, slave, slaveorbit, phi, lambda, height, MAXITER, CRITERTIM);
		s_l = TempRefObject12.argvalue;
		s_p = TempRefObject13.argvalue;
		masterDEMline(i,j) = m_l;
		masterDEMpixel(i,j) = m_p;
			  deltaDEMline(i,j) = s_l-m_l;
			  deltaDEMpixel(i,j) = s_p-m_p;

			  lambda += DEMdeltalon;
			} // loop DEM pixels

		  // ______ update latitude of next line ______
		  phi -= DEMdeltalat; // upper left is max. value
		  } // loop DEM lines


		// Write results to output files 
		PROGRESS << "Writing radar coded DEM to file: " << buffer+1;
		PROGRESS.print();

		demofile << DEM;
		masterdemlineoutfile << masterDEMline;
		masterdempixeloutfile << masterDEMpixel;
		deltademlineoutfile << deltaDEMline;
		deltadempixeloutfile << deltaDEMpixel;

		masterDEMline.resize(1, 1); //deallocate
		masterDEMpixel.resize(1, 1); //deallocate
		deltaDEMline.resize(1, 1); //deallocate
		deltaDEMpixel.resize(1, 1); //deallocate
		DEM.resize(1, 1); //deallocate
		} // buffer loop

	  demofile.close();
	  masterdemlineoutfile.close();
	  masterdempixeloutfile.close();
	  deltademlineoutfile.close();
	  deltadempixeloutfile.close();


	  //===================================================================
	  //============ End first loop: radarcode DEM ========================
	  //===================================================================


	  //===================================================================
	  //============ Second loop: interpolation               =============
	  //============ (radar geometry)                         =============
	  //===================================================================

	  INFO << "Start interpolation....";
	  INFO.print();

	  // ______ Line/pixel of first point in original master coordinates ______
	  final int32 Nlinesml = master.currentwindow.lines();
	  final int32 Npixelsml = master.currentwindow.pixels();

	  final real8 veryfirstline = new real8(master.currentwindow.linelo);
	  final real8 verylastline = new real8(master.currentwindow.linehi);
	  final real8 firstpixel = new real8(master.currentwindow.pixlo);
	  final real8 lastpixel = new real8(master.currentwindow.pixhi);


	  //Determine range-azimuth spacing ratio, needed for proper triangulation
	  cn P1;
	  cn P2;
	  cn P3;
	  cn P4;
	  RefObject<cn> TempRefObject14 = new RefObject<cn>(P1);
	  lp2xyz(veryfirstline, firstpixel, ellips, master, masterorbit, TempRefObject14, MAXITER, CRITERPOS);
	  P1 = TempRefObject14.argvalue;
	  RefObject<cn> TempRefObject15 = new RefObject<cn>(P2);
	  lp2xyz(veryfirstline, lastpixel, ellips, master, masterorbit, TempRefObject15, MAXITER, CRITERPOS);
	  P2 = TempRefObject15.argvalue;
	  RefObject<cn> TempRefObject16 = new RefObject<cn>(P3);
	  lp2xyz(verylastline, firstpixel, ellips, master, masterorbit, TempRefObject16, MAXITER, CRITERPOS);
	  P3 = TempRefObject16.argvalue;
	  RefObject<cn> TempRefObject17 = new RefObject<cn>(P4);
	  lp2xyz(verylastline, lastpixel, ellips, master, masterorbit, TempRefObject17, MAXITER, CRITERPOS);
	  P4 = TempRefObject17.argvalue;

	  final real8 r_spacing = (P1.min(P2).norm() + P3.min(P4).norm()) / 2 /(lastpixel - firstpixel);
	  final real8 az_spacing = (P1.min(P3).norm() + P2.min(P4).norm()) /2 /(verylastline - veryfirstline);
	  final real8 r_az_ratio = r_spacing/az_spacing;

	  INFO << "Master azimuth spacing: " << az_spacing;
	  INFO.print();
	  INFO << "Master range spacing: " << r_spacing;
	  INFO.print();
	  INFO << "Range-azimuth spacing ratio: " << r_az_ratio;
	  INFO.print();

	  // ______ Compute buffer size interpolation______
	  final real8 Nlinesml_possible = BUFFERMEMSIZE / (6 *8 *Npixelsml);
	  bufferlines = int32(Math.ceil(Nlinesml_possible));
	  if (bufferlines > Nlinesml)
		  bufferlines =Nlinesml;
	  numfullbuffers = Nlinesml / bufferlines;
	  restlines = Nlinesml % bufferlines;
	  extrabuffer = (restlines == 0) ? 0 : 1;

	  // ______ Extra info ______
	  INFO << "Interpolation in: " << numfullbuffers << " buffers of " << bufferlines << " lines and " << extrabuffer << " extra buffer of " << restlines << " lines.";
	  INFO.print();


	  // ______ Open output files ______
	  ofstream deltalineofile = new ofstream("delta_line.raw", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(deltalineofile, "delta_line.raw", __FILE__, __LINE__);

	  ofstream deltapixelofile = new ofstream("delta_pixel.raw", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(deltapixelofile, "delta_pixel.raw", __FILE__, __LINE__);

	  // if request for height in radar coordinates l,p
	  ofstream refdemheiofile;
	  if (outputrefdemhei ==true)
		{
		RefObject<ofstream> TempRefObject18 = new RefObject<ofstream>(refdemheiofile);
		openfstream(TempRefObject18, demassistinput.forefdemhei, generalinput.overwrit);
		refdemheiofile = TempRefObject18.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		bk_assert(refdemheiofile, demassistinput.forefdemhei, __FILE__, __LINE__);
		}

	  // ______ interpolation loop per buffer ______
	  for (register int32 buffer = 0; buffer < numfullbuffers + extrabuffer; ++buffer)
		{

		// Determine indices for buffer
		final int32 blines = (buffer == numfullbuffers) ? restlines : bufferlines;
		final real8 firstline_buffer = veryfirstline+buffer *bufferlines;
		final real8 lastline_buffer = firstline_buffer+blines-1;

		// ______ Extra info ______
		PROGRESS << "Interpolation buffer# [l0:lN, p0:pN]: " << buffer+1 << " [" << firstline_buffer << ": " << lastline_buffer << ", " << firstpixel << ": " << lastpixel << "]";
		PROGRESS.print();

		// Get corners of buffer
		real8 phimin_az;
		real8 phimax_az;
		real8 lambdamin_az;
		real8 lambdamax_az;
		RefObject<real8> TempRefObject19 = new RefObject<real8>(phimin_az);
		RefObject<real8> TempRefObject20 = new RefObject<real8>(phimax_az);
		RefObject<real8> TempRefObject21 = new RefObject<real8>(lambdamin_az);
		RefObject<real8> TempRefObject22 = new RefObject<real8>(lambdamax);
		RefObject<int32> TempRefObject23 = new RefObject<int32>(indexphi0DEM);
		RefObject<int32> TempRefObject24 = new RefObject<int32>(indexphiNDEM);
		RefObject<int32> TempRefObject25 = new RefObject<int32>(indexlambda0DEM);
		RefObject<int32> TempRefObject26 = new RefObject<int32>(indexlambdaNDEM);
		getcorners(firstline_buffer, lastline_buffer, firstpixel, lastpixel, extralat, extralong, phimax, lambdamin, DEMdeltalat, DEMdeltalon, NrowsDEM, NcolsDEM, ellips, master, masterorbit, TempRefObject19, TempRefObject20, TempRefObject21, TempRefObject22, TempRefObject23, TempRefObject24, TempRefObject25, TempRefObject26);
		phimin_az = TempRefObject19.argvalue;
		phimax_az = TempRefObject20.argvalue;
		lambdamin_az = TempRefObject21.argvalue;
		lambdamax = TempRefObject22.argvalue;
		indexphi0DEM = TempRefObject23.argvalue;
		indexphiNDEM = TempRefObject24.argvalue;
		indexlambda0DEM = TempRefObject25.argvalue;
		indexlambdaNDEM = TempRefObject26.argvalue;

		window zerooffset = new window(0,0,0,0);
		window winfromfile = new window(indexphi0DEM,indexphiNDEM, indexlambda0DEM,indexlambdaNDEM);
		final int32 NrowsDEM_buffer = indexphiNDEM-indexphi0DEM+1;
		final int32 NcolsDEM_buffer = indexlambdaNDEM-indexlambda0DEM+1;

		PROGRESS << "Reading input for interpolation buffer: " << buffer+1;
		PROGRESS.print();

		// read x,y
		matrix<real8> DEMline_buffer = new matrix(NrowsDEM_buffer, NcolsDEM_buffer);
		matrix<real8> DEMpixel_buffer = new matrix(NrowsDEM_buffer, NcolsDEM_buffer);

		readfile(DEMline_buffer,"master_demline.temp",NrowsDEM,winfromfile,zerooffset);
		readfile(DEMpixel_buffer,"master_dempixel.temp",NrowsDEM,winfromfile,zerooffset);

		// read z (multiple, number can easily be increased, e.g. simulated intensity)
		int32 Nz = 2; //number of z
		matrix<real8> input_buffer = new matrix(NrowsDEM_buffer *Nz, NcolsDEM_buffer);
		matrix<real8> temp_input_buffer = new matrix(NrowsDEM_buffer, NcolsDEM_buffer);
		if (outputrefdemhei ==true)
		  {
		  Nz += 1;
		  input_buffer.resize(NrowsDEM_buffer *Nz, NcolsDEM_buffer);
		  }

		readfile(temp_input_buffer,"delta_demline.temp",NrowsDEM,winfromfile,zerooffset);
		input_buffer.setdata(0, 0, temp_input_buffer);
		readfile(temp_input_buffer,"delta_dempixel.temp",NrowsDEM,winfromfile,zerooffset);
		input_buffer.setdata(NrowsDEM_buffer, 0, temp_input_buffer);
		Nz = 2;
		if (outputrefdemhei ==true)
		  {
			Nz += 1;
			/// i would like to use real4, test later on
			matrix<real4> dem_input = new matrix(NrowsDEM_buffer, NcolsDEM_buffer);
			readfile(dem_input,demassistinput.fodem,NrowsDEM,winfromfile,zerooffset);
			for (register int32 i =0 ; i < NrowsDEM_buffer ; i ++)
			  for(register int32 j = 0; j < NcolsDEM_buffer; j++)
				temp_input_buffer(i,j) = real8(dem_input(i,j));
			input_buffer.setdata(NrowsDEM_buffer * (Nz-1), 0, temp_input_buffer);
		  }

		// initialize output array
		Nz = 2;
		matrix<real8> output_buffer = new matrix(blines * Nz, Npixelsml);
		if (outputrefdemhei ==true)
		  {
			Nz += 1;
			output_buffer.resize(blines * Nz, Npixelsml);
		  }

		// interpolation
		RefObject<matrix<real8>> TempRefObject27 = new RefObject<matrix<real8>>(output_buffer);
		griddatalinear(DEMline_buffer, DEMpixel_buffer, input_buffer, firstline_buffer, lastline_buffer, firstpixel, lastpixel, 1, 1, r_az_ratio, 0, NODATA, TempRefObject27);
		output_buffer = TempRefObject27.argvalue;

		deltalineofile << output_buffer(window(0, blines - 1, 0, Npixelsml -1));
		deltapixelofile << output_buffer(window(blines, 2 * blines - 1, 0, Npixelsml -1));
		Nz = 2;
		if (outputrefdemhei ==true)
		  {
			Nz += 1;
			refdemheiofile << output_buffer(window((Nz-1) * blines,Nz * blines - 1, 0, Npixelsml -1));
		  }

		DEMline_buffer.resize(1, 1);
		DEMpixel_buffer.resize(1, 1);
		input_buffer.resize(1, 1);
		temp_input_buffer.resize(1, 1);
		output_buffer.resize(1, 1);

	  } // end loop azimuth direction

	  INFO << "Closing output files";
	  INFO.print();

	  deltalineofile.close();
	  deltapixelofile.close();
	  if (outputrefdemhei ==true) // For Zbigniew Perski
		refdemheiofile.close();

	  //===================================================================
	  //============ End second loop: interpolation           =============
	  //============ (radar geometry)                         =============
	  //===================================================================


	  //===================================================================
	  //============ Determine inverse transformation         =============
	  //============ (slave corners only, needed for overlap) =============
	  //===================================================================

	  real8 line;
	  real8 pixel;
	  real8 deltaline_slave00;
	  real8 deltapixel_slave00;
	  real8 deltaline_slave0N;
	  real8 deltapixel_slave0N;
	  real8 deltaline_slaveN0;
	  real8 deltapixel_slaveN0;
	  real8 deltaline_slaveNN;
	  real8 deltapixel_slaveNN;
	  real8 phimin_az;
	  real8 phimax_az;
	  real8 lambdamin_az;
	  real8 lambdamax_az;


	  for (register int16 corner = 0 ; corner < 4 ; corner ++)
		{

		  PROGRESS << "Radarcoding slave corner: " << corner+1;
		  PROGRESS.print();

		  switch (corner)
			{
			case 0:
			  {
				line =slave.currentwindow.linelo;
				pixel =slave.currentwindow.pixlo;
				break;
			  }
			case 1:
			  {
				line =slave.currentwindow.linelo;
				pixel =slave.currentwindow.pixhi;
				break;
			  }
			case 2:
			  {
				line =slave.currentwindow.linehi;
				pixel =slave.currentwindow.pixlo;
				break;
			  }
			case 3:
			  {
				line =slave.currentwindow.linehi;
				pixel =slave.currentwindow.pixhi;
				break;
			  }
			default:
			  PRINT_ERROR("totally impossible, checked input.");

			}

		  //use getcorners with line,line,pixel,pixel for single point
		  RefObject<real8> TempRefObject28 = new RefObject<real8>(phimin_az);
		  RefObject<real8> TempRefObject29 = new RefObject<real8>(phimax_az);
		  RefObject<real8> TempRefObject30 = new RefObject<real8>(lambdamin_az);
		  RefObject<real8> TempRefObject31 = new RefObject<real8>(lambdamax);
		  RefObject<int32> TempRefObject32 = new RefObject<int32>(indexphi0DEM);
		  RefObject<int32> TempRefObject33 = new RefObject<int32>(indexphiNDEM);
		  RefObject<int32> TempRefObject34 = new RefObject<int32>(indexlambda0DEM);
		  RefObject<int32> TempRefObject35 = new RefObject<int32>(indexlambdaNDEM);
		  getcorners(line, line, pixel, pixel, extralat, extralong, lat0file, lon0file, DEMdeltalat, DEMdeltalon, numberoflatpixels, numberoflonpixels, ellips, slave, slaveorbit, TempRefObject28, TempRefObject29, TempRefObject30, TempRefObject31, TempRefObject32, TempRefObject33, TempRefObject34, TempRefObject35);
		  phimin_az = TempRefObject28.argvalue;
		  phimax_az = TempRefObject29.argvalue;
		  lambdamin_az = TempRefObject30.argvalue;
		  lambdamax = TempRefObject31.argvalue;
		  indexphi0DEM = TempRefObject32.argvalue;
		  indexphiNDEM = TempRefObject33.argvalue;
		  indexlambda0DEM = TempRefObject34.argvalue;
		  indexlambdaNDEM = TempRefObject35.argvalue;


		  NcolsDEM = indexlambdaNDEM-indexlambda0DEM+1;
		  NrowsDEM = indexphiNDEM-indexphi0DEM+1;
		  final real8 upperleftphi = lat0file-indexphi0DEM *DEMdeltalat;
		  final real8 upperleftlambda = lon0file+indexlambda0DEM *DEMdeltalon;

		  window zerooffset = new window(0,0,0,0);
		  window winfromfile = new window(indexphi0DEM,indexphiNDEM, indexlambda0DEM,indexlambdaNDEM);


		  // ______ Read in DEM in matrix R4 (raw data, no header) _______

		  matrix<real4> DEM = new matrix(NrowsDEM, NcolsDEM);

		  switch (demassistinput.iformatflag)
			{
			  // ______ Read as short BE, then convert to host order ______
			case FORMATI2_BIGENDIAN:
			  {
				matrix<int16> DEMi2 = new matrix(NrowsDEM, NcolsDEM);
				readfile(DEMi2,demassistinput.firefdem,numberoflatpixels,winfromfile,zerooffset);
				for (int32 iii =0; iii<DEM.lines(); ++iii)
				  for (int32 jjj =0; jjj<DEM.pixels(); ++jjj)
					DEM(iii,jjj) = real4(ntohs(DEMi2(iii,jjj))); // cast to real4
				DEMi2.resize(1, 1); // dealloc...
				INFO.print("Read crop of input DEM: format: SHORT SIGNED INT.");
				break;
			  }

			case FORMATI2:
			  {
				matrix<int16> DEMi2 = new matrix(NrowsDEM, NcolsDEM);
				readfile(DEMi2,demassistinput.firefdem,numberoflatpixels,winfromfile,zerooffset);
				for (int32 iii =0; iii<DEM.lines(); ++iii)
				  for (int32 jjj =0; jjj<DEM.pixels(); ++jjj)
					DEM(iii,jjj) = DEMi2(iii,jjj); // cast to real4
				DEMi2.resize(1, 1); // dealloc...
				INFO.print("Read crop of input DEM: format: SHORT SIGNED INT.");
				break;
			  }

			case FORMATR4:
			  readfile(DEM,demassistinput.firefdem,numberoflatpixels,winfromfile,zerooffset);
			  INFO.print("Read crop of input DEM: format: REAL4.");
			  break;
			case FORMATR8:
			  {
				matrix<real8> DEMr8 = new matrix(NrowsDEM, NcolsDEM);
				readfile(DEMr8,demassistinput.firefdem,numberoflatpixels,winfromfile,zerooffset);
				for (int32 iii =0; iii<DEM.lines(); ++iii)
				  for (int32 jjj =0; jjj<DEM.pixels(); ++jjj)
					DEM(iii,jjj) = DEMr8(iii,jjj); // cast to real4
				DEMr8.resize(1, 1); // dealloc...
				INFO.print("Read crop of input DEM: format: REAL8.");
				break;
			  }
			default:
			  PRINT_ERROR("totally impossible, checked input.");
			  //throw(unhandled_case_error);
			}


		  // radarcode dem
		  matrix<real8> slaveDEMline = new matrix(DEM.lines(), DEM.pixels());
		  matrix<real8> slaveDEMpixel = new matrix(DEM.lines(), DEM.pixels());
		  matrix<real8> deltaDEMline = new matrix(DEM.lines(), DEM.pixels());
		  matrix<real8> deltaDEMpixel = new matrix(DEM.lines(), DEM.pixels());

		  // --- Loop DEM ---
		  real8 phi;
		  real8 lambda;
		  real8 height;
		  real8 m_l;
		  real8 m_p;
		  real8 s_l;
		  real8 s_p;


		  phi = upperleftphi;
		  for (i =0; i<DEM.lines(); ++i)
			{

			  lambda = upperleftlambda;
			  for (j =0; j<DEM.pixels(); ++j)
				{
				  height = DEM(i,j);
				  RefObject<real8> TempRefObject36 = new RefObject<real8>(m_l);
				  RefObject<real8> TempRefObject37 = new RefObject<real8>(m_p);
				  ell2lp(TempRefObject36, TempRefObject37, ellips, master, masterorbit, phi, lambda, height, MAXITER, CRITERTIM);
				  m_l = TempRefObject36.argvalue;
				  m_p = TempRefObject37.argvalue;
				  RefObject<real8> TempRefObject38 = new RefObject<real8>(s_l);
				  RefObject<real8> TempRefObject39 = new RefObject<real8>(s_p);
				  ell2lp(TempRefObject38, TempRefObject39, ellips, slave, slaveorbit, phi, lambda, height, MAXITER, CRITERTIM);
				  s_l = TempRefObject38.argvalue;
				  s_p = TempRefObject39.argvalue;
				  slaveDEMline(i,j) = s_l;
				  slaveDEMpixel(i,j) = s_p;
				  deltaDEMline(i,j) = m_l-s_l;
				  deltaDEMpixel(i,j) = m_p-s_p;

				  lambda += DEMdeltalon;
				} // loop DEM pixels

			  // ______ update latitude of next line ______
			  phi -= DEMdeltalat; // upper left is max. value
			} // loop DEM lines



		  // interpolate to slave corner
		  matrix<real8> input_buffer = new matrix(DEM.lines()*2, DEM.pixels());
		  input_buffer.setdata(0, 0, deltaDEMline);
		  input_buffer.setdata(DEM.lines(), 0, deltaDEMpixel);

		  matrix<real8> output_buffer = new matrix(2, 1);

		  RefObject<matrix<real8>> TempRefObject40 = new RefObject<matrix<real8>>(output_buffer);
		  griddatalinear(slaveDEMline, slaveDEMpixel, input_buffer, line, line, pixel, pixel, 1, 1, r_az_ratio, 0, NODATA, TempRefObject40);
		  output_buffer = TempRefObject40.argvalue;

		  switch (corner)
			{
			case 0:
			  {
				deltaline_slave00 = output_buffer(0,0);
				deltapixel_slave00 = output_buffer(1,0);
				INFO << "Deltaline_slave00: " << deltaline_slave00;
				INFO.print();
				INFO << "Deltapixel_slave00: " << deltapixel_slave00;
				INFO.print();
				break;
			  }
			case 1:
			  {
				deltaline_slave0N = output_buffer(0,0);
				deltapixel_slave0N = output_buffer(1,0);
				INFO << "Deltaline_slave0N: " << deltaline_slave0N;
				INFO.print();
				INFO << "Deltapixel_slave0N: " << deltapixel_slave0N;
				INFO.print();
				break;
			  }
			case 2:
			  {
				deltaline_slaveN0 = output_buffer(0,0);
				deltapixel_slaveN0 = output_buffer(1,0);
				INFO << "Deltaline_slaveN0: " << deltaline_slaveN0;
				INFO.print();
				INFO << "Deltapixel_slaveN0: " << deltapixel_slaveN0;
				INFO.print();
				break;
			  }
			case 3:
			  {
				deltaline_slaveNN = output_buffer(0,0);
				deltapixel_slaveNN = output_buffer(1,0);
				INFO << "Deltaline_slaveNN: " << deltaline_slaveNN;
				INFO.print();
				INFO << "Deltapixel_slaveNN: " << deltapixel_slaveNN;
				INFO.print();
				break;
			  }
			default:
			  PRINT_ERROR("totally impossible, checked input.");

			}

		}

	  //===================================================================
	  //============ End determine inverse transformation     =============
	  //============ (slave corners only, needed for overlap) =============
	  //===================================================================


	  // ====== Write output information ======
	  String croppeddemi = new String(new char[ONE27]);
	  croppeddemi = "NO output requested";
	  if (outputdemi)
		  croppeddemi = demassistinput.fodemi;
	  INFO << "Min. value of input DEM covering master: " << min_input_dem;
	  INFO.print();
	  INFO << "Max. value of input DEM covering master: " << max_input_dem;
	  INFO.print();

	  ofstream scratchlogfile = new ofstream("scratchlogdemassist", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchlogfile, "demassist: scratchlogdemassist", __FILE__, __LINE__);
	  scratchlogfile << "\n*******************************************************************" << "\n* " << processcontrol[pr_i_demassist] << "\n*******************************************************************" << "\n1) DEM source file:                   \t" << demassistinput.firefdem << "\nFormat:                               \t";
		switch (demassistinput.iformatflag)
		  {
		  case FORMATI2:
			{
			scratchlogfile << "SHORT SIGNED INTEGER (HOST ENDIANNESS)";
			break;
			}
		  case FORMATI2_BIGENDIAN:
			{
			scratchlogfile << "SHORT SIGNED INTEGER, BIG ENDIAN";
			break;
			}
		  case FORMATR4:
			{
			scratchlogfile << "REAL4 SIGNED FLOAT";
			break;
			}
		  case FORMATR8:
			{
			scratchlogfile << "REAL8 SIGNED DOUBLE";
			break;
			}
		  default:
			{
			scratchlogfile << "UNKNOWN? IMPOSSIBLE...";
			break;
			}
		  }
	//    << "\nMean value:                           \t" <<  meancroppedDEM
	  scratchlogfile << "\nByte order:                           \t" << "check yourself..." << "\nNumber of lines:                      \t" << numberoflatpixels << "\nNumber of pixels:                     \t" << numberoflonpixels << "\nResolution latitude:                  \t" << rad2deg(DEMdeltalat) << " [deg]" << "\nResolution longitude:                 \t" << rad2deg(DEMdeltalon) << " [deg]" << "\nMost West point in input DEM:         \t" << rad2deg(lon0file) << "\nMost East point in input DEM:         \t" << rad2deg(lonNfile) << "\nMost South point in input DEM:        \t" << rad2deg(latNfile) << "\nMost North point in input DEM:        \t" << rad2deg(lat0file) << "\nMin. value of input DEM covering master: " << min_input_dem << "\nMax. value of input DEM covering master: " << max_input_dem << "\n2) Output file cropped DEM:           \t" << demassistinput.fodem << "\nFormat:                               \t" << "REAL4" << "\nByte order:                           \t" << "(same as host)" << "\n3) Output file interpolated crop DEM: \t" << croppeddemi << "\nFormat:                               \t" << "REAL4" << "\nByte order:                           \t" << "(same as host)" << "\nNumber of lines (multilooked):        \t" << Nlinesml << "\nNumber of pixels (multilooked):       \t" << Npixelsml << "\nDeltaline_slave00_dem:                    \t" << deltaline_slave00 << "\nDeltapixel_slave00_dem:                   \t" << deltapixel_slave00 << "\nDeltaline_slave0N_dem:                    \t" << deltaline_slave0N << "\nDeltapixel_slave0N_dem:                   \t" << deltapixel_slave0N << "\nDeltaline_slaveN0_dem:                    \t" << deltaline_slaveN0 << "\nDeltapixel_slaveN0_dem:                   \t" << deltapixel_slaveN0 << "\nDeltaline_slaveNN_dem:                    \t" << deltaline_slaveNN << "\nDeltapixel_slaveNN_dem:                   \t" << deltapixel_slaveNN << "\n*******************************************************************\n\n";
	  scratchlogfile.close();


	  ofstream scratchresfile = new ofstream("scratchresdemassist", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile, "demassist: scratchresdemassist", __FILE__, __LINE__);
	  scratchresfile << "\n\n*******************************************************************" << "\n*_Start_" << processcontrol[pr_i_demassist] << "\n*******************************************************************";
	  scratchresfile << "\nDEM source file:                      \t" << demassistinput.firefdem << "\nMin. of input DEM:                    \t" << min_input_dem << "\nMax. of input DEM:                    \t" << max_input_dem << "\nFirst_line (w.r.t. original_master):   \t" << master.currentwindow.linelo << "\nLast_line (w.r.t. original_master):   \t" << master.currentwindow.linehi << "\nFirst_pixel (w.r.t. original_master): \t" << master.currentwindow.pixlo << "\nLast_pixel (w.r.t. original_master):  \t" << master.currentwindow.pixhi << "\nNumber of lines:        \t" << Nlinesml << "\nNumber of pixels:       \t" << Npixelsml << "\nDeltaline_slave00_dem:      \t" << deltaline_slave00 << "\nDeltapixel_slave00_dem:     \t" << deltapixel_slave00 << "\nDeltaline_slave0N_dem:      \t" << deltaline_slave0N << "\nDeltapixel_slave0N_dem:     \t" << deltapixel_slave0N << "\nDeltaline_slaveN0_dem:      \t" << deltaline_slaveN0 << "\nDeltapixel_slaveN0_dem:     \t" << deltapixel_slaveN0 << "\nDeltaline_slaveNN_dem:      \t" << deltaline_slaveNN << "\nDeltapixel_slaveNN_dem:     \t" << deltapixel_slaveNN << "\n*******************************************************************" << "\n* End_" << processcontrol[pr_i_demassist] << "_NORMAL" << "\n*******************************************************************\n";
	  scratchresfile.close();


	  } // END demassist


	//***************************************************************
	// *    radarcodedem                                              *
	// *                                                              *
	// * Compute reference phase based on DEM (SRTM)                  *
	// * DEM on equiangular grid (lat/lon) assumed                    *
	// * DEM seems stored from North to South                         *
	// *                                                              *
	// * Freek van Leijen, 26-Sep-2007                                *
	// ***************************************************************
	public static void radarcodedem(input_gen generalinput, input_ell ellips, input_comprefdem refdeminput, slcimage master, slcimage slave, productinfo interferogram, RefObject<orbit> masterorbit, RefObject<orbit> slaveorbit)
	  {
	  TRACE_FUNCTION("radarcodedem (FvL 26-Sep-2007)")

	  final int32 MAXITER = 10;
	  final real8 CRITERPOS = 1e-6;
	  final real8 CRITERTIM = 1e-10;

	  final real8 lat0file = refdeminput.demlatleftupper; // first pix on disk w02090
	  final real8 lon0file = refdeminput.demlonleftupper; // first pix on disk
	  final real8 DEMdeltalat = refdeminput.demdeltalat; // in radians
	  final real8 DEMdeltalon = refdeminput.demdeltalon; // in radians
	  final int32 numberoflonpixels = refdeminput.demcols; // NCOLS on file
	  final int32 numberoflatpixels = refdeminput.demrows; // NROWS on file
	  final real8 NODATA = refdeminput.demnodata; // (BK 4 may 2001)
	  boolean onlyrefphasetopo = !refdeminput.includerefpha; // true: phase DEM w.r.t. ellipsoid
	  final boolean outputdemi = specified(refdeminput.fodemi); // if spec. then output
	  final boolean outputh2ph = specified(refdeminput.foh2ph); // if spec. then output, added by FvL
	  final boolean outputrefdemhei = specified(refdeminput.forefdemhei);

	  // _____ start added by MA _____
	  boolean mlookedIFG = false; // true: ifg is multilooked
	  int32 mlL = interferogram.multilookL; // initialize multilookfactor
	  int32 mlP = interferogram.multilookP;
	  final int32 &ifgmlL = interferogram.multilookL; // multilookfactor of interferogram
	  final int32 &ifgmlP = interferogram.multilookP; // multilookfactor of interferogram
	  if (ifgmlL != 1 || ifgmlP != 1) // [MA] additional entry for Coherence comptation using refdem.
		{ // always do computation without multilooking
		  mlL = 1; // set multilookfactor for interpolation
		  mlP = 1; // set multilookfactor for interpolation
		  mlookedIFG = true; // dealing with mlooked ifg.
		}
	  // _____ end added by MA _____

	  final real8 m_min4picdivlam = (-4.0 *PI *SOL)/master.wavelength;
	  final real8 s_min4picdivlam = (-4.0 *PI *SOL)/slave.wavelength;
	  DEBUG << "master wavelength = " << master.wavelength;
	  DEBUG.print();
	  DEBUG << "slave  wavelength = " << slave.wavelength;
	  DEBUG.print();

	  final real8 latNfile = lat0file-DEMdeltalat*(numberoflatpixels-1); // upper=max. lat value
	  final real8 lonNfile = lon0file+DEMdeltalon*(numberoflonpixels-1); // left=min. lon value

	  // ______ Extra info ______
	  INFO << "DEM input: w/e/s/n:          \t" << rad2deg(lon0file) << "/" << rad2deg(lonNfile) << "/" << rad2deg(latNfile) << "/" << rad2deg(lat0file);
	  INFO.print();

	  // ______ Get corners of interferogram (approx) to select DEM ______
	  // ______ in radians (if height were zero)______
	  real8 extralat = (1.5 *DEMdeltalat + deg2rad(4.0/25.0));
	  real8 extralong = (1.5 *DEMdeltalon + deg2rad(4.0/25.0));

	  real8 phimin;
	  real8 phimax;
	  real8 lambdamin;
	  real8 lambdamax;
	  int32 indexphi0DEM;
	  int32 indexphiNDEM;
	  int32 indexlambda0DEM;
	  int32 indexlambdaNDEM;
	  final uint ifglinelo = interferogram.win.linelo; // [MA] win no-mlooked master coords
	  final uint ifglinehi = interferogram.win.linehi;
	  final uint ifgpixlo = interferogram.win.pixlo;
	  final uint ifgpixhi = interferogram.win.pixhi;

	  RefObject<real8> TempRefObject = new RefObject<real8>(phimin);
	  RefObject<real8> TempRefObject2 = new RefObject<real8>(phimax);
	  RefObject<real8> TempRefObject3 = new RefObject<real8>(lambdamin);
	  RefObject<real8> TempRefObject4 = new RefObject<real8>(lambdamax);
	  RefObject<int32> TempRefObject5 = new RefObject<int32>(indexphi0DEM);
	  RefObject<int32> TempRefObject6 = new RefObject<int32>(indexphiNDEM);
	  RefObject<int32> TempRefObject7 = new RefObject<int32>(indexlambda0DEM);
	  RefObject<int32> TempRefObject8 = new RefObject<int32>(indexlambdaNDEM);
	  getcorners(ifglinelo, ifglinehi, ifgpixlo, ifgpixhi, extralat, extralong, lat0file, lon0file, DEMdeltalat, DEMdeltalon, numberoflatpixels, numberoflonpixels, ellips, master, masterorbit, TempRefObject, TempRefObject2, TempRefObject3, TempRefObject4, TempRefObject5, TempRefObject6, TempRefObject7, TempRefObject8);
	  phimin = TempRefObject.argvalue;
	  phimax = TempRefObject2.argvalue;
	  lambdamin = TempRefObject3.argvalue;
	  lambdamax = TempRefObject4.argvalue;
	  indexphi0DEM = TempRefObject5.argvalue;
	  indexphiNDEM = TempRefObject6.argvalue;
	  indexlambda0DEM = TempRefObject7.argvalue;
	  indexlambdaNDEM = TempRefObject8.argvalue;

	  // ______ Extra info ______
	  INFO << "DEM input required: w/e/s/n: \t" << rad2deg(lambdamin) << "/" << rad2deg(lambdamax) << "/" << rad2deg(phimin) << "/" << rad2deg(phimax);
	  INFO.print();
	  INFO << "For window (l0,lN,p0,pN):    \t" << ifglinelo << " " << ifglinehi << " " << ifgpixlo << " " << ifgpixhi;
	  INFO.print();


	  // ______ Check corners of DEM ______
	  // check if DEM is appropriate for interferogram
	  // DEM should at least partially cover IFG
	  // note: phi is [90:-90]
	  if (phimax <= latNfile) // DEM is more north than IFG
		{
		ERROR << "IFG outside DEM: most South latitude: " << rad2deg(latNfile) << " [deg]; IFG requires: " << rad2deg(phimax) << " [deg]";
		PRINT_ERROR(ERROR.get_str())
		throw(some_error);
		}
	  // DEM is more south than IFG
	  if (phimin >= lat0file) // largest latitude at first line of file
		{
		ERROR << "IFG outside DEM: most North latitude: " << rad2deg(lat0file) << " [deg]; IFG requires: " << rad2deg(phimax) << " [deg]";
		PRINT_ERROR(ERROR.get_str())
		throw(some_error);
		}
	  if (lambdamax <= lon0file)
		{
		ERROR << "IFG outside DEM: most West longitude: " << rad2deg(lon0file) << " [deg]; IFG window requires: " << rad2deg(lambdamax) << " [deg]";
		PRINT_ERROR(ERROR.get_str())
		throw(some_error);
		}
	  if (lambdamin >= lonNfile)
		{
		ERROR << "IFG outside DEM: most East longitude: " << rad2deg(lonNfile) << " [deg]; IFG window requires: " << rad2deg(lambdamin) << " [deg]";
		PRINT_ERROR(ERROR.get_str())
		throw(some_error);
		}


	  //===================================================================
	  //============ First loop: radarcode DEM ============================
	  //============ (DEM geometry)            ============================
	  //===================================================================

	  int32 numvalid = 0; // number of good values, not NODATA in buffer
	  int32 numNODATA = 0; // number of NODATA values in buffer
	  real8 meancroppedDEM = 0.0; // to detect byte order problems, formats
	  real8 min_input_dem = 100000.0; // stats
	  real8 max_input_dem = -100000.0; // stats

	  // ______ Compute buffer size radarcoding DEM______
	  final real8 BUFFERMEMSIZE = generalinput.memory; // Bytes
	  final int32 NcolsDEM = indexlambdaNDEM-indexlambda0DEM+1;
	  final int32 NrowsDEM = indexphiNDEM-indexphi0DEM+1;
	  final real8 Nrows_possible_DEM = BUFFERMEMSIZE / (5 *8 *NcolsDEM);
	  int32 bufferlines = Math.ceil(Nrows_possible_DEM);
	  if (bufferlines>NrowsDEM)
		  bufferlines =NrowsDEM;
	  int32 numfullbuffers = NrowsDEM / bufferlines;
	  int32 restlines = NrowsDEM % bufferlines;
	  int32 extrabuffer = (restlines == 0) ? 0 : 1;

	  // ______ Extra info ______
	  INFO << "Radar coding of DEM in: " << numfullbuffers << " buffers of " << bufferlines << " lines and " << extrabuffer << " extra buffer of " << restlines << " lines.";
	  INFO.print();
	  INFO << "NcolsDEM: " << NcolsDEM;
	  INFO.print();
	  INFO << "NrowsDEM: " << NrowsDEM;
	  INFO.print();



	  // ______ Open (temporary) output files ______
	  // DEM heights 
	  ofstream demofile;
	  RefObject<ofstream> TempRefObject9 = new RefObject<ofstream>(demofile);
	  openfstream(TempRefObject9, refdeminput.fodem, generalinput.overwrit);
	  demofile = TempRefObject9.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(demofile, refdeminput.fodem, __FILE__, __LINE__);

	  // master line coordinates of DEM
	  ofstream masterdemlineoutfile = new ofstream("master_demline.temp2", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(masterdemlineoutfile, "master_demline.temp2", __FILE__, __LINE__);

	  // master pixel coordinates of DEM
	  ofstream masterdempixeloutfile = new ofstream("master_dempixel.temp2", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(masterdempixeloutfile, "master_dempixel.temp2", __FILE__, __LINE__);

	  // ref phase in DEM geometry
	  ofstream demrefphaseoutfile = new ofstream("dem_refphase.temp", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(demrefphaseoutfile, "dem_refphase.temp", __FILE__, __LINE__);

	  // h2ph factor in DEM geometry
	  ofstream demh2phoutfile;
	  if (outputh2ph ==true)
		{
		  RefObject<ofstream> TempRefObject10 = new RefObject<ofstream>(demh2phoutfile);
		  openfstream(TempRefObject10, "dem_h2ph.temp", generalinput.overwrit);
		  demh2phoutfile = TempRefObject10.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  bk_assert(demh2phoutfile, "dem_h2ph.temp", __FILE__, __LINE__);
		}

	  // ______ DEM loop per buffer ______
//C++ TO JAVA CONVERTER NOTE: 'register' variable declarations are not supported in Java:
//ORIGINAL LINE: register int32 j,i;
	  int32 j; // DEM index grid counter, register j first to ensure allocation
	  int32 i;
	  for (register int32 buffer =0; buffer<numfullbuffers+extrabuffer; ++buffer)
		{

		 // Determine indices for buffer
		final int32 indexphi0BUFFER = indexphi0DEM+buffer *bufferlines;
		final int32 blines = (buffer == numfullbuffers) ? restlines : bufferlines;
		final int32 indexphiNBUFFER = indexphi0BUFFER+(blines-1);
		matrix<real4> DEM = new matrix(blines, NcolsDEM);

		// ______ Extra info ______
		PROGRESS << "Buffer# [l0:lN, p0:pN]: " << buffer+1 << " [" << indexphi0BUFFER << ": " << indexphiNBUFFER << ", " << indexlambda0DEM << ": " << indexlambdaNDEM << "]";
		PROGRESS.print();

		// ______ lat/lon for first pixel in matrix read from file ______
		// ______ upper is max. latitude, left is min. longitude ______
		final real8 upperleftphi = lat0file-indexphi0BUFFER *DEMdeltalat;
		final real8 upperleftlambda = lon0file+indexlambda0DEM *DEMdeltalon;

		window zerooffset = new window(0,0,0,0);
		window winfromfile = new window(indexphi0BUFFER,indexphiNBUFFER, indexlambda0DEM,indexlambdaNDEM);

		// ______ Read in grdfile of DEM in matrix R4 (raw data, no header) _______
		// ______ added formats (BK 4-May-2001) ______
		PROGRESS << "Reading crop of DEM for buffer: " << buffer+1;
		PROGRESS.print();
		DEBUG.print("Reading input DEM into real4 matrix (buffer).");
		switch (refdeminput.iformatflag)
		  {
		  // ______ Read as short BE, then convert to host order ______
		  case FORMATI2_BIGENDIAN:
			{
			matrix<int16> DEMi2 = new matrix(blines, NcolsDEM);
			readfile(DEMi2,refdeminput.firefdem,numberoflatpixels,winfromfile,zerooffset);
			for (int32 iii =0; iii<DEM.lines(); ++iii)
			  for (int32 jjj =0; jjj<DEM.pixels(); ++jjj)
				DEM(iii,jjj) = real4(ntohs(DEMi2(iii,jjj))); // cast to real4
			DEMi2.resize(1, 1); // dealloc...
			INFO.print("Read crop of input DEM: format: SHORT SIGNED INT.");
			break;
			}

		  case FORMATI2:
			{
			matrix<int16> DEMi2 = new matrix(blines, NcolsDEM);
			readfile(DEMi2,refdeminput.firefdem,numberoflatpixels,winfromfile,zerooffset);
			for (int32 iii =0; iii<DEM.lines(); ++iii)
			  for (int32 jjj =0; jjj<DEM.pixels(); ++jjj)
				DEM(iii,jjj) = DEMi2(iii,jjj); // cast to real4
			DEMi2.resize(1, 1); // dealloc...
			INFO.print("Read crop of input DEM: format: SHORT SIGNED INT.");
			break;
			}

		  case FORMATR4:
			readfile(DEM,refdeminput.firefdem,numberoflatpixels,winfromfile,zerooffset);
			INFO.print("Read crop of input DEM: format: REAL4.");
			break;
		  case FORMATR8:
			{
			matrix<real8> DEMr8 = new matrix(blines, NcolsDEM);
			readfile(DEMr8,refdeminput.firefdem,numberoflatpixels,winfromfile,zerooffset);
			for (int32 iii =0; iii<DEM.lines(); ++iii)
			  for (int32 jjj =0; jjj<DEM.pixels(); ++jjj)
				DEM(iii,jjj) = DEMr8(iii,jjj); // cast to real4
			DEMr8.resize(1, 1); // dealloc...
			INFO.print("Read crop of input DEM: format: REAL8.");
			break;
			}
		  default:
			PRINT_ERROR("totally impossible, checked input.")
			//throw(unhandled_case_error);
		  }


		// ----- Loop over DEM for stats ------------------------
		real8 min_dem_buffer = 100000.0;
		real8 max_dem_buffer = -100000.0;
		for (i =0; i<DEM.lines(); ++i)
		  {
		  // ----- Loop over oversampled matrix in x ------
		  for (j =0; j<DEM.pixels(); ++j)
			{
			if(DEM(i,j)!=NODATA)
			  {
			  numvalid++;
			  meancroppedDEM += DEM(i,j); // divide by numvalid later
			  if (DEM(i,j)<min_dem_buffer) //buffer
				  min_dem_buffer =DEM(i,j);
			  if (DEM(i,j)>max_dem_buffer) // stats
				  max_dem_buffer =DEM(i,j);
			  }
			else
			  {
			  numNODATA++;
			  }
			} //loop dem for stats
		  } //loop dem for stats
		min_input_dem = min(min_input_dem,min_dem_buffer); //global stats
		max_input_dem = max(max_input_dem,max_dem_buffer); //global stats


		// ====== Radarcoding DEM ==============================
		// ______ DEM contains values from leftupper with ______
		// ______ spacing (DEMdeltalat,DEMdeltalon) ______
		// ______ Transform DEM to l,p,refphase ______
		PROGRESS.print("Converting DEM to radar system for this buffer.");
		final int32 NpointsDEM = DEM.size();
		final int32 NpixelsDEM = DEM.pixels();
		// ______ Extra info ______
		INFO << "Number of points in DEM: " << NpointsDEM;
		INFO.print();

		matrix<real8> masterDEMline = new matrix(DEM.lines(), DEM.pixels());
		matrix<real8> masterDEMpixel = new matrix(DEM.lines(), DEM.pixels());
		matrix<real8> ref_phase_array = new matrix(DEM.lines(), DEM.pixels());
		matrix<real8> h2ph_array = new matrix(DEM.lines(), DEM.pixels());

		// --- Loop DEM ---
		cn P;
		real8 phi;
		real8 lambda;
		real8 height;
		real8 l;
		real8 p;
		real8 ref_phase;


		phi = upperleftphi;
		for (i =0; i<DEM.lines(); ++i)
		  {
		  if ((i%100)==0)
			{
			// ______ Extra info ______
			PROGRESS << "Radarcoding DEM line: " << i << " (" << Math.floor(.5+(100.*real8(i)/real8DEM.lines())) << "%)";
			PROGRESS.print();
			}

		  lambda = upperleftlambda;
		  for (j =0; j<DEM.pixels(); ++j)
			{
		height = DEM(i,j);
		RefObject<real8> TempRefObject11 = new RefObject<real8>(l);
		RefObject<real8> TempRefObject12 = new RefObject<real8>(p);
		ell2lp(TempRefObject11, TempRefObject12, ellips, master, masterorbit, phi, lambda, height, MAXITER, CRITERTIM);
		l = TempRefObject11.argvalue;
		p = TempRefObject12.argvalue;
		masterDEMline(i,j) = l;
		masterDEMpixel(i,j) = p;
		P = ellips.ell2xyz(phi,lambda,height); // returns P(x,y,z)

		real8 t_range_master;
		real8 t_azi_master;
		RefObject<real8> TempRefObject13 = new RefObject<real8>(t_azi_master);
		RefObject<real8> TempRefObject14 = new RefObject<real8>(t_range_master);
		xyz2t(TempRefObject13, TempRefObject14, master, masterorbit, P, MAXITER, CRITERTIM);
		t_azi_master = TempRefObject13.argvalue;
		t_range_master = TempRefObject14.argvalue;
		real8 t_azi_slave;
		real8 t_range_slave;
		RefObject<real8> TempRefObject15 = new RefObject<real8>(t_azi_slave);
		RefObject<real8> TempRefObject16 = new RefObject<real8>(t_range_slave);
		xyz2t(TempRefObject15, TempRefObject16, slave, slaveorbit, P, MAXITER, CRITERTIM);
		t_azi_slave = TempRefObject15.argvalue;
		t_range_slave = TempRefObject16.argvalue;

	  if (outputh2ph ==true)
	{
		// compute h2ph factor
		cn Psat_master = masterorbit.argvalue.getxyz(t_azi_master);
		cn Psat_slave = slaveorbit.argvalue.getxyz(t_azi_slave);
		real8 B = Psat_master.dist(Psat_slave);
		real8 Bpar = SOL*(t_range_master-t_range_slave);
		cn r1 = Psat_master.min(P);
		cn r2 = Psat_slave.min(P);
		// real8 theta = Psat_master.angle(r1);  // look angle
		real8 theta = P.angle(r1); // incidence angle
		real8 theta_slave = P.angle(r2); // incidence angle slave
		real8 Bperp = (theta > theta_slave) ? Math.sqrt(sqr(B)-sqr(Bpar)) : -Math.sqrt(sqr(B)-sqr(Bpar));

		h2ph_array(i,j) = Bperp/(t_range_master *SOL *Math.sin(theta));
	}


		if (onlyrefphasetopo) // do not include flat earth phase
		  {
			RefObject<cn> TempRefObject17 = new RefObject<cn>(P);
			lp2xyz(l, p, ellips, master, masterorbit, TempRefObject17, MAXITER, CRITERPOS); // P returned
			P = TempRefObject17.argvalue;

			real8 t_range_flatearth;
			real8 t_azi_dummy;

			RefObject<real8> TempRefObject18 = new RefObject<real8>(t_azi_dummy);
			RefObject<real8> TempRefObject19 = new RefObject<real8>(t_range_flatearth);
			xyz2t(TempRefObject18, TempRefObject19, slave, slaveorbit, P, MAXITER, CRITERTIM); // P on h=0
			t_azi_dummy = TempRefObject18.argvalue;
			t_range_flatearth = TempRefObject19.argvalue;
			ref_phase = s_min4picdivlam *t_range_flatearth- s_min4picdivlam *t_range_slave;
		  }
		else // include flatearth, ref.pha = phi_topo+phi_flatearth
		  {
			ref_phase = m_min4picdivlam *master.pix2tr(p)- s_min4picdivlam *t_range_slave;
		  }

			ref_phase_array(i,j) = ref_phase;

			lambda += DEMdeltalon;
			} // loop DEM pixels

		  // ______ update latitude of next line ______
		  phi -= DEMdeltalat; // upper left is max. value
		  } // loop DEM lines


		// Write results to output files 
		PROGRESS << "Writing radar coded DEM to file: " << buffer+1;
		PROGRESS.print();
		demofile << DEM;
		masterdemlineoutfile << masterDEMline;
		masterdempixeloutfile << masterDEMpixel;
		demrefphaseoutfile << ref_phase_array;
	if (outputh2ph ==true)
	  demh2phoutfile << h2ph_array;

		masterDEMline.resize(1, 1); //deallocate
		masterDEMpixel.resize(1, 1); //deallocate
		DEM.resize(1, 1); //deallocate
		ref_phase_array(1,1); //deallocate
		h2ph_array(1,1); //deallocate
		} // buffer loop

	  demofile.close();
	  masterdemlineoutfile.close();
	  masterdempixeloutfile.close();
	  demrefphaseoutfile.close();
	if (outputh2ph ==true)
	  demh2phoutfile.close();


	  //===================================================================
	  //============ End first loop: radarcode DEM ========================
	  //============ (DEM geometry)            ============================
	  //===================================================================


	  //===================================================================
	  //============ Second loop: interpolation               =============
	  //============ (radar geometry)                         =============
	  //===================================================================

	  INFO << "Start interpolation....";
	  INFO.print();

	  // ______ Line/pixel of first point in original master coordinates ______
	  // ______ maybe this should be changed to be x+(ml/2) ?? but depends on
	  // ______ definition of range_to_first_bin is to center or start..
	  // Bert Kampes, 08-Apr-2005: chose center by adding ml/2
	  final int32 Nlinesml = interferogram.win.lines() / mlL; // ifg lines when mlL = 1 (no multilooking)
	  final int32 Npixelsml = interferogram.win.pixels() / mlP;
	  final int32 ifgNlinesml = interferogram.win.lines() / ifgmlL; // for the result file when mlL != 1
	  final int32 ifgNpixelsml = interferogram.win.pixels() / ifgmlP;
	  final real8 offset = 0;

	//cerr << "xNFO:   linesnoml:    " << Nlinesml <<  " pixnoml:    " <<  Npixelsml << endl;
	//cerr << "xNFO:   ifglinesml: " << ifgNlinesml <<  " ifgpixml: " <<  ifgNpixelsml << endl;

	  final real8 veryfirstline = real8(ifglinelo) + (real8(mlL)-1.0)/2.0;
	  final real8 verylastline = veryfirstline + real8((Nlinesml-1)*mlL);
	  final real8 firstpixel = real8(ifgpixlo) + (real8(mlP)-1.0)/2.0;
	  final real8 lastpixel = firstpixel + real8((Npixelsml-1)*mlP);

	//cerr << "xNFO:   vl0:vlN " << veryfirstline << ":"  << verylastline << " p1:pN " <<  firstpixel << ":" << lastpixel << endl;

	  //Determine range-azimuth spacing ratio, needed for proper triangulation
	  cn P1;
	  cn P2;
	  cn P3;
	  cn P4;
	  RefObject<cn> TempRefObject20 = new RefObject<cn>(P1);
	  lp2xyz(veryfirstline, firstpixel, ellips, master, masterorbit, TempRefObject20, MAXITER, CRITERPOS);
	  P1 = TempRefObject20.argvalue;
	  RefObject<cn> TempRefObject21 = new RefObject<cn>(P2);
	  lp2xyz(veryfirstline, lastpixel, ellips, master, masterorbit, TempRefObject21, MAXITER, CRITERPOS);
	  P2 = TempRefObject21.argvalue;
	  RefObject<cn> TempRefObject22 = new RefObject<cn>(P3);
	  lp2xyz(verylastline, firstpixel, ellips, master, masterorbit, TempRefObject22, MAXITER, CRITERPOS);
	  P3 = TempRefObject22.argvalue;
	  RefObject<cn> TempRefObject23 = new RefObject<cn>(P4);
	  lp2xyz(verylastline, lastpixel, ellips, master, masterorbit, TempRefObject23, MAXITER, CRITERPOS);
	  P4 = TempRefObject23.argvalue;

	  final real8 r_spacing = (P1.min(P2).norm() + P3.min(P4).norm()) / 2 /(lastpixel - firstpixel);
	  final real8 az_spacing = (P1.min(P3).norm() + P2.min(P4).norm()) /2 /(verylastline - veryfirstline);
	  final real8 r_az_ratio = r_spacing/az_spacing;

	  INFO << "Interferogram azimuth spacing: " << az_spacing;
	  INFO.print();
	  INFO << "Interferogram range spacing: " << r_spacing;
	  INFO.print();
	  INFO << "Range-azimuth spacing ratio: " << r_az_ratio;
	  INFO.print();

	  // ______ Compute buffer size interpolation______
	  final real8 Nlinesml_possible = BUFFERMEMSIZE / (6 *8 *Npixelsml);
	  bufferlines = int32(Math.ceil(Nlinesml_possible)); // initialized

	  if (mlookedIFG == true) // if ifg is multilooked by a factor
		{
		  bufferlines = int32(Math.floor(Nlinesml_possible/ifgmlL) * ifgmlL); // [HB] Herrmann provided the fix:
																		   // Successive bufferlines must have a size which is multiple of multilooking factor
																		   // unless data can fit completely to an initial single buffer. 
																		   // Extra buffer will scale correctly and rounding due to multilooking 
																		   // will yield correct number lines for the output file.
																		   // Ex: floor(2097/25)*2 + floor(1578/25)    =  229 lines (wrong)  (bufferlines/mlL)*Nbuffers+extrabufferlines/mlL 
		} // Ex: floor(2075/25)*25*2 + floor(1622/25) = 230 lines (correct)
	  else // no-multilooking
		{
		  bufferlines = int32(Math.floor(Nlinesml_possible)); // [MA] instead of ceil, prefered floor to use less mem
		}
	  if (bufferlines > Nlinesml) // if bufferlines > Nlines then shrink bufferlines to Nlines, no extra buffer requested.
		  bufferlines =Nlinesml;
	  numfullbuffers = Nlinesml / bufferlines;
	  restlines = Nlinesml % bufferlines; // the number of lines in extra buffer
	  extrabuffer = (restlines == 0) ? 0 : 1;

	  // ______ Extra info ______
	  INFO << "Interpolation in: " << numfullbuffers << " buffers of " << bufferlines << " lines and " << extrabuffer << " extra buffer of " << restlines << " lines.";
	  INFO.print();


	  // ______ Open output files ______
	  ofstream refdemofile; // refdem phase
	  RefObject<ofstream> TempRefObject24 = new RefObject<ofstream>(refdemofile);
	  openfstream(TempRefObject24, refdeminput.forefdem, generalinput.overwrit);
	  refdemofile = TempRefObject24.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(refdemofile, refdeminput.forefdem, __FILE__, __LINE__);

	  // _____ start added by MA _____
	  ofstream refdemofilenoML; // [MA] refdem phase no-multilooked
	  { // local scope practice
		String fname = String(refdeminput.forefdem) + ".noML"; // new name as m_s_refdemphase.raw.noML
		if (mlookedIFG == true) // if ifg is multilooked by a factor
		  {
			RefObject<ofstream> TempRefObject25 = new RefObject<ofstream>(refdemofilenoML);
			openfstream(TempRefObject25, fname.c_str(), generalinput.overwrit);
			refdemofilenoML = TempRefObject25.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
			bk_assert(refdemofilenoML, fname.c_str(), __FILE__, __LINE__);
		  }
		else // no-multilooking
		  {
			if(!remove(fname.c_str())) // when success report removed.
			  {
			  WARNING << "Removed existing " << fname << "file.";
			  WARNING.print();
			  }
		  }
	  }
	  // _____ end added by MA _____

	  // if request for height in radar coordinates l,p
	  ofstream refdemheiofile; // for Zbigniew Perski
	  if (outputrefdemhei ==true)
		{
		RefObject<ofstream> TempRefObject26 = new RefObject<ofstream>(refdemheiofile);
		openfstream(TempRefObject26, refdeminput.forefdemhei, generalinput.overwrit);
		refdemheiofile = TempRefObject26.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		bk_assert(refdemheiofile, refdeminput.forefdemhei, __FILE__, __LINE__);
		}

	  // if request for h2ph in radar coordinates l,p
	  ofstream h2phofile;
	  if (outputh2ph ==true)
		{
		RefObject<ofstream> TempRefObject27 = new RefObject<ofstream>(h2phofile);
		openfstream(TempRefObject27, refdeminput.foh2ph, generalinput.overwrit);
		h2phofile = TempRefObject27.argvalue;
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		bk_assert(h2phofile, refdeminput.foh2ph, __FILE__, __LINE__);
		}

	  // ______ interpolation loop per buffer ______
	  for (register int32 buffer = 0; buffer < numfullbuffers + extrabuffer; ++buffer)
		{

		// Determine indices for buffer
		final int32 blines = (buffer == numfullbuffers) ? restlines : bufferlines;
		final real8 firstline_buffer = veryfirstline+buffer *bufferlines *mlL;
		final real8 lastline_buffer = firstline_buffer+(blines-1)*mlL;

		// ______ Extra info ______
		PROGRESS << "Interpolation buffer# [l0:lN, p0:pN]: " << buffer+1 << " [" << firstline_buffer << ": " << lastline_buffer << ", " << firstpixel << ": " << lastpixel << "]";
		PROGRESS.print();

		// Get corners of buffer
		real8 phimin_az;
		real8 phimax_az;
		real8 lambdamin_az;
		real8 lambdamax_az;
		RefObject<real8> TempRefObject28 = new RefObject<real8>(phimin_az);
		RefObject<real8> TempRefObject29 = new RefObject<real8>(phimax_az);
		RefObject<real8> TempRefObject30 = new RefObject<real8>(lambdamin_az);
		RefObject<real8> TempRefObject31 = new RefObject<real8>(lambdamax);
		RefObject<int32> TempRefObject32 = new RefObject<int32>(indexphi0DEM);
		RefObject<int32> TempRefObject33 = new RefObject<int32>(indexphiNDEM);
		RefObject<int32> TempRefObject34 = new RefObject<int32>(indexlambda0DEM);
		RefObject<int32> TempRefObject35 = new RefObject<int32>(indexlambdaNDEM);
		getcorners(firstline_buffer+offset, lastline_buffer+offset, firstpixel+offset, lastpixel+offset, extralat, extralong, phimax, lambdamin, DEMdeltalat, DEMdeltalon, NrowsDEM, NcolsDEM, ellips, master, masterorbit, TempRefObject28, TempRefObject29, TempRefObject30, TempRefObject31, TempRefObject32, TempRefObject33, TempRefObject34, TempRefObject35);
		phimin_az = TempRefObject28.argvalue;
		phimax_az = TempRefObject29.argvalue;
		lambdamin_az = TempRefObject30.argvalue;
		lambdamax = TempRefObject31.argvalue;
		indexphi0DEM = TempRefObject32.argvalue;
		indexphiNDEM = TempRefObject33.argvalue;
		indexlambda0DEM = TempRefObject34.argvalue;
		indexlambdaNDEM = TempRefObject35.argvalue;

		window zerooffset = new window(0,0,0,0);
		window winfromfile = new window(indexphi0DEM,indexphiNDEM, indexlambda0DEM,indexlambdaNDEM);
		final int32 NrowsDEM_buffer = indexphiNDEM-indexphi0DEM+1;
		final int32 NcolsDEM_buffer = indexlambdaNDEM-indexlambda0DEM+1;

		PROGRESS << "Reading input for interpolation buffer: " << buffer+1;
		PROGRESS.print();

		// read x,y
		matrix<real8> DEMline_buffer = new matrix(NrowsDEM_buffer, NcolsDEM_buffer);
		matrix<real8> DEMpixel_buffer = new matrix(NrowsDEM_buffer, NcolsDEM_buffer);

		readfile(DEMline_buffer,"master_demline.temp2",NrowsDEM,winfromfile,zerooffset);
		readfile(DEMpixel_buffer,"master_dempixel.temp2",NrowsDEM,winfromfile,zerooffset);

		// read z (multiple, number can easily be increased, e.g. simulated intensity)
		int32 Nz = 1; //number of z
		matrix<real8> input_buffer = new matrix(NrowsDEM_buffer *Nz, NcolsDEM_buffer);
		matrix<real8> temp_input_buffer = new matrix(NrowsDEM_buffer, NcolsDEM_buffer);
		if (outputrefdemhei ==true)
		  {
		  Nz += 1;
		  input_buffer.resize(NrowsDEM_buffer *Nz, NcolsDEM_buffer);
		  }
		if (outputh2ph ==true)
		  {
		  Nz += 1;
		  input_buffer.resize(NrowsDEM_buffer *Nz, NcolsDEM_buffer);
		  }

		readfile(temp_input_buffer,"dem_refphase.temp",NrowsDEM,winfromfile,zerooffset);
		input_buffer.setdata(0, 0, temp_input_buffer);
		Nz = 1;
		if (outputrefdemhei ==true)
		  {
			Nz += 1;
			/// i would like to use real4, test later on
			matrix<real4> dem_input = new matrix(NrowsDEM_buffer, NcolsDEM_buffer);
			readfile(dem_input,refdeminput.fodem,NrowsDEM,winfromfile,zerooffset);
			for (register int32 i =0 ; i < NrowsDEM_buffer ; i ++)
			  for(register int32 j = 0; j < NcolsDEM_buffer; j++)
				temp_input_buffer(i,j) = real8(dem_input(i,j));
			input_buffer.setdata(NrowsDEM_buffer * (Nz-1), 0, temp_input_buffer);
		  }
		if (outputh2ph ==true)
		  {
			Nz += 1;
			readfile(temp_input_buffer,"dem_h2ph.temp",NrowsDEM,winfromfile,zerooffset);
			input_buffer.setdata(NrowsDEM_buffer * (Nz-1), 0, temp_input_buffer);
		  }

		// initialize output array
		Nz = 1;
		matrix<real8> output_buffer = new matrix(blines * Nz, Npixelsml);

		if (outputrefdemhei ==true)
		  {
			Nz += 1;
			output_buffer.resize(blines * Nz, Npixelsml);
		  }
		if (outputh2ph ==true)
		  {
			Nz += 1;
			output_buffer.resize(blines * Nz, Npixelsml);
		  }


		// interpolation
		RefObject<matrix<real8>> TempRefObject36 = new RefObject<matrix<real8>>(output_buffer);
		griddatalinear(DEMline_buffer, DEMpixel_buffer, input_buffer, firstline_buffer, lastline_buffer, firstpixel, lastpixel, mlL, mlP, r_az_ratio, offset, NODATA, TempRefObject36);
		output_buffer = TempRefObject36.argvalue;


		//MA multilooking will start here.
		//MA cast all output files to type real4 
		matrix<real4> output_layer = new matrix(blines, Npixelsml);

		Nz = 1; // matrix 3rd dimension counter, 1 --> PHASE

		for (register int32 i =0 ; i < blines ; i++)
		  for(register int32 j = 0; j < Npixelsml; j++)
			output_layer(i,j) = real4(output_buffer(i,j)); // real8 --> real4
	   //     convert_type(output_buffer,output_layer); // MA

	//cerr << "refphase: blines: " << blines << " Npixelsml " << Npixelsml << endl;

	  // _____ start added by MA _____
		if (mlookedIFG == true) // [MA] if ifg is multilooked by a factor
		  {
			refdemofile << multilook(output_layer, ifgmlL, ifgmlP); // multilook to output
			refdemofilenoML << output_layer; // default output: PHASE
		  }
		else
		  {
			refdemofile << output_layer;
			//refdemofile << output_buffer(window(0, blines - 1, 0, Npixelsml -1 ));
		  }
	  // _____ end added by MA _____

		if (outputrefdemhei ==true)
		  {
			Nz += 1;
			for (register int32 i = blines * (Nz-1) ; i < blines * Nz ; i++)
			  for(register int32 j = 0; j < Npixelsml; j++)
				output_layer(i-(blines * (Nz-1)),j) = real4(output_buffer(i,j)); // real8 --> real4
			//refdemheiofile << output_layer;            // output reference dem heights
			(mlookedIFG == true) ? refdemheiofile << multilook(output_layer, ifgmlL, ifgmlP) : refdemheiofile << output_layer;
			//refdemheiofile << output_buffer(window((Nz-1) * blines,Nz * blines - 1, 0, Npixelsml -1 ));
		  }
		if (outputh2ph ==true)
		  {
			Nz += 1;
			for (register int32 i = blines * (Nz-1) ; i < blines * Nz ; i++)
			  for(register int32 j = 0; j < Npixelsml; j++)
				output_layer(i-(blines * (Nz-1)),j) = real4(output_buffer(i,j)); // real8 --> real4
			//h2phofile  << output_layer;            // output h2ph matrix
			(mlookedIFG == true) ? h2phofile << multilook(output_layer, ifgmlL, ifgmlP) : h2phofile << output_layer;
			//h2phofile << output_buffer(window((Nz-1) * blines,Nz * blines - 1, 0, Npixelsml -1 ));
		  }

		DEMline_buffer.resize(1, 1); // deallocate
		DEMpixel_buffer.resize(1, 1);
		input_buffer.resize(1, 1);
		temp_input_buffer.resize(1, 1);
		output_buffer.resize(1, 1);

	  } // end loop azimuth direction

	  INFO << "Closing output files";
	  INFO.print();

	  refdemofile.close(); // has the same multilook as interferrogram
	  if (mlookedIFG ==true)
		refdemofilenoML.close(); // [MA] if interferogram is mlooked then this is
									// generated as non-multilooked for coherence estimation.
	  if (outputrefdemhei ==true) // For Zbigniew Perski
		refdemheiofile.close();
	  if (outputh2ph ==true)
		h2phofile.close();

	  //===================================================================
	  //============ End second loop: interpolation           =============
	  //============ (radar geometry)                         =============
	  //===================================================================

	  // === Clean up temporary outfiles === [MA]

	  if (remove("master_demline.temp2")) // remove files
		WARNING.print("code 101: could not remove file: master_demline.temp2.");
	  if (remove("master_dempixel.temp2"))
		WARNING.print("code 101: could not remove file: master_dempixel.temp2.");
	  if (remove("dem_refphase.temp"))
		WARNING.print("code 101: could not remove file: dem_refphase.temp.");
	  if (outputh2ph ==true && remove("dem_h2ph.temp")) // output available
		WARNING.print("code 101: could not remove file: dem_h2ph.temp.");

	  // === Clean up Done. ===

	  // ====== Write output information ======
	  String croppeddemi = new String(new char[ONE27]);
	  croppeddemi = "NO output requested";
	  if (outputdemi)
		  croppeddemi = refdeminput.fodemi;
	  INFO << "Min. value of input DEM covering interferogram: " << min_input_dem;
	  INFO.print();
	  INFO << "Max. value of input DEM covering interferogram: " << max_input_dem;
	  INFO.print();

	  ofstream scratchlogfile = new ofstream("scratchlogcomprefdem", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchlogfile, "comprefdem: scratchlogcomprefdem", __FILE__, __LINE__);
	  scratchlogfile << "\n*******************************************************************" << "\n* " << processcontrol[pr_i_comprefdem] << "\n*******************************************************************" << "\n1) DEM source file:                   \t" << refdeminput.firefdem << "\nFormat:                               \t";
		switch (refdeminput.iformatflag)
		  {
		  case FORMATI2:
			{
			scratchlogfile << "SHORT SIGNED INTEGER (HOST ENDIANNESS)";
			break;
			}
		  case FORMATI2_BIGENDIAN:
			{
			scratchlogfile << "SHORT SIGNED INTEGER, BIG ENDIAN";
			break;
			}
		  case FORMATR4:
			{
			scratchlogfile << "REAL4 SIGNED FLOAT";
			break;
			}
		  case FORMATR8:
			{
			scratchlogfile << "REAL8 SIGNED DOUBLE";
			break;
			}
		  default:
			{
			scratchlogfile << "UNKNOWN? IMPOSSIBLE...";
			break;
			}
		  }
	//    << "\nMean value:                           \t" <<  meancroppedDEM
	// this is not correct, only stats per buffer...
	//    << "\n\n----- Other STATS -----"
	//    << "\nTotal points in cropped DEM:          \t" << numpoints
	//    << "\nNumber of valid points in DEM:        \t" << numvalid
	//    << " (" << 100*numvalid/numpoints << "%)"
	//    << "\nNumber of NODATA points in DEM:       \t" << numNODATA
	//    << " (" << 100*numNODATA/numpoints << "%)"
	//    << "\nMean height in meters at valid points:\t" << meancroppedDEM
	  scratchlogfile << "\nByte order:                           \t" << "check yourself..." << "\nNumber of lines:                      \t" << numberoflatpixels << "\nNumber of pixels:                     \t" << numberoflonpixels << "\nResolution latitude:                  \t" << rad2deg(DEMdeltalat) << " [deg]" << "\nResolution longitude:                 \t" << rad2deg(DEMdeltalon) << " [deg]" << "\nMost West point in input DEM:         \t" << rad2deg(lon0file) << "\nMost East point in input DEM:         \t" << rad2deg(lonNfile) << "\nMost South point in input DEM:        \t" << rad2deg(latNfile) << "\nMost North point in input DEM:        \t" << rad2deg(lat0file) << "\nMin. value of input DEM covering interferogram: " << min_input_dem << "\nMax. value of input DEM covering interferogram: " << max_input_dem << "\n2) Output file cropped DEM:           \t" << refdeminput.fodem << "\nFormat:                               \t" << "REAL4" << "\nByte order:                           \t" << "(same as host)" << "\n3) Output file interpolated crop DEM: \t" << croppeddemi << "\nFormat:                               \t" << "REAL4" << "\nByte order:                           \t" << "(same as host)" << "\n4) Output file synthetic phase:       \t" << refdeminput.forefdem << "\nFormat:                               \t" << "REAL4" << "\nByte order:                           \t" << "(same as host)" << "\nNumber of lines (multilooked):        \t" << ifgNlinesml << "\nNumber of pixels (multilooked):       \t" << ifgNpixelsml << "\n*******************************************************************\n\n";
	  scratchlogfile.close();


	  ofstream scratchresfile = new ofstream("scratchrescomprefdem", ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  bk_assert(scratchresfile, "comprefdem: scratchrescomprefdem", __FILE__, __LINE__);
	  scratchresfile << "\n\n*******************************************************************" << "\n*_Start_" << processcontrol[pr_i_comprefdem] << "\n*******************************************************************";
	  if (onlyrefphasetopo ==true)
		  scratchresfile << "\nInclude_flatearth:                 \tNo";
	  else
		  scratchresfile << "\nInclude_flatearth:                 \tYes";
	  scratchresfile << "\nDEM source file:                      \t" << refdeminput.firefdem << "\nMin. of input DEM:                    \t" << min_input_dem << "\nMax. of input DEM:                    \t" << max_input_dem << "\nData_output_file:                     \t" << refdeminput.forefdem << "\nData_output_format:                   \t" << "real4" << "\nFirst_line (w.r.t. original_master):  \t" << interferogram.win.linelo << "\nLast_line (w.r.t. original_master):   \t" << interferogram.win.linehi << "\nFirst_pixel (w.r.t. original_master): \t" << interferogram.win.pixlo << "\nLast_pixel (w.r.t. original_master):  \t" << interferogram.win.pixhi << "\nMultilookfactor_azimuth_direction:    \t" << ifgmlL << "\nMultilookfactor_range_direction:      \t" << ifgmlP << "\nNumber of lines (multilooked):        \t" << ifgNlinesml << "\nNumber of pixels (multilooked):       \t" << ifgNpixelsml << "\n*******************************************************************" << "\n* End_" << processcontrol[pr_i_comprefdem] << "_NORMAL" << "\n*******************************************************************\n";
	  scratchresfile.close();


	  } // END radarcodedem


	//***************************************************************
	// *    getcorners                                                *
	// *                                                              *
	// * Get corners of window (approx) to select DEM in radians (if  *
	// * height were zero)                                            *
	// *                                                              *
	// * Implementation:                                              *
	// * 1) calculate phi, lambda of corners                          *
	// * 2) select the extreme values                                 *
	// * 3) add extra overlap                                         *
	// * 4) determine the indices in the file                         *
	// *                                                              *
	// *    Freek van Leijen, 07-AUG-2006                             *
	// *                                                              *
	// ***************************************************************
	public static void getcorners(real8 l0, real8 lN, real8 p0, real8 pN, real8 extralat, real8 extralong, real8 lat0, real8 long0, real8 DEMdeltalat, real8 DEMdeltalong, int32 Nlatpixels, int32 Nlongpixels, input_ell ellips, slcimage master, RefObject<orbit> masterorbit, RefObject<real8> phimin, RefObject<real8> phimax, RefObject<real8> lambdamin, RefObject<real8> lambdamax, RefObject<int32> indexphi0DEM, RefObject<int32> indexphiNDEM, RefObject<int32> indexlambda0DEM, RefObject<int32> indexlambdaNDEM)
	  {
	  TRACE_FUNCTION("getcorners (FvL 07-AUG-2006)")

	  DEBUG << "l0 :" << l0;
	  DEBUG.print();
	  DEBUG << "lN :" << lN;
	  DEBUG.print();
	  DEBUG << "p0 :" << p0;
	  DEBUG.print();
	  DEBUG << "pN :" << pN;
	  DEBUG.print();
	  DEBUG << "extralat :" << extralat;
	  DEBUG.print();
	  DEBUG << "extralong :" << extralong;
	  DEBUG.print();
	  DEBUG << "lat0 :" << lat0;
	  DEBUG.print();
	  DEBUG << "long0 :" << long0;
	  DEBUG.print();
	  DEBUG << "DEMdeltalat :" << DEMdeltalat;
	  DEBUG.print();
	  DEBUG << "DEMdeltalong :" << DEMdeltalong;
	  DEBUG.print();
	  DEBUG << "Nlatpixels :" << Nlatpixels;
	  DEBUG.print();
	  DEBUG << "Nlongpixels :" << Nlongpixels;
	  DEBUG.print();

	  real8 phi;
	  real8 lambda;
	  real8 height;
	  lp2ell(l0,p0, ellips, master, masterorbit.argvalue, phi, lambda, height); // returned
	  real8 phil0p0 = phi;
	  real8 lambdal0p0 = lambda;

	  lp2ell(lN,p0, ellips, master, masterorbit.argvalue, phi, lambda, height); // returned
	  real8 philNp0 = phi;
	  real8 lambdalNp0 = lambda;

	  lp2ell(lN,pN, ellips, master, masterorbit.argvalue, phi, lambda, height); // returned
	  real8 philNpN = phi;
	  real8 lambdalNpN = lambda;

	  lp2ell(l0,pN, ellips, master, masterorbit.argvalue, phi, lambda, height); // returned
	  real8 phil0pN = phi;
	  real8 lambdal0pN = lambda;

	  // ______ Select DEM values based on rectangle outside l,p border ______
	  phimin.argvalue = min(min(min(phil0p0,philNp0),philNpN),phil0pN);
	  phimax.argvalue = max(max(max(phil0p0,philNp0),philNpN),phil0pN);
	  lambdamin.argvalue = min(min(min(lambdal0p0,lambdalNp0),lambdalNpN),lambdal0pN);
	  lambdamax.argvalue = max(max(max(lambdal0p0,lambdalNp0),lambdalNpN),lambdal0pN);

	  // ______ a little bit extra at edges to be sure ______ 
	  phimin.argvalue -= extralat;
	  phimax.argvalue += extralat;
	  lambdamax.argvalue += extralong;
	  lambdamin.argvalue -= extralong;

	  DEBUG << "phimin :" << phimin.argvalue;
	  DEBUG.print();
	  DEBUG << "phimax :" << phimax.argvalue;
	  DEBUG.print();
	  DEBUG << "lambdamin :" << lambdamin.argvalue;
	  DEBUG.print();
	  DEBUG << "lambdamax :" << lambdamax.argvalue;
	  DEBUG.print();

	  // ______ Get indices of DEM needed ______
	  // ______ Index boundary: [0:numberofx-1] ______

	  indexphi0DEM.argvalue = int32(Math.floor((lat0-phimax.argvalue)/DEMdeltalat));
	  if (indexphi0DEM.argvalue < 0)
		{
		WARNING << "indexphi0DEM: " << indexphi0DEM.argvalue;
		WARNING.print();
		indexphi0DEM.argvalue =0; // default start at first
		WARNING.print("DEM does not cover entire interferogram.");
		WARNING.print("input DEM should be extended to the North.");
		}
	  indexphiNDEM.argvalue = int32(Math.ceil((lat0-phimin.argvalue)/DEMdeltalat));
	  if (indexphiNDEM.argvalue > Nlatpixels-1)
		{
		  WARNING << "indexphiNDEM: " << indexphiNDEM.argvalue;
		  WARNING.print();
		indexphiNDEM.argvalue =Nlatpixels-1;
		WARNING.print("DEM does not cover entire interferogram.");
		WARNING.print("input DEM should be extended to the South.");
		}
	  indexlambda0DEM.argvalue = int32(Math.floor((lambdamin.argvalue-long0)/DEMdeltalong));
	  if (indexlambda0DEM.argvalue < 0)
		{
		  WARNING << "indexlambda0DEM: " << indexlambda0DEM.argvalue;
		  WARNING.print();
		indexlambda0DEM.argvalue =0; // default start at first
		WARNING.print("DEM does not cover entire interferogram.");
		WARNING.print("input DEM should be extended to the West.");
		}
	  indexlambdaNDEM.argvalue = int32(Math.ceil((lambdamax.argvalue-long0)/DEMdeltalong));
	  if (indexlambdaNDEM.argvalue > Nlongpixels-1)
		{
		  WARNING << "indexlambdaNDEM: " << indexlambdaNDEM.argvalue;
		  WARNING.print();
		indexlambdaNDEM.argvalue =Nlongpixels-1;
		WARNING.print("DEM does not cover entire interferogram.");
		WARNING.print("input DEM should be extended to the East.");
		}

		DEBUG << "indexphi0DEM :" << indexphi0DEM.argvalue;
		DEBUG.print();
		DEBUG << "indexphiNDEM :" << indexphiNDEM.argvalue;
		DEBUG.print();
		DEBUG << "indexlambda0DEM :" << indexlambda0DEM.argvalue;
		DEBUG.print();
		DEBUG << "indexlambdaNDEM :" << indexlambdaNDEM.argvalue;
		DEBUG.print();


	  } // END getcorners


	//***************************************************************
	// *    griddatalinear (naming after Matlab function)             *
	// *                                                              *
	// *    Implementation after GMT function triangulate.c           *
	// ***************************************************************
	  public static void griddatalinear(matrix<real8> x_in, matrix<real8> y_in, matrix<real8> z_in, real8 x_min, real8 x_max, real8 y_min, real8 y_max, int32 x_inc, int32 y_inc, real8 r_az_ratio, real8 offset, real8 NODATA, RefObject<matrix<real8>> grd)
	{
	  TRACE_FUNCTION("griddatalinear (LG&FvL 13-AUG-2006)")

	  INFO << "griddataLinear interpolation.";
	  INFO.print();

	  int32 i;
	  int32 j;
	  int32 k;
	  int32 ij;
	  int32 p;
	  int32 i_min;
	  int32 i_max;
	  int32 j_min;
	  int32 j_max;
	  int32 n;
	  int32 nx;
	  int32 ny;
	  int32 zLoops;
	  int32 zLoop;
	  int32 zBlockSize;
	  int32 indexFirstPoint;
	  int32 zInterpolateBlockSize;
	  real8[] vx = new real8[4];
	  real8[] vy = new real8[4];
	  real8 xkj;
	  real8 xlj;
	  real8 ykj;
	  real8 ylj;
	  real8 zj;
	  real8 zk;
	  real8 zl;
	  real8 zlj;
	  real8 zkj;
	  real8 xp;
	  real8 yp;
	  real8 f; // linear interpolation parameters
	  real8 a = null;
	  real8 b = null;
	  real8 c = null;
	  triangulateio In;
	  triangulateio Out;
	  triangulateio vorOut;


	  // Initialize variables

	  n = zBlockSize = x_in.size(); // block size of x and y coordination

	  // How many groups of z value should be interpolated 
	  if ((z_in.size() % zBlockSize) != 0)
		{
		  INFO << "The input of the DEM buffer and z is not the same...";
		  INFO.print();
		  return;
		}
	  else
		zLoops = z_in.size()/x_in.size();

	  a = new real8[zLoops];
	  b = new real8[zLoops];
	  c = new real8[zLoops];
	  if (a == null || b == null || c == null)
		{
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
		  ERROR << "Memory ERROR in source file: " << __FILE__ << " at line: " << __LINE__;
		  PRINT_ERROR(ERROR.get_str());
		  throw(memory_error);
		}

	  nx = grd.argvalue.lines()/zLoops;
	  ny = grd.argvalue.pixels();
	  zInterpolateBlockSize = grd.argvalue.size()/zLoops;

	  // Set everything to 0 and NULL 
//C++ TO JAVA CONVERTER TODO TASK: The memory management function 'memset' has no equivalent in Java:
	  memset ((Object) In, 0, sizeof (triangulateio));
//C++ TO JAVA CONVERTER TODO TASK: The memory management function 'memset' has no equivalent in Java:
	  memset ((Object) Out, 0, sizeof (triangulateio));
//C++ TO JAVA CONVERTER TODO TASK: The memory management function 'memset' has no equivalent in Java:
	  memset ((Object) vorOut, 0, sizeof (triangulateio));

	  // Allocate memory for input points 
	  In.numberofpoints = n;
	  In.pointlist = new real8 [2 * n];

	  // Copy x,y points to In structure array 

	  for (i = j = 0; i < n; i++)
		{
		  In.pointlist[j++] = *(x_in[0] + i);
		  In.pointlist[j++] = (*(y_in[0] + i)) * r_az_ratio;
		  // to eliminate the effect of difference in range and azimuth spacing;
		}

	//   Call Jonathan Shewchuk's triangulate algorithm.  This is 64-bit safe since
	//   * all the structures use 4-byte ints (longs are used internally). 

	  triangulate ("zIQB", In, Out, vorOut);

	  int32 link = Out.trianglelist; // List of node numbers to return via link
	  int32 np = Out.numberoftriangles;

	  for (k = ij = 0; k < np; k++)
		{
		  DEBUG << "k of np, ij: " << k << " of " << np << ", :" << ij;
		  DEBUG.print();
		  //Store the Index of the first Point of this triangle.
		  indexFirstPoint = ij;

		  vx[0] = vx[3] = *(x_in[0] + link[ij]);
		  vy[0] = vy[3] = *(y_in[0] + link[ij]);
		  ij++;
		  vx[1] = *(x_in[0] + link[ij]);
		  vy[1] = *(y_in[0]+link[ij]);
		  ij++;
		  vx[2] = *(x_in[0] + link[ij]);
		  vy[2] = *(y_in[0]+link[ij]);
		  ij++;

		  if (vx[0] == NODATA || vx[1] == NODATA || vx[2] == NODATA)
			  continue;
		  if (vy[0] == NODATA || vy[1] == NODATA || vy[2] == NODATA)
			  continue;

		  // Compute grid indices the current triangle may cover.
		  xp = min (min (vx[0], vx[1]), vx[2]);
		  i_min = x_to_i (xp, x_min, x_inc, offset, nx);
		  //INFO << "xp: " << xp;
		  //INFO.print();
		  xp = max (max (vx[0], vx[1]), vx[2]);
		  i_max = x_to_i (xp, x_min, x_inc, offset, nx);
		  //INFO << "xp: " << xp;
		  //INFO.print();
		  yp = min (min (vy[0], vy[1]), vy[2]);
		  j_min = y_to_j (yp, y_min, y_inc, offset, ny);
		  //INFO << "yp: " << yp;
		  //INFO.print();
		  yp = max (max (vy[0], vy[1]), vy[2]);
		  j_max = y_to_j (yp, y_min, y_inc, offset, ny);
		  //INFO << "yp: " << yp;
		  //INFO.print();
		  // Adjustments for triangles outside -R region. 
		  // Triangle to the left or right. 
		  if ((i_max < 0) || (i_min >= nx))
			  continue;
		  // Triangle Above or below 
		  if ((j_max < 0) || (j_min >= ny))
			  continue;
		  // Triangle covers boundary, left or right. 
		  if (i_min < 0)
			  i_min = 0;
		  if (i_max >= nx)
			  i_max = nx - 1;
		  // Triangle covers boundary, top or bottom. 
		  if (j_min < 0)
			  j_min = 0;
		  if (j_max >= ny)
			  j_max = ny - 1;
		  // for (kk = 0; kk<npar;kk++) {  //do for each parameter LIUG
		  // read zj, zk, zl (instead of above) LIUG
		  // Find equation for the plane as z = ax + by + c 
		  xkj = vx[1] - vx[0];
		  ykj = vy[1] - vy[0];
		  xlj = vx[2] - vx[0];
		  ylj = vy[2] - vy[0];

		  f = 1.0 / (xkj * ylj - ykj * xlj);

		  for(zLoop = 0 ; zLoop < zLoops; zLoop++)
	  {
		zj = *(z_in[0] + zLoop * zBlockSize + link[indexFirstPoint]);
		zk = *(z_in[0] + zLoop * zBlockSize + link[indexFirstPoint + 1]);
		zl = *(z_in[0] + zLoop * zBlockSize + link[indexFirstPoint + 2]);
		zkj = zk - zj;
		zlj = zl - zj;
		a[zLoop] = -f * (ykj * zlj - zkj * ylj);
		b[zLoop] = -f * (zkj * xlj - xkj * zlj);
		c[zLoop] = -a[zLoop] * vx[1] - b[zLoop] * vy[1] + zk;
	  }

		  for (i = i_min; i <= i_max; i++)
			{
		xp = i_to_x (i, x_min, x_max, x_inc, offset, nx);
		p = i * ny + j_min;
		for (j = j_min; j <= j_max; j++, p++)
				{
				  yp = j_to_y (j, y_min, y_max, y_inc, offset, ny);
				  if (pointintriangle(vx, vy, xp, yp) == 0) // Outside
					  continue;

			for(zLoop = 0 ; zLoop < zLoops; zLoop++)
		*(grd.argvalue[0] + zLoop * zInterpolateBlockSize + p) = a[zLoop] * xp + b[zLoop] * yp + c[zLoop];
				} //LIUG
			}
		}

	  if (a != null)
		  a = null;
	  if (b != null)
		  b = null;
	  if (c != null)
		  c = null;
	  if(In.pointlist) // only this use delete
		  In.pointlist = null;
	  if(Out.pointlist)
//C++ TO JAVA CONVERTER TODO TASK: The memory management function 'free' has no equivalent in Java:
		  free(Out.pointlist);
	  if(Out.trianglelist)
//C++ TO JAVA CONVERTER TODO TASK: The memory management function 'free' has no equivalent in Java:
		  free(Out.trianglelist);


	} //END griddatalinear


	//***************************************************************
	// *    pointintriangle                                           *
	// *                                                              *
	// *    Liu, Guang and Freek van Leijen, 16-AUG-2006              *
	// *                                                              *
	// *                                                              *
	// ***************************************************************
	public static int pointintriangle(RefObject<_double> xt, RefObject<_double> yt, _double x, _double y)
	{

	  int iRet0 = ((xt.argvalue[2] - xt.argvalue[0]) * (y - yt.argvalue[0])) > ((x - xt.argvalue[0]) * (yt.argvalue[2] - yt.argvalue[0])) ? 1:-1;
	  int iRet1 = ((xt.argvalue[0] - xt.argvalue[1]) * (y - yt.argvalue[1])) > ((x - xt.argvalue[1]) * (yt.argvalue[0] - yt.argvalue[1])) ? 1:-1;
	  int iRet2 = ((xt.argvalue[1] - xt.argvalue[2]) * (y - yt.argvalue[2])) > ((x - xt.argvalue[2]) * (yt.argvalue[1] - yt.argvalue[2])) ? 1:-1;

	  if ((iRet0 >0 && iRet1 > 0 && iRet2 > 0) || (iRet0 <0 && iRet1 < 0 && iRet2 < 0))
		return 1;
	  else
		return 0;

	} //END pointintriangle
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
// * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/referencephase.cc,v $ *
// * $Revision: 3.22 $                                            *
// * $Date: 2006/05/18 11:09:20 $                                 *
// * $Author: kampes $                                            *
// *                                                              *
// * -computation flat earth correction.                          *
// * -computation radarcoding dem + interpolation to (1,1) grid   *
// ***************************************************************


//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
//#if WIN32
  // Jia defined this.
  // Bert Kampes, 24-Aug-2005
//C++ TO JAVA CONVERTER WARNING: The following #include directive was ignored:
//  #include "winsock2.h"
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
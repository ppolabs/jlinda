package org.jdoris.core.todo_classes;

/**
 * User: pmar@ppolabs.com
 * Date: 4/8/11
 * Time: 5:45 PM
 */
public class todo_classes {

    public class productinfo {
/*        public:
          char          file[EIGHTY];                   // current filename
          // ______ window / multilook factors ______
          window        win;                            // current window, line(1:N) etc
          uint          multilookL;                     // multilookfactor in line (azi) dir.
          uint          multilookP;                     // multilookfactor in pixel (ra) dir.
          // ______ file format ______
          int16         formatflag;                     // current read formatflag


          // ______ Public function in struct ______
          // ______ constructor ______
          productinfo()
            {
            formatflag = -1;// undefined
            multilookL =  1;
            multilookP =  1;
            } // rest ==0

          // ______ fill it from info in resultfiles ______
          void fillproductinfo(const char *file, const char *iden);

          // ______ assignment operator ______
          productinfo& operator = (productinfo X)
            {
            if (this != &X)
              {
              strcpy(file,X.file);
              win        = X.win;
              multilookL = X.multilookL;
              multilookP = X.multilookP;
              formatflag = X.formatflag;
              }
            return *this;
            };

          // ______ show content ______
          inline void showdata() const                  // show content
            {DEBUG << "\ncurrent file: \t" << file
                   << "\nformatflag:   \t" << formatflag
                   << "\nmultilook:    \t" << multilookL << " " << multilookP
                   << "\nwindow:       \t" << win.linelo << " " << win.linehi
                                    << " " << win.pixlo  << " " << win.pixhi;
             DEBUG.print();
            }

          // ______ read data from file ______
          matrix<real4> readphase(window win) const;

          // ______ read data from file ______
          matrix<complr4> readdata(window win) const;
          matrix<real4> readdatar4(window win) const; // [MA]*/
    }

    public class inputgeneral {
/*
        char          logfile[4*ONE27];
        char          m_resfile[4*ONE27];
        char          s_resfile[4*ONE27];
        char          i_resfile[4*ONE27];
        uint          memory;                 // available mem. in Bytes
        bool          process[NUMPROCESSES];  // if .[i] != 0 => process step_(i+1)
        bool          interactive;            // if true, pause
        bool          overwrit;               // 0: don't overwrite existing data output files
        int16         orb_interp;             // method for orbit interpolation
        int32         dumpbaselineL;          // #lines to dump baseline param.
        int32         dumpbaselineP;          // #lines to dump baseline param.
        int32         preview;                // generate sunraster preview file
                                              // 0: no; 1: sh files; 2: sh sh_files.
        real4         terrain_height;         // mean terrain height, or of a point.
*/
    }

    public class input_ell {            // ellips a,b
/*
  {
  private:
  real8         e2;                     // squared first  eccentricity (derived)
  real8         e2b;                    // squared second eccentricity (derived)
  // ______ Helpers should be private ______
  inline void set_ecc1st_sqr()          // first ecc.
    {e2=1.0-sqr(b/a);}//  faster than e2=(sqr(a)-sqr(b))/sqr(a);
  inline void set_ecc2nd_sqr()          // second ecc.
    {e2b=sqr(a/b)-1.0;}// faster than e2b=(sqr(a)-sqr(b))/sqr(b);

  public:
  real8         a;                      // semi major
  real8         b;                      // semi minor
  char          name[EIGHTY];
  inline void set_name(const char *s)    {strcpy(name,s);}
  // ______ Default constructor ______
  input_ell()
    {a   = WGS84_A;
     b   = WGS84_B;
     e2  = 0.00669438003551279091;
     e2b = 0.00673949678826153145;
     //set_ecc1st_sqr();// compute e2
     //set_ecc2nd_sqr();// compute e2b
     set_name("WGS84");
    }
  // ______ constructor ellips(a,b) ______
  input_ell(const real8 &semimajor, const real8 &semiminor)
    {a = semimajor;
     b = semiminor;
     set_ecc1st_sqr();// compute e2 (not required for zero-doppler iter.)
     set_ecc2nd_sqr();// compute e2b (not required for zero-doppler iter.)
     //set_name("unknown");
    }
  // ______ Copy constructor ______
  input_ell(const input_ell& ell)
    {a=ell.a; b=ell.b; e2=ell.e2; e2b=ell.e2b; strcpy(name,ell.name);}
  // ______ Destructor ______
  ~input_ell()
    {;}// nothing to destruct that isn't destructed automatically
  // ______ Public function in struct ______
  inline input_ell& operator = (const input_ell &ell)// assignment operator
    {
    if (this != &ell)
      {a=ell.a; b=ell.b; e2=ell.e2; e2b=ell.e2b; strcpy(name,ell.name);}
    return *this;
    }
  inline void showdata() const
    {
    INFO << "ELLIPSOID: \tEllipsoid used (orbit, output): " << name << ".";
    INFO.print();
    INFO << "ELLIPSOID: a   = " << setw(15) << setprecision(13) << a;
    INFO.print();
    INFO << "ELLIPSOID: b   = " << setw(15) << setprecision(13) << b;
    INFO.print();
    INFO << "ELLIPSOID: e2  = " << e2;
    INFO.print();
    INFO << "ELLIPSOID: e2' = " << e2b;
    INFO.print();
    INFO.reset();
    }

  // ______ Convert xyz cartesian coordinates to ______
  // ______ Geodetic ellipsoid coordinates latlonh ______
  */
/****************************************************************
 *    xyz2ell                                                   *
 *                                                              *
 * Converts geocentric cartesian coordinates in the XXXX        *
 *  reference frame to geodetic coordinates.                    *
 *  method of bowring see globale en locale geodetische systemen*
 * input:                                                       *
 *  - ellipsinfo, xyz, (phi,lam,hei)                            *
 * output:                                                      *
 *  - void (updated lam<-pi,pi>, phi<-pi,pi>, hei)              *
 *                                                              *
 *    Bert Kampes, 05-Jan-1999                                  *
 ****************************************************************//*

  inline void xyz2ell(const cn &xyz, real8 &phi, real8 &lambda, real8 &height) const
    {
    TRACE_FUNCTION("xyz2ell (BK 05-Jan-1999)");
    const real8 r    = sqrt(sqr(xyz.x)+sqr(xyz.y));
    const real8 nu   = atan2((xyz.z*a),(r*b));
    const real8 sin3 = pow(sin(nu),3);
    const real8 cos3 = pow(cos(nu),3);
    phi              = atan2((xyz.z+e2b*b*sin3),(r-e2*a*cos3));
    lambda           = atan2(xyz.y,xyz.x);
    const real8 N    = a / sqrt(1.0-e2*sqr(sin(phi)));
    height           = (r/cos(phi)) - N;
    } // END xyz2ell
  // --- Same but without height ---
  inline void xyz2ell(const cn &xyz, real8 &phi, real8 &lambda) const
    {
    TRACE_FUNCTION("xyz2ell (BK 05-Jan-1999)");
    const real8 r    = sqrt(sqr(xyz.x)+sqr(xyz.y));
    const real8 nu   = atan2((xyz.z*a),(r*b));
    const real8 sin3 = pow(sin(nu),3);
    const real8 cos3 = pow(cos(nu),3);
    phi              = atan2((xyz.z+e2b*b*sin3),(r-e2*a*cos3));
    lambda           = atan2(xyz.y,xyz.x);
    } // END xyz2ell


  */
/****************************************************************
 *    ell2xyz                                                   *
 *                                                              *
 * Converts wgs84 ellipsoid cn to geocentric cartesian coord.   *
 * input:                                                       *
 *  - phi,lam,hei (geodetic co-latitude, longitude, [rad] h [m] *
 * output:                                                      *
 *  - cn XYZ                                                    *
 *                                                              *
 *    Bert Kampes, 05-Jan-1999                                  *
 ****************************************************************//*

  cn ell2xyz(const real8 &phi, const real8 &lambda, const real8 &height) const
    {
    const real8 N      = a / sqrt(1.0-e2*sqr(sin(phi)));
    const real8 Nph    = N + height;
    return cn(Nph * cos(phi) * cos(lambda),
              Nph * cos(phi) * sin(lambda),
             (Nph - e2*N)    * sin(phi));
    } // END ell2xyz


  // --- Test program for ellips class ----------------------------
  // --- This test is executed in inittest() --------------------------
  inline void test()
    {
    // constructors
    input_ell X;
    DEBUG << "ELLIPS: "  << X.name
          << ": X.a="  << X.a  << "; X.b="   << X.b
          << "; X.e2=" << X.e2 << "; X.e2b=" << X.e2b;
    DEBUG.print();
    }
*/
    }

    public class input_filtphase {
/*        // ______ ______
struct input_filtphase                  // arguments for phase filter
  {
  int16         method;                 // method selector
  char          fofiltphase[4*ONE27];   // output filename
  char          fifiltphase[4*ONE27];   // input filename
  uint          finumlines;             // number of lines input
  // ______ method goldstein ______
  real8         alpha;                  // weighting
  int32         blocksize;              // blocksize filtered blocks
  int32         overlap;                // half overlap
  // ______ method goldstein, spatial conv. and spectral ______
  matrix<real4> kernel;                 // e.g. [1 1 1]
  // ______ method spatial conv. and spectral ______
  char          fikernel2d[4*ONE27];    // input filename
  };*/
    }

    public class input_filtazi {
/*
        int16         method;                 // method selector
        int32         fftlength;              // length per buffer
        int32         overlap;                // 0.5overlap each buffer
        real8         hammingalpha;           // alpha for hamming, 1 is no
        char          foname[4*ONE27];        // output filename passed to routine
        char          fomaster[4*ONE27];      // output filename
        char          foslave[4*ONE27];       // output filename
        int16         oformatflag;            // output format [cr4] ci16
*/
    }

    public class input_filtrange {

        /*
              int16         method;                 // method selector
              int32         oversample;             // factor
              bool          doweightcorrel;         // weighting of correlation values
              int32         nlmean;                 // number of lines to take mean of
              int32         fftlength;              // length for adaptive
              int32         overlap;                // half overlap between blocks of fftlength
              real8         hammingalpha;           // alpha for hamming
              real8         SNRthreshold;           // spectral peak estimation
              real8         terrainslope;           // [rad] porbits method only
              char          fomaster[4*ONE27];      // output filename
              char          foslave[4*ONE27];       // output filename
              int16         oformatflag;            // output format [cr4] ci16
        */


    }

    public class input_comprefpha {                // arguments for flatearth correction.
        /*
          //char                idflatearth[EIGHTY];
          char          ifpositions[4*ONE27];   // input file name for positions
          int16         method;                 // method selector
          int32         degree;                 // degree of polynomial
          int32         Npoints;                // number of observations
        */
    }

    public class input_resample {                   // arguments for resampling slave
        /*
          int16         method;                 // method selector (interpolator) (%100 == Npoints)
          char          fileout[4*ONE27];
          int16         oformatflag;            // output format [cr4] ci16
          window        dbow_geo;               // cut out of original master.geo
          window        dbow;                   // cut out of original master.radar
          bool          shiftazi;               // [true] shift spectrum to 0
        */
    }

    public class input_interfero {                  // arguments for computation interferogram
        /*
        int16         method;                 // method selector
        char          focint[4*ONE27];                // optional output filename complex interferogram.
        char          foint[4*ONE27];         //  ~ of interferogram (phase).
      //  char        foflatearth[EIGHTY];    //  ~ of correction (flatearth) model (phase)
                                              //  these are flags as well as arguments.
                                              //  one is man (else no output)
        uint          multilookL;             // multilookfactor in line dir.
        uint          multilookP;             // multilookfactor in pixel dir.
        */
    }

    public class input_coherence {                  // arguments for computation coherence
        /*
          int16         method;                 // method selector
          char          focoh[4*ONE27];         // opt output filename of real coherence image.
          char          foccoh[4*ONE27];                //  ~ of complex coherence image.
                                                //  these are flags as well as arguments.
          uint          multilookL;             // multilookfactor in line dir.
          uint          multilookP;             // multilookfactor in pixel dir.
          uint          cohsizeL;               // size of estimation window coherence
          uint          cohsizeP;               // size of estimation window coherence
        */
    }

    public class input_subtrrefpha {               // arguments for subtract 'flat earth'
        /*
          int16         method;                 // method selector
          uint          multilookL;             // multilookfactor in line dir.
          uint          multilookP;             // multilookfactor in pixel dir.
          char          focint[4*ONE27];                // output filename complex interferogram
          char          forefpha[4*ONE27];      // output filename complex refpha
          char          foh2ph[4*ONE27];                // output filename h2ph, added by FvL
          bool          dumponlyrefpha;         // do nothing except dump refpha
        */
    }

    public class input_comprefdem {                // arguments for reference phase from DEM

        /*
        //  int16       method;                 // method selector
          char          firefdem[4*ONE27];      // input filename reference dem
          int16         iformatflag;            // input format [signed short]
          uint          demrows;                // number of
          uint          demcols;                // number of
          real8         demdeltalat;            // radians
          real8         demdeltalon;            // radians
          real8         demlatleftupper;        // radians
          real8         demlonleftupper;        // radians
          real8         demnodata;              // identifier/flag
        //  real8               extradense;             // extra interpolation factor (4)
          char          forefdem[4*ONE27];      // output filename reference phase
          char          foh2ph[4*ONE27];                // output perp. baseline, added by FvL
          char          forefdemhei[4*ONE27];   // output filename DEM in radarcoord.
          bool          includerefpha;          // flag to include_flatearth correction
          char          fodem[4*ONE27];         // flag+name output of cropped dem
          char          fodemi[4*ONE27];                // flag+name output of interpolated dem
        */
    }

    public class input_subtrrefdem {               // arguments for subtract reference DEM
        /*
          int16         method;                 // method selector
          int32         offsetL;                // offset applied before subtraction
          int32         offsetP;                // offset applied before subtraction
          char          focint[4*ONE27];                // output filename complex interferogram
        */
    }

}

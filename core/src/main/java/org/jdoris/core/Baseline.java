package org.jdoris.core;

import org.apache.log4j.Logger;
import org.esa.nest.util.Constants;
import org.jblas.Decompose;
import org.jblas.DoubleMatrix;
import org.jblas.Solve;

import static org.jblas.MatrixFunctions.abs;
import static org.jdoris.core.MathUtilities.normalize;

/**
 * User: pmar@ppolabs.com
 * Date: 3/18/11
 * Time: 9:28 PM
 */

// BASELINE is a new type (class) that is initialized using orbits
// and then either models the baseline parameters such as Bperp
// or can give them exact.
// Usage in the programs is something like in main do:
// BASELINE baseline;
// baseline.init(orbit1,orbit2,product?master?);
// and pass baseline to subprograms, there use baseline.get_bperp(x,y) etc.
// Bert Kampes, 31-Mar-2005
// For stability of normalmatrix, internally data are normalized
// using line/1024  pixel/1024  height/1024
// probably it is better to only do this in model_param
// but I also did it in eval_param because no time to check for correctness.
// Bert Kampes, 02-Aug-2005

public class Baseline {

    static Logger logger = Logger.getLogger(Baseline.class.getName());

    private boolean initialized;        //
    private double master_wavelength;   // tmp for now used for h_amb
    private double nearrange;           // range=nearrange+drange_dp*pixel
    private double drange_dp;           // range=nearrange+drange_dp*pixel
    private double orbit_convergence;   //    tmp for now constant
    private double orbit_heading;       //    tmp for now NOT USED
    private double L_min;               // for normalization
    private double L_max;               // for normalization
    private double P_min;               // for normalization
    private double P_max;               // for normalization
    private double H_min;               // height at which parameters are modeled
    private double H_max;               // to model phase=f(line,pix,hei) and hei=g(line,pix,phase)

    private long N_coeffs;              // ==10th degree of 3D-poly to model B(l,p,h)
    // --- B(l,p,h) = a000 +
    //                a100*l   + a010*p   + a001*h   +
    //                a110*l*p + a101*l*h + a011*p*h +
    //                a200*l^2 + a020*p^2 + a002*h^2

    // TODO: refactor to SmallDoubleMatrix class
    // --- Coefficients ---
    private DoubleMatrix BPERP_cf;          // perpendicular baseline
    private DoubleMatrix BPAR_cf;           // parallel baseline
    private DoubleMatrix THETA_cf;          // viewing angle
    private DoubleMatrix THETA_INC_cf;      // incidence angle to satellite
    //double avg_height_ambiguity;  //center height ambiguity

    // ellipsoid axes
    private static double ell_a = Constants.semiMajorAxis;
    private static double ell_b = Constants.semiMinorAxis;

    public Baseline() {

        logger.trace("Baseline class initialized");

        initialized = false;
        N_coeffs = 10;
        L_min = 0.0;        // for normalization
        L_max = 25000.0;    // for normalization
        P_min = 0.0;        // for normalization
        P_max = 5000.0;     // for normalization
        H_min = 0.0;        // height at which baseline is computed.
        H_max = 5000.0;     // height at which baseline is computed.
        master_wavelength = 0.0;
        orbit_convergence = 0.0; // tmp for now constant
        orbit_heading = 0.0;     // tmp for now NOT USED

    }


    /**************************************************************
     * --- B(l,p,h) = a000 +
     *                a100*l   + a010*p   + a001*h   +
     *                a110*l*p + a101*l*h + a011*p*h +
     *                a200*l^2 + a020*p^2 + a002*h^2
     **************************************************************/
    public double polyval(final DoubleMatrix C,
                          final double line,
                          final double pixel,
                          final double height) throws Exception {

        logger.debug("BASELINE CLASS::polyval");

        if (C.length != 10) {
            throw new Exception();
        } else {

            return
                    C.get(0, 0) +
                            C.get(1, 0) * line + C.get(2, 0) * pixel + C.get(3, 0) * height +
                            C.get(4, 0) * line * pixel + C.get(5, 0) * line * height + C.get(6, 0) * pixel * height +
                            C.get(7, 0) * sqr(line) + C.get(8, 0) * sqr(pixel) + C.get(9, 0) * sqr(height);

        }
    }

    private double sqr(double value) {
        return Math.pow(value, 2);
    }

    private double sqrt(double value) {
        return Math.sqrt(value);
    }

    /****************************************************************
     * Return baselineparameters                                    *
     ****************************************************************/
    public void BBparBperpTheta(double B, double Bpar, double Bperp, double theta,
                                final Point master, final Point point, final Point slave) {

        logger.debug("BASELINE CLASS::BBparBperpTheta");
        B = master.distance(slave); // baseline. abs. value (in plane master,point,slave)
        final double range1 = master.distance(point);
        final double range2 = slave.distance(point);
        Bpar = range1 - range2;  // parallel baseline, sign ok
        final Point r1 = master.min(point);// points from P to M
        final Point r2 = slave.min(point);
        theta = master.angle(r1);// viewing angle
        Bperp = sqr(B) - sqr(Bpar);
        if (Bperp < 0.0) {
            Bperp = 0.0;
        } else if (theta > master.angle(r2)) { // perpendicular baseline, sign ok
            Bperp = Math.sqrt(Bperp);
        } else {
            Bperp = -Math.sqrt(Bperp);
        }
    }

    /****************************************************************
     * returns incidence angle in radians based on coordinate of    *
     * point P on ellips and point M in orbit                       *
     ****************************************************************/
    public double IncidenceAngle(final Point master, final Point point) {

        logger.debug("BASELINE CLASS: IncidenceAngle");
        final Point r1 = master.min(point);// points from P to M
        final double inc = point.angle(r1);// incidence angle (assume P on ellips)

        return inc;

    } // END BBparBperpTheta

    private void model_parameters(final SLCImage master, final SLCImage slave, Orbit masterorbit, Orbit slaveorbit) {

        // TODO: with throwables!
        if (masterorbit.is_initialized() == false) {
            logger.debug("Baseline cannot be computed, master orbit not initialized.");
            return;
        } else if (slaveorbit.is_initialized() == false) {
            logger.debug("Baseline cannot be computed, slave orbit not initialized.");
            return;
        }

        // --- Get on with it ---------------------------------------
        if (initialized == true) {
            logger.warn("baseline already initialized??? (returning)");
            return;
        }
        initialized = true;//
        master_wavelength = master.getRadarWavelength();//

        // ______ Model r=nearrange+drange_dp*p, p starts at 1 ______
        nearrange = master.pix2range(1.0);
        drange_dp = master.pix2range(2.0) - master.pix2range(1.0);
        nearrange -= drange_dp;// (p starts at 1)

        // ______ Set min/max for normalization ______
        L_min = master.currentWindow.linelo;// also used during polyval
        L_max = master.currentWindow.linehi;// also used during polyval
        P_min = master.currentWindow.pixlo;// also used during polyval
        P_max = master.currentWindow.pixhi;// also used during polyval
        H_min = 0.0;// also used during polyval
        H_max = 5000.0;// also used during polyval

        // ______ Loop counters ______
        int cnt = 0;// matrix index
        final int N_pointsL = 10;// every 10km in azimuth
        final int N_pointsP = 10;// every 10km ground range
        final int N_heights = 4;// one more level than required for poly
        final double deltapixels = master.currentWindow.pixels() / N_pointsP;
        final double deltalines = master.currentWindow.lines() / N_pointsL;
        final double deltaheight = (H_max - H_min) / N_heights;

        // ______ Matrices for modeling Bperp (BK 21-mar-01) ______
        // --- For stability of normalmatrix, fill AMATRIX with normalized line, etc.
        DoubleMatrix BPERP = new DoubleMatrix(N_pointsL * N_pointsP * N_heights, 1);// perpendicular baseline
        DoubleMatrix BPAR = new DoubleMatrix(N_pointsL * N_pointsP * N_heights, 1);// parallel baseline
        DoubleMatrix THETA = new DoubleMatrix(N_pointsL * N_pointsP * N_heights, 1);// viewing angle
        DoubleMatrix THETA_INC = new DoubleMatrix(N_pointsL * N_pointsP * N_heights, 1);// inc. angle
        DoubleMatrix AMATRIX = new DoubleMatrix(N_pointsL * N_pointsP * N_heights, (int) N_coeffs);// design matrix


        // ______ Loop over heights, lines, pixels to compute baseline param. ______
        for (long k = 0; k < N_heights; ++k)  // height levels
        {
            final double HEIGHT = H_min + k * deltaheight;
//        input_ell ELLIPS(ellips.a + HEIGHT, ellips.b + HEIGHT);
            for (long i = 0; i < N_pointsL; ++i) // azimuthlines
            {
                final double line = master.currentWindow.linelo + i * deltalines;
                Point P;          // point, returned by lp2xyz
                double s_tazi;  // returned by xyz2t
                double s_trange;// returned by xyz2t
                final long MAXITER = 10;
                final double CRITERPOS = 1e-6;
                final double CRITERTIM = 1e-10;

                // ______ Azimuth time for this line ______
                final double m_tazi = master.line2ta(line);
                // ______ xyz for master satellite from time ______
                final Point M = masterorbit.getXYZ(m_tazi);
                // ______ Loop over a pixels to compute baseline param. ______
                for (long j = 0; j < N_pointsP; ++j) // rangepixels
                {
                    final double pixel = master.currentWindow.pixlo + j * deltapixels;

                    // ______ Range time for this pixel ______
                    //final double m_trange = master.pix2tr(pixel);
//                lp2xyz(line, pixel, ELLIPS, master, masterorbit, P, MAXITER, CRITERPOS);
                    // TODO: check this!
                    P = masterorbit.lp2xyz(line, pixel, master, masterorbit);

                    // ______ Compute xyz for slave satellite from P ______
//                xyz2t(s_tazi, s_trange, slave, slaveorbit, P, MAXITER, CRITERTIM);

                    // TODO: check and refactor this!
                    Point temp = slaveorbit.xyz2t(P, slave, slaveorbit);
                    s_tazi = temp.y;
                    s_trange = temp.x;

                    // ______ Slave position ______
                    final Point S = slaveorbit.getXYZ(s_tazi);

                    // ______ Compute angle between near parallel orbits ______
                    final Point Mdot = masterorbit.getXYZDot(m_tazi);
                    final Point Sdot = slaveorbit.getXYZDot(s_tazi);
                    final double angleorbits = Mdot.angle(Sdot);
                    logger.trace("Angle between orbits master-slave (at l,p= " + line + "," + pixel + ") = " +
                            rad2deg(angleorbits) + " [deg]");

                    orbit_convergence = angleorbits;// assume constant; store in member

                    //const real8 heading = angle(Mdot,[1 0 0])?
                    //orbit_heading = 0.0;// not yet used

                    // ====== The baseline parameters, derived from the positions (x,y,z) ======
                    // ______ alpha is angle counterclockwize(B, vlak met normal=rho1=rho2)
                    // ______ theta is angle counterclockwize(rho1=M, r1=M-P, r2=S-P)
                    double B = 0;
                    double Bpar = 0;
                    double Bperp = 0;
                    double theta = 0;
                    BBparBperpTheta(B, Bpar, Bperp, theta, M, P, S);// return B etc.
                    final double theta_inc = IncidenceAngle(M, P);// [rad]

                    // ______ Modelling of Bperp(l,p) = a00 + a10*l + a01*p ______
                    BPERP.put(cnt, 0, Bperp);
                    BPAR.put(cnt, 0, Bpar);
                    THETA.put(cnt, 0, theta);
                    THETA_INC.put(cnt, 0, theta_inc);

                    // --- B(l,p,h) = a000 +
                    //                a100*l   + a010*p   + a001*h   +
                    //                a110*l*p + a101*l*h + a011*p*h +
                    //                a200*l^2 + a020*p^2 + a002*h^2

                    AMATRIX.put(cnt, 0, 1.0);
                    AMATRIX.put(cnt, 1, normalize(line, L_min, L_max));

                    AMATRIX.put(cnt, 2, normalize(pixel, P_min, P_max));
                    AMATRIX.put(cnt, 3, MathUtilities.normalize(HEIGHT, H_min, H_max));
                    AMATRIX.put(cnt, 4, normalize(line, L_min, L_max) * normalize(pixel, P_min, P_max));
                    AMATRIX.put(cnt, 5, normalize(line, L_min, L_max) * normalize(HEIGHT, H_min, H_max));
                    AMATRIX.put(cnt, 6, normalize(pixel, P_min, P_max) * normalize(HEIGHT, H_min, H_max));
                    AMATRIX.put(cnt, 7, sqr(normalize(line, L_min, L_max)));
                    AMATRIX.put(cnt, 8, sqr(normalize(pixel, P_min, P_max)));
                    AMATRIX.put(cnt, 9, sqr(normalize(HEIGHT, H_min, H_max)));
                    cnt++;

                    // ______ B/alpha representation of baseline ______
                    final double alpha = (Bpar == 0 && Bperp == 0) ? Double.NaN : theta - Math.atan2(Bpar, Bperp);            // sign ok atan2

                    // ______ hor/vert representation of baseline ______
                    final double Bh = B * Math.cos(alpha);                        // sign ok
                    final double Bv = B * Math.sin(alpha);                        // sign ok
                    // ______ Height ambiguity: [h] = -lambda/4pi * (r1sin(theta)/Bperp) * phi==2pi ______
                    // TODO: check sign of infinity!!!
                    final double hambiguity = (Bperp == 0) ? Double.POSITIVE_INFINITY : -master.getRadarWavelength() * (M.min(P)).norm() * Math.sin(theta) / (2.0 * Bperp);


                    // ______ Some extra info if in debug mode ______

                    logger.debug("The baseline parameters for (l,p,h) = " + " line " + ", " + pixel + ", " + HEIGHT);
                    logger.debug("\talpha (deg), BASELINE: \t" + rad2deg(alpha) + " \t" + B);

                    logger.debug("\tBpar, Bperp:      \t" + Bpar + " \t" + Bperp);
                    logger.debug("\tBh, Bv:           \t" + Bh + " \t" + Bv);
                    logger.debug("\tHeight ambiguity: \t" + hambiguity);
                    logger.debug("\ttheta (deg):      \t" + rad2deg(theta));
                    logger.debug("\ttheta_inc (deg):  \t" + rad2deg(theta_inc));
                    logger.debug("\tM (x,y,z) = " + M.toString());
                    logger.debug("\tS (x,y,z) = " + S.toString());
                    logger.debug("\tP (x,y,z) = " + P.toString());
                } // loop pixels
            } // loop lines
        } // loop heights

//          // ====== Model the Bperp as 2d polynomial of degree 1 ======
        DoubleMatrix N = matTxmat(AMATRIX, AMATRIX);
//        DoubleMatrix rhsBperp = matTxmat(AMATRIX,BPERP);
//        DoubleMatrix rhsBpar  = matTxmat(AMATRIX,BPAR);
//        DoubleMatrix rhsT     = matTxmat(AMATRIX,THETA);
//        DoubleMatrix rhsT_INC = matTxmat(AMATRIX,THETA_INC);
//        DoubleMatrix Qx_hat   = N;
//
        // TODO: check this!
        final DoubleMatrix Qx_hat = Decompose.cholesky(N);
////        choles(Qx_hat);               // Cholesky factorisation normalmatrix

        final DoubleMatrix rhsBperp = Solve.solvePositive(AMATRIX, BPERP);
        final DoubleMatrix rhsBpar = Solve.solvePositive(AMATRIX, BPAR);
        final DoubleMatrix rhsT = Solve.solvePositive(AMATRIX, THETA);
        final DoubleMatrix rhsT_INC = Solve.solvePositive(AMATRIX, THETA_INC);

//        solvechol(Qx_hat,rhsBperp);   // Solution Bperp coefficients in rhsB
//        solvechol(Qx_hat,rhsBpar);    // Solution Theta coefficients in rhsT
//        solvechol(Qx_hat,rhsT);       // Solution Theta coefficients in rhsT
//        solvechol(Qx_hat,rhsT_INC);   // Solution Theta_inc coefficients in rhsT_INC
//        invertchol(Qx_hat);           // Covariance matrix of normalized unknowns

        final DoubleMatrix y_hatBperp = AMATRIX.mmul(rhsBperp);
        final DoubleMatrix e_hatBperp = BPERP.sub(y_hatBperp);

//      // ______Some other stuff, normalization is ok______
//      //matrix<real8> Qy_hat = AMATRIX * (matxmatT(Qx_hat,AMATRIX));
//      matrix<real8> y_hatBperp  = AMATRIX * rhsBperp;
//      matrix<real8> e_hatBperp  = BPERP - y_hatBperp;
//      //matrix<real8> Qe_hat  = Qy - Qy_hat;
//      //matrix<real8> y_hatT  = AMATRIX * rhsT;
//      //matrix<real8> e_hatT  = THETA - y_hatT;
//
        // === Copy estimated coefficients to private members ===
        BPERP_cf = rhsBperp;//
        BPAR_cf = rhsBpar;//
        THETA_cf = rhsT;//
        THETA_INC_cf = rhsT_INC;//

        // ______Test inverse______
        for (int i = 0; i < Qx_hat.rows; i++) {
            for (int j = 0; j < i; j++) {
                Qx_hat.put(j, i, Qx_hat.get(i, j)); // repair matrix
            }
        }
        final double maxdev = abs((N.mmul(Qx_hat).sub(DoubleMatrix.eye(Qx_hat.rows)))).max();

//        final double maxdev = max(abs(N*Qx_hat-eye(real8(Qx_hat.lines()))));

        logger.debug("BASELINE: max(abs(N*inv(N)-I)) = " + maxdev);

        if (maxdev > .01) {
            logger.warn("BASELINE: max. deviation N*inv(N) from unity = " + maxdev + ". This is larger than .01: do not use this!");
        } else if (maxdev > .001) {
            logger.warn("BASELINE: max. deviation N*inv(N) from unity = " + maxdev + ". This is between 0.01 and 0.001 (maybe not use it)");
        }


        // ______ Output solution and give max error ______
        // --- B(l,p,h) = a000 +
        //                a100*l   + a010*p   + a001*h   +
        //                a110*l*p + a101*l*h + a011*p*h +
        //                a200*l^2 + a020*p^2 + a002*h^2
        logger.debug("--------------------");
        logger.debug("Result of modeling: Bperp(l,p) = a000 + a100*l + a010*p + a001*h + ");
        logger.debug(" a110*l*p + a101*l*h + a011*p*h + a200*l^2 + a020*p^2 + a002*h^2");
        logger.debug("l,p,h in normalized coordinates [-2:2].");
        logger.debug("Bperp_a000 = " + rhsBperp.get(0, 0));
        logger.debug("Bperp_a100 = " + rhsBperp.get(1, 0));
        logger.debug("Bperp_a010 = " + rhsBperp.get(2, 0));
        logger.debug("Bperp_a001 = " + rhsBperp.get(3, 0));
        logger.debug("Bperp_a110 = " + rhsBperp.get(4, 0));
        logger.debug("Bperp_a101 = " + rhsBperp.get(5, 0));
        logger.debug("Bperp_a011 = " + rhsBperp.get(6, 0));
        logger.debug("Bperp_a200 = " + rhsBperp.get(7, 0));
        logger.debug("Bperp_a020 = " + rhsBperp.get(8, 0));
        logger.debug("Bperp_a002 = " + rhsBperp.get(9, 0));
        double maxerr = (abs(e_hatBperp)).max();
        if (maxerr > 2.00)//
        {
            logger.warn("Max. error bperp modeling at 3D datapoints: " + maxerr + "m");

        } else {
            logger.info("Max. error bperp modeling at 3D datapoints: " + maxerr + "m");
        }
        logger.debug("--------------------");
        logger.debug("Range: r(p) = r0 + dr*p");
        logger.debug("l and p in unnormalized, absolute, coordinates (1:N).");
        final double range1 = master.pix2range(1.0);
        final double range5000 = master.pix2range(5000.0);
        final double drange = (range5000 - range1) / 5000.0;
        logger.debug("range = " + (range1 - drange) + " + " + drange + "*p");

    } // END model_parameters()

// ___ Return RANGE to user ___

    public double get_range(final double pixel) {
        return nearrange + drange_dp * pixel;
    }
    //        inline real8 get_range(const real8 pixel) const
    //          {
    //          return nearrange + drange_dp*pixel;
    //          };// END get_bperp()

// === Polyval modeled quantities ===
// --- B(l,p,h) = a000 +
//                a100*l   + a010*p   + a001*h   +
//                a110*l*p + a101*l*h + a011*p*h +
//                a200*l^2 + a020*p^2 + a002*h^2
//
// l,p,h coefficients take normalized input
// Bert Kampes, 25-Aug-2005

    // ___ Return BPERP to user ___
    public double get_bperp(final double line, final double pixel, final double height) throws Exception {
        return polyval(BPERP_cf,
                normalize(line, L_min, L_max),
                normalize(pixel, P_min, P_max),
                normalize(height, H_min, H_max));
        //      return polyval(BPERP_cf,
        //        normalize(line,L_min,L_max),
        //        normalize(pixel,P_min,P_max),
        //        normalize(height,H_min,H_max));
    }

    public double get_bperp(final double line, final double pixel) throws Exception {
        return get_bperp(line, pixel, 0);
    }

    // ___ Return BPAR to user ___
    public double get_bpar(final double line, final double pixel, final double height) throws Exception {
        return polyval(BPAR_cf,
                normalize(line, L_min, L_max),
                normalize(pixel, P_min, P_max),
                normalize(height, H_min, H_max));
    }

    // ___ Return BPAR to user ___
    public double get_bpar(final double line, final double pixel) throws Exception {
        return get_bpar(line, pixel, 0);
    }

    // ___ Return THETA to user ___
    public double get_theta(final double line, final double pixel, final double height) throws Exception {
        return polyval(THETA_cf,
                normalize(line, L_min, L_max),
                normalize(pixel, P_min, P_max),
                normalize(height, H_min, H_max));
    }

    // ___ Return THETA_INC to user ___
    public double get_theta_inc(final double line, final double pixel, final double height) throws Exception {
        return polyval(THETA_INC_cf,
                normalize(line, L_min, L_max),
                normalize(pixel, P_min, P_max),
                normalize(height, H_min, H_max));
    }


    // === Derived quantities: do not normalize these!!! ===

    // ___ Return B to user ___
    public double get_b(final double line, final double pixel, final double height) throws Exception {
        return Math.sqrt(sqr(get_bpar(line, pixel, height)) + sqr(get_bperp(line, pixel, height)));
    }// END get_b()

    // ___ Return alpha baseline orientation to user ___
    public double get_alpha(final double line, final double pixel, final double height) throws Exception {
        final double Bperp = get_bperp(line, pixel, height);
        final double Bpar = get_bpar(line, pixel, height);
        final double theta = get_theta(line, pixel, height);
        final double alpha = (Bpar == 0 && Bperp == 0) ? Double.NaN : theta - Math.atan2(Bpar, Bperp);            // sign ok atan2
        return alpha;// sign ok
    }// END get_bhor()

    // ___ Return Bh to user ___
    public double get_bhor(final double line, final double pixel, final double height) throws Exception {
        final double B = get_b(line, pixel, height);
        final double alpha = get_alpha(line, pixel, height);
        return B * Math.cos(alpha);// sign ok
    }// END get_bhor()

    // ___ Return Bv to user ___
    public double get_bvert(final double line, final double pixel, final double height) throws Exception {
        final double B = get_b(line, pixel, height);
        final double alpha = get_alpha(line, pixel, height);
        return B * Math.sin(alpha);// sign ok
    }// END get_bvert()

    // ___ Return Height ambiguity to user ___
    public double get_hamb(final double line, final double pixel, final double height) throws Exception {
        //final double theta     =  get_theta(line,pixel,height);
        final double theta_inc = get_theta_inc(line, pixel, height);
        final double Bperp = get_bperp(line, pixel, height);
        final double range_MP = get_range(pixel);// >

        final double h_amb = (Bperp == 0) ? Double.POSITIVE_INFINITY : // inf
                -master_wavelength * range_MP * Math.sin(theta_inc) / (2.0 * Bperp);// this is wrt local
        //-master_wavelength*range_MP*sin(theta)/(2.0*Bperp);// this is wrt
        return h_amb;
    }// END get_hamb()

    // ___ Return orbit convergence to user ___
    //    public double get_orb_conv(final double line, final double pixel, final double height=0.0) const
    //      {
    //      // do not use l,p..
    //      return orbit_convergence;
    //      };// END get_orb_conv()

    // --- Dump overview of all ---
    void dump(final double line, final double pixel, final double height) throws Exception {
        if (initialized == false) {
            logger.debug("Exiting dumpbaseline, not initialized.");
            return;
        }
        // ______ Modeled quantities ______
        final double Bperp = get_bperp(line, pixel, height);
        final double Bpar = get_bpar(line, pixel, height);
        final double theta = get_theta(line, pixel, height);
        final double theta_inc = get_theta_inc(line, pixel, height);
        // ______ Derived quantities ______
        final double B = get_b(line, pixel, height);
        final double alpha = get_alpha(line, pixel, height);
        final double Bh = get_bhor(line, pixel, height);
        final double Bv = get_bvert(line, pixel, height);
        final double h_amb = get_hamb(line, pixel, height);
        // ______ Height ambiguity: [h] = -lambda/4pi * (r1sin(theta)/Bperp) * phi==2pi ______
        // ====== Write output to screen as INFO ======
        logger.info("The baseline parameters for (l,p,h) = "
                + line + ", " + pixel + ", " + height);

        logger.info("\tBpar, Bperp:      \t" + Bpar + " \t" + Bperp);

        logger.debug("\tB, alpha (deg):  \t" + B + " \t" + rad2deg(alpha));
        logger.debug("\tBh, Bv:          \t" + Bh + " \t" + Bv);
        logger.info("\tHeight ambiguity: \t" + h_amb);
        logger.info("\tLook angle (deg): \t" + rad2deg(theta));
        logger.debug("\tIncidence angle (deg): \t" + rad2deg(theta_inc));

    }// END dump()


    // HELPER FUNCTIONS!
    private DoubleMatrix matTxmat(DoubleMatrix matrix1, DoubleMatrix matrix2) {
        return matrix1.transpose().mmul(matrix2);
    }

    private double rad2deg(double rad) {
        return Math.toDegrees(rad);
    }

}
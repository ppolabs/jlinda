package org.jdoris.core;

import org.apache.log4j.Logger;
import org.esa.beam.framework.datamodel.MetadataElement;
import org.esa.nest.datamodel.AbstractMetadata;
import org.jblas.DoubleMatrix;
import org.jdoris.core.io.ResFile;

import java.io.File;


public final class Orbit {

    static Logger logger = Logger.getLogger(Orbit.class.getName());


    public static String interpMethod;
    public static int numStateVectors;
    public static double[] data_X;
    public static double[] data_Y;
    public static double[] data_Z;
    public static double[] coeff_X;
    public static double[] coeff_Y;
    public static double[] coeff_Z;
    public static int poly_degree;
    static int MAXITER = 10;

    // for SPLINE interpolator
    public static long klo; // index in timevector correct
    public static long khi; // +part picewize polynomial

    public static double[] time;
    static double CRITERPOS = Math.pow(10, -6);
    static double CRITERTIM = Math.pow(10, -10);
    static int refHeight = 0;

    // ellipsoid axes
    private double ell_a = Constants.WGS84_A;
    private double ell_b = Constants.WGS84_B;
    private double SOL = Constants.SOL;

    public Orbit() {
    }

    public Orbit(double[] timeVector, double[] xVector, double[] yVector, double[] zVector) {

        time = timeVector;
        data_X = xVector;
        data_Y = yVector;
        data_Z = zVector;

        numStateVectors = time.length;

    }

    public Orbit(double[][] stateVectors) {
        setOrbit(stateVectors);
    }

    public void parseOrbit(File file) throws Exception {
        ResFile resFile = new ResFile(file);
        setOrbit(resFile.parseOrbit());
    }

    public void setOrbit(double[][] stateVectors) {

        numStateVectors = stateVectors.length;

        for (int i = 0; i < stateVectors.length; i++) {
            time[i]   = stateVectors[i][0];
            data_X[i] = stateVectors[i][1];
            data_Y[i] = stateVectors[i][2];
            data_Z[i] = stateVectors[i][3];
        }
    }


    // TODO: refactor this one, split in definition and interpolation
    public Orbit(MetadataElement nestMetadataElement) {

        final AbstractMetadata.OrbitStateVector[] orbitStateVectors = AbstractMetadata.getOrbitStateVectors(nestMetadataElement);

        numStateVectors = orbitStateVectors.length;

//        time = new double[numStateVectors];
//        data_X = new double[numStateVectors];
//        data_Y = new double[numStateVectors];
//        data_Z = new double[numStateVectors];

        for (int i = 0; i < numStateVectors; i++) {
            // convert time to seconds of the acquisition day
            time[i] = (orbitStateVectors[i].time_mjd -
                    (int) orbitStateVectors[i].time_mjd) * (24 * 3600); // Modified Julian Day 2000 (MJD2000)
            data_X[i] = orbitStateVectors[i].x_pos;
            data_Y[i] = orbitStateVectors[i].y_pos;
            data_Z[i] = orbitStateVectors[i].z_pos;
        }

    }

    private void computeCoefficients(Orbit orbit, int degree) throws Exception {

        // TODO: switch on interpolation method, either spline method or degree of polynomial
        logger.info("Computing coefficients for orbit polyfit degree: " + interpMethod);
        // compute coefficients of orbit interpolator
        // final int getPolyDegree = 3; //HARDCODED because of how AbstractedMetadata Handles orbits
        // poly_degree = MathUtilities.getPolyDegree();
        poly_degree = degree;
        coeff_X = MathUtilities.polyFit(new DoubleMatrix(time), new DoubleMatrix(data_X), poly_degree);
        coeff_Y = MathUtilities.polyFit(new DoubleMatrix(time), new DoubleMatrix(data_Y), poly_degree);
        coeff_Z = MathUtilities.polyFit(new DoubleMatrix(time), new DoubleMatrix(data_Z), poly_degree);

    }

    public static int getNumStateVectors() {
        return numStateVectors;
    }

    // TODO: for splines!
    private void getKloKhi() {
    }

    // TODO: make generic so it can work with arrays of lines as well: see matlab implementation
    public Point lp2xyz(double line, double pixel, SLCImage slcimage) throws Exception {

        Point satellitePosition;
        Point satelliteVelocity;
        Point ellipsoidPosition;

        // TODO: check notations and difference between NEST and DORIS
        double azTime = slcimage.line2ta(line);
        double rgTime = slcimage.pix2tr(pixel);

        satellitePosition = getXYZ(azTime);
        satelliteVelocity = getXYZDot(azTime);

        // initial value
        ellipsoidPosition = slcimage.getApproxXYZCentreOriginal();

        // allocate matrices
        double[] equationSet = new double[3];
        double[][] partialsXYZ = new double[3][3];

        // iterate
        for (int iter = 0; iter <= MAXITER; iter++) {
            //   update equations and slove system
            double dsat_Px = ellipsoidPosition.x - satellitePosition.x;   // vector of 'satellite to P on ellipsoid'
            double dsat_Py = ellipsoidPosition.y - satellitePosition.y;   // vector of 'satellite to P on ellipsoid'
            double dsat_Pz = ellipsoidPosition.z - satellitePosition.z;   // vector of 'satellite to P on ellipsoid'

            equationSet[0] =
                    -(satelliteVelocity.x * dsat_Px +
                            satelliteVelocity.y * dsat_Py +
                            satelliteVelocity.z * dsat_Pz);

            equationSet[1] =
                    -(dsat_Px * dsat_Px +
                            dsat_Py * dsat_Py +
                            dsat_Pz * dsat_Pz - Math.pow(SOL * rgTime, 2));

            equationSet[2] =
                    -((ellipsoidPosition.x * ellipsoidPosition.x + ellipsoidPosition.y * ellipsoidPosition.y) / (Math.pow(ell_a + refHeight, 2)) +
                            Math.pow(ellipsoidPosition.z / (ell_b + refHeight), 2) - 1.0);

            partialsXYZ[0][0] = satelliteVelocity.x;
            partialsXYZ[0][1] = satelliteVelocity.y;
            partialsXYZ[0][2] = satelliteVelocity.z;
            partialsXYZ[1][0] = 2 * dsat_Px;
            partialsXYZ[1][1] = 2 * dsat_Py;
            partialsXYZ[1][2] = 2 * dsat_Pz;
            partialsXYZ[2][0] = (2 * ellipsoidPosition.x) / (Math.pow(ell_a + refHeight, 2));
            partialsXYZ[2][1] = (2 * ellipsoidPosition.y) / (Math.pow(ell_a + refHeight, 2));
            partialsXYZ[2][2] = (2 * ellipsoidPosition.z) / (Math.pow(ell_a + refHeight, 2));

            // solve system [NOTE!] orbit has to be normalized, otherwise close to singular
            // DoubleMatrix ellipsoidPositionSolution = Solve.solve(partialsXYZ, equationSet);
            double[] ellipsoidPositionSolution = MathUtilities.solve33(partialsXYZ, equationSet);

            // update solution
            ellipsoidPosition.x = ellipsoidPosition.x + ellipsoidPositionSolution[0];
            ellipsoidPosition.y = ellipsoidPosition.y + ellipsoidPositionSolution[1];
            ellipsoidPosition.z = ellipsoidPosition.z + ellipsoidPositionSolution[2];

            // check convergence
            if (Math.abs(ellipsoidPositionSolution[0]) < CRITERPOS &&
                    Math.abs(ellipsoidPositionSolution[1]) < CRITERPOS &&
                    Math.abs(ellipsoidPositionSolution[2]) < CRITERPOS) {
//                System.out.println("INFO: ellipsoidPosition (converged) = " + ellipsoidPosition);
                break;
            } else if (iter >= MAXITER) {
                MAXITER = MAXITER + 1;
                System.out.println("WARNING: line, pix -> x,y,z: maximum iterations (" + MAXITER + ") reached. " + "Criterium (m): " + CRITERPOS +
                        "dx,dy,dz=" + ellipsoidPositionSolution[0] + ", " + ellipsoidPositionSolution[1] + ", " + ellipsoidPositionSolution[2]);
            }
        }

        return ellipsoidPosition;
    }

    public Point lp2xyz(Point sarPixel, SLCImage slcimage) throws Exception {
        return lp2xyz(sarPixel.y, sarPixel.x, slcimage);
    }

    public Point xyz2orb(Point pointOnEllips, SLCImage slcimage) {

        // Initial value azimuth time
        Point posSat;
        double solution = 0.0;
        long iter;
        double tAzi = slcimage.line2ta(0.5 * slcimage.getCurrentWindow().linehi - slcimage.getCurrentWindow().linelo);

        for (iter = 0; iter < MAXITER; ++iter) {

            // update equations
            posSat = getXYZ(tAzi);
            Point velSat = getXYZDot(tAzi);
            Point accSat = getXYZDotDot(tAzi);
            Point delta = pointOnEllips.min(posSat);

            solution = -eq1_Doppler(velSat, delta) / eq1_Doppler_dt(delta, velSat, accSat);
            tAzi += solution;

            // ______ Check convergence ______
            if (Math.abs(solution) < CRITERTIM)                   // dta
                break;

        }

        // ______ Check for number of iterations ______
        if (iter >= MAXITER) {
            logger.warn("x,y,z -> line, pix: maximum iterations (" + MAXITER + ") reached. " +
                    "Criterium (s):" + CRITERTIM + "dta (s)=" + solution);
        }

        // Compute range time
        // ____ Update equations _____
        return posSat = getXYZ(tAzi);

    }

    public Point xyz2t(Point position, SLCImage slcimage) {

        Point delta;
        Point returnVector = new Point();

        // inital value
        double timeAzimuth = slcimage.line2ta(0.5 * slcimage.getApproxRadarCentreOriginal().y);

        int iter;
        double solution = 0;
        for (iter = 0; iter <= MAXITER; ++iter) {
            Point satellitePosition = getXYZ(timeAzimuth);
            Point satelliteVelocity = getXYZDot(timeAzimuth);
            Point satelliteAcceleration = getXYZDotDot(timeAzimuth);
            delta = position.min(satellitePosition);

            // update solution
            solution = -1 * (satelliteVelocity.x * delta.x + satelliteVelocity.y * delta.y + satelliteVelocity.z * delta.z) /
                    (satelliteAcceleration.x * delta.x + satelliteAcceleration.y * delta.y + satelliteAcceleration.z * delta.z -
                            Math.pow(satelliteVelocity.x, 2) - Math.pow(satelliteVelocity.y, 2) - Math.pow(satelliteVelocity.z, 2));

            timeAzimuth += solution;

            if (Math.abs(solution) < CRITERTIM) {
                break;
            }

        }
        // ______ Check number of iterations _____
        if (iter >= MAXITER) {
            logger.warn("x,y,z -> line, pix: maximum iterations (" + MAXITER + ") reached. " + "Criterium (s):" + CRITERTIM + "dta (s)=" + solution);
        }

        // ====== Compute range time ======
        // ______ Update equations ______
        final Point satellitePosition = getXYZ(timeAzimuth);

        delta = position.min(satellitePosition);

        double timeRange = Math.sqrt(Math.pow(delta.x, 2) + Math.pow(delta.y, 2) + Math.pow(delta.z, 2)) / SOL;

        returnVector.y = timeAzimuth;
        returnVector.x = timeRange;

        return returnVector;

    }

    public Point xyz2lp(Point position, SLCImage slcimage) {

        Point returnPixel = new Point();

        // ______ Compute tazi, tran ______
        Point time = xyz2t(position, slcimage);

        // ______ Convert time to pixel ______
        // ______ (Converged) Result is in returnline/pixel ______
        returnPixel.y = slcimage.ta2line(time.y);
        returnPixel.x = slcimage.tr2pix(time.x);

        return returnPixel;
    }

    public Point ell2lp(double[] phi_lam_height, SLCImage slcimage) throws Exception {
        Point xyz = Ellipsoid.ell2xyz(phi_lam_height);
        return xyz2lp(xyz, slcimage);
    }

    public double[] lp2ell(Point position, SLCImage slcimage) throws Exception {
        Point xyz = lp2xyz(position, slcimage);
        return Ellipsoid.xyz2ell(xyz);
    }

    // TODO: legacy support, implementation from baseline class
    @Deprecated
    public void computeBaseline() {

    }

    public Point getXYZ(double azTime) {

        //TODO: sanity check!
        Point satelliteXYZPosition = new Point();

        // normalize time
        azTime = (azTime - time[time.length / 2]) / 10;

        satelliteXYZPosition.x = MathUtilities.polyVal1d(azTime, coeff_X);
        satelliteXYZPosition.y = MathUtilities.polyVal1d(azTime, coeff_Y);
        satelliteXYZPosition.z = MathUtilities.polyVal1d(azTime, coeff_Z);

        return satelliteXYZPosition;  //To change body of created methods use File | Settings | File Templates.
    }

    public Point getXYZDot(double azTime) {

        //TODO: sanity check
        Point satelliteVelocity = new Point();

        // normalize time
        azTime = (azTime - time[time.length / 2]) / 10;

        // NOTE: orbit interpolator is simple polynomial
        satelliteVelocity.x = coeff_X[1];
        satelliteVelocity.y = coeff_Y[1];
        satelliteVelocity.z = coeff_Z[1];

        for (int i = 2; i <= poly_degree; ++i) {
            double powT = (double) i * Math.pow(azTime, (double) (i - 1));
            satelliteVelocity.x += coeff_X[i] * powT;
            satelliteVelocity.y += coeff_Y[i] * powT;
            satelliteVelocity.z += coeff_Z[i] * powT;
        }

        return satelliteVelocity.divByScalar(10.0d);

    }

    public Point getXYZDotDot(double azTime) {

        //TODO: sanity check
        Point satelliteAcceleration = new Point();

        // normalize time
        azTime = (azTime - time[time.length / 2]) / 10.0d;

        // NOTE: orbit interpolator is simple polynomial
        // 2a_2 + 2*3a_3*t^1 + 3*4a_4*t^2...

        for (int i = 2; i <= poly_degree; ++i) {
            double powT = (double) ((i - 1) * i) * Math.pow(azTime, (double) (i - 2));
            satelliteAcceleration.x += coeff_X[i] * powT;
            satelliteAcceleration.y += coeff_Y[i] * powT;
            satelliteAcceleration.z += coeff_Z[i] * powT;
        }
        return satelliteAcceleration.divByScalar(100.0d);

    }

    public double eq1_Doppler(Point velocity, Point position) {
        return velocity.in(position);
    }

    private double eq1_Doppler_dt(Point dsatP, Point velocity, Point accerelation) {
        return accerelation.in(dsatP) - Math.pow(velocity.x, 2) - Math.pow(velocity.y, 2) - Math.pow(velocity.z, 2);
    }

    public double eq2_Range(Point dSatP, double rangeTime) {
        return dSatP.in(dSatP) - Math.pow(SOL * rangeTime, 2);
    }

    public double eq3_Ellipsoid(Point P) {
        return ((Math.pow(P.x, 2) + Math.pow(P.y, 2)) / Math.pow(ell_a, 2)) +
                Math.pow(P.z / ell_b, 2) - 1.0;
    }

    public double eq3_Ellipsoid(Point P, double semimajor_a, double semiminor_b) {
        return ((Math.pow(P.x, 2) + Math.pow(P.y, 2)) / Math.pow(semimajor_a, 2)) +
                Math.pow(P.z / semiminor_b, 2) - 1.0;
    }

    // TODO: implement check on initialization of orbit class: JAVA should have better facility for this
    public boolean is_initialized() {
        return true;
    }

    public Point[][] dumpOrbit() {

        if (numStateVectors == 0) {
            System.out.println("Exiting Orbit.dumporbit(), no orbit data available.");
        }

        double dt = 0.;

        logger.info("dumporbits: MAXITER: " + MAXITER + "; " + "\n" +
                    "          CRITERPOS: " + CRITERPOS + " m; " + "\n" +
                    "          CRITERTIM: " + CRITERTIM + " s");

        //  ______ Evaluate polynomial orbit for t1:dt:tN ______
        int outputlines = 1 + (int) ((time[numStateVectors - 1] - time[0]) / dt);

        double tAzi = time[0];
        Point[][] dumpedOrbit = new Point[outputlines][3];

        for (int i = 0; i < outputlines; ++i) {
            dumpedOrbit[i][0] = getXYZ(tAzi);
            dumpedOrbit[i][1] = getXYZDot(tAzi);
            dumpedOrbit[i][2] = getXYZDotDot(tAzi);
            tAzi += dt;
        }


        // ______ dump coeff. as well for testing ... ______
        //  #ifdef __DEBUG
        //  if (ID==MASTERID)
        //    {
        //    DEBUG.print("dumping files m_t, m_x, m_y, m_z, m_cx, m_cy, m_cz for spline interpolation.");
        //    dumpasc("m_t",time);
        //    dumpasc("m_x",data_x);  dumpasc("m_y",data_y);  dumpasc("m_z",data_z);
        //    dumpasc("m_cx",coef_x); dumpasc("m_cy",coef_y); dumpasc("m_cz",coef_z);
        //    }
        //  else
        //    {
        //    DEBUG.print("dumping files s_t, s_x, s_y, s_z, s_cx, s_cy, s_cz for spline interpolation.");
        //    dumpasc("s_t",time);
        //    dumpasc("s_x",data_x);  dumpasc("s_y",data_y);  dumpasc("s_z",data_z);
        //    dumpasc("s_cx",coef_x); dumpasc("s_cy",coef_y); dumpasc("s_cz",coef_z);
        //    }
        //  #endif

        return dumpedOrbit;

    }

    public void showOrbit() {
        // TODO: refactor to new matrix class!
        logger.info("Time of orbit ephemerides: " + time.toString());
        logger.info("Orbit ephemerides x:" + data_X.toString());
        logger.info("Orbit ephemerides y:" + data_Y.toString());
        logger.info("Orbit ephemerides z:" + data_Z.toString());
        logger.info("Estimated coefficients x(t):" + coeff_X.toString());
        logger.info("Estimated coefficients y(t):" + coeff_Y.toString());
        logger.info("Estimated coefficients z(t):" + coeff_Z.toString());
    }


}
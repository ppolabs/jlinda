package org.jdoris.core;

import org.apache.log4j.Logger;
import org.esa.beam.framework.datamodel.MetadataElement;
import org.esa.nest.datamodel.AbstractMetadata;
import org.jblas.DoubleMatrix;
import org.jdoris.core.io.ResFile;

import java.io.File;


public final class Orbit {

    static Logger logger = Logger.getLogger(Orbit.class.getName());

    private String interpMethod;

    private boolean isInterpolated = false;
    private boolean isPropagated = false;

    private int numStateVectors;

    private double[] time;
    private double[] data_X;
    private double[] data_Y;
    private double[] data_Z;
    private double[] coeff_X;
    private double[] coeff_Y;
    private double[] coeff_Z;
    private int poly_degree;

    private int MAXITER = 10;
    private final static double CRITERPOS = Math.pow(10, -6);
    private final static double CRITERTIM = Math.pow(10, -10);
    private final static int refHeight = 0;

//    // for SPLINE interpolator
//    private long klo; // index in timevector correct
//    private long khi; // +part picewize polynomial

    // ellipsoid axes
    private final static double ell_a = Constants.WGS84_A;
    private final static double ell_b = Constants.WGS84_B;
    private final static double SOL = Constants.SOL;

    public Orbit() {
    }

    public Orbit(double[] timeVector, double[] xVector, double[] yVector, double[] zVector, int degree) throws Exception {

        numStateVectors = time.length;

        // state vectors
        time = timeVector;
        data_X = xVector;
        data_Y = yVector;
        data_Z = zVector;

        // polynomial coefficients
        poly_degree = degree;
        computeCoefficients();
    }

    public Orbit(double[][] stateVectors, int degree) throws Exception {
        setOrbit(stateVectors);

        this.poly_degree = degree;
        computeCoefficients();
    }

    public void parseOrbit(File file) throws Exception {
        ResFile resFile = new ResFile(file);
        setOrbit(resFile.parseOrbit());
    }

    // TODO: refactor this one, split in definition and interpolation
    public Orbit(MetadataElement nestMetadataElement, int degree) throws Exception {

        final AbstractMetadata.OrbitStateVector[] orbitStateVectors = AbstractMetadata.getOrbitStateVectors(nestMetadataElement);

        numStateVectors = orbitStateVectors.length;

        for (int i = 0; i < numStateVectors; i++) {
            // convert time to seconds of the acquisition day
            time[i] = (orbitStateVectors[i].time_mjd -
                    (int) orbitStateVectors[i].time_mjd) * (24 * 3600); // Modified Julian Day 2000 (MJD2000)
            data_X[i] = orbitStateVectors[i].x_pos;
            data_Y[i] = orbitStateVectors[i].y_pos;
            data_Z[i] = orbitStateVectors[i].z_pos;
        }

        poly_degree = degree;
        computeCoefficients();
    }

    public void setOrbit(double[][] stateVectors) {

        numStateVectors = stateVectors.length;

        time = new double[numStateVectors];
        data_X = new double[numStateVectors];
        data_Y = new double[numStateVectors];
        data_Z = new double[numStateVectors];

        for (int i = 0; i < stateVectors.length; i++) {
            time[i] = stateVectors[i][0];
            data_X[i] = stateVectors[i][1];
            data_Y[i] = stateVectors[i][2];
            data_Z[i] = stateVectors[i][3];
        }
    }

    // TODO: switch on interpolation method, either spline method or degree of polynomial
    private void computeCoefficients() throws Exception {

        logger.info("Computing coefficients for orbit polyfit degree: " + poly_degree);
//        poly_degree = degree;
        this.coeff_X = MathUtilities.polyFit(new DoubleMatrix(time), new DoubleMatrix(data_X), poly_degree);
        this.coeff_Y = MathUtilities.polyFit(new DoubleMatrix(time), new DoubleMatrix(data_Y), poly_degree);
        this.coeff_Z = MathUtilities.polyFit(new DoubleMatrix(time), new DoubleMatrix(data_Z), poly_degree);

        isInterpolated = true;

    }

    public void computeCoefficients(int degree) throws Exception {
        poly_degree = degree;
        computeCoefficients();
    }


    // TODO: for splines!
    private void getKloKhi() {
    }

    // TODO: make generic so it can work with arrays of lines as well: see matlab implementation
    public Point lph2xyz(double line, double pixel, double height, SLCImage slcimage) throws Exception {

        Point satellitePosition;
        Point satelliteVelocity;
        Point ellipsoidPosition; // returned

        // allocate matrices
        double[] equationSet = new double[3];
        double[][] partialsXYZ = new double[3][3];

        double azTime = slcimage.line2ta(line);
        double rgTime = slcimage.pix2tr(pixel);

        satellitePosition = getXYZ(azTime);
        satelliteVelocity = getXYZDot(azTime);

        // initial value
        ellipsoidPosition = slcimage.getApproxXYZCentreOriginal();

        // iterate for the solution
        for (int iter = 0; iter <= MAXITER; iter++) {

            // update equations and solve system

            Point dsat_P = ellipsoidPosition.min(satellitePosition);

            equationSet[0] = -eq1_Doppler(satelliteVelocity, dsat_P);
            equationSet[1] = -eq2_Range(dsat_P, rgTime);
            equationSet[2] = -eq3_Ellipsoid(ellipsoidPosition, height);

            partialsXYZ[0][0] = satelliteVelocity.x;
            partialsXYZ[0][1] = satelliteVelocity.y;
            partialsXYZ[0][2] = satelliteVelocity.z;
            partialsXYZ[1][0] = 2 * dsat_P.x;
            partialsXYZ[1][1] = 2 * dsat_P.y;
            partialsXYZ[1][2] = 2 * dsat_P.z;
            partialsXYZ[2][0] = (2 * ellipsoidPosition.x) / (Math.pow(ell_a + height, 2));
            partialsXYZ[2][1] = (2 * ellipsoidPosition.y) / (Math.pow(ell_a + height, 2));
            partialsXYZ[2][2] = (2 * ellipsoidPosition.z) / (Math.pow(ell_a + height, 2));

            // solve system [NOTE!] orbit has to be normalized, otherwise close to singular
            // DoubleMatrix ellipsoidPositionSolution = Solve.solve(partialsXYZ, equationSet);
//            double[] ellipsoidPositionSolution = Solve.solve(new DoubleMatrix(partialsXYZ), new DoubleMatrix(equationSet)).toArray();
            double[] ellipsoidPositionSolution = MathUtilities.solve33(partialsXYZ, equationSet);

            // update solution
            ellipsoidPosition.x += ellipsoidPositionSolution[0];
            ellipsoidPosition.y += ellipsoidPositionSolution[1];
            ellipsoidPosition.z += ellipsoidPositionSolution[2];

            logger.debug("ellipsoidPosition.x = " + ellipsoidPosition.x);
            logger.debug("ellipsoidPosition.y = " + ellipsoidPosition.y);
            logger.debug("ellipsoidPosition.z = " + ellipsoidPosition.z);

            // check convergence
            if (Math.abs(ellipsoidPositionSolution[0]) < CRITERPOS &&
                    Math.abs(ellipsoidPositionSolution[1]) < CRITERPOS &&
                    Math.abs(ellipsoidPositionSolution[2]) < CRITERPOS) {
                logger.info("INFO: ellipsoidPosition (converged) = " + ellipsoidPosition);
                break;

            } else if (iter >= MAXITER) {
                MAXITER = MAXITER + 1;
                logger.warn("line, pix -> x,y,z: maximum iterations (" + MAXITER + ") reached. ");
                logger.warn("Criterium (m): " + CRITERPOS + " dx,dy,dz = "
                        + ellipsoidPositionSolution[0] + ", " + ellipsoidPositionSolution[1] + ", " + ellipsoidPositionSolution[2]);

                if (MAXITER > 11) {
                    logger.error("lp2xyz : MAXITER limit reached! lp2xyz() estimation is diverging?!");
                    throw new Exception("Orbit.lp2xyz : MAXITER limit reached! lp2xyz() estimation is diverging?!");
                }

            }
        }

        return ellipsoidPosition;
    }

    public Point lp2xyz(Point sarPixel, SLCImage slcimage) throws Exception {
        return lph2xyz(sarPixel.y, sarPixel.x, 0, slcimage);
    }

    public Point lp2xyz(double line, double pixel, SLCImage slcimage) throws Exception {
        return lph2xyz(line, pixel, 0, slcimage);
    }

    public Point xyz2orb(Point pointOnEllips, SLCImage slcimage) {
        // return satellite position
        // Point pointTime = xyz2t(pointOnEllips,slcimage);
        return getXYZ(xyz2t(pointOnEllips,slcimage).y); // inlined
    }

    public Point xyz2t(Point pointOnEllips, SLCImage slcimage) {

        Point delta;

        // inital value
        double timeAzimuth = slcimage.line2ta(0.5 * slcimage.getApproxRadarCentreOriginal().y);

        int iter;
        double solution = 0;
        for (iter = 0; iter <= MAXITER; ++iter) {
            Point satellitePosition = getXYZ(timeAzimuth);
            Point satelliteVelocity = getXYZDot(timeAzimuth);
            Point satelliteAcceleration = getXYZDotDot(timeAzimuth);
            delta = pointOnEllips.min(satellitePosition);

            // update solution
            solution = -eq1_Doppler(satelliteVelocity, delta) / eq1_Doppler_dt(delta, satelliteVelocity, satelliteAcceleration);
            timeAzimuth += solution;

            if (Math.abs(solution) < CRITERTIM) {
                break;
            }
        }

        // Check number of iterations
        if (iter >= MAXITER) {
            logger.warn("x,y,z -> line, pix: maximum iterations (" + MAXITER + ") reached. " + "Criterium (s):" + CRITERTIM + "dta (s)=" + solution);
        }

        // Compute range time

        // Update equations
        Point satellitePosition = getXYZ(timeAzimuth);
        delta = pointOnEllips.min(satellitePosition);
        double timeRange = delta.norm() / SOL;

        return new Point(timeRange, timeAzimuth);
    }


    public Point xyz2lp(Point pointOnEllips, SLCImage slcimage) {

        Point sarPixel = new Point();

        // Compute tazi, tran
        Point time = xyz2t(pointOnEllips, slcimage);

        // Convert time to pixel
        sarPixel.x = slcimage.tr2pix(time.x);
        sarPixel.y = slcimage.ta2line(time.y);

        return sarPixel;
    }

    public Point ell2lp(double[] phi_lam_height, SLCImage slcimage) throws Exception {
        Point pointOnEllips = Ellipsoid.ell2xyz(phi_lam_height);
        return xyz2lp(pointOnEllips, slcimage);
    }

    public double[] lp2ell(Point sarPixel, SLCImage slcimage) throws Exception {
        Point xyz = lp2xyz(sarPixel, slcimage);
        return Ellipsoid.xyz2ell(xyz);
    }

    // TODO: legacy support, implementation from baseline class
    @Deprecated
    public void computeBaseline() {

    }

    public Point getXYZ(double azTime) {

        if (azTime < time[0] || azTime > time[numStateVectors - 1]) {
            logger.warn("getXYZ() interpolation at: " + azTime + " is outside interval time axis: ("
                + time[0] + ", " + time[numStateVectors-1] + ").");
        }

        Point satelliteXYZPosition = new Point();

        // normalize time
        azTime = (azTime - time[time.length / 2]) / 10;

        satelliteXYZPosition.x = MathUtilities.polyVal1d(azTime, coeff_X);
        satelliteXYZPosition.y = MathUtilities.polyVal1d(azTime, coeff_Y);
        satelliteXYZPosition.z = MathUtilities.polyVal1d(azTime, coeff_Z);

        return satelliteXYZPosition;  //To change body of created methods use File | Settings | File Templates.
    }

    public Point getXYZDot(double azTime) {

        if (azTime < time[0] || azTime > time[numStateVectors - 1]) {
            logger.warn("getXYZDot() interpolation at: " + azTime + " is outside interval time axis: ("
                + time[0] + ", " + time[numStateVectors-1] + ").");
        }

        Point satelliteVelocity = new Point();

        //TODO: spline support!

        // normalize time
        azTime = (azTime - time[numStateVectors / 2]) / 10;

        int DEGREE = coeff_X.length - 1;

        satelliteVelocity.x = coeff_X[1];
        satelliteVelocity.y = coeff_Y[1];
        satelliteVelocity.z = coeff_Z[1];

        for (int i = 2; i <= DEGREE; ++i) {
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

    public double eq1_Doppler(Point satVelocity, Point pointOnEllips) {
        return satVelocity.in(pointOnEllips);
    }

    private double eq1_Doppler_dt(Point pointEllipsSat, Point satVelocity, Point satAcceleration) {
        return satAcceleration.in(pointEllipsSat) - Math.pow(satVelocity.x, 2) - Math.pow(satVelocity.y, 2) - Math.pow(satVelocity.z, 2);
    }

    public double eq2_Range(Point pointEllipsSat, double rgTime) {
        return pointEllipsSat.in(pointEllipsSat) - Math.pow(SOL * rgTime, 2);
    }

    public double eq3_Ellipsoid(Point pointOnEllips, double height) {
        return ((Math.pow(pointOnEllips.x, 2) + Math.pow(pointOnEllips.y, 2)) / Math.pow(ell_a + height, 2)) +
                Math.pow(pointOnEllips.z / (ell_b + height), 2) - 1.0;
    }

    public double eq3_Ellipsoid(Point pointOnEllips) {
        return eq3_Ellipsoid(pointOnEllips, 0);
    }

    public double eq3_Ellipsoid(Point pointOnEllips, double semiMajorA, double semiMinorB, double height) {
        return ((Math.pow(pointOnEllips.x, 2) + Math.pow(pointOnEllips.y, 2)) / Math.pow(semiMajorA+height, 2)) +
                Math.pow(pointOnEllips.z / (semiMinorB+height), 2) - 1.0;
    }

    // TODO: sanity checks
    public Point[][] dumpOrbit() {

        if (numStateVectors == 0) {
            System.out.println("Exiting Orbit.dumporbit(), no orbit data available.");
        }

        final double dt = 0.25;

        logger.info("dumporbits: MAXITER: " + MAXITER + "; " +
                "CRITERPOS: " + CRITERPOS + "m; " +
                "CRITERTIM: " + CRITERTIM + "s");

        //  ______ Evaluate polynomial orbit for t1:dt:tN ______
        int outputLines = 1 + (int) ((time[numStateVectors - 1] - time[0]) / dt);

        double azTime = time[0];
        Point[][] dumpedOrbit = new Point[outputLines][3];

        for (int i = 0; i < outputLines; ++i) {
            dumpedOrbit[i][0] = getXYZ(azTime);
            dumpedOrbit[i][1] = getXYZDot(azTime);
            dumpedOrbit[i][2] = getXYZDotDot(azTime);
            azTime += dt;
        }

        return dumpedOrbit;

    }

    public void showOrbit() {

        logger.info("Time of orbit ephemerides: " + time.toString());
        logger.info("Orbit ephemerides x:" + data_X.toString());
        logger.info("Orbit ephemerides y:" + data_Y.toString());
        logger.info("Orbit ephemerides z:" + data_Z.toString());

        if (isInterpolated) {
            logger.info("Estimated coefficients x(t):" + coeff_X.toString());
            logger.info("Estimated coefficients y(t):" + coeff_Y.toString());
            logger.info("Estimated coefficients z(t):" + coeff_Z.toString());
        }
    }

    public int getNumStateVectors() {
        return numStateVectors;
    }

    public double[] getTime() {
        return time;
    }

    public double[] getData_X() {
        return data_X;
    }

    public double[] getData_Y() {
        return data_Y;
    }

    public double[] getData_Z() {
        return data_Z;
    }

    public double[] getCoeff_X() {
        return coeff_X;
    }

    public double[] getCoeff_Y() {
        return coeff_Y;
    }

    public double[] getCoeff_Z() {
        return coeff_Z;
    }

    public int getPoly_degree() {
        return poly_degree;
    }

    public void setPoly_degree(int degree) {
        poly_degree = degree;
    }

    public boolean isInterpolated() {
        return isInterpolated;
    }

}
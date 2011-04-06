package org.jdoris.core;

import org.esa.beam.framework.datamodel.GeoPos;
import org.esa.beam.framework.datamodel.MetadataElement;
import org.esa.nest.datamodel.AbstractMetadata;
import org.esa.nest.util.Constants;
import org.esa.nest.util.GeoUtils;
import org.jblas.DoubleMatrix;
import org.jblas.Solve;


public final class Orbit {

    public static String interpMethod;
    public static int numStateVectors;
    public static DoubleMatrix data_X;
    public static DoubleMatrix data_Y;
    public static DoubleMatrix data_Z;
    public static DoubleMatrix coeff_X;
    public static DoubleMatrix coeff_Y;
    public static DoubleMatrix coeff_Z;
    public static int poly_degree;
    static int MAXITER = 10;

    // for SPLINE interpolator
    public static long klo; // index in timevector correct
    public static long khi; // +part picewize polynomial

    public static DoubleMatrix time;
    static double CRITERPOS = Math.pow(10, -6);
    static double CRITERTIM = Math.pow(10, -10);
    static double SOL = Constants.lightSpeed;
    static int refHeight = 0;

    // ellipsoid axes
    private static double ell_a = Constants.semiMajorAxis;
    private static double ell_b = Constants.semiMinorAxis;

    public Orbit() {
    }

    public Orbit(DoubleMatrix timeVector, DoubleMatrix xVector, DoubleMatrix yVector, DoubleMatrix zVector) {

        time = timeVector;
        data_X = xVector;
        data_Y = yVector;
        data_Z = zVector;

        numStateVectors = time.rows;

    }

    // TODO: refactor this one, split in definition and interpolation
    public Orbit(MetadataElement nestMetadataElement) {

        final AbstractMetadata.OrbitStateVector[] orbitStateVectors = AbstractMetadata.getOrbitStateVectors(nestMetadataElement);

        numStateVectors = orbitStateVectors.length;

        time = new DoubleMatrix(numStateVectors);
        data_X = new DoubleMatrix(numStateVectors);
        data_Y = new DoubleMatrix(numStateVectors);
        data_Z = new DoubleMatrix(numStateVectors);

        for (int i = 0; i < numStateVectors; i++) {
            // convert time to seconds of the acquisition day
            time.put(i, (orbitStateVectors[i].time_mjd -
                    (int) orbitStateVectors[i].time_mjd) * (24 * 3600)); // Modified Julian Day 2000 (MJD2000)
            data_X.put(i, orbitStateVectors[i].x_pos);
            data_Y.put(i, orbitStateVectors[i].y_pos);
            data_Z.put(i, orbitStateVectors[i].z_pos);
        }

    }

    private void computeCoefficients(Orbit orbit, int degree) {

        // TODO: switch on interpolation method, either spline method or degree of polynomial
        //      INFO << "Computing coefficients for orbit polyfit degree: "
        //           << interp_method;
        // compute coefficients of orbit interpolator
        // final int getPolyDegree = 3; //HARDCODED because of how AbstractedMetadata Handles orbits
        // poly_degree = MathUtilities.getPolyDegree();
        poly_degree = degree;
        coeff_X = MathUtilities.polyFit(time, data_X, poly_degree);
        coeff_Y = MathUtilities.polyFit(time, data_Y, poly_degree);
        coeff_Z = MathUtilities.polyFit(time, data_Z, poly_degree);

    }

    public static int getNumStateVectors() {
        return numStateVectors;
    }

    // TODO:
    private void getKloKhi() {
    }

    // TODO: make generic so it can work with arrays of lines as well: see matlab implementation
    // TODO: switch to Cn class
    public Point lp2xyz(double line, double pixel, SLCImage slcimage, Orbit orbits) {

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
        DoubleMatrix equationSet = DoubleMatrix.zeros(3);
        DoubleMatrix partialsXYZ = DoubleMatrix.zeros(3, 3);

        // iterate
        for (int iter = 0; iter <= MAXITER; iter++) {
            //   update equations and slove system
            double dsat_Px = ellipsoidPosition.x - satellitePosition.x;   // vector of 'satellite to P on ellipsoid'
            double dsat_Py = ellipsoidPosition.y - satellitePosition.y;   // vector of 'satellite to P on ellipsoid'
            double dsat_Pz = ellipsoidPosition.z - satellitePosition.z;   // vector of 'satellite to P on ellipsoid'

            equationSet.put(0,
                    -(satelliteVelocity.x * dsat_Px +
                            satelliteVelocity.y * dsat_Py +
                            satelliteVelocity.z * dsat_Pz));

            equationSet.put(1,
                    -(dsat_Px * dsat_Px +
                            dsat_Py * dsat_Py +
                            dsat_Pz * dsat_Pz - Math.pow(SOL * rgTime, 2)));

            equationSet.put(2,
                    -((ellipsoidPosition.x * ellipsoidPosition.x + ellipsoidPosition.y * ellipsoidPosition.y) / (Math.pow(ell_a + refHeight, 2)) +
                            Math.pow(ellipsoidPosition.z / (ell_b + refHeight), 2) - 1.0));

            partialsXYZ.put(0, 0, satelliteVelocity.x);
            partialsXYZ.put(0, 1, satelliteVelocity.y);
            partialsXYZ.put(0, 2, satelliteVelocity.z);
            partialsXYZ.put(1, 0, 2 * dsat_Px);
            partialsXYZ.put(1, 1, 2 * dsat_Py);
            partialsXYZ.put(1, 2, 2 * dsat_Pz);
            partialsXYZ.put(2, 0, (2 * ellipsoidPosition.x) / (Math.pow(ell_a + refHeight, 2)));
            partialsXYZ.put(2, 1, (2 * ellipsoidPosition.y) / (Math.pow(ell_a + refHeight, 2)));
            partialsXYZ.put(2, 2, (2 * ellipsoidPosition.z) / (Math.pow(ell_a + refHeight, 2)));

            // solve system [NOTE!] orbit has to be normalized, otherwise close to singular
            DoubleMatrix ellipsoidPositionSolution = Solve.solve(partialsXYZ, equationSet);
            // DoubleMatrix ellipsoidPositionSolution = solve33(partialsXYZ, equationSet);

            // update solution
            ellipsoidPosition.x = ellipsoidPosition.x + ellipsoidPositionSolution.get(0);
            ellipsoidPosition.y = ellipsoidPosition.y + ellipsoidPositionSolution.get(1);
            ellipsoidPosition.z = ellipsoidPosition.z + ellipsoidPositionSolution.get(2);

            // check convergence
            if (Math.abs(ellipsoidPositionSolution.get(0)) < CRITERPOS &&
                    Math.abs(ellipsoidPositionSolution.get(1)) < CRITERPOS &&
                    Math.abs(ellipsoidPositionSolution.get(2)) < CRITERPOS) {
//                System.out.println("INFO: ellipsoidPosition (converged) = " + ellipsoidPosition);
                break;
            } else if (iter >= MAXITER) {
                MAXITER = MAXITER + 1;
                System.out.println("WARNING: line, pix -> x,y,z: maximum iterations (" + MAXITER + ") reached. " + "Criterium (m): " + CRITERPOS +
                        "dx,dy,dz=" + ellipsoidPositionSolution.get(0) + ", " + ellipsoidPositionSolution.get(1) + ", " + ellipsoidPositionSolution.get(2, 0));
            }
        }

        return ellipsoidPosition;
    }

    public Point lp2xyz(Point sarPixel, SLCImage slcImage, Orbit orbit) {

        Point satellitePosition;
        Point satelliteVelocity;
        Point ellipsoidPosition = new Point();

        // TODO: check notations and difference between NEST and DORIS
        double azTime = slcImage.line2ta(sarPixel.x);
        double rgTime = slcImage.pix2tr(sarPixel.y);

        satellitePosition = getXYZ(azTime);
        satelliteVelocity = getXYZDot(azTime);

        // initial value
//        ellipsoidPosition = slcImage.getApproxXYZCentreOriginal();
//        ellipsoidPosition = slcImage.getApproxXYZCentreOriginal();

        // TODO: switch to ARRAYs
        // allocate matrices
        DoubleMatrix equationSet = DoubleMatrix.zeros(3);
        DoubleMatrix partialsXYZ = DoubleMatrix.zeros(3, 3);

        // iterate
        for (int iter = 0; iter <= MAXITER; iter++) {
            //   update equations and slove system

            double dsat_Px = ellipsoidPosition.x - satellitePosition.x;   // vector of 'satellite to P on ellipsoid'
            double dsat_Py = ellipsoidPosition.y - satellitePosition.y;   // vector of 'satellite to P on ellipsoid'
            double dsat_Pz = ellipsoidPosition.z - satellitePosition.z;   // vector of 'satellite to P on ellipsoid'

            equationSet.put(0,
                    -(satelliteVelocity.x * dsat_Px +
                            satelliteVelocity.y * dsat_Py +
                            satelliteVelocity.z * dsat_Pz));

            equationSet.put(1,
                    -(dsat_Px * dsat_Px +
                            dsat_Py * dsat_Py +
                            dsat_Pz * dsat_Pz - Math.pow(SOL * rgTime, 2)));

            equationSet.put(2,
                    -((ellipsoidPosition.x * ellipsoidPosition.x + ellipsoidPosition.y * ellipsoidPosition.y) / (Math.pow(ell_a + refHeight, 2)) +
                            Math.pow(ellipsoidPosition.z / (ell_b + refHeight), 2) - 1.0));

            partialsXYZ.put(0, 0, satelliteVelocity.x);
            partialsXYZ.put(0, 1, satelliteVelocity.y);
            partialsXYZ.put(0, 2, satelliteVelocity.z);
            partialsXYZ.put(1, 0, 2 * dsat_Px);
            partialsXYZ.put(1, 1, 2 * dsat_Py);
            partialsXYZ.put(1, 2, 2 * dsat_Pz);
            partialsXYZ.put(2, 0, (2 * ellipsoidPosition.x) / (Math.pow(ell_a + refHeight, 2)));
            partialsXYZ.put(2, 1, (2 * ellipsoidPosition.y) / (Math.pow(ell_a + refHeight, 2)));
            partialsXYZ.put(2, 2, (2 * ellipsoidPosition.z) / (Math.pow(ell_a + refHeight, 2)));

            // solve system [NOTE!] orbit has to be normalized, otherwise close to singular
            DoubleMatrix ellipsoidPositionSolution = Solve.solve(partialsXYZ, equationSet); // use JBlas call
            // DoubleMatrix ellipsoidPositionSolution = solve33(partialsXYZ, equationSet);

            // update solution
            ellipsoidPosition.x = ellipsoidPosition.x + ellipsoidPositionSolution.get(0);
            ellipsoidPosition.y = ellipsoidPosition.y + ellipsoidPositionSolution.get(1);
            ellipsoidPosition.z = ellipsoidPosition.z + ellipsoidPositionSolution.get(2);

            // check convergence
            if (Math.abs(ellipsoidPositionSolution.get(0)) < CRITERPOS &&
                    Math.abs(ellipsoidPositionSolution.get(1)) < CRITERPOS &&
                    Math.abs(ellipsoidPositionSolution.get(2)) < CRITERPOS) {
//                System.out.println("INFO: ellipsoidPosition (converged) = " + ellipsoidPosition);
                break;
            } else if (iter >= MAXITER) {
                MAXITER = MAXITER + 1;
                System.out.println("WARNING: line, pix -> x,y,z: maximum iterations (" + MAXITER + ") reached. " + "Criterium (m): " + CRITERPOS +
                        "dx,dy,dz=" + ellipsoidPositionSolution.get(0) + ", " + ellipsoidPositionSolution.get(1) + ", " + ellipsoidPositionSolution.get(2, 0));
            }
        }

        return ellipsoidPosition;
    }

    public Point lp2xyz(Point sarPixel, SLCImage slcImage) {

        Point satellitePosition;
        Point satelliteVelocity;
        Point ellipsoidPosition = new Point();

        // TODO: check notations and difference between NEST and DORIS
        double azTime = slcImage.line2ta(sarPixel.x);
        double rgTime = slcImage.pix2tr(sarPixel.y);

        satellitePosition = getXYZ(azTime);
        satelliteVelocity = getXYZDot(azTime);

        // initial value
//        ellipsoidPosition = slcImage.getApproxXYZCentreOriginal();
//        ellipsoidPosition = slcImage.getApproxXYZCentreOriginal();

        // TODO: switch to ARRAYs
        // allocate matrices
        DoubleMatrix equationSet = DoubleMatrix.zeros(3);
        DoubleMatrix partialsXYZ = DoubleMatrix.zeros(3, 3);

        // iterate
        for (int iter = 0; iter <= MAXITER; iter++) {
            //   update equations and slove system

            double dsat_Px = ellipsoidPosition.x - satellitePosition.x;   // vector of 'satellite to P on ellipsoid'
            double dsat_Py = ellipsoidPosition.y - satellitePosition.y;   // vector of 'satellite to P on ellipsoid'
            double dsat_Pz = ellipsoidPosition.z - satellitePosition.z;   // vector of 'satellite to P on ellipsoid'

            equationSet.put(0,
                    -(satelliteVelocity.x * dsat_Px +
                            satelliteVelocity.y * dsat_Py +
                            satelliteVelocity.z * dsat_Pz));

            equationSet.put(1,
                    -(dsat_Px * dsat_Px +
                            dsat_Py * dsat_Py +
                            dsat_Pz * dsat_Pz - Math.pow(SOL * rgTime, 2)));

            equationSet.put(2,
                    -((ellipsoidPosition.x * ellipsoidPosition.x + ellipsoidPosition.y * ellipsoidPosition.y) / (Math.pow(ell_a + refHeight, 2)) +
                            Math.pow(ellipsoidPosition.z / (ell_b + refHeight), 2) - 1.0));

            partialsXYZ.put(0, 0, satelliteVelocity.x);
            partialsXYZ.put(0, 1, satelliteVelocity.y);
            partialsXYZ.put(0, 2, satelliteVelocity.z);
            partialsXYZ.put(1, 0, 2 * dsat_Px);
            partialsXYZ.put(1, 1, 2 * dsat_Py);
            partialsXYZ.put(1, 2, 2 * dsat_Pz);
            partialsXYZ.put(2, 0, (2 * ellipsoidPosition.x) / (Math.pow(ell_a + refHeight, 2)));
            partialsXYZ.put(2, 1, (2 * ellipsoidPosition.y) / (Math.pow(ell_a + refHeight, 2)));
            partialsXYZ.put(2, 2, (2 * ellipsoidPosition.z) / (Math.pow(ell_a + refHeight, 2)));

            // solve system [NOTE!] orbit has to be normalized, otherwise close to singular
            DoubleMatrix ellipsoidPositionSolution = Solve.solve(partialsXYZ, equationSet);
            // DoubleMatrix ellipsoidPositionSolution = solve33(partialsXYZ, equationSet);

            // update solution
            ellipsoidPosition.x = ellipsoidPosition.x + ellipsoidPositionSolution.get(0);
            ellipsoidPosition.y = ellipsoidPosition.y + ellipsoidPositionSolution.get(1);
            ellipsoidPosition.z = ellipsoidPosition.z + ellipsoidPositionSolution.get(2);

            // check convergence
            if (Math.abs(ellipsoidPositionSolution.get(0)) < CRITERPOS &&
                    Math.abs(ellipsoidPositionSolution.get(1)) < CRITERPOS &&
                    Math.abs(ellipsoidPositionSolution.get(2)) < CRITERPOS) {
//                System.out.println("INFO: ellipsoidPosition (converged) = " + ellipsoidPosition);
                break;
            } else if (iter >= MAXITER) {
                MAXITER = MAXITER + 1;
                System.out.println("WARNING: line, pix -> x,y,z: maximum iterations (" + MAXITER + ") reached. " + "Criterium (m): " + CRITERPOS +
                        "dx,dy,dz=" + ellipsoidPositionSolution.get(0) + ", " + ellipsoidPositionSolution.get(1) + ", " + ellipsoidPositionSolution.get(2, 0));
            }
        }

        return ellipsoidPosition;
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
            System.out.println("WARNING: x,y,z -> line, pix: maximum iterations (" + MAXITER + ") reached. " +
                    "Criterium (s):" + CRITERTIM + "dta (s)=" + solution);
        }

        // Compute range time
        // ____ Update equations _____
        return posSat = getXYZ(tAzi);

    }

    public Point xyz2t(Point position, SLCImage slcimage, Orbit orbit) {

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
/*
            delta.x = position.x - satellitePosition.x;
            delta.y = position.y - satellitePosition.y;
            delta.z = position.z - satellitePosition.z;
*/

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
            System.out.println("WARNING: x,y,z -> line, pix: maximum iterations (" + MAXITER + ") reached. " + "Criterium (s):" + CRITERTIM + "dta (s)=" + solution);
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
/*
            delta.x = position.x - satellitePosition.x;
            delta.y = position.y - satellitePosition.y;
            delta.z = position.z - satellitePosition.z;
*/

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
            System.out.println("WARNING: x,y,z -> line, pix: maximum iterations (" + MAXITER + ") reached. " + "Criterium (s):" + CRITERTIM + "dta (s)=" + solution);
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

    // TODO
    public Point ell2lp(GeoPos philamheight, SLCImage slcimage) {
        double[] xyz = new double[3];
        GeoUtils.geo2xyz(philamheight, xyz);
        return xyz2lp(new Point(xyz), slcimage);
    }

    // TODO
    public GeoPos lp2ell(Point position, SLCImage slcimage) {
        GeoPos returnPos = new GeoPos();
        Point xyz = lp2xyz(position, slcimage);
        GeoUtils.xyz2geo(xyz.toArray(), returnPos);
        return returnPos;
    }

    public Point[][] dumpOrbit() {

        if (numStateVectors == 0) {
            System.out.println("Exiting Orbit.dumporbit(), no orbit data available.");
        }

        double dt = 0.;

        //  INFO << "dumporbits: MAXITER: "   << MAXITER   << "; "
        //                   << "CRITERPOS: " << CRITERPOS << " m; "
        //                   << "CRITERTIM: " << CRITERTIM << " s";
        //  INFO.print();

        //  ______ Evaluate polynomial orbit for t1:dt:tN ______
        int outputlines = 1 + (int) ((time.get(numStateVectors - 1, 0) - time.get(0, 0)) / dt);
        double tAzi = time.get(0, 0);

        Point[][] dumpedOrbit = new Point[(int) outputlines][3];

        for (int i = 0; i < outputlines; ++i) {
            dumpedOrbit[i][0] = getXYZ(tAzi);
            dumpedOrbit[i][1] = getXYZDot(tAzi);
            dumpedOrbit[i][2] = getXYZDotDot(tAzi);
            tAzi += dt;
        }


        //  // ______ dump coeff. as well for testing ... ______
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
        System.out.println("Time of orbit ephemerides: " + time.toString());
        System.out.println("Orbit ephemerides x:" + data_X.toString());
        System.out.println("Orbit ephemerides y:" + data_Y.toString());
        System.out.println("Orbit ephemerides z:" + data_Z.toString());
        System.out.println("Estimated coefficients x(t):" + coeff_X.toString());
        System.out.println("Estimated coefficients y(t):" + coeff_Y.toString());
        System.out.println("Estimated coefficients z(t):" + coeff_Z.toString());
    }

    // TODO
    public void computeBaseline() {

    }

    public Point getXYZ(double azTime) {

        //TODO: sanity check!
        Point satelliteXYZPosition = new Point();

        // normalize time
        azTime = (azTime - time.get(time.length / 2)) / 10;

        satelliteXYZPosition.x = MathUtilities.polyVal1d(azTime, coeff_X);
        satelliteXYZPosition.y = MathUtilities.polyVal1d(azTime, coeff_Y);
        satelliteXYZPosition.z = MathUtilities.polyVal1d(azTime, coeff_Z);

        return satelliteXYZPosition;  //To change body of created methods use File | Settings | File Templates.
    }

    public Point getXYZDot(double azTime) {

        //TODO: sanity check
        Point satelliteVelocity = new Point();

        // normalize time
        azTime = (azTime - time.get(time.length / 2)) / 10;

        // NOTE: orbit interpolator is simple polynomial
        satelliteVelocity.x = coeff_X.get(1);
        satelliteVelocity.y = coeff_Y.get(1);
        satelliteVelocity.z = coeff_Z.get(1);

        for (int i = 2; i <= poly_degree; ++i) {
            double powT = (double) i * Math.pow(azTime, (double) (i - 1));
            satelliteVelocity.x += coeff_X.get(i) * powT;
            satelliteVelocity.y += coeff_Y.get(i) * powT;
            satelliteVelocity.z += coeff_Z.get(i) * powT;
        }

        return satelliteVelocity.divByScalar(10.0d);

    }

    public Point getXYZDotDot(double azTime) {

        //TODO: sanity check
        Point satelliteAcceleration = new Point();

        // normalize time
        azTime = (azTime - time.get(time.length / 2)) / 10.0d;

        // NOTE: orbit interpolator is simple polynomial
        // 2a_2 + 2*3a_3*t^1 + 3*4a_4*t^2...

        for (int i = 2; i <= poly_degree; ++i) {
            double powT = (double) ((i - 1) * i) * Math.pow(azTime, (double) (i - 2));
            satelliteAcceleration.x += coeff_X.get(i) * powT;
            satelliteAcceleration.y += coeff_Y.get(i) * powT;
            satelliteAcceleration.z += coeff_Z.get(i) * powT;
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

}
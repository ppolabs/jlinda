package org.jdoris.core;

import org.apache.log4j.Logger;
import org.jblas.Decompose;
import org.jblas.DoubleMatrix;
import org.jblas.Solve;

import java.awt.*;

import static org.jblas.MatrixFunctions.abs;

public class MathUtilities {

    static Logger logger = Logger.getLogger(MathUtilities.class.getName());

    private static int polyDegree;

    public MathUtilities() {
    }

    public static int getPolyDegree() {
        return polyDegree;
    }

    /* Solve33
       Solves setof 3 equations by straightforward (no pivotting) LU
       y=Ax (unknown x)

        input:
          - matrix righthandside 3x1 (y)
          - matrix partials 3x3 (A)

        output/return:
         - matrix result 3x1 unknown
    */
    public static double[] solve33(double[][] A, double[] rhs) throws Exception {

        double[] result = new double[3];

        if (A[0].length != 3 || A.length != 3) {
            throw new Exception("solve33: input: size of A not 33.");
        }
        if (rhs.length != 3) {
            throw new Exception("solve33: input: size rhs not 3x1.");
        }

        // ______  real8 L10, L20, L21: used lower matrix elements
        // ______  real8 U11, U12, U22: used upper matrix elements
        // ______  real8 b0,  b1,  b2:  used Ux=b
        final double L10 = A[1][0] / A[0][0];
        final double L20 = A[2][0] / A[0][0];
        final double U11 = A[1][1] - L10 * A[0][1];
        final double L21 = (A[2][1] - (A[0][1] * L20)) / U11;
        final double U12 = A[1][2] - L10 * A[0][2];
        final double U22 = A[2][2] - L20 * A[0][2] - L21 * U12;

        // ______ Solution: forward substitution ______
        final double b0 = rhs[0];
        final double b1 = rhs[1] - b0 * L10;
        final double b2 = rhs[2] - b0 * L20 - b1 * L21;

        // ______ Solution: backwards substitution ______
        result[2] = b2 / U22;
        result[1] = (b1 - U12 * result[2]) / U11;
        result[0] = (b0 - A[0][1] * result[1] - A[0][2] * result[2]) / A[0][0];

        return result;

    }

    /**
     * solve22                                                   *
     * *
     * Solves setof 2 equations by straightforward substitution     *
     * y=Ax (unknown x)                                             *
     * *
     * input:                                                       *
     * - matrix<real8> righthandside 2x1 (y)                       *
     * - matrix<real8> partials 2x2 (A)                            *
     * output:                                                      *
     * - matrix<real8> result 2x1 unknown dx,dy,dz                 *
     * *
     */
    public static double[] solve22(double[][] A, double[] y) throws Exception {

        double[] result = new double[2];

        if (A[0].length != 2 || A.length != 2) {
            throw new Exception("solve33: input: size of A not 33.");
        }
        if (y.length != 2) {
            throw new Exception("solve33: input: size y not 3x1.");
        }

        // Direct Solution
        result[1] = (y[0] - ((A[0][0] / A[1][0]) * y[1])) / (A[0][1] - ((A[0][0] * A[1][1]) / A[1][0]));
        result[0] = (y[0] - A[0][1] * result[1]) / A[0][0];

        return result;

    }

    public double[][] distributePoints(int numOfPoints, Rectangle rectangle) {

        double lines = (rectangle.getMaxY() - rectangle.getMinY() + 1);
        double pixels = (rectangle.getMaxX() - rectangle.getMinX() + 1);

        double[][] result = new double[numOfPoints][2];

        // ______ Distribution for dl=dp ______
        double wp = Math.sqrt(numOfPoints / (lines / pixels));   // wl: #windows in line direction
        double wl = numOfPoints / wp;                   // wp: #windows in pixel direction
        if (wl < wp) {
            // switch wl,wp : later back
            wl = wp;
        }

        double wlint = Math.ceil(wl); // round largest
        double deltal = (lines - 1) / (wlint - 1);
        double totp = Math.ceil(pixels * wlint);
        double deltap = (totp - 1) / (numOfPoints - 1);
        double p = -deltap;
        double l = 0.;
        double lcnt = 0;
        int i;
        for (i = 0; i < numOfPoints; i++) {
            p += deltap;
            while (Math.ceil(p) >= pixels) // ceil
            {
                p -= pixels;
                lcnt++;
            }
            l = lcnt * deltal;

            result[i][0] = (int) Math.ceil(l);
            result[i][1] = (int) Math.ceil(p);
        }

//        // ______ Correct distribution to window ______
//        for (i=0; i<numOfPoints; i++){
//            result[i][0] += (int)rectangle.getMinY();
//            result[i][1] += (int)rectangle.getMinX();
//        }
        return result;
    }

    public static DoubleMatrix absMatrix(DoubleMatrix matrix) {
        for (int i = 0; i < matrix.rows; i++) {
            for (int j = 0; j < matrix.columns; j++) {
                matrix.put(i, j, Math.abs(matrix.get(i, j)));
            }
        }
        return matrix;
    }

    public static DoubleMatrix matrixPower(DoubleMatrix data, double scalar) {
        for (int i = 0; i < data.rows; ++i) {
            for (int j = 0; j < data.columns; ++j) {
                data.put(i, j, Math.pow(data.get(i, j), scalar));
            }
        }
        return data;
    }

    public static DoubleMatrix invertCholesky(DoubleMatrix matrix) {
        int numOfRows = matrix.rows;
        double sum;
        int i, j, k;
// ______ Compute inv(L) store in lower of A ______
        for (i = 0; i < numOfRows; ++i) {
            matrix.put(i, i, 1. / matrix.get(i, i));
            for (j = i + 1; j < numOfRows; ++j) {
                sum = 0.;
                for (k = i; k < j; ++k) {
                    sum -= matrix.get(j, k) * matrix.get(k, i);
                }
                matrix.put(j, i, sum / matrix.get(j, j));
            }
        }
// ______ Compute inv(A)=inv(LtL) store in lower of A ______
        for (i = 0; i < numOfRows; ++i) {
            for (j = i; j < numOfRows; ++j) {
                sum = 0.;
                for (k = j; k < numOfRows; ++k) {
                    sum += matrix.get(k, i) * matrix.get(k, j);                 // transpose
                }
                matrix.put(j, i, sum);
            }
        }
        return matrix;
    }

    public static double normalize(double data, int min, int max) {
        data -= (float) (0.5 * (min + max));
        data /= (float) (0.25 * (max - min));
        return data;
    }

    public static double normalize(double data, double min, double max) {
        data -= (0.5 * (min + max));
        data /= (0.25 * (max - min));
        return data;
    }


    public int numberOfCoefficients(int degree) {
//        return 0;  //To change body of created methods use File | Settings | File Templates.
        return (int) (0.5 * (Math.pow(degree + 1, 2) + degree + 1));
    }

    public int degreeFromCoefficients(int numOfCoefficients) {
//        return 0;  //To change body of created methods use File | Settings | File Templates.
        // TODO: validate this?!
        return (int) (0.5 * (-1 + (int) (Math.sqrt(1 + 8 * numOfCoefficients)))) - 1;
    }

    /****************************************************************
     *    polyfit                                                   *
     *                                                              *
     * Compute coefficients of x=a0+a1*t+a2*t^2+a3*t3 polynomial    *
     * for orbit interpolation.  Do this to facilitate a method     *
     * in case only a few datapoints are given.                     *
     * Data t is normalized approximately [-x,x], then polynomial   *
     * coefficients are computed.  For poly_val this is repeated    *
     * see getxyz, etc.                                             *
     *                                                              *
     * input:                                                       *
     *  - matrix by getdata with time and position info             *
     * output:                                                      *
     *  - matrix with coeff.                                        *
     *    (input for interp. routines)                              */
    public static double[] polyFit(DoubleMatrix time, DoubleMatrix y, int DEGREE) throws Exception {

        if (time.length != y.length) {
            logger.error("polyfit: require same size vectors.");
            throw new Exception("polyfit: require same size vectors.");
        }

        // Normalize _posting_ for numerical reasons
        final int numOfPoints = time.length;
        logger.debug("Normalizing t axis for least squares fit");
        DoubleMatrix normPosting = time.sub(time.get(numOfPoints / 2)).div(10.0);

        // Check redundancy
        final int numOfUnknowns = DEGREE + 1;
        logger.debug("Degree of orbit interpolating polynomial: " + DEGREE);
        logger.debug("Number of unknowns: " + numOfUnknowns);
        logger.debug("Number of data points (orbit): " + numOfPoints);

        if (numOfPoints < numOfUnknowns) {
            logger.error("Number of points is smaller than parameters solved for.");
            throw new Exception("Number of points is smaller than parameters solved for.");
        }

        // Set up system of equations to solve coeff
        logger.debug("Setting up linear system of equations");
        DoubleMatrix A = new DoubleMatrix(numOfPoints, numOfUnknowns);// designmatrix
        for (int j = 0; j <= DEGREE; j++) {
            DoubleMatrix normPostingTemp = normPosting.dup();
            normPostingTemp = matrixPower(normPostingTemp, (double) j);
            A.putColumn(j, normPostingTemp);
        }

        logger.debug("Solving lin. system of equations with Cholesky.");
        // Fit polynomial through computed vector of phases
//        DoubleMatrix y = y.dup();
        DoubleMatrix N = A.transpose().mmul(A);
        DoubleMatrix rhs = A.transpose().mmul(y);

        // solution seems to be OK up to 10^-09!
        DoubleMatrix x = Solve.solve(N, rhs);

        // JBLAS returns UPPER triangular, while we work with the LOWER triangular
        DoubleMatrix Qx_hat = Decompose.cholesky(N).transpose();

        // get covarinace matrix of normalized unknowns
        Qx_hat = invertCholesky(Qx_hat); // this could be more efficient

        // ______Test inverse______
        // repair matrix!! (see doris.core code for numerical example)
        for (int i = 0; i < Qx_hat.rows; i++) {
            for (int j = 0; j < i; j++) {
                Qx_hat.put(j, i, Qx_hat.get(i, j));
            }
        }


        double maxDeviation = abs(N.mmul(Qx_hat).sub(DoubleMatrix.eye(Qx_hat.rows))).max();
        logger.debug("polyfit orbit: max(abs(N*inv(N)-I)) = " + maxDeviation);

        // ___ report max error... (seems sometimes this can be extremely large) ___
        if (maxDeviation > 1e-6) {
            logger.warn("polyfit orbit: max(abs(N*inv(N)-I)) = " + maxDeviation);
            logger.warn("polyfit orbit interpolation unstable!");
        }

        // work out residuals
        DoubleMatrix y_hat = A.mmul(x);
        DoubleMatrix e_hat = y.sub(y_hat);

        DoubleMatrix e_hat_abs = abs(e_hat);

        // TODO: absMatrix(e_hat_abs).max() there is a simpleBlas function that implements this!
        // 0.05 is already 1 wavelength! (?)
        if (absMatrix(e_hat_abs).max() > 0.02) {
            logger.warn("WARNING: Max. approximation error at datapoints (x,y,or z?): " + abs(e_hat).max() + " m");

        }
        else {
            logger.info("Max. approximation error at datapoints (x,y,or z?): " + abs(e_hat).max() + " m");
        }

        logger.debug("REPORTING POLYFIT LEAST SQUARES ERRORS");
        logger.debug(" time \t\t\t y \t\t\t yhat  \t\t\t ehat");
        for (int i = 0; i < numOfPoints; i++) {
            logger.debug(" " + time.get(i) + "\t" + y.get(i) + "\t" + y_hat.get(i) + "\t" + e_hat.get(i));
        }

        for (int i = 0; i < numOfPoints - 1; i++) {
            // ___ check if dt is constant, not necessary for me, but may ___
            // ___ signal error in header data of SLC image ___
            double dt = time.get(i + 1) - time.get(i);
            logger.debug("Time step between point " + i + 1 + " and " + i + "= " + dt);

            if (Math.abs(dt - (time.get(1) - time.get(0))) > 0.001)// 1ms of difference we allow...
                logger.warn("WARNING: Orbit: data does not have equidistant time interval?");
        }

        return x.toArray();
    }



    public static double polyVal1d(double x, double[] coefficients) {
        double sum = 0.0;
        for (int d = coefficients.length - 1; d >= 0; --d) {
            sum *= x;
            sum += coefficients[d];
        }
        return sum;
    }

    public DoubleMatrix polyValOnGrid(DoubleMatrix x, DoubleMatrix y, final DoubleMatrix coeff, int degreee) {
        if (x.length != x.rows) {
            System.out.println("WARNING: polyal functions require (x) standing data vectors!");
        }

        if (y.length != y.rows) {
            System.out.println("WARNING: polyal functions require (y) standing data vectors!");
        }

        if (coeff.length != coeff.rows) {
            System.out.println("WARNING: polyal functions require (coeff) standing data vectors!");
        }

        if (degreee < -1) {
            System.out.println("WARNING: polyal degree < -1 ????");
        }

//        if (x.length > y.length) {
//            System.out.println("WARNING: polyal function, x larger than y, while optimized for y larger x");
//        }

        if (degreee == -1) {
            degreee = degreeFromCoefficients(coeff.length);
        }

        // evaluate polynomial
        DoubleMatrix result = new DoubleMatrix(new double[x.length][y.length]);
        int i;
        int j;

//        double sum = coeff.get(0,0);
        switch (degreee) {
            case 0:
                result.put(0, 0, coeff.get(0, 0));
                break;
            case 1:
                double c00 = coeff.get(0, 0);
                double c10 = coeff.get(1, 0);
                double c01 = coeff.get(2, 0);
                for (j = 0; j < result.columns; j++) {
                    double c00pc01y1 = c00 + c01 * y.get(j, 0);
                    for (i = 0; i < result.rows; i++) {
                        result.put(i, j, c00pc01y1 + c10 * x.get(i, 0));
                    }
                }
                break;
            case 2:
                c00 = coeff.get(0, 0);
                c10 = coeff.get(1, 0);
                c01 = coeff.get(2, 0);
                double c20 = coeff.get(3, 0);
                double c11 = coeff.get(4, 0);
                double c02 = coeff.get(5, 0);
                for (j = 0; j < result.columns; j++) {
                    double y1 = y.get(j, 0);
                    double c00pc01y1 = c00 + c01 * y1;
                    double c02y2 = c02 * Math.pow(y1, 2);
                    double c11y1 = c11 * y1;
                    for (i = 0; i < result.rows; i++) {
                        double x1 = x.get(i, 0);
                        result.put(i, j, c00pc01y1
                                + c10 * x1
                                + c20 * Math.pow(x1, 2)
                                + c11y1 * x1
                                + c02y2);
                    }
                }
                break;
            case 3:
                c00 = coeff.get(0, 0);
                c10 = coeff.get(1, 0);
                c01 = coeff.get(2, 0);
                c20 = coeff.get(3, 0);
                c11 = coeff.get(4, 0);
                c02 = coeff.get(5, 0);
                double c30 = coeff.get(6, 0);
                double c21 = coeff.get(7, 0);
                double c12 = coeff.get(8, 0);
                double c03 = coeff.get(9, 0);
                for (j = 0; j < result.columns; j++) {
                    double y1 = y.get(j, 0);
                    double y2 = Math.pow(y1, 2);
                    double c00pc01y1 = c00 + c01 * y1;
                    double c02y2 = c02 * y2;
                    double c11y1 = c11 * y1;
                    double c21y1 = c21 * y1;
                    double c12y2 = c12 * y2;
                    double c03y3 = c03 * y1 * y2;
                    for (i = 0; i < result.rows; i++) {
                        double x1 = x.get(i, 0);
                        double x2 = Math.pow(x1, 2);
                        result.put(i, j, c00pc01y1
                                + c10 * x1
                                + c20 * x2
                                + c11y1 * x1
                                + c02y2
                                + c30 * x1 * x2
                                + c21y1 * x2
                                + c12y2 * x1
                                + c03y3);
                    }
                }
                break;

            case 4:
                c00 = coeff.get(0, 0);
                c10 = coeff.get(1, 0);
                c01 = coeff.get(2, 0);
                c20 = coeff.get(3, 0);
                c11 = coeff.get(4, 0);
                c02 = coeff.get(5, 0);
                c30 = coeff.get(6, 0);
                c21 = coeff.get(7, 0);
                c12 = coeff.get(8, 0);
                c03 = coeff.get(9, 0);
                double c40 = coeff.get(10, 0);
                double c31 = coeff.get(11, 0);
                double c22 = coeff.get(12, 0);
                double c13 = coeff.get(13, 0);
                double c04 = coeff.get(14, 0);
                for (j = 0; j < result.columns; j++) {
                    double y1 = y.get(j, 0);
                    double y2 = Math.pow(y1, 2);
                    double c00pc01y1 = c00 + c01 * y1;
                    double c02y2 = c02 * y2;
                    double c11y1 = c11 * y1;
                    double c21y1 = c21 * y1;
                    double c12y2 = c12 * y2;
                    double c03y3 = c03 * y1 * y2;
                    double c31y1 = c31 * y1;
                    double c22y2 = c22 * y2;
                    double c13y3 = c13 * y2 * y1;
                    double c04y4 = c04 * y2 * y2;
                    for (i = 0; i < result.rows; i++) {
                        double x1 = x.get(i, 0);
                        double x2 = Math.pow(x1, 2);
                        result.put(i, j, c00pc01y1
                                + c10 * x1
                                + c20 * x2
                                + c11y1 * x1
                                + c02y2
                                + c30 * x1 * x2
                                + c21y1 * x2
                                + c12y2 * x1
                                + c03y3
                                + c40 * x2 * x2
                                + c31y1 * x2 * x1
                                + c22y2 * x2
                                + c13y3 * x1
                                + c04y4);
                    }
                }
                break;
            case 5:
                c00 = coeff.get(0, 0);
                c10 = coeff.get(1, 0);
                c01 = coeff.get(2, 0);
                c20 = coeff.get(3, 0);
                c11 = coeff.get(4, 0);
                c02 = coeff.get(5, 0);
                c30 = coeff.get(6, 0);
                c21 = coeff.get(7, 0);
                c12 = coeff.get(8, 0);
                c03 = coeff.get(9, 0);
                c40 = coeff.get(10, 0);
                c31 = coeff.get(11, 0);
                c22 = coeff.get(12, 0);
                c13 = coeff.get(13, 0);
                c04 = coeff.get(14, 0);
                double c50 = coeff.get(15, 0);
                double c41 = coeff.get(16, 0);
                double c32 = coeff.get(17, 0);
                double c23 = coeff.get(18, 0);
                double c14 = coeff.get(19, 0);
                double c05 = coeff.get(20, 0);
                for (j = 0; j < result.columns; j++) {
                    double y1 = y.get(j, 0);
                    double y2 = Math.pow(y1, 2);
                    double y3 = y2 * y1;
                    double c00pc01y1 = c00 + c01 * y1;
                    double c02y2 = c02 * y2;
                    double c11y1 = c11 * y1;
                    double c21y1 = c21 * y1;
                    double c12y2 = c12 * y2;
                    double c03y3 = c03 * y3;
                    double c31y1 = c31 * y1;
                    double c22y2 = c22 * y2;
                    double c13y3 = c13 * y3;
                    double c04y4 = c04 * y2 * y2;
                    double c41y1 = c41 * y1;
                    double c32y2 = c32 * y2;
                    double c23y3 = c23 * y3;
                    double c14y4 = c14 * y2 * y2;
                    double c05y5 = c05 * y3 * y2;
                    for (i = 0; i < result.rows; i++) {
                        double x1 = x.get(i, 0);
                        double x2 = Math.pow(x1, 2);
                        double x3 = x1 * x2;
                        result.put(i, j, c00pc01y1
                                + c10 * x1
                                + c20 * x2
                                + c11y1 * x1
                                + c02y2
                                + c30 * x3
                                + c21y1 * x2
                                + c12y2 * x1
                                + c03y3
                                + c40 * x2 * x2
                                + c31y1 * x3
                                + c22y2 * x2
                                + c13y3 * x1
                                + c04y4
                                + c50 * x3 * x2
                                + c41y1 * x2 * x2
                                + c32y2 * x3
                                + c23y3 * x2
                                + c14y4 * x1
                                + c05y5);
                    }
                }
                break;

            // ______ solve up to 5 efficiently, do rest in loop ______
            default:
                c00 = coeff.get(0, 0);
                c10 = coeff.get(1, 0);
                c01 = coeff.get(2, 0);
                c20 = coeff.get(3, 0);
                c11 = coeff.get(4, 0);
                c02 = coeff.get(5, 0);
                c30 = coeff.get(6, 0);
                c21 = coeff.get(7, 0);
                c12 = coeff.get(8, 0);
                c03 = coeff.get(9, 0);
                c40 = coeff.get(10, 0);
                c31 = coeff.get(11, 0);
                c22 = coeff.get(12, 0);
                c13 = coeff.get(13, 0);
                c04 = coeff.get(14, 0);
                c50 = coeff.get(15, 0);
                c41 = coeff.get(16, 0);
                c32 = coeff.get(17, 0);
                c23 = coeff.get(18, 0);
                c14 = coeff.get(19, 0);
                c05 = coeff.get(20, 0);
                for (j = 0; j < result.columns; j++) {
                    double y1 = y.get(j, 0);
                    double y2 = Math.pow(y1, 2);
                    double y3 = y2 * y1;
                    double c00pc01y1 = c00 + c01 * y1;
                    double c02y2 = c02 * y2;
                    double c11y1 = c11 * y1;
                    double c21y1 = c21 * y1;
                    double c12y2 = c12 * y2;
                    double c03y3 = c03 * y3;
                    double c31y1 = c31 * y1;
                    double c22y2 = c22 * y2;
                    double c13y3 = c13 * y3;
                    double c04y4 = c04 * y2 * y2;
                    double c41y1 = c41 * y1;
                    double c32y2 = c32 * y2;
                    double c23y3 = c23 * y3;
                    double c14y4 = c14 * y2 * y2;
                    double c05y5 = c05 * y3 * y2;
                    for (i = 0; i < result.rows; i++) {
                        double x1 = x.get(i, 0);
                        double x2 = Math.pow(x1, 2);
                        double x3 = x1 * x2;
                        result.put(i, j, c00pc01y1
                                + c10 * x1
                                + c20 * x2
                                + c11y1 * x1
                                + c02y2
                                + c30 * x3
                                + c21y1 * x2
                                + c12y2 * x1
                                + c03y3
                                + c40 * x2 * x2
                                + c31y1 * x3
                                + c22y2 * x2
                                + c13y3 * x1
                                + c04y4
                                + c50 * x3 * x2
                                + c41y1 * x2 * x2
                                + c32y2 * x3
                                + c23y3 * x2
                                + c14y4 * x1
                                + c05y5);
                    }
                }

                final int STARTDEGREE = 6;
                final int STARTCOEFF = degreeFromCoefficients(STARTDEGREE - 1);   // 5-> 21 6->28 7->36 etc.
                for (j = 0; j < result.columns; j++) {
                    double yy = y.get(j, 0);
                    for (i = 0; i < result.rows; i++) {
                        double xx = x.get(i, 0);        // ??? this seems to be wrong (BK 9-feb-00)
                        double sum = 0.;
                        int coeffindex = STARTCOEFF;
                        for (int l = STARTDEGREE; l <= degreee; l++) {
                            for (int k = 0; k <= l; k++) {
                                sum += coeff.get(coeffindex, 0) * Math.pow(xx, (double) (l - k)) * Math.pow(yy, (double) (k));
                                coeffindex++;
                            }
                        }
                        result.put(i, j, result.get(i, j) + sum);
                    }
                }
        } // switch degree

        return result;
    }


    // TODO: refactor and better integrate helper functions for Baseline.class
    // HELPER FUNCTIONS == from Baseline class
    public static DoubleMatrix matTxmat(DoubleMatrix matrix1, DoubleMatrix matrix2) {
        return matrix1.transpose().mmul(matrix2);
    }

    public static double rad2deg(double rad) {
        return Math.toDegrees(rad);
    }


    public static double sqr(double value) {
        return Math.pow(value, 2);
    }

    public static double sqrt(double value) {
        return Math.sqrt(value);
    }


}
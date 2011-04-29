package org.jdoris.core.utils;

import org.apache.log4j.Logger;
import org.jblas.Decompose;
import org.jblas.DoubleMatrix;
import org.jblas.Solve;

import static org.jblas.MatrixFunctions.abs;

public class PolyUtils {
    public static Logger logger = Logger.getLogger(PolyUtils.class.getName());

    public static double normalize(double data, final int min, final int max) {
        data -= (float) (0.5 * (min + max));
        data /= (float) (0.25 * (max - min));
        return data;
    }

    public static double normalize(double data, final double min, final double max) {
        data -= (0.5 * (min + max));
        data /= (0.25 * (max - min));
        return data;
    }

    public static int degreeFromCoefficients(int numOfCoefficients) {
//        return 0;  //To change body of created methods use File | Settings | File Templates.
        // TODO: validate this?!
        return (int) (0.5 * (-1 + (int) (Math.sqrt(1 + 8 * numOfCoefficients)))) - 1;
    }

    /**
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
     *    (input for interp. routines)
     */
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
            normPostingTemp = LinearAlgebraUtils.matrixPower(normPostingTemp, (double) j);
            A.putColumn(j, normPostingTemp);
        }

        // Fit polynomial through computed vector of phases
        logger.debug("Solving lin. system of equations with Cholesky.");

//        DoubleMatrix y = y.dup();
        DoubleMatrix N = A.transpose().mmul(A);
        DoubleMatrix rhs = A.transpose().mmul(y);

        // solution seems to be OK up to 10^-09!
        DoubleMatrix x = Solve.solve(N, rhs);

        // JBLAS returns UPPER triangular, while we work(!!!!) with the LOWER triangular
        DoubleMatrix Qx_hat = Decompose.cholesky(N).transpose();

        // get covarinace matrix of normalized unknowns
        Qx_hat = LinearAlgebraUtils.invertCholesky(Qx_hat); // this could be more efficient

        // ______Test inverse______
        // repair matrix!! (see doris.core code for numerical example)
        logger.debug("reparing cholesky decomposed matrix");
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
        if (LinearAlgebraUtils.absMatrix(e_hat_abs).max() > 0.02) {
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

    public static DoubleMatrix polyValOnGrid(DoubleMatrix x, DoubleMatrix y, final DoubleMatrix coeff, int degreee) {
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

    public static int numberOfCoefficients(final int degree) {
//        return 0;  //To change body of created methods use File | Settings | File Templates.
        return (int) (0.5 * (Math.pow(degree + 1, 2) + degree + 1));
    }
}

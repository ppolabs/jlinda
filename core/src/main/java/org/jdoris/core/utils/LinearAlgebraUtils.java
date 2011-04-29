package org.jdoris.core.utils;

import org.apache.log4j.Logger;
import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jdoris.core.Window;

public class LinearAlgebraUtils {

    static Logger logger = Logger.getLogger(LinearAlgebraUtils.class.getName());


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

        // real8 L10, L20, L21: used lower matrix elements
        // real8 U11, U12, U22: used upper matrix elements
        // real8 b0,  b1,  b2:  used Ux=b
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
     * solve22
     * Solves setof 2 equations by straightforward substitution
     * y=Ax (unknown x)
     * input:
     * - matrix<real8> righthandside 2x1 (y)
     * - matrix<real8> partials 2x2 (A)
     * output:
     * - matrix<real8> result 2x1 unknown dx,dy,dz
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
        final int numOfRows = matrix.rows;
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

    /**
     * setdata(B, winB, A, winA):
     * set winB of B to winA of A
     * if winB==0 defaults to totalB, winA==0 defaults to totalA
     * first line matrix =0 (?)
     */
    public static void setdata(DoubleMatrix B, Window winB, DoubleMatrix A, Window winA) {

        // Check default request
        if (winB.linehi == 0 && winB.pixhi == 0) {
            winB.linehi = B.rows - 1;
            winB.pixhi = B.columns - 1;
        }
        if (winA.linehi == 0 && winA.pixhi == 0) {
            winA.linehi = A.rows - 1;
            winA.pixhi = A.columns - 1;
        }

        // TODO: for now errors only logged, introduce exceptions
        // More sanity checks
        if (((winB.linehi - winB.linelo) != (winA.linehi - winA.linelo)) ||
                ((winB.pixhi - winB.pixlo) != (winA.pixhi - winA.pixlo)))
            logger.error("setdata: wrong input.");

        if (winB.linehi < winB.linelo || winB.pixhi < winB.pixlo)
            logger.error("setdata: wrong input.1");

        if ((winB.linehi > B.rows - 1) ||
                (winB.pixhi > B.columns - 1))
            logger.error("setdata: wrong input.2");

        if ((winA.linehi > A.rows - 1) ||
                (winA.pixhi > A.columns - 1))
            logger.error("setdata: wrong input.3");

        // Fill data
        int sizeLin = (int) winA.pixels();
        for (int i = (int) winB.linelo; i <= winB.linehi; i++) {

            int startA = (int) (i * A.columns + winA.pixlo);
            int startB = (int) (i * B.columns + winB.pixlo);

            System.arraycopy(A.data, startA, B.data, startB, sizeLin);

        }
    }

    /**
     * setdata(B, winB, A, winA):
     * set winB of B to winA of A
     * if winB==0 defaults to totalB, winA==0 defaults to totalA
     * first line matrix =0 (?)
     */
    public static void setdata(ComplexDoubleMatrix B, Window winB, ComplexDoubleMatrix A, Window winA) {

        // Check default request
        if (winB.linehi == 0 && winB.pixhi == 0) {
            winB.linehi = B.rows - 1;
            winB.pixhi = B.columns - 1;
        }
        if (winA.linehi == 0 && winA.pixhi == 0) {
            winA.linehi = A.rows - 1;
            winA.pixhi = A.columns - 1;
        }

        // TODO: for now errors only logged, introduce exceptions
        // More sanity checks
        if (((winB.linehi - winB.linelo) != (winA.linehi - winA.linelo)) ||
                ((winB.pixhi - winB.pixlo) != (winA.pixhi - winA.pixlo)))
            logger.error("setdata: wrong input.");

        if (winB.linehi < winB.linelo || winB.pixhi < winB.pixlo)
            logger.error("setdata: wrong input.1");

        if ((winB.linehi > B.rows - 1) ||
                (winB.pixhi > B.columns - 1))
            logger.error("setdata: wrong input.2");

        if ((winA.linehi > A.rows - 1) ||
                (winA.pixhi > A.columns - 1))
            logger.error("setdata: wrong input.3");

        // Fill data
        int sizeLin = (int) winA.pixels() * 2;
        for (int i = (int) winB.linelo; i <= winB.linehi; i++) {

            int startA = (int) (i * A.columns * 2 + winA.pixlo * 2);
            int startB = (int) (i * B.columns * 2 + winB.pixlo * 2);

            System.arraycopy(A.data, startA, B.data, startB, sizeLin);
        }
    }

    public static ComplexDoubleMatrix dotmult(ComplexDoubleMatrix A, ComplexDoubleMatrix B) {
        return A.mmul(B);
    }

    public static void dotmultIn(ComplexDoubleMatrix A, ComplexDoubleMatrix B) {
        A.mmul(B);
    }

    /**
    * B.fliplr()
    * Mirror in center vertical (flip left right).
    */
    public static void fliplr(DoubleMatrix A) {

        int nrows = A.rows;
        int ncols = A.columns;

        if (A.rows == 1) {
            double tmp;
            for (int i=0; i< (ncols/2); ++i) {
                tmp = A.get(1, i);
                A.put(1,i,A.get(1, nrows - i));
                A.put(1, nrows - 1, tmp);
            }
        } else {
            for (int i = 0; i < (ncols / 2); ++i)     // floor
            {
                DoubleMatrix tmp1 = A.getColumn(i);
                DoubleMatrix tmp2 = A.getColumn(ncols - i - 1);
                A.putColumn(i, tmp2);
                A.putColumn(ncols - i - 1, tmp1);
            }
        }
    }

    // TODO: refactor and better integrate helper functions for Baseline.class
    // HELPER FUNCTIONS == from Baseline class
    public static DoubleMatrix matTxmat(DoubleMatrix matrix1, DoubleMatrix matrix2) {
        return matrix1.transpose().mmul(matrix2);
    }

    public static ComplexDoubleMatrix matTxmat(ComplexDoubleMatrix matrix1, ComplexDoubleMatrix matrix2) {
        return matrix1.transpose().mmul(matrix2);
    }

    public static void setdata(ComplexDoubleMatrix B, ComplexDoubleMatrix A, Window winA) {
        setdata(B, new Window(0, B.rows, 0, B.columns), A, winA);
    }

    /**
     * wshift(A,n)                                                  *
     * circular shift of vector A by n pixels. positive n for    *
     * right to left shift.                                      *
     * implementation: WSHIFT(A,n) == WSHIFT(A,n-sizeA);         *
     * A is changed itself!                                      *
     */

    public static void wshift(DoubleMatrix A, int n) throws Exception {

        if (n >= A.length) {
            System.err.println("wshift: shift larger than matrix not implemented.");
            throw new Exception();
        }

        if (!A.isVector()) {
            System.err.println("wshift: only vectors");
            throw new Exception();
        }

        // positive only, use rem!  n = n%A.nsize;
        if (n == 0) return;
        if (n < 0) n += A.length;

        DoubleMatrix Res = new DoubleMatrix(A.rows, A.columns);

        //  n always >0 here
        System.arraycopy(A.data, n, Res.data, 0, A.length - n);
        System.arraycopy(A.data, 0, Res.data, A.length - n, n);

        A.copy(Res);

    }
}

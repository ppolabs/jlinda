package org.jdoris.core.utils;

import org.apache.log4j.Logger;
import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jdoris.core.Window;

import static org.jblas.MatrixFunctions.abs;
import static org.jblas.MatrixFunctions.pow;

public class LinearAlgebraUtils {

    static Logger logger = Logger.getLogger(LinearAlgebraUtils.class.getName());

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
            throw new Exception("solve22: input: size of A not 22.");
        }
        if (y.length != 2) {
            throw new Exception("solve22: input: size y not 2x1.");
        }

        // Direct Solution
        result[1] = (y[0] - ((A[0][0] / A[1][0]) * y[1])) / (A[0][1] - ((A[0][0] * A[1][1]) / A[1][0]));
        result[0] = (y[0] - A[0][1] * result[1]) / A[0][0];

        return result;

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

    public static DoubleMatrix absMatrix(DoubleMatrix inMatrix) {
        return abs(inMatrix);
    }

    public static double[][] absMatrix(double[][] inMatrix) {
        for (int i = 0; i < inMatrix.length; i++) {
            for (int j = 0; j < inMatrix[0].length; j++) {
                inMatrix[i][j] = Math.abs(inMatrix[i][j]);
            }
        }
        return inMatrix;
    }

    public static DoubleMatrix matrixPower(DoubleMatrix data, double scalar) {
        return pow(data, scalar);
    }

    public static double[][] matrixPower(double[][] inMatrix, double scalar) {
        for (int i = 0; i < inMatrix.length; ++i) {
            for (int j = 0; j < inMatrix[0].length; ++j) {
                inMatrix[i][j] = Math.pow(inMatrix[i][j], scalar);
            }
        }
        return inMatrix;
    }

    public static void invertChol_inplace(double[][] inMatrix) {
        final int numOfRows = inMatrix.length;
        double sum;
        int i, j, k;
        // Compute inv(L) store in lower of inMatrix
        for (i = 0; i < numOfRows; ++i) {
            inMatrix[i][i] = 1. / inMatrix[i][i];
            for (j = i + 1; j < numOfRows; ++j) {
                sum = 0.;
                for (k = i; k < j; ++k) {
                    sum -= inMatrix[j][k] * inMatrix[k][i];
                }
                inMatrix[j][i] = sum / inMatrix[j][j];
            }
        }
        // Compute inv(inMatrix)=inv(LtL) store in lower of inMatrix
        for (i = 0; i < numOfRows; ++i) {
            for (j = i; j < numOfRows; ++j) {
                sum = 0.;
                for (k = j; k < numOfRows; ++k) {
                    sum += inMatrix[k][i] * inMatrix[k][j];
                }
                inMatrix[j][i] = sum;
            }
        }
    }

    public static double[][] invertChol(double[][] inMatrix) {
        double[][] outMatrix = inMatrix.clone();
        invertChol_inplace(outMatrix);
        return outMatrix;
    }

    public static void invertChol_inplace(DoubleMatrix inMatrix) {
        invertChol_inplace(inMatrix.toArray2());
    }

    public static DoubleMatrix invertChol(DoubleMatrix inMatrix) {
        DoubleMatrix outMatrix = inMatrix.dup();
        invertChol_inplace(outMatrix);
        return outMatrix;
    }

    public static ComplexDoubleMatrix dotmult(ComplexDoubleMatrix A, ComplexDoubleMatrix B) {
        return A.mmul(B);
    }

    public static void dotmult_inplace(ComplexDoubleMatrix A, ComplexDoubleMatrix B) {
        A.mmuli(B);
    }


    /**
     * B.fliplr()
     * Mirror in center vertical (flip left right).
     */
    public static void fliplr(DoubleMatrix A) {

        final int nRows = A.rows;
        final int nCols = A.columns;

        if (nRows == 1) {
            double tmp;
            for (int i = 0; i < (nCols / 2); ++i) {
                tmp = A.get(1, i);
                A.put(1, i, A.get(1, nRows - i));
                A.put(1, nRows - 1, tmp);
            }
        } else {
            for (int i = 0; i < (nCols / 2); ++i)     // floor
            {
                DoubleMatrix tmp1 = A.getColumn(i);
                DoubleMatrix tmp2 = A.getColumn(nCols - i - 1);
                A.putColumn(i, tmp2);
                A.putColumn(nCols - i - 1, tmp1);
            }
        }
    }

    public static DoubleMatrix matTxmat(DoubleMatrix matrix1, DoubleMatrix matrix2) {
        return matrix1.transpose().mmul(matrix2);
    }

    public static ComplexDoubleMatrix matTxmat(ComplexDoubleMatrix matrix1, ComplexDoubleMatrix matrix2) {
        return matrix1.transpose().mmul(matrix2);
    }

    /**
     * wshift(A,n)                                                  *
     * circular shift of vector A by n pixels. positive n for    *
     * right to left shift.                                      *
     * implementation: WSHIFT(A,n) == WSHIFT(A,n-sizeA);         *
     * A is changed itself!                                      *
     */
    public static void wshift_inplace(DoubleMatrix A, int n) throws Exception {

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

    public static DoubleMatrix wshift(DoubleMatrix inMatrix, int n) throws Exception {
        DoubleMatrix outMatrix = inMatrix.dup();
        wshift_inplace(outMatrix, n);
        return outMatrix;
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

    public static void setdata(ComplexDoubleMatrix B, ComplexDoubleMatrix A, Window winA) {
        setdata(B, new Window(0, B.rows, 0, B.columns), A, winA);
    }


    /**
     * *************************************************************
     * choles(A);   cholesky factorisation internal implementation  *
     * lower triangle of A is changed on output                     *
     * upper reamins un referenced                                  *
     * this one is a lot slower then veclib and there may be more   *
     * efficient implementations.                                   *
     * **************************************************************
     */
    public static void choles_inplace(double[][] A) {
        final int N = A.length;
        double sum;
        for (int i = 0; i < N; ++i) {
            for (int j = i; j < N; ++j) {
                sum = A[i][j];
                for (int k = i - 1; k >= 0; --k) {
                    sum -= A[i][k] * A[j][k];
                }
                if (i == j) {
                    if (sum <= 0.) {
                        logger.error("choles: internal: A not pos. def.");
                    }
                    A[i][i] = Math.sqrt(sum);
                } else {
                    A[j][i] = sum / A[i][i];
                }
            }
        }
    }

}

package org.jdoris.core.utils;

import org.jblas.DoubleMatrix;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

public class LinearAlgebraUtilsTest {

//    final static private int NUM_ROWS = 4;
//    final static private int NUM_COLS = 4;

    private static DoubleMatrix A_PASCAL_22 = new DoubleMatrix(2, 2);
    private static DoubleMatrix A_PASCAL_33 = new DoubleMatrix(3, 3);
    private static double[][] A_PASCAL_33_SQUARED = new double[3][3];
    private static final double[] X_22 = new double[]{3, 1};

    private static final double[] X_33 = new double[]{3, 1, 4};
    private static final double[] SOL_22_EXPECTED = new double[]{5, -2};

    private static final double[] SOL_33_EXPECTED = new double[]{10,-12, 5};
    private static final double[] SOL_33_EXPECTED_ABS = new double[]{10, 12, 5};
    private static final double DELTA = 1e-06;

    private static DoubleMatrix A_PASCAL_33_CHOL_EXPECTED = new DoubleMatrix(3, 3);
    private static DoubleMatrix A_PASCAL_33_CHOL_INV_EXPECTED = new DoubleMatrix(3, 3);

    @BeforeClass
    public static void setUpTestMatrices() throws Exception {


        A_PASCAL_22 = DoubleMatrix.ones(2,2);
        A_PASCAL_22.put(1, 1, 2);

        A_PASCAL_33 = DoubleMatrix.ones(3, 3);
        A_PASCAL_33.put(1, 1, 2);
        A_PASCAL_33.put(1, 2, 3);
        A_PASCAL_33.put(2, 1, 3);
        A_PASCAL_33.put(2, 2, 6);

        A_PASCAL_33_SQUARED = new double[][]{{1, 1, 1}, {1, 8, 27}, {1, 27, 216}};

        A_PASCAL_33_CHOL_EXPECTED = A_PASCAL_33.dup();
        // define lower triangular block
        A_PASCAL_33_CHOL_EXPECTED.put(1, 0, 1);
        A_PASCAL_33_CHOL_EXPECTED.put(1, 1, 1);
        A_PASCAL_33_CHOL_EXPECTED.put(2, 0, 1);
        A_PASCAL_33_CHOL_EXPECTED.put(2, 1, 2);
        A_PASCAL_33_CHOL_EXPECTED.put(2, 2, 1);

        // define inverted matrix
        A_PASCAL_33_CHOL_INV_EXPECTED.put(0, 0, 3);
        A_PASCAL_33_CHOL_INV_EXPECTED.put(0, 1, -3);
        A_PASCAL_33_CHOL_INV_EXPECTED.put(0, 2, 1);
        A_PASCAL_33_CHOL_INV_EXPECTED.put(1, 0, -3);
        A_PASCAL_33_CHOL_INV_EXPECTED.put(1, 1, 5);
        A_PASCAL_33_CHOL_INV_EXPECTED.put(1, 2, -2);
        A_PASCAL_33_CHOL_INV_EXPECTED.put(2, 0, 1);
        A_PASCAL_33_CHOL_INV_EXPECTED.put(2, 1, -2);
        A_PASCAL_33_CHOL_INV_EXPECTED.put(2, 2, 1);

    }

    @Test
    public void testSolve22() throws Exception {
        double[] SOL_22_ACTUAL = LinearAlgebraUtils.solve22(A_PASCAL_22.toArray2(), X_22);
        Assert.assertArrayEquals(SOL_22_EXPECTED, SOL_22_ACTUAL, DELTA);
    }

    @Test
    public void testSolve33() throws Exception {
        double[] SOL_33_ACTUAL = LinearAlgebraUtils.solve33(A_PASCAL_33.toArray2(), X_33);
        Assert.assertArrayEquals(SOL_33_EXPECTED, SOL_33_ACTUAL, DELTA);

    }

    @Test
    public void testAbsMatrixJBLAS() throws Exception {
        DoubleMatrix SOL_33_ACTUAL_ABS = LinearAlgebraUtils.absMatrix(new DoubleMatrix(SOL_33_EXPECTED));
        Assert.assertEquals(new DoubleMatrix(SOL_33_EXPECTED_ABS), SOL_33_ACTUAL_ABS);
    }

    @Test
    public void testAbsMatrixArrays() throws Exception {
        double[][] SOL_33_EXPECTED_TEMP = new double[][]{SOL_33_EXPECTED,SOL_33_EXPECTED};
        double[][] SOL_33_EXPECTED_ABS_TEMP = new double[][]{SOL_33_EXPECTED_ABS,SOL_33_EXPECTED_ABS};

        double[][] SOL_33_ACTUAL_ABS = LinearAlgebraUtils.absMatrix(SOL_33_EXPECTED_TEMP);
        Assert.assertArrayEquals(SOL_33_EXPECTED_ABS_TEMP, SOL_33_ACTUAL_ABS);
    }

    @Test
    public void testMatrixPowerJBLAS() throws Exception {
        Assert.assertEquals(new DoubleMatrix(A_PASCAL_33_SQUARED),
                LinearAlgebraUtils.matrixPower(A_PASCAL_33, 3));
    }

    @Test
    public void testMatrixPowerArrays() throws Exception {
        Assert.assertEquals(A_PASCAL_33_SQUARED,
                LinearAlgebraUtils.matrixPower(A_PASCAL_33.toArray2(), 3));
    }


    @Test
    public void testCholesky() throws Exception {
        double[][] A_PASCAL_33_CHOL_ARRAY_ACTUAL = A_PASCAL_33.toArray2();
        LinearAlgebraUtils.choles_inplace(A_PASCAL_33_CHOL_ARRAY_ACTUAL);
        Assert.assertEquals(A_PASCAL_33_CHOL_EXPECTED.toArray2(),A_PASCAL_33_CHOL_ARRAY_ACTUAL);
    }

    @Test
    public void testInvertCholesky() throws Exception {

        double[][] A_PASCAL_33_CHOL_ARRAY_ACTUAL = A_PASCAL_33.toArray2();
        LinearAlgebraUtils.choles_inplace(A_PASCAL_33_CHOL_ARRAY_ACTUAL);
        LinearAlgebraUtils.invertChol_inplace(A_PASCAL_33_CHOL_ARRAY_ACTUAL);

        // assume squared
        for (int i = 0; i < A_PASCAL_33_CHOL_ARRAY_ACTUAL.length; i++) {
            for (int j = 0; j < i; j++) {
                A_PASCAL_33_CHOL_ARRAY_ACTUAL[j][i] = A_PASCAL_33_CHOL_ARRAY_ACTUAL[i][j];
                A_PASCAL_33_CHOL_ARRAY_ACTUAL[j][i] = A_PASCAL_33_CHOL_ARRAY_ACTUAL[i][j];
            }
        }

        Assert.assertEquals(A_PASCAL_33_CHOL_INV_EXPECTED.toArray2(), A_PASCAL_33_CHOL_ARRAY_ACTUAL);


    }

//    @Test
//    public void testSetdata() throws Exception {
//
//    }
//
//    @Test
//    public void testSetdata() throws Exception {
//
//    }
//
//    @Test
//    public void testDotmult() throws Exception {
//
//    }
//
//    @Test
//    public void testDotmultIn() throws Exception {
//
//    }
//
//    @Test
//    public void testFliplr() throws Exception {
//
//    }
//
//    @Test
//    public void testMatTxmat() throws Exception {
//
//    }
//
//    @Test
//    public void testMatTxmat() throws Exception {
//
//    }
//
//    @Test
//    public void testSetdata() throws Exception {
//
//    }
//
//    @Test
//    public void testWshift() throws Exception {
//
//    }

}

package org.jdoris.core.utils;

import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

public class SpectralUtilsTest {

    private static DoubleMatrix realMatrix_EXPECTED;
    private static ComplexDoubleMatrix complexMatrix_EXPECTED;

    private static ComplexDoubleMatrix fftVector_EXPECTED;


    private static ComplexDoubleMatrix fftMatrix_dim1_EXPECTED;
    private static ComplexDoubleMatrix fftMatrix_dim2_EXPECTED;


    @BeforeClass
    public static void setUpTestData() {

        realMatrix_EXPECTED = new DoubleMatrix(new double[][]{{2, 4, 8, 16}, {4, 4, 8, 16}, {8, 8, 8, 16}, {16, 16, 16, 16}});
        complexMatrix_EXPECTED = new ComplexDoubleMatrix(realMatrix_EXPECTED, realMatrix_EXPECTED);

        fftMatrix_dim1_EXPECTED = new ComplexDoubleMatrix(
                new DoubleMatrix(new double[][]{{30, 32, 40, 64}, {-18, -16, -8, 0}, {-10, -8, -8, 0}, {6, 8, 8, 0}}),
                new DoubleMatrix(new double[][]{{30, 32, 40, 64}, {6, 8, 8, 0}, {-10, -8, -8, 0}, {-18, -16, -8, 0}}));

        fftMatrix_dim2_EXPECTED = new ComplexDoubleMatrix(
                new DoubleMatrix(new double[][]{{30, -18, -10, 6}, {32, -16, -8, 8}, {40, -8, -8, 8}, {64, 0, 0, 0}}),
                new DoubleMatrix(new double[][]{{30, 6, -10, -18}, {32, 8, -8, -16}, {40, 8, 8, -8}, {64, 0, 0, 0}}));

//        DoubleMatrix real = new DoubleMatrix(new double[]{30, -18, -10, 6});
//        DoubleMatrix complex = new DoubleMatrix(new double[]{30, 6, -10, -18});
//        fftVector_EXPECTED = new ComplexDoubleMatrix(real, complex);

    }

    @Test
    public void testFourier1D_inplace() throws Exception {

        ComplexDoubleMatrix tempVectorMatrix = complexMatrix_EXPECTED.getColumn(0);
        SpectralUtils.fft1D_inplace(tempVectorMatrix, tempVectorMatrix.length);

        Assert.assertEquals(fftMatrix_dim1_EXPECTED.getColumn(0), tempVectorMatrix);

       SpectralUtils.invfft1D_inplace(tempVectorMatrix, tempVectorMatrix.length);
        Assert.assertEquals(complexMatrix_EXPECTED.getColumn(0), tempVectorMatrix);

    }

    @Test
    public void testFourier1D() throws Exception {

        ComplexDoubleMatrix tempVectorMatrix = complexMatrix_EXPECTED.getColumn(0);
        ComplexDoubleMatrix fftVector_EXPECTED = SpectralUtils.fft1D(tempVectorMatrix, tempVectorMatrix.length);

        Assert.assertEquals(fftMatrix_dim1_EXPECTED.getColumn(0), fftVector_EXPECTED);

        ComplexDoubleMatrix complexVector_EXPECTED = SpectralUtils.invfft1D(fftVector_EXPECTED, tempVectorMatrix.length);
        Assert.assertEquals(complexMatrix_EXPECTED.getColumn(0), complexVector_EXPECTED);

    }

    @Test
    public void testFft_inplace() throws Exception {

        ComplexDoubleMatrix tempMatrix = complexMatrix_EXPECTED;
        SpectralUtils.fft_inplace(tempMatrix, 1);

        Assert.assertEquals(fftMatrix_dim1_EXPECTED, tempMatrix);

        SpectralUtils.invfft_inplace(tempMatrix, 1);
        Assert.assertEquals(complexMatrix_EXPECTED, tempMatrix);

    }

    @Test
    public void testFft() throws Exception {

        ComplexDoubleMatrix fftMatrix_dim1_ACTUAL = SpectralUtils.fft(complexMatrix_EXPECTED, 1);
        Assert.assertEquals(fftMatrix_dim1_EXPECTED, fftMatrix_dim1_ACTUAL);

        ComplexDoubleMatrix ifftMatrix_dim_ACTUAL = SpectralUtils.invfft(fftMatrix_dim1_ACTUAL, 1);
        Assert.assertEquals(complexMatrix_EXPECTED, ifftMatrix_dim_ACTUAL);

    }

//    @Test
//    public void testIfftshift() throws Exception {
//
//    }
//
//    @Test
//    public void testFft2d() throws Exception {
//
//    }
//
//    @Test
//    public void testFft2d() throws Exception {
//
//    }
//
//    @Test
//    public void testIfft2d() throws Exception {
//
//    }
}

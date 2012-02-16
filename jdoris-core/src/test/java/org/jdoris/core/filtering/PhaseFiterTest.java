package org.jdoris.core.filtering;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;
import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.FloatMatrix;
import org.jdoris.core.utils.SpectralUtils;
import org.jdoris.core.utils.WeightWindows;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.jdoris.core.io.DataReader.readCplxFloatData;
import static org.jdoris.core.io.DataReader.readFloatData;

/**
 * User: pmar@ppolabs.com
 * Date: 6/8/11
 * Time: 12:04 PM
 */
public class PhaseFiterTest {

    private static final String testDirectoryGoldstein = "/d2/etna_test/phaseFiltTest/goldstein/";
    private static final String testDirectorySpectral = "/d2/etna_test/phaseFiltTest/spectral/";
    private static final String testDirectorySpatial = "/d2/etna_test/phaseFiltTest/spatialconv/";
    private static final double DELTA_01 = 1e-01;
    private static final double DELTA_03 = 1e-03;
    private static final double DELTA_005 = 5e-01;
    private DoubleMatrix kernel2d;

    private static Logger initLog() {
        String filePathToLog4JProperties = "log4j.properties";
        Logger logger = Logger.getLogger(PhaseFiter.class);
        PropertyConfigurator.configure(filePathToLog4JProperties);
        return logger;
    }

    @BeforeClass
    public static void setUp() {
        initLog();
    }


    @Test
    public void testGoldstein_no_smoothing() throws Exception {

        // load input data
        // master data
        String cplxDataFileName = testDirectoryGoldstein + "data_input.cr4.swap";
        ComplexDoubleMatrix cplxData = readCplxFloatData(cplxDataFileName, 32, 512);

        String cplxDataFilteredFileName = testDirectoryGoldstein + "data_filtered.cr4.swap";
        ComplexDoubleMatrix cplxDataFilt_EXPECTED = readCplxFloatData(cplxDataFilteredFileName, 32, 512);

        double[] smoothKernel = null;
        float alpha = (float) 0.5;
        int overlap = 12;
        ComplexDoubleMatrix cplxDataFilt_ACTUAL = PhaseFiter.goldstein(cplxData, alpha, overlap, smoothKernel);

        Assert.assertArrayEquals(cplxDataFilt_EXPECTED.real().toArray(), cplxDataFilt_ACTUAL.real().toArray(), DELTA_01);
        Assert.assertArrayEquals(cplxDataFilt_EXPECTED.imag().toArray(), cplxDataFilt_ACTUAL.imag().toArray(), DELTA_01);

    }

    @Test
    public void testGoldstein_with_smoothing() throws Exception {

        // load input data
        // master data
        String cplxDataFileName = testDirectoryGoldstein + "data_input.cr4.swap";
        ComplexDoubleMatrix cplxData = readCplxFloatData(cplxDataFileName, 32, 512);

        String cplxDataFilteredFileName = testDirectoryGoldstein + "data_filtered_smooth.cr4.swap";
        ComplexDoubleMatrix cplxDataFilt_EXPECTED = readCplxFloatData(cplxDataFilteredFileName, 32, 512);

        double[] smoothKernel = {0.2, 0.2, 0.2, 0.2, 0.2};
        float alpha = (float) 0.5;
        int overlap = 12;
        ComplexDoubleMatrix cplxDataFilt_ACTUAL = PhaseFiter.goldstein(cplxData, alpha, overlap, smoothKernel);

        Assert.assertArrayEquals(cplxDataFilt_EXPECTED.real().toArray(), cplxDataFilt_ACTUAL.real().toArray(), DELTA_01);
        Assert.assertArrayEquals(cplxDataFilt_EXPECTED.imag().toArray(), cplxDataFilt_ACTUAL.imag().toArray(), DELTA_01);

    }

    @Before
    public void setUpFilterKernel2d() {

        final int blockSize = 32;
        final DoubleMatrix ones = DoubleMatrix.ones(1, blockSize);
        final DoubleMatrix hammingWind = new DoubleMatrix(WeightWindows.hamming(blockSize));
        kernel2d = hammingWind.mmul(ones);
        kernel2d = kernel2d.mul(kernel2d.transpose());
    }

    @Test
    public void TestArrangingKernel() throws Exception {

        final int scalefactor = 1;
        DoubleMatrix arrangedKernel2d_ACTUAL = PhaseFiter.arrangeKernel2d(kernel2d, scalefactor);

        String kernel2dFileName = testDirectorySpectral + "kernel2d.r4.swap";
        FloatMatrix arrangedKernel2d_EXPECTED = readFloatData(kernel2dFileName, 32, 32);

        Assert.assertArrayEquals(arrangedKernel2d_EXPECTED.toArray(), arrangedKernel2d_ACTUAL.toFloat().toArray(), (float) DELTA_01);

    }

    @Test
    public void testSpectral() throws Exception {

        // load input data
        String cplxDataFileName = testDirectorySpectral + "data_input.cr4.swap";
        ComplexDoubleMatrix cplxData = readCplxFloatData(cplxDataFileName, 32, 512);

        String cplxDataFilteredFileName = testDirectorySpectral + "data_filtered.cr4.swap";
        ComplexDoubleMatrix cplxDataFilt_EXPECTED = readCplxFloatData(cplxDataFilteredFileName, 32, 512);

        int overlap = 12;
        ComplexDoubleMatrix cplxDataFilt_ACTUAL = PhaseFiter.spectralfilt(cplxData, kernel2d, overlap);

        Assert.assertArrayEquals(cplxDataFilt_EXPECTED.real().toArray(), cplxDataFilt_ACTUAL.real().toArray(), DELTA_005);
        Assert.assertArrayEquals(cplxDataFilt_EXPECTED.imag().toArray(), cplxDataFilt_ACTUAL.imag().toArray(), DELTA_005);

    }


    @Test
    public void testSpatialConvFilter() throws Exception {

        String kernelFileName = testDirectorySpatial + "kernel2d.cr4.swap";
        ComplexDoubleMatrix kernel2d_EXPECTED = readCplxFloatData(kernelFileName, 128, 128);

        final double[] kernel = {0.2, 0.2, 0.2, 0.2, 0.2};
        DoubleMatrix kernel2d = PhaseFiter.defineRectKernel2d(128, kernel);

        ComplexDoubleMatrix kernel2d_ACTUAL = new ComplexDoubleMatrix(kernel2d);
        SpectralUtils.fft2D_inplace(kernel2d_ACTUAL);
        kernel2d_ACTUAL.conji();

        Assert.assertArrayEquals(kernel2d_EXPECTED.real().toArray(), kernel2d_ACTUAL.real().toArray(), DELTA_03);
        Assert.assertArrayEquals(kernel2d_EXPECTED.imag().toArray(), kernel2d_ACTUAL.imag().toArray(), DELTA_03);

    }

    @Test
    public void testSpatialConv() throws Exception {

        // load input data
        String cplxDataFileName = testDirectorySpatial + "data_input.cr4.swap";
        ComplexDoubleMatrix cplxData = readCplxFloatData(cplxDataFileName, 128, 512);

        String cplxDataFilteredFileName = testDirectorySpatial + "data_filtered.cr4.swap";
        ComplexDoubleMatrix cplxDataFilt_EXPECTED = readCplxFloatData(cplxDataFilteredFileName, 128, 512);

        final double[] kernel = {0.2, 0.2, 0.2, 0.2, 0.2};
        DoubleMatrix kernel2d = PhaseFiter.defineRectKernel2d(128, kernel);

        ComplexDoubleMatrix kernel2d_ACTUAL = new ComplexDoubleMatrix(kernel2d);
        SpectralUtils.fft2D_inplace(kernel2d_ACTUAL);
        kernel2d_ACTUAL.conji();

        int overlap = 2;

        ComplexDoubleMatrix cplxDataFilt_ACTUAL = PhaseFiter.convbuffer(cplxData, kernel2d_ACTUAL, overlap);

        Assert.assertArrayEquals(cplxDataFilt_EXPECTED.real().toArray(), cplxDataFilt_ACTUAL.real().toArray(), DELTA_01);
        Assert.assertArrayEquals(cplxDataFilt_EXPECTED.imag().toArray(), cplxDataFilt_ACTUAL.imag().toArray(), DELTA_01);

    }



}

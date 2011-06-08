package org.jdoris.core.filtering;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;
import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.jdoris.core.io.DataReader.readCplxFloatData;

/**
 * User: pmar@ppolabs.com
 * Date: 6/8/11
 * Time: 12:04 PM
 */
public class PhaseFiterTest {

    private static final String testDirectory = "/d2/etna_test/phaseFiltTest/";
    private static final double DELTA = 1e-01;

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
    public void testGoldstein() throws Exception {

        // load input data
        // master data
        String cplxDataFileName = testDirectory + "data_input.cr4.swap";
        ComplexDoubleMatrix cplxData = readCplxFloatData(cplxDataFileName, 32, 512);

        String cplxDataFilteredFileName = testDirectory + "data_filtered.cr4.swap";
        ComplexDoubleMatrix cplxDataFilt_EXPECTED = readCplxFloatData(cplxDataFilteredFileName, 32, 512);

        DoubleMatrix smoothKernel = null;
        float alpha = (float) 0.5;
        int overlap = 12;
        ComplexDoubleMatrix cplxDataFilt_ACTUAL = PhaseFiter.goldstein(cplxData, alpha, overlap, smoothKernel);

        Assert.assertArrayEquals(cplxDataFilt_EXPECTED.real().toArray(), cplxDataFilt_ACTUAL.real().toArray(), DELTA);
        Assert.assertArrayEquals(cplxDataFilt_EXPECTED.imag().toArray(), cplxDataFilt_ACTUAL.imag().toArray(), DELTA);

    }
}

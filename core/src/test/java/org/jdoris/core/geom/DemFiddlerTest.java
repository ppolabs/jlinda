package org.jdoris.core.geom;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;
import org.jdoris.core.Orbit;
import org.jdoris.core.SLCImage;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.File;

import static org.jdoris.core.io.DataReader.readDoubleData;

public class DemFiddlerTest {

    static Logger logger = Logger.getLogger(DemFiddler.class.getName());
    private static final File resFile = new File("/d2/etna_test/master_temp.res");

    static DemFiddler dem;

    static SLCImage metadata;
    static Orbit orbit;

    private static final double DELTA_08 = 1e-08;
    private static final double DELTA_06 = 1e-06;

    public static Logger initLog() {
        String filePathToLog4JProperties = "log4j.properties";
        Logger logger = Logger.getLogger(DemFiddler.class);
        PropertyConfigurator.configure(filePathToLog4JProperties);
        return logger;
    }

    @Before
    public void setUp() throws Exception {

        initLog();

        // initialize
        dem = new DemFiddler();

        //// Tile corners : radar coordinates
        dem.l0 = 10000;
        dem.lN = 10127;
        dem.p0 = 1500;
        dem.pN = 2011;

        // extras [rad]
        dem.extraLat = 0.0028143434188408478;
        dem.extraLon = 0.0028143434188408478;

        //// dem params
        // spacing in radians
        dem.demDeltaLat = 1.4544410433280261e-05;
        dem.demDeltaLon = 1.4544410433280261e-05;
        // dem size in pixels
        dem.nLatPixels = 3601;
        dem.nLonPixels = 3601;
        dem.lat0 = 0.68067840827778847;
        dem.lon0 = 0.24434609527920614;

        // initialize metadata
        metadata = new SLCImage();
        metadata.parseResFile(resFile);

        orbit = new Orbit();
        orbit.parseOrbit(resFile);
        orbit.computeCoefficients(3);

    }

    @Test
    public void testGetDEMCorners() throws Exception {

        dem.getDEMCorners(metadata, orbit);

        double phiMin_EXPECTED = 0.65198531114095126;
        Assert.assertEquals(phiMin_EXPECTED, dem.phiMin, DELTA_08);

        double phiMax_EXPECTED = 0.65803878531122906;
        Assert.assertEquals(phiMax_EXPECTED, dem.phiMax, DELTA_08);

        double lambdaMin_EXPECTED = 0.2584835617746124;
        Assert.assertEquals(lambdaMin_EXPECTED, dem.lambdaMin, DELTA_08);

        double lambdaMax_EXPECTED = 0.26619645946965714;
        Assert.assertEquals(lambdaMax_EXPECTED, dem.lambdaMax, DELTA_08);

        int indexPhi0DEM_EXPECTED = 1556;
        Assert.assertEquals(indexPhi0DEM_EXPECTED, dem.indexPhi0DEM, DELTA_08);

        int indexPhiNDEM_EXPECTED = 1973;
        Assert.assertEquals(indexPhiNDEM_EXPECTED, dem.indexPhiNDEM, DELTA_08);

        int indexLambda0DEM_EXPECTED = 972;
        Assert.assertEquals(indexLambda0DEM_EXPECTED, dem.indexLambda0DEM, DELTA_08);

        int indexLambdaNDEM_EXPECTED = 1503;
        Assert.assertEquals(indexLambdaNDEM_EXPECTED, dem.indexLambdaNDEM, DELTA_08);

    }

    @Test
    public void testGridData() throws Exception {

        // define params
        double firstline_buffer = 10000;
        double lastline_buffer = 10127;
        double first_pixel = 1500;
        double lastpixel = 2011;
        double mlL = 1;
        double mlP = 1;
        double r_az_ratio = 5.2487532186594095;
        float offset = 0;
        final float NODATA = -32768;

        final int nRows = 418;
        final int nCols = 532;
        double[][] DEMline_buffer = new double[nRows][nCols];  // dem_line_buffer.r4
        double[][] DEMpixel_buffer = new double[nRows][nCols]; // dem_pixel_buffer.r4
        double[][] input_buffer = new double[nRows][nCols];  // dem_buffer.r4

        double[][] grd_EXPECTED; // output_buffer.r4

        /* grid input tile */
        DemFiddler demFiddler = new DemFiddler();
        demFiddler.setGrd(new double[128][512]);

        // load test data
        String testDataDir = "/d2/etna_test/demTest/";
        String bufferFileName;

        bufferFileName = testDataDir + "dem_line_buffer.r8.swap";
        DEMline_buffer = readDoubleData(bufferFileName, nRows, nCols).toArray2();

        bufferFileName = testDataDir + "dem_pixel_buffer.r8.swap";
        DEMpixel_buffer = readDoubleData(bufferFileName, nRows, nCols).toArray2();

        bufferFileName = testDataDir + "dem_buffer.r8.swap";
        input_buffer = readDoubleData(bufferFileName, nRows, nCols).toArray2();

        bufferFileName = testDataDir + "output_buffer.r8.swap";
        grd_EXPECTED = readDoubleData(bufferFileName, 128, 512).toArray2();

        /* computation */


        long t0 = System.currentTimeMillis();
        demFiddler.gridData(DEMline_buffer, DEMpixel_buffer, input_buffer,
                firstline_buffer, first_pixel,
                mlL, mlP, r_az_ratio,
                offset, NODATA);
        long t1 = System.currentTimeMillis();
        logger.info("Data set gridded in " + (0.001 * (t1 - t0)) + " sec");


        /* assert result */
        for (int i = 0; i < grd_EXPECTED.length; i++) {
            double[] doubles = grd_EXPECTED[i];
            Assert.assertArrayEquals(doubles, demFiddler.grd[i], DELTA_06);
        }

    }
}

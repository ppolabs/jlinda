package org.jdoris.core.geom;

import org.jdoris.core.Orbit;
import org.jdoris.core.SLCImage;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.File;

/**
 * User: pmar@ppolabs.com
 * Date: 6/14/11
 * Time: 12:47 PM
 */
public class DemFiddlerTest {

    private static final File resFile = new File("/d2/etna_test/master_temp.res");

    static DemFiddler dem;

    static SLCImage metadata;
    static Orbit orbit;

    private static double DELTA_08 = 1e-08;

    @Before
    public void setUp() throws Exception {

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
        Assert.assertEquals(phiMax_EXPECTED,dem.phiMax, DELTA_08);

        double lambdaMin_EXPECTED = 0.2584835617746124;
        Assert.assertEquals(lambdaMin_EXPECTED,dem.lambdaMin, DELTA_08);

        double lambdaMax_EXPECTED = 0.26619645946965714;
        Assert.assertEquals(lambdaMax_EXPECTED,dem.lambdaMax, DELTA_08);

        int indexPhi0DEM_EXPECTED = 1556;
        Assert.assertEquals(indexPhi0DEM_EXPECTED,dem.indexPhi0DEM, DELTA_08);

        int indexPhiNDEM_EXPECTED = 1973;
        Assert.assertEquals(indexPhiNDEM_EXPECTED,dem.indexPhiNDEM, DELTA_08);

        int indexLambda0DEM_EXPECTED = 972;
        Assert.assertEquals(indexLambda0DEM_EXPECTED,dem.indexLambda0DEM, DELTA_08);

        int indexLambdaNDEM_EXPECTED = 1503;
        Assert.assertEquals(indexLambdaNDEM_EXPECTED,dem.indexLambdaNDEM, DELTA_08);

    }
}

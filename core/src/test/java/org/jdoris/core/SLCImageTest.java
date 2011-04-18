package org.jdoris.core;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.text.ParseException;

public class SLCImageTest {

    private static File resFile;
    private static SLCImage master = new SLCImage();
    private static final double delta = 1E-05;
    private static final double delta_1 = 1E-04;
    private static final int lineExpected = 16475;
    private static final double timeAzExpected = 38199.5916316067;
    private static final int pixelExpected = 3151;
    private static final double timeRgExpected = 0.00286525698652786;
    private static final double dopplerExpected = 933.516721989199;
    @BeforeClass
    public static void setUp() throws ParseException {
        resFile = new File("/d2/delft_cr_asar.res");
        master.parseResFile(resFile);
    }


//    @Test
//    public void testParseResFile() throws Exception {
//
//
//
//    }

    @Test
    public void testPix2tr() throws Exception {
        double timeRgActual = master.pix2tr(pixelExpected);
        Assert.assertEquals(timeRgExpected,timeRgActual,delta);
    }

    @Test
    public void testTr2pix() throws Exception {
        double pixelActual = master.tr2pix(timeRgExpected);
        Assert.assertEquals(pixelExpected,pixelActual,delta);

    }

    @Test
    public void testPix2fdc() throws Exception {
        double dopplerActual = master.pix2fdc(pixelExpected);
        System.out.println("dopplerActual = " + dopplerActual);
        Assert.assertEquals(dopplerExpected,dopplerActual,delta);
    }

    @Test
    public void testLine2ta() throws Exception {
        double timeAzActual = master.line2ta(lineExpected);
        Assert.assertEquals(timeAzExpected, timeAzActual, delta);

    }


    @Test
    public void testTa2line() throws Exception {
        double lineActual = master.ta2line(timeAzExpected);
        Assert.assertEquals(lineExpected, lineActual, delta_1);

    }
}

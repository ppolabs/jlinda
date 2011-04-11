package org.jdoris.core;

import org.esa.beam.framework.datamodel.GeoPos;
import org.esa.beam.framework.datamodel.MetadataElement;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.nest.datamodel.AbstractMetadata;
import org.esa.nest.util.Constants;
import org.esa.nest.util.GeoUtils;

/**
 * User: pmar@ppolabs.com
 * Date: 2/17/11
 * Time: 12:40 PM
 */
public final class SLCImage {

    // file & format
    private static String fileName;
    private static int formatFlag; // not used

    // sensor
    private static String sensor;
    private static String sarProcessor;
    private static double radar_wavelength; // TODO: close this modifier

    // geo & orientation
    private static Point approxRadarCentreOriginal; // use PixelPos as double!
    private static GeoPos approxGeoCentreOriginal;
    private static Point approxXYZCentreOriginal;

    private static double averageHeight;

    // azimuth annotations
    private static double PRF;
    private static double azimuthBandwidth;
    private static double tAzi1;

    // range annotations
    private static double rsr2x;
    private static double rangeBandwidth;
    private static double tRange1;

    //    // doppler
//    private static double[] f_DC; // TODO
    private static double f_DC_a0;                // constant term Hz
    private static double f_DC_a1;                // linear term Hz/s
    private static double f_DC_a2;                // quadratic term Hz/s/s

    // ______ offset = X(l,p) - X(L,P) ______
    // ______ Where l,p are in the local slave coordinate system and ______
    // ______ where L,P are in the local master coordinate system ______
    // ______ These variables are stored in the slaveinfo variable only ______
    private static int coarseOrbitOffsetL;     // orbit offset in line (azimuth) direction
    private static int coarseOrbitOffsetP;     // orbit offset in pixel (range) direction
    private static int coarseOffsetL;          // offset in line (azimuth) direction
    private static int coarseOffsetP;          // offset in pixel (range) direction

    // oversampling factors
    private static int ovsAz;                 // oversampling of SLC
    private static int ovsRg;                 // oversampling of SLC

    // timing errors
    private static int azTimingError;        // timing error in azimuth direction
    // relative to master geometry, or
    // absolute timing error of master
    // units: lines

    private static int rgTimingError;        // timing error in range direction
    // relative to master geometry, or
    // absolute timing error of master
    // units: pixels

    private static boolean absTimingErrorFlag;   // FALSE if master time is NOT updated,
    // true if it is

    //    private static Rectangle originalWindow;       // position and size of the full scene
    Window originalWindow;       // position and size of the full scene
    Window currentWindow;        // position and size of the subset
    Window slaveMasterOffsets;   // overlapping slave window in master coordinates


    public SLCImage() {

        sensor = "SLC_ERS";                    // default (vs. SLC_ASAR, JERS, RSAT)
        sarProcessor = "SARPR_VMP";            // (VMP (esa paf) or ATLANTIS or TUDELFT) // TODO PGS update?
        formatFlag = 0;                        // format of file on disk

        approxXYZCentreOriginal.x = 0.0;
        approxXYZCentreOriginal.y = 0.0;
        approxXYZCentreOriginal.z = 0.0;

        radar_wavelength = 0.0565646;          // [m] default ERS2
        tAzi1 = 0.0;                           // [s] sec of day
        tRange1 = 5.5458330 / 2.0e3;           // [s] one way, default ERS2
        PRF = 1679.902;                        // [Hz] default ERS2
        azimuthBandwidth = 1378.0;             // [Hz] default ERS2
        f_DC_a0 = 0.0;                         // [Hz] default ERS2
        f_DC_a1 = 0.0;
        f_DC_a2 = 0.0;
        rsr2x = 18.9624680 * 2.0e6;            // [Hz] default ERS2
        rangeBandwidth = 15.55e6;              // [Hz] default ERS2

        coarseOffsetL = 0;                     // by default
        coarseOffsetP = 0;                     // by default
        coarseOrbitOffsetL = 0;                // by default
        coarseOrbitOffsetP = 0;                // by default

        ovsRg = 1;                             // by default
        ovsAz = 1;                             // by default

        absTimingErrorFlag = false;
        azTimingError = 0;                     // by default, unit lines
        rgTimingError = 0;                     // by default, unit pixels

        currentWindow  = new Window(1, 25000, 1, 5000);
        originalWindow = new Window(1, 25000, 1, 5000);
//        slavemasteroffsets.l00  = 0;               // window in master coordinates
//        slavemasteroffsets.p00  = 0;
//        slavemasteroffsets.l0N  = 0;
//        slavemasteroffsets.p0N  = 0;
//        slavemasteroffsets.lN0  = 0;
//        slavemasteroffsets.pN0  = 0;
//        slavemasteroffsets.lNN  = 0;
//        slavemasteroffsets.pNN  = 0;
    }


    public SLCImage(MetadataElement element) {

        // units [meters]
        radar_wavelength = (Constants.lightSpeed / Math.pow(10, 6)) / element.getAttributeDouble(AbstractMetadata.radar_frequency);

        // units [Hz]
        PRF = element.getAttributeDouble(AbstractMetadata.pulse_repetition_frequency);

        // work with seconds of the day!
        ProductData.UTC t_azi1_UTC = element.getAttributeUTC(AbstractMetadata.first_line_time);
        tAzi1 = (t_azi1_UTC.getMJD() - (int) t_azi1_UTC.getMJD()) * 24 * 3600;

        // 2 times range sampling rate [HZ]
        rsr2x = (element.getAttributeDouble(AbstractMetadata.range_sampling_rate) * Math.pow(10, 6) * 2);

        // one way (!!!) time to first range pixels [sec]
        tRange1 = element.getAttributeDouble(AbstractMetadata.slant_range_to_first_pixel) / Constants.lightSpeed;

        approxRadarCentreOriginal.x = element.getAttributeDouble(AbstractMetadata.num_samples_per_line) / 2.0d;  // x direction is range!
        approxRadarCentreOriginal.y = element.getAttributeDouble(AbstractMetadata.num_output_lines) / 2.0d;  // y direction is azimuth

        // TODO: replace computation of the centre using getGeoPos()
        // simple averaging of the corners : as approximation accurate enough
        approxGeoCentreOriginal.lat = (float) ((element.getAttributeDouble(AbstractMetadata.first_near_lat) +
                element.getAttributeDouble(AbstractMetadata.first_far_lat) +
                element.getAttributeDouble(AbstractMetadata.last_near_lat) +
                element.getAttributeDouble(AbstractMetadata.last_far_lat)) / 4);

        approxGeoCentreOriginal.lon = (float) ((element.getAttributeDouble(AbstractMetadata.first_near_long) +
                element.getAttributeDouble(AbstractMetadata.first_far_long) +
                element.getAttributeDouble(AbstractMetadata.last_near_long) +
                element.getAttributeDouble(AbstractMetadata.last_far_long)) / 4);

        double[] xyz = new double[3];
        GeoUtils.geo2xyz(getApproxGeoCentreOriginal(), xyz);

        approxXYZCentreOriginal.x = xyz[0];
        approxXYZCentreOriginal.y = xyz[1];
        approxXYZCentreOriginal.z = xyz[2];

        // set dopplers
        final AbstractMetadata.DopplerCentroidCoefficientList[] dopplersArray = AbstractMetadata.getDopplerCentroidCoefficients(element);

        // TODO: check correctness of this!!
        f_DC_a0 = dopplersArray[1].coefficients[0];
        f_DC_a1 = dopplersArray[1].coefficients[1];
        f_DC_a2 = dopplersArray[1].coefficients[2];

    }

    /*---  RANGE CONVERSIONS ----*/

    // Convert pixel number to range time (1 is first pixel)
    public double pix2tr(double pixel) {
        return tRange1 + ((pixel - 1.0) / rsr2x);
    }

    // Convert pixel number to range (1 is first pixel)
    public double pix2range(double pixel) {
        return Constants.lightSpeed * pix2tr(pixel);
    }

    // Convert range time to pixel number (1 is first pixel)
    public double tr2pix(double rangeTime) {
        return 1.0 + (rsr2x * (rangeTime - tRange1));
    }

    // Convert range pixel to fDC (1 is first pixel, can be ovs)
    public double pix2fdc(double pixel) {
        double tau = (pixel - 1.0) / (rsr2x / 2.0);// two-way time
        return f_DC_a0 + (f_DC_a1 * tau) + (f_DC_a2 * Math.sqrt(tau));
    }

    /*---  AZIMUTH CONVERSIONS ----*/

    // Convert line number to azimuth time (1 is first line)
    public double line2ta(double line) {
        return tAzi1 + ((line - 1.0) / PRF);
    }

    // Convert azimuth time to line number (1 is first line)
    public double ta2line(double azitime) {
        return 1.0 + PRF * (azitime - tAzi1);
    }


    /* Getters and setters for Encapsulation */

    public double getRadarWavelength() {
        return radar_wavelength;
    }

    public Point getApproxRadarCentreOriginal() {
        return approxRadarCentreOriginal;
    }

    public GeoPos getApproxGeoCentreOriginal() {
        return approxGeoCentreOriginal;
    }

    public Point getApproxXYZCentreOriginal() {
        return approxXYZCentreOriginal;
    }

    public Window getCurrentWindow() {
        return currentWindow;
    }

    public double getRsr2x() {
        return rsr2x;
    }

    public double getPRF() {
        return PRF;
    }

    public double getAzimuthBandwidth() {
        return azimuthBandwidth;
    }

    public double getF_DC_a0() {
        return f_DC_a0;
    }

    public double getF_DC_a1() {
        return f_DC_a1;
    }

    public double getF_DC_a2() {
        return f_DC_a2;
    }

    public int getCoarseOffsetP() {
        return coarseOffsetP;
    }

}

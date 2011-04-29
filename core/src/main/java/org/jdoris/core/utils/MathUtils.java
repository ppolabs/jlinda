package org.jdoris.core.utils;

import java.awt.*;

public class MathUtils {

//    static Logger logger = Logger.getLogger(MathUtils.class.getName());

    public static double[][] distributePoints(int numOfPoints, Rectangle rectangle) {

        double lines = (rectangle.getMaxY() - rectangle.getMinY() + 1);
        double pixels = (rectangle.getMaxX() - rectangle.getMinX() + 1);

        double[][] result = new double[numOfPoints][2];

        // ______ Distribution for dl=dp ______
        double wp = Math.sqrt(numOfPoints / (lines / pixels));   // wl: #windows in line direction
        double wl = numOfPoints / wp;                   // wp: #windows in pixel direction
        if (wl < wp) {
            // switch wl,wp : later back
            wl = wp;
        }

        double wlint = Math.ceil(wl); // round largest
        double deltal = (lines - 1) / (wlint - 1);
        double totp = Math.ceil(pixels * wlint);
        double deltap = (totp - 1) / (numOfPoints - 1);
        double p = -deltap;
        double l = 0.;
        double lcnt = 0;
        int i;
        for (i = 0; i < numOfPoints; i++) {
            p += deltap;
            while (Math.ceil(p) >= pixels) // ceil
            {
                p -= pixels;
                lcnt++;
            }
            l = lcnt * deltal;

            result[i][0] = (int) Math.ceil(l);
            result[i][1] = (int) Math.ceil(p);
        }

//        // ______ Correct distribution to window ______
//        for (i=0; i<numOfPoints; i++){
//            result[i][0] += (int)rectangle.getMinY();
//            result[i][1] += (int)rectangle.getMinX();
//        }
        return result;
    }


    public static boolean isodd(long value) {
        return value % 2 == 0;
    }

    public static boolean ispower2(long value) {
        return value == 1 || value == 2 || value == 4 || value == 8 || value == 16 ||
                value == 32 || value == 64 || value == 128 || value == 256 ||
                value == 512 || value == 1024 || value == 2048 || value == 4096;
    }


    public static double rad2deg(double rad) {
        return Math.toDegrees(rad);
    }

    public static double sqr(double value) {
        return Math.pow(value, 2);
    }

    public static double sqrt(double value) {
        return Math.sqrt(value);
    }

}


package org.jdoris.core.todo_classes;

import org.apache.log4j.Logger;
import org.jdoris.core.Point;

public class input_ell {

    Logger logger = Logger.getLogger(input_ell.class.getName());

    private double e2;  // squared first  eccentricity (derived)
    private double e2b; // squared second eccentricity (derived)

    public double a; // semi major
    public double b; // semi minor
    public String name;

    // first ecc.
    private void set_ecc1st_sqr() {
        //  faster than e2=(sqr(a)-sqr(b))/sqr(a)
        e2 = 1.0 - Math.pow(b / a, 2);
    }

    // second ecc.
    private void set_ecc2nd_sqr() {
        // faster than e2b=(sqr(a)-sqr(b))/sqr(b);
        e2b = Math.pow(a / b, 2) - 1.0;
    }

    public input_ell() {
        a = Constants.WGS84_A;
        b = Constants.WGS84_B;
        e2 = 0.00669438003551279091;
        e2b = 0.00673949678826153145;
        //set_ecc1st_sqr();// compute e2
        //set_ecc2nd_sqr();// compute e2b
        name = "WGS84";
    }


    public input_ell(final double semimajor, final double semiminor) {
        a = semimajor;
        b = semiminor;
        set_ecc1st_sqr();// compute e2 (not required for zero-doppler iter.)
        set_ecc2nd_sqr();// compute e2b (not required for zero-doppler iter.)
        //set_name("unknown");

    }

    public input_ell(input_ell ell) {
        a = ell.a;
        b = ell.b;
        e2 = ell.e2;
        e2b = ell.e2b;
        name = ell.name;
    }

    public void showdata() {
        logger.info("ELLIPSOID: \tEllipsoid used (orbit, output): " + name + ".");
        logger.info("ELLIPSOID: a   = " + a);
        logger.info("ELLIPSOID: b   = " + b);
        logger.info("ELLIPSOID: e2  = " + e2);
        logger.info("ELLIPSOID: e2' = " + e2b);
    }

    /*
    *  Convert xyz cartesian coordinates to
    *  Geodetic ellipsoid coordinates latlonh
    *    xyz2ell                                                   *
    *                                                              *
    * Converts geocentric cartesian coordinates in the XXXX        *
    *  reference frame to geodetic coordinates.                    *
    *  method of bowring see globale en locale geodetische systemen*
    * input:                                                       *
    *  - ellipsinfo, xyz, (phi,lam,hei)                            *
    * output:                                                      *
    *  - void (updated lam<-pi,pi>, phi<-pi,pi>, hei)              *
    *                                                              *
    ****************************************************************/
    public void xyz2ell(final Point xyz, double phi, double lambda, double height) {
        double r = Math.sqrt(Math.pow(xyz.x, 2) + Math.pow(xyz.y, 2));
        double nu = Math.atan2((xyz.z * a), (r * b));
        double sin3 = Math.pow(Math.sin(nu), 3);
        double cos3 = Math.pow(Math.cos(nu), 3);
        phi = Math.atan2((xyz.z + e2b * b * sin3), (r - e2 * a * cos3));
        lambda = Math.atan2(xyz.y, xyz.x);
        double N = a / Math.sqrt(1.0 - e2 * Math.pow(Math.sin(phi), 2));
        height = (r / Math.cos(phi)) - N;
    }

    public void xyz2ell(final Point xyz, double phi, double lambda) {
        double r = Math.sqrt(Math.pow(xyz.x, 2) + Math.pow(xyz.y, 2));
        double nu = Math.atan2((xyz.z * a), (r * b));
        double sin3 = Math.pow(Math.sin(nu), 3);
        double cos3 = Math.pow(Math.cos(nu), 3);
        phi = Math.atan2((xyz.z + e2b * b * sin3), (r - e2 * a * cos3));
        lambda = Math.atan2(xyz.y, xyz.x);
    }


    /**
     * *************************************************************
     * ell2xyz                                                   *
     * *
     * Converts wgs84 ellipsoid cn to geocentric cartesian coord.   *
     * input:                                                       *
     * - phi,lam,hei (geodetic co-latitude, longitude, [rad] h [m] *
     * output:                                                      *
     * - cn XYZ                                                    *
     * *
     * Bert Kampes, 05-Jan-1999                                  *
     * ***************************************************************
     */
    public Point ell2xyz(final double phi, final double lambda, final double height) {
        final double N = a / Math.sqrt(1.0 - e2 * Math.pow(Math.sin(phi), 2));
        final double Nph = N + height;
        return new Point(Nph * Math.cos(phi) * Math.cos(lambda),
                Nph * Math.cos(phi) * Math.sin(lambda),
                (Nph - e2 * N) * Math.sin(phi));
    }

}

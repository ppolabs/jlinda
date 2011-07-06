package org.jdoris.core.utils;

import org.esa.beam.framework.datamodel.GeoPos;
import org.jdoris.core.*;

public class GeoUtils {

    @Deprecated
    public static GeoPos[] computeCorners(final SLCImage meta, final Orbit orbit, final Window tile,
                                                       final double phiExtra, final double lambdaExtra) throws Exception {

        GeoPos[] corners = new GeoPos[2];
        double[] phiAndLambda;

        final double l0 = tile.linelo;
        final double lN = tile.linehi;
        final double p0 = tile.pixlo;
        final double pN = tile.pixhi;

        // compute Phi, Lambda for Tile corners
        phiAndLambda = orbit.lp2ell(new Point(p0, l0), meta);
        final double phi_l0p0 = phiAndLambda[0];
        final double lambda_l0p0 = phiAndLambda[1];

        phiAndLambda = orbit.lp2ell(new Point(p0, lN), meta);
        final double phi_lNp0 = phiAndLambda[0];
        final double lambda_lNp0 = phiAndLambda[1];

        phiAndLambda = orbit.lp2ell(new Point(pN, lN), meta);
        final double phi_lNpN = phiAndLambda[0];
        final double lambda_lNpN = phiAndLambda[1];

        phiAndLambda = orbit.lp2ell(new Point(pN, l0), meta);
        final double phi_l0pN = phiAndLambda[0];
        final double lambda_l0pN = phiAndLambda[1];

        //// Select DEM values based on rectangle outside l,p border ////
        // phi
        double phiMin = Math.min(Math.min(Math.min(phi_l0p0, phi_lNp0), phi_lNpN), phi_l0pN);
        double phiMax = Math.max(Math.max(Math.max(phi_l0p0, phi_lNp0), phi_lNpN), phi_l0pN);
        // lambda
        double lambdaMin = Math.min(Math.min(Math.min(lambda_l0p0, lambda_lNp0), lambda_lNpN), lambda_l0pN);
        double lambdaMax = Math.max(Math.max(Math.max(lambda_l0p0, lambda_lNp0), lambda_lNpN), lambda_l0pN);

        // a little bit extra at edges to be sure
        phiMin -= phiExtra;
        phiMax += phiExtra;
        // lambda
        lambdaMax += lambdaExtra;
        lambdaMin -= lambdaExtra;

        corners[0] = new GeoPos((float) (phiMax * Constants.RTOD), (float) (lambdaMin * Constants.RTOD));
        corners[1] = new GeoPos((float) (phiMin * Constants.RTOD), (float) (lambdaMax * Constants.RTOD));

        return corners;
    }

    public static GeoPos[] computeCorners(final SLCImage meta, final Orbit orbit, final Window tile) throws Exception {

        GeoPos[] corners = new GeoPos[2];
        double[] phiAndLambda;

        final double l0 = tile.linelo;
        final double lN = tile.linehi;
        final double p0 = tile.pixlo;
        final double pN = tile.pixhi;

        // compute Phi, Lambda for Tile corners
        phiAndLambda = orbit.lp2ell(new Point(p0, l0), meta);
        final double phi_l0p0 = phiAndLambda[0];
        final double lambda_l0p0 = phiAndLambda[1];

        phiAndLambda = orbit.lp2ell(new Point(p0, lN), meta);
        final double phi_lNp0 = phiAndLambda[0];
        final double lambda_lNp0 = phiAndLambda[1];

        phiAndLambda = orbit.lp2ell(new Point(pN, lN), meta);
        final double phi_lNpN = phiAndLambda[0];
        final double lambda_lNpN = phiAndLambda[1];

        phiAndLambda = orbit.lp2ell(new Point(pN, l0), meta);
        final double phi_l0pN = phiAndLambda[0];
        final double lambda_l0pN = phiAndLambda[1];

        //// Select DEM values based on rectangle outside l,p border ////
        // phi
        double phiMin = Math.min(Math.min(Math.min(phi_l0p0, phi_lNp0), phi_lNpN), phi_l0pN);
        double phiMax = Math.max(Math.max(Math.max(phi_l0p0, phi_lNp0), phi_lNpN), phi_l0pN);
        // lambda
        double lambdaMin = Math.min(Math.min(Math.min(lambda_l0p0, lambda_lNp0), lambda_lNpN), lambda_l0pN);
        double lambdaMax = Math.max(Math.max(Math.max(lambda_l0p0, lambda_lNp0), lambda_lNpN), lambda_l0pN);

        corners[0] = new GeoPos((float) (phiMax * Constants.RTOD), (float) (lambdaMin * Constants.RTOD));
        corners[1] = new GeoPos((float) (phiMin * Constants.RTOD), (float) (lambdaMax * Constants.RTOD));

        return corners;
    }

    public static GeoPos defineExtraPhiLam_smart() {

        // for maximum height of the tile estimate the overlapping are

        return new GeoPos();
    }


    public static GeoPos defineExtraPhiLam(double latDelta, double lonDelta) {
        // TODO: introduce methods for dynamic scaling of extra lambda/phi depending on average tile Height!
//        lambdaExtra = (1.5 * latDelta + (4.0 / 25.0) * Constants.DTOR); // for himalayas!
//        phiExtra = (1.5 * lonDelta + (4.0 / 25.0) * Constants.DTOR);
        double latExtra = (1.5 * lonDelta + (0.1 / 25.0) * Constants.DTOR);
        double lonExtra = (1.5 * latDelta + (0.1 / 25.0) * Constants.DTOR); // for Etna

        return new GeoPos((float) latExtra, (float) lonExtra);
    }

    public static void extendCorners(GeoPos extra, GeoPos[] coordinates) {

        if (coordinates.length > 2) {
            throw new IllegalArgumentException("More then two corner points in GeoPos array!");
        }

        coordinates[0].lat += (extra.lat * Constants.RTOD);
        coordinates[0].lon -= (extra.lon * Constants.RTOD);

        coordinates[1].lat -= (extra.lat * Constants.RTOD);
        coordinates[1].lon += (extra.lon * Constants.RTOD);
    }

    public static GeoPos defineExtraPhiLam(final double height, final Window window, final SLCImage meta, final Orbit orbit) throws Exception {

        // compute Phi, Lambda for Tile corners
        double[] latLon_ONELLIPS = orbit.lp2ell(new Point(window.pixlo, window.pixhi), meta);
        double[] latLon_HEIGHT = orbit.lph2ell(new Point(window.pixlo, window.pixhi, height), meta);

        float latExtra = (float) Math.abs(latLon_HEIGHT[0] - latLon_ONELLIPS[0]);
        float lonExtra = (float) Math.abs(latLon_HEIGHT[1] - latLon_ONELLIPS[1]);

        return new GeoPos(latExtra, lonExtra);
    }
}

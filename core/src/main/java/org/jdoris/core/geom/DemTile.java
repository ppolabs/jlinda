package org.jdoris.core.geom;

import org.jdoris.core.*;

public class DemTile {

    //// logger
    static org.apache.log4j.Logger logger = org.apache.log4j.Logger.getLogger(DemTile.class.getName());

     //// topoPhase global params
    double lat0;
    double lon0;
    long nLatPixels;
    long nLonPixels;
    double latitudeDelta;
    double longitudeDelta;
    // no data value
    double nodata;

    //// actual demTileData
    double[][] data;

    /// extent : coordinates in radians
    double lambdaExtra;
    double phiExtra;

    /// tile coordinates in phi,lam
    double phiMin;
    double phiMax;
    double lambdaMin;
    double lambdaMax;

    /// tile index
    int indexPhi0DEM;
    int indexPhiNDEM;
    int indexLambda0DEM;
    int indexLambdaNDEM;

    //// tile stats
    boolean statsComputed = false;
    long totalNumPoints;
    long numNodata = 0;
    long numValid;
    double meanValue;
    double minValue;
    double maxValue;


    public DemTile() {
    }

    public DemTile(double lat0, double lon0, long nLatPixels, long nLonPixels,
                   double latitudeDelta, double longitudeDelta, long nodata) {
        this.lat0 = lat0;
        this.lon0 = lon0;
        this.nLatPixels = nLatPixels;
        this.nLonPixels = nLonPixels;
        this.latitudeDelta = latitudeDelta;
        this.longitudeDelta = longitudeDelta;
        this.nodata = nodata;
    }

    public double[][] getData() {
        return data;
    }

    public void setData(double[][] data) {
        this.data = data;
    }

    // get corners of tile (approx) to select DEM
    //	in radians (if height were zero)
    private void setExtraPhiLam() {
        lambdaExtra = (1.5 * latitudeDelta + (4.0 / 25.0) * Constants.DTOR);
        phiExtra = (1.5 * longitudeDelta + (4.0 / 25.0) * Constants.DTOR);
    }

    // TODO: stub for computing statistics of DEM tile
    // ----- Loop over DEM for stats ------------------------
    public void stats() throws Exception{

        // inital values
        double min_dem_buffer = 100000.0;
        double max_dem_buffer = -100000.0;

        try {
            totalNumPoints = data.length * data[0].length;
            for (double[] aData : data) {
                for (int j = 0; j < data[0].length; j++) {
                    if (aData[j] != nodata) {
                        numValid++;
                        meanValue += aData[j];           // divide by numValid later
                        if (aData[j] < min_dem_buffer)
                            min_dem_buffer = aData[j];  //buffer
                        if (aData[j] > max_dem_buffer)
                            max_dem_buffer = aData[j];  // stats
                    } else {
                        numNodata++;
                        System.out.println("dataValue = " + aData[j]);
                    }
                }
            }
        } catch (Exception e) {
            logger.error("Something went wrong when computing DEM tile stats");
            logger.error("Is DEM tile declared?");
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        //global stats
        minValue = Math.min(minValue, min_dem_buffer);
        maxValue = Math.max(maxValue, max_dem_buffer);
        meanValue /= numValid;

        statsComputed = true;
        showStats();

    }

    private void showStats() {

        if (statsComputed) {
            logger.info("DEM Tile Stats");
            logger.info("------------------------------------------------");
            logger.info("Total number of points: " + totalNumPoints);
            logger.info("Number of valid points: " + numValid);
            logger.info("Number of NODATA points: " + numNodata);
            logger.info("Max height in meters at valid points: " + maxValue);
            logger.info("Min height in meters at valid points: " + minValue);
            logger.info("Mean height in meters at valid points: " + meanValue);
        } else {
            logger.warn("DEM Tile Stats");
            logger.warn("------------------------------------------------");
            logger.warn("DemTile.stats() method not invoked!");
        }

    }

    public void computeDemCorners(SLCImage meta, Orbit orbit, Window tile) throws Exception {

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
        phiMin = Math.min(Math.min(Math.min(phi_l0p0, phi_lNp0), phi_lNpN), phi_l0pN);
        phiMax = Math.max(Math.max(Math.max(phi_l0p0, phi_lNp0), phi_lNpN), phi_l0pN);
        // lambda
        lambdaMin = Math.min(Math.min(Math.min(lambda_l0p0, lambda_lNp0), lambda_lNpN), lambda_l0pN);
        lambdaMax = Math.max(Math.max(Math.max(lambda_l0p0, lambda_lNp0), lambda_lNpN), lambda_l0pN);

        // a little bit extra at edges to be sure

        // redefine it: no checks whether there are previous declarations
        setExtraPhiLam();

        // phi
        phiMin -= phiExtra;
        phiMax += phiExtra;
        // lambda
        lambdaMax += lambdaExtra;
        lambdaMin -= lambdaExtra;

        indexPhi0DEM = (int) (Math.floor((lat0 - phiMax) / latitudeDelta));
        indexPhiNDEM = (int) (Math.ceil((lat0 - phiMin) / latitudeDelta));
        indexLambda0DEM = (int) (Math.floor((lambdaMin - lon0) / longitudeDelta));
        indexLambdaNDEM = (int) (Math.ceil((lambdaMax - lon0) / longitudeDelta));

        //// sanity checks ////
        if (indexPhi0DEM < 0) {
            TopoPhase.logger.warn("indexPhi0DEM: " + indexPhi0DEM);
            indexPhi0DEM = 0;   // reset to default start at first
            TopoPhase.logger.warn("DEM does not cover entire interferogram/tile.");
            TopoPhase.logger.warn("input DEM should be extended to the North.");
        }

        if (indexPhiNDEM > nLatPixels - 1) {
            TopoPhase.logger.warn("indexPhiNDEM: " + indexPhi0DEM);
            indexPhiNDEM = (int) (nLatPixels - 1);
            TopoPhase.logger.warn("DEM does not cover entire interferogram/tile.");
            TopoPhase.logger.warn("input DEM should be extended to the South.");
        }

        if (indexLambda0DEM < 0) {
            TopoPhase.logger.warn("indexLambda0DEM: " + indexLambda0DEM);
            indexLambda0DEM = 0;    // default start at first
            TopoPhase.logger.warn("DEM does not cover entire interferogram/tile.");
            TopoPhase.logger.warn("input DEM should be extended to the West.");
        }

        if (indexLambdaNDEM > nLonPixels - 1) {
            TopoPhase.logger.warn("indexLambdaNDEM: " + indexLambdaNDEM);
            indexLambdaNDEM = (int) (nLonPixels - 1);
            TopoPhase.logger.warn("DEM does not cover entire interferogram/tile.");
            TopoPhase.logger.warn("input DEM should be extended to the East.");
        }
    }

}

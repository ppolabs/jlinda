package org.jdoris.core.geom;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import org.apache.log4j.Logger;
import org.jdoris.core.Orbit;
import org.jdoris.core.Point;
import org.jdoris.core.SLCImage;
import org.jdoris.core.delaunay.FastDelaunayTriangulator;
import org.jdoris.core.delaunay.Triangle;
import org.jdoris.core.delaunay.TriangulationException;

import java.util.ArrayList;
import java.util.List;

/**
 * User: pmar@ppolabs.com
 * Date: 6/14/11
 * Time: 12:16 PM
 */
public class DemFiddler {

    //// logger
    static Logger logger = Logger.getLogger(DemFiddler.class.getName());

    //// fields
    // dem tile parameters!!
    // here it is assumed that these are virgin!
    double l0;
    double lN;
    double p0;
    double pN;

    double lat0;
    double lon0;
    double demDeltaLat;
    double demDeltaLon;
    long nLatPixels;
    long nLonPixels;

    double extraLat;
    double extraLon;

    //// output values
    // [i,j] index within DEM
    int indexPhi0DEM;
    int indexPhiNDEM;
    int indexLambda0DEM;
    int indexLambdaNDEM;

    // [phi,lambda] : of corners in radians
    double phiMin;
    double phiMax;
    double lambdaMin;
    double lambdaMax;
    //    private double[] x_in;
//    private double[] y_in;
//    private double[][] z_in;
    double[][] grd;
//    private double r_az_ratio;
//    private double NODATA = 99999;

    public DemFiddler() {
    }

    public void getDEMCorners(final SLCImage meta,
                              final Orbit orbit) throws Exception {

        double[] phiAndLambda;

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
        // phi
        phiMin -= extraLat;
        phiMax += extraLat;
        // lambda
        lambdaMax += extraLon;
        lambdaMin -= extraLon;

        sanityCheck();

    }

    public void sanityCheck() {

        indexPhi0DEM = (int) (Math.floor((lat0 - phiMax) / demDeltaLat));
        if (indexPhi0DEM < 0) {
            logger.warn("indexPhi0DEM: " + indexPhi0DEM);
            indexPhi0DEM = 0;   // reset to default start at first
            logger.warn("DEM does not cover entire interferogram/tile.");
            logger.warn("input DEM should be extended to the North.");
        }

        indexPhiNDEM = (int) (Math.ceil((lat0 - phiMin) / demDeltaLat));
        if (indexPhiNDEM > nLatPixels - 1) {
            logger.warn("indexPhiNDEM: " + indexPhi0DEM);
            indexPhiNDEM = (int) (nLatPixels - 1);
            logger.warn("DEM does not cover entire interferogram/tile.");
            logger.warn("input DEM should be extended to the South.");
        }

        indexLambda0DEM = (int) (Math.floor((lambdaMin - lon0) / demDeltaLon));
        if (indexLambda0DEM < 0) {
            logger.warn("indexLambda0DEM: " + indexLambda0DEM);
            indexLambda0DEM = 0;    // default start at first
            logger.warn("DEM does not cover entire interferogram/tile.");
            logger.warn("input DEM should be extended to the West.");
        }

        indexLambdaNDEM = (int) (Math.ceil((lambdaMax - lon0) / demDeltaLon));
        if (indexLambdaNDEM > nLonPixels - 1) {
            logger.warn("indexLambdaNDEM: " + indexLambdaNDEM);
            indexLambdaNDEM = (int) (nLonPixels - 1);
            logger.warn("DEM does not cover entire interferogram/tile.");
            logger.warn("input DEM should be extended to the East.");
        }
    }

    public double[][] gridData(final double[][] x_in, final double[][] y_in, final double[][] z_in,
                               double x_min, final double x_max, double y_min, final double y_max,
                               double x_inc, double y_inc, final double r_az_ratio, double offset,
                               final double NODATA) {


        int i, j, k, ij;
        long p;
        long i_min, i_max, j_min, j_max;
        int n, nx, ny;
        int zLoops, zLoop, zBlockSize, zInterpolateBlockSize;
        int indexFirstPoint;

        double[] vx = new double[4];
        double[] vy = new double[4];
        grd = new double[128][512];

        double xkj, xlj;
        double ykj, ylj;
        double zj, zk;
        double zl, zlj;
        double zkj;
        double xp, yp;
        double f; // linear interpolation parameters

        Coordinate[] In;

        // Initialize variables
        final int x_in_dim = x_in.length * x_in[0].length;
        final int z_in_dim = z_in.length * z_in[0].length;
        zBlockSize = x_in_dim; // block size of x and y coordination
        n = zBlockSize;

        // How many groups of z value should be interpolated
        if ((z_in_dim % zBlockSize) != 0) {
            logger.warn("The input of the DEM buffer and z is not the same...");
            return null;
        } else {
            zLoops = z_in.length / x_in.length;
        }

        // containers
        double[] a = new double[zLoops];
        double[] b = new double[zLoops];
        double[] c = new double[zLoops];

        nx = grd.length / zLoops;
        ny = grd[0].length;
//        zInterpolateBlockSize = grd.length * ny / zLoops;

//        In = new Coordinate[n];
//        int tmpCounter = 0;
//        // Copy x,y points to In structure array
//        for (i = 0; i < x_in.length; i++) {
//            for (j = 0; j < x_in[0].length; j++) {
//                In[tmpCounter] = new Coordinate(x_in[i][j], y_in[i][j] * r_az_ratio, z_in[i][j]);
//                tmpCounter++;
//            }
//        }

        // TODO: integrate initialization of In object and GeometryFactory
        // organize input data
        logger.trace("DelaunayTriangulator with " + n + " points");
        long t0 = System.currentTimeMillis();
        List<Geometry> list = new ArrayList<Geometry>();
        GeometryFactory gf = new GeometryFactory();
        for (i = 0; i < x_in.length; i++) {
            for (j = 0; j < x_in[0].length; j++) {
                list.add(gf.createPoint(new Coordinate(x_in[i][j], y_in[i][j] * r_az_ratio, z_in[i][j])));
            }
        }
        long t1 = System.currentTimeMillis();
        logger.info("Input set constructed in " + (0.001 * (t1 - t0)) + " sec");

        // triangulate input data
        long t2 = System.currentTimeMillis();
        FastDelaunayTriangulator FDT = new FastDelaunayTriangulator();
        try {
            FDT.triangulate(list.iterator());
        } catch (TriangulationException te) {
            te.printStackTrace();
        }
        long t3 = System.currentTimeMillis();
        logger.info("Data set triangulated in " + (0.001 * (t3 - t2)) + " sec");

        long t4 = System.currentTimeMillis();
        //// loop over triangles
        for (Triangle triangle : FDT.triangles) {

            // store triangle coordinates in local variables
            vx[0] = vx[3] = triangle.getA().x;
            vy[0] = vy[3] = triangle.getA().y / r_az_ratio;

            vx[1] = triangle.getB().x;
            vy[1] = triangle.getB().y / r_az_ratio;

            vx[2] = triangle.getC().x;
            vy[2] = triangle.getC().y / r_az_ratio;

            // check whether something is no-data
            if (vx[0] == NODATA || vx[1] == NODATA || vx[2] == NODATA)
                continue;
            if (vy[0] == NODATA || vy[1] == NODATA || vy[2] == NODATA)
                continue;

            /* Compute grid indices the current triangle may cover.*/
            xp = Math.min(Math.min(vx[0], vx[1]), vx[2]);
            i_min = coordToIndex(xp, x_min, x_inc, offset);

            xp = Math.max(Math.max(vx[0], vx[1]), vx[2]);
            i_max = coordToIndex(xp, x_min, x_inc, offset);

            yp = Math.min(Math.min(vy[0], vy[1]), vy[2]);
            j_min = coordToIndex(yp, y_min, y_inc, offset);

            yp = Math.max(Math.max(vy[0], vy[1]), vy[2]);
            j_max = coordToIndex(yp, y_min, y_inc, offset);

            /* Adjustments for triangles outside -R region. */
            /* Triangle to the left or right. */
            if ((i_max < 0) || (i_min >= nx))
                continue;
            /* Triangle Above or below */
            if ((j_max < 0) || (j_min >= ny))
                continue;
            /* Triangle covers boundary, left or right. */
            if (i_min < 0)
                i_min = 0;
            if (i_max >= nx)
                i_max = nx - 1;
            /* Triangle covers boundary, top or bottom. */
            if (j_min < 0)
                j_min = 0;
            if (j_max >= ny)
                j_max = ny - 1;

            /* Find equation for the plane as z = ax + by + c */
            xkj = vx[1] - vx[0];
            ykj = vy[1] - vy[0];
            xlj = vx[2] - vx[0];
            ylj = vy[2] - vy[0];

            f = 1.0 / (xkj * ylj - ykj * xlj);

            for (zLoop = 0; zLoop < zLoops; zLoop++) {
                zj = triangle.getA().z;
                zk = triangle.getB().z;
                zl = triangle.getC().z;
                zkj = zk - zj;
                zlj = zl - zj;
                a[zLoop] = -f * (ykj * zlj - zkj * ylj);
                b[zLoop] = -f * (zkj * xlj - xkj * zlj);
                c[zLoop] = -a[zLoop] * vx[1] - b[zLoop] * vy[1] + zk;
            }

            for (i = (int)i_min; i <= i_max; i++) {

                xp = indexToCoord(i, x_min, x_inc, offset);

                for (j = (int)j_min; j <= j_max; j++) {

                    yp = indexToCoord(j, y_min, y_inc, offset);

                    if (!pointInTriangle(vx, vy, xp, yp))
                        continue; /* Outside */

                    for (zLoop = 0; zLoop < zLoops; zLoop++) {
                        grd[i][j] = a[zLoop] * xp + b[zLoop] * yp + c[zLoop];
                    }
                }
            }
        }
        long t5 = System.currentTimeMillis();
        logger.info("Data set interpolated in " + (0.001 * (t5 - t4)) + " sec");

        return grd;

    }

    private static boolean pointInTriangle(double[] xt, double[] yt, double x, double y) {
        int iRet0 = ((xt[2] - xt[0]) * (y - yt[0])) > ((x - xt[0]) * (yt[2] - yt[0])) ? 1 : -1;
        int iRet1 = ((xt[0] - xt[1]) * (y - yt[1])) > ((x - xt[1]) * (yt[0] - yt[1])) ? 1 : -1;
        int iRet2 = ((xt[1] - xt[2]) * (y - yt[2])) > ((x - xt[2]) * (yt[1] - yt[2])) ? 1 : -1;

        return (iRet0 > 0 && iRet1 > 0 && iRet2 > 0) || (iRet0 < 0 && iRet1 < 0 && iRet2 < 0);
    }

    private long coordToIndex(final double coord, final double coord0, final double deltaCoord, final double offset) {
        return irint((((coord - coord0) / (deltaCoord)) - offset));
    }

    private double indexToCoord(final long idx, final double coord0, final double deltaCoord, final double offset) {
        return (coord0 + idx * deltaCoord + offset);
    }


    private long irint(final double coord) {
        return ((long ) rint(coord));
    }

    private double rint(final double coord) {
        return Math.floor(coord + 0.5);
    }


}

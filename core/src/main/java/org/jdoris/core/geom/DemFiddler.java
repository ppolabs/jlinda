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
    private double[] x_in;
    private double[] y_in;
    private double[][] z_in;
    private double[][] grd;
    private double r_az_ratio;
    private double NODATA = 99999;

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

    public void GridData() {

        int i;
        int j;
        int k;
        int ij;
        int p;
        int i_min;
        int i_max;
        int j_min;
        int j_max;
        int n;
        int nx;
        int ny;
        int zLoops;
        int zLoop;
        int zBlockSize;
        int indexFirstPoint;
        int zInterpolateBlockSize;
        double[] vx = new double[4];
        double[] vy = new double[4];
        double xkj;
        double xlj;
        double ykj;
        double ylj;
        double zj;
        double zk;
        double zl;
        double zlj;
        double zkj;
        double xp;
        double yp;
        double f; // linear interpolation parameters

        Coordinate[] In = null;
        Coordinate[] Out = null;

        // Initialize variables
        zBlockSize = x_in.length; // block size of x and y coordination
        n = zBlockSize;


        // How many groups of z value should be interpolated
        if ((z_in.length % zBlockSize) != 0) {
            System.out.println("The input of the DEM buffer and z is not the same...");
            return;
        } else {
            zLoops = z_in.length / x_in.length;
        }

        double[] a = new double[zLoops];
        double[] b = new double[zLoops];
        double[] c = new double[zLoops];

//        DoubleMatrix a = new DoubleMatrix(zLoops);
//        DoubleMatrix b = new DoubleMatrix(zLoops);
//        DoubleMatrix c = new DoubleMatrix(zLoops);

/*
        if (a == null || b == null || c == null)
		{
		  ERROR << "Memory ERROR in source file: " << __FILE__ << " at line: " << __LINE__;
		  PRINT_ERROR(ERROR.get_str());
		  throw(memory_error);
		}
*/
        // TODO: check the orientation
        // if JBLASS used
//        nx = grd.rows / zLoops;
//        ny = grd.columns;
//        zInterpolateBlockSize = grd.length / zLoops;

        // if ARRAYS used
        nx = grd.length / zLoops;
        ny = grd[0].length;
//        zInterpolateBlockSize = grd.length * ny / zLoops;

//        In = new Coordinate[2 * n];
        In = new Coordinate[n];

        // Copy x,y points to In structure array
        for (i = 0; i < n; i++) {
            In[i] = new Coordinate(x_in[i], y_in[i] * r_az_ratio);
        }

        logger.trace("testFastDelaunayTriangulator with " + In.length + " points");
        long t0 = System.currentTimeMillis();
        List<Geometry> list = new ArrayList<Geometry>();
        GeometryFactory gf = new GeometryFactory();
        for (Coordinate coord : In) {
            list.add(gf.createPoint(coord));
        }
        long t1 = System.currentTimeMillis();
        logger.info("Input set constructed in %10.3f sec\n" + (0.001 * (t1 - t0)));

        long t2 = System.currentTimeMillis();
        FastDelaunayTriangulator FDT = new FastDelaunayTriangulator();
        try {
            FDT.triangulate(list.iterator());
        } catch (TriangulationException te) {
            te.printStackTrace();
        }
        long t3 = System.currentTimeMillis();

        logger.info("   triangulated in %10.3f sec\n" + (0.001 * (t3 - t2)));

        // here it loops through triangles!
        for (Triangle triangle : FDT.triangles) {

            // store the index of the first Point of this triangle
            indexFirstPoint = triangle.getIndex(triangle.getA());

            // get coordinates from trianglelist
            vx[0] = vx[3] = triangle.getA().x;
            vy[0] = vy[3] = triangle.getA().y;

            vx[1] = triangle.getB().x;
            vy[1] = triangle.getB().y;

            vx[2] = triangle.getB().x;
            vy[2] = triangle.getB().y;

            // check whether something is no-data
            if (vx[0] == NODATA || vx[1] == NODATA || vx[2] == NODATA)
                continue;
            if (vy[0] == NODATA || vy[1] == NODATA || vy[2] == NODATA)
                continue;

            /* Compute grid indices the current triangle may cover.*/
            xp = Math.min(Math.min(vx[0], vx[1]), vx[2]);
            double x_min = 0;
            long x_inc = 0;
            double offset = 0;
            i_min = (int) coordToIndex(xp, x_min, x_inc, offset);
            //INFO << "xp: " << xp;
            //INFO.print();

            xp = Math.max(Math.max(vx[0], vx[1]), vx[2]);
            i_max = (int) coordToIndex(xp, x_min, x_inc, offset);
            //INFO << "xp: " << xp;
            //INFO.print();

            yp = Math.min(Math.min(vy[0], vy[1]), vy[2]);
            double y_min = 0;
            long y_inc = 0;
            j_min = (int) coordToIndex(yp, y_min, y_inc, offset);

            yp = Math.max(Math.max(vy[0], vy[1]), vy[2]);
            j_max = (int) coordToIndex(yp, y_min, y_inc, offset);

            /* Adjustments for triangles outside -R region. */
            /* Triangle to the left or right. */
            if ((i_max < 0) || (i_min >= nx)) continue;
            /* Triangle Above or below */
            if ((j_max < 0) || (j_min >= ny)) continue;
            /* Triangle covers boundary, left or right. */
            if (i_min < 0) i_min = 0;
            if (i_max >= nx) i_max = nx - 1;
            /* Triangle covers boundary, top or bottom. */
            if (j_min < 0) j_min = 0;
            if (j_max >= ny) j_max = ny - 1;

            /* Find equation for the plane as z = ax + by + c */
            xkj = vx[1] - vx[0];
            ykj = vy[1] - vy[0];
            xlj = vx[2] - vx[0];
            ylj = vy[2] - vy[0];

            f = 1.0 / (xkj * ylj - ykj * xlj);

            for (zLoop = 0; zLoop < zLoops; zLoop++) {
                zj = z_in[zLoop][triangle.getIndex(triangle.getA())];
                zk = z_in[zLoop][triangle.getIndex(triangle.getB())];
                zl = z_in[zLoop][zBlockSize + triangle.getIndex(triangle.getC())];
                zkj = zk - zj;
                zlj = zl - zj;
                a[zLoop] = -f * (ykj * zlj - zkj * ylj);
                b[zLoop] = -f * (zkj * xlj - xkj * zlj);
                c[zLoop] = -a[zLoop] * vx[1] - b[zLoop] * vy[1] + zk;
            }

            for (i = i_min; i <= i_max; i++) {

                xp = indexToCoord(i, x_min, x_inc, offset);

                for (j = j_min; j <= j_max; j++) {

                    yp = indexToCoord(j, y_min, y_inc, offset);

                    if (!pointInTriangle(vx, vy, xp, yp))
                        continue; /* Outside */

                    for (zLoop = 0; zLoop < zLoops; zLoop++) {
                        grd[i][j] = a[zLoop] + xp + b[zLoop] * yp + c[zLoop];
                    }
                }
            }
        }


    }

    private static boolean pointInTriangle(double[] xt, double[] yt, double x, double y) {
        int iRet0 = ((xt[2] - xt[0]) * (y - yt[0])) > ((x - xt[0]) * (yt[2] - yt[0])) ? 1 : -1;
        int iRet1 = ((xt[0] - xt[1]) * (y - yt[1])) > ((x - xt[1]) * (yt[0] - yt[1])) ? 1 : -1;
        int iRet2 = ((xt[1] - xt[2]) * (y - yt[2])) > ((x - xt[2]) * (yt[1] - yt[2])) ? 1 : -1;

//        if ((iRet0 > 0 && iRet1 > 0 && iRet2 > 0) || (iRet0 < 0 && iRet1 < 0 && iRet2 < 0))
//            return true;
//        else
//            return false;

        return (iRet0 > 0 && iRet1 > 0 && iRet2 > 0) || (iRet0 < 0 && iRet1 < 0 && iRet2 < 0);

    }

    private long coordToIndex(final double coord, final double coord0, final long deltaCoord, final double offset) {
        return (irint(((((coord) - (coord0)) / (deltaCoord)) - (offset))));
    }

    private double indexToCoord(final long idx, final double coord0, final long deltaCoord, final double offset) {
        return (coord0 + idx * deltaCoord + offset);
    }


    private long irint(final double coord) {
        return ((long) rint(coord));
    }

    private double rint(final double coord) {
        return Math.floor(coord + 0.5);
    }


}

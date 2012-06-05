package org.jdoris.core.geocode;

import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;
import org.jdoris.core.Baseline;
import org.jdoris.core.Orbit;
import org.jdoris.core.SLCImage;
import org.jdoris.core.Window;
import org.jdoris.core.utils.MathUtils;
import org.jdoris.core.utils.PolyUtils;

/**
 * DInSAR prototype class
 * <p/>
 * Differential insar with an unwrapped topo interferogram (hgt or real4 format) and a wrapped(!) defo interf.
 * <p/>
 * The topography is removed from the deformation interferogram by the formula (prime ' denotes defo pair):
 * <p/>
 * dr = lambda\4pi * [phi' - phi(Bperp'/Bperp)]
 * phi_diff = phi(Bperp'/Bperp) - phi'
 * <p/>
 * where Bperp is the perpendicular baseline for points on the ellipsoid (and not the true one)!
 * <p/>
 * Implementation details: First the baseline for a number of points for topo and defo interferograms computed,
 * and then the ratio between baselines is modeled by a 2D polynomial of degree 1. Then this polynomial is
 * evaluated to compute the new (defo only) phase according to the equation above.
 * It is assumed that the interferogram's are coregistered on each other and have the same dimensions.
 * <p/>
 * Input:
 * -input parameters
 * -orbits
 * -info on input files
 * <p/>
 * Output:
 * -complex float file with differential phase.
 * (set to (0,0) for not ok unwrapped parts)
 */

public class DInSAR {

    final private SLCImage masterMeta;
    final private SLCImage slaveDefoMeta;
    final private SLCImage topoSlaveMeta;
    final private Orbit masterOrbit;
    final private Orbit slaveDefoOrbit;
    final private Orbit slaveTopoOrbit;

    private volatile Window dataWindow;
    private volatile Window tileWindow;

    private volatile DoubleMatrix topoData;
    private volatile ComplexDoubleMatrix defoData;

    public DInSAR(SLCImage masterMeta, SLCImage slaveDefoMeta, SLCImage topoSlaveMeta, Orbit masterOrbit, Orbit slaveDefoOrbit, Orbit slaveTopoOrbit) {
        this.masterMeta = masterMeta;
        this.slaveDefoMeta = slaveDefoMeta;
        this.topoSlaveMeta = topoSlaveMeta;
        this.masterOrbit = masterOrbit;
        this.slaveDefoOrbit = slaveDefoOrbit;
        this.slaveTopoOrbit = slaveTopoOrbit;

        dataWindow = masterMeta.getCurrentWindow();
    }

    public void setDataWindow(Window dataWindow) {
        this.dataWindow = dataWindow;
    }

    public void setTileWindow(Window tileWindow) {
        this.tileWindow = tileWindow;
    }

    public void setTopoData(DoubleMatrix topoData) {
        this.topoData = topoData;
    }

    public void setDefoData(ComplexDoubleMatrix defoData) {
        this.defoData = defoData;
    }

    public ComplexDoubleMatrix getDefoData() {
        return defoData;
    }

    public void dinsar() throws Exception {

        // Normalization factors for polynomial
        double minL = dataWindow.linelo;
        double maxL = dataWindow.linehi;
        double minP = dataWindow.pixlo;
        double maxP = dataWindow.pixhi;

        // Model perpendicular baseline for master and defo
        // .....compute B on grid every 500 lines, 100 pixels
        // .....in window for topoData/defoData
        final int nPoints = 200;
        final int[][] positionArray = MathUtils.distributePoints(nPoints, masterMeta.getCurrentWindow());
        DoubleMatrix position = new DoubleMatrix(nPoints, 2);
        for (int i = 0; i < nPoints; i++) {
            position.put(i, 0, positionArray[i][0]);
            position.put(i, 1, positionArray[i][1]);
        }

        // model baselines
        Baseline topoBaseline = new Baseline();
        topoBaseline.model(masterMeta, topoSlaveMeta, masterOrbit, slaveTopoOrbit);

        Baseline defoBaseline = new Baseline();
        defoBaseline.model(masterMeta, slaveDefoMeta, masterOrbit, slaveDefoOrbit);

        double lastline = -1.0;

        DoubleMatrix bPerpTopo = new DoubleMatrix(nPoints);
        DoubleMatrix bPerpDefo = new DoubleMatrix(nPoints);

        for (int i = 0; i < nPoints; ++i) {
            bPerpTopo.put(i, topoBaseline.getBperp(position.get(i, 0), position.get(i, 1), 0));
            bPerpDefo.put(i, defoBaseline.getBperp(position.get(i, 0), position.get(i, 1), 0));
        }

        // Now model ratio bPerpDefo/bPerpTopo as linear ______
        //   ...r(l,p) = a00 + a10*l + a01*p ______
        //   ...give stats on max. error ______

        DoubleMatrix baselineRatio = bPerpDefo.div(bPerpTopo);

        DoubleMatrix rhs = new DoubleMatrix(PolyUtils.polyFit2D(normalize(position.getColumn(0), minL, maxL),
                normalize(position.getColumn(1), minP, maxP), baselineRatio, 1));

//        /** per Tile: scale TOPO and subtract from DEFO : FOR TILE!!! */
//        int numlines = (int) (tileWindow.lines() / masterMeta.getMlAz());
//        int numpixels = (int) (tileWindow.pixels() / masterMeta.getMlRg());
//        float firstline = (float) (tileWindow.linelo + (masterMeta.getMlAz() - 1.) / 2.);
//        float firstpixel = (float) (tileWindow.pixlo + (masterMeta.getMlAz() - 1.) / 2.);

        DoubleMatrix azimuthAxisNormalize = DoubleMatrix.linspace((int) tileWindow.linelo, (int) tileWindow.linehi, defoData.rows);
        normalize_inplace(azimuthAxisNormalize, minL, maxL);

        DoubleMatrix rangeAxisNormalize = DoubleMatrix.linspace((int) tileWindow.pixlo, (int) tileWindow.pixhi, defoData.columns);
        normalize_inplace(rangeAxisNormalize, minP, maxP);

        DoubleMatrix ratio = PolyUtils.polyval(azimuthAxisNormalize, rangeAxisNormalize, rhs, PolyUtils.degreeFromCoefficients(rhs.length));

        DoubleMatrix scaledTopo = topoData.mul(ratio);
        ComplexDoubleMatrix ratioBaselinesCplx = new ComplexDoubleMatrix(MatrixFunctions.cos(scaledTopo), MatrixFunctions.sin(scaledTopo).neg());

        for (int i = 0; i < defoData.length; i++) {
            if (defoData.data[i] == Double.NaN) {
                defoData.data[i] = 0.0d;
            }
        }
        defoData.muli(ratioBaselinesCplx);

////        DoubleMatrix ratioline = new DoubleMatrix(linspace((int) firstpixel, (int) (firstpixel + numpixels), (int) numpixels));
//        DoubleMatrix ratioline = DoubleMatrix.linspace((int) firstpixel, (int) (firstpixel + numpixels - 1), numpixels);
//        normalize_inplace(ratioline, minP, maxP);
//        ratioline.muli(rhs[2]);  //     a01*p
//        ratioline.addi(rhs[0]);  // a00+a01*p
//
//        // read in matrices line by line, correct phase ______
////        ComplexDoubleMatrix DEFO = new ComplexDoubleMatrix(1, numpixels);    // buffer
//
//        for (int i = 0; i < numlines - 1; ++i) {
//
//            // ratio = a00 + a10*LINE + a01*PIXEL
//            double line = firstline + i * masterMeta.getMlAz();
//            DoubleMatrix ratio = ratioline.add(rhs[1] * normalize2(line, minL, maxL));
//
//            // read from file, correct, write to file ______
////            const window filewin(i + 1, i + 1, 1, numpixels);
//            DoubleMatrix topoRow = topoData.getRow(i); //readphase(filewin);
//            ComplexDoubleMatrix DEFO = defoData.getRow(i);
//
//            // seems faster, but how to check for unwrapping?
//            //TOPO *= ratio;    // scaled topoData to defoData baseline
//            //DEFO *= complr4(cos(TOPO),-sin(TOPO));
//            // better matrix <int32> index = TOPO.find(NaN); later reset??
//            // and topoData=topoData*ratio; and defoData(index)=(0,0);
//            // but how to implement this best in matrixclass?
//            // BK 24-Oct-2000
//            for (int j = 0; j < numpixels - 1; ++j) {
//                if (topoRow.get(0, j) == Double.NaN) {
//                    DEFO.put(0, j, new ComplexDouble(0, 0));
//                } else {
//                    double realPart = ratio.get(j) * topoRow.get(j);
//                    double imagPart = ratio.get(j) * topoRow.get(j);
//                    DEFO.put(j, DEFO.get(j).mul(new ComplexDouble(Math.cos(realPart), Math.sin(imagPart))));
////                    DEFO.put(0, j, DEFO.get(0, j).mul(new ComplexDouble(Math.cos(realPart), Math.sin(imagPart))));
//                }
////                (TOPO(0, j) == NaN) ?                // if unwrapping ok, then subtract scaled phase
////                        DEFO(0, j) = complr4(0.0, 0.0) :
////                        DEFO(0, j) *= complr4(fast_cos(ratio(0, j) * TOPO(0, j)), fast_min_sin(ratio(0, j) * TOPO(0, j)));
//            }
//
//            tempData.putRow(i, DEFO);
//
//        }

    }

    public static double[] linspace(final int lower, final int upper, final int size) {
        double[] result = new double[size];
        for (int i = 0; i < size; i++) {
            double t = (double) i / (size - 1);
            result[i] = lower * (1 - t) + t * upper;
        }
        return result;
    }

    public static double[] linspace(final int lower, final int upper) {
        return linspace(lower, upper, 100);
    }

    private void normalize_inplace(DoubleMatrix data, double min, double max) {
        data.subi(.5 * (min + max));
        data.divi(.25 * (max - min));
    }

    private DoubleMatrix normalize(DoubleMatrix in, double min, double max) {
        DoubleMatrix out = in.dup();
        out.subi(.5 * (min + max));
        out.divi(.25 * (max - min));
        return out;
    }

}

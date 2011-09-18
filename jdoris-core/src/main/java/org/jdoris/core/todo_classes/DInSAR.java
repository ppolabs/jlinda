package org.jdoris.core.todo_classes;


// integration class for DInSAR -- to be removed after initialization -- not to be commited!

import org.jblas.DoubleMatrix;
import org.jblas.Solve;
import org.jdoris.core.Baseline;
import org.jdoris.core.Ellipsoid;
import org.jdoris.core.Orbit;
import org.jdoris.core.SLCImage;
import org.jdoris.core.utils.LinearAlgebraUtils;

import static org.jdoris.core.utils.PolyUtils.normalize2;

/**
 * *************************************************************
 * dinsar                                                    *
 * Differential insar with an unwrapped topo interferogram      *
 * (hgt or real4 format) and a wrapped(!) defo interf.          *
 * if r4 then NaN==-999 is problem with unwrapping, else hgt    *
 * if ampl. =0 then problem flagged with unwrapping.            *
 * The topography is removed from the deformation interferogram *
 * by the formula (prime ' denotes defo pair):                  *
 * dr       = lambda\4pi * [phi' - phi(Bperp'/Bperp)]           *
 * phi_diff = phi(Bperp'/Bperp) - phi'                          *
 * where Bperp is the perpendicular baseline for points on the  *
 * ellipsoid (and not the true one)!                            *
 * I implemented this by computing the baseline for a number    *
 * of points for topo and defo, and then modeling the ratio     *
 * as a 2D polynomial of degree 1 for the image.                *
 * Then evaluating this to compute the new phase (defo only).   *
 * I assume the interferogram files are coregistered on each    *
 * other and have the same dimensions.                          *
 * *
 * If TOPOMASTER file is empty (" "), then use current master   *
 * res file for master orbit (== 3pass), else get orbit         *
 * (==4pass).                                                   *
 * *
 * Input:                                                       *
 * -input parameters                                           *
 * -orbits                                                     *
 * -info on input files                                        *
 * -                                                           *
 * Output:                                                      *
 * -complex float file with differential phase.                *
 * (set to (0,0) for not ok unwrapped parts)                  *
 * *
 * See also Zebker, 1994.                                       *
 * #%// BK 22-Sep-2000                                            *
 * **************************************************************
 */

public class DInSAR {

    public static void dinsar(
            todo_classes.inputgeneral input_general,
            todo_classes.input_dinsar dinsarinput,
            Ellipsoid ellips,
            SLCImage master,
            Orbit masterorbit,
            SLCImage defoslave,
            Orbit defoorbit,
            todo_classes.productinfo defointerferogram
    ) throws Exception {

        SLCImage toposlave = null;                                 // info on slave image
        todo_classes.productinfo topounwrappedinterf = null;       // interferogram
        Orbit toposlaveorbit = null;                               // always fill

        // ______ Check 4 pass if topomaster is specified (diff than m_res) ______
        boolean FOURPASS = false;                 // assume 4pass

        // ______ Normalization factors for polynomial ______
        double minL = master.getCurrentWindow().linelo;
        double maxL = master.getCurrentWindow().linehi;
        double minP = master.getCurrentWindow().pixlo;
        double maxP = master.getCurrentWindow().pixhi;

        // ====== Model perpendicular baseline for master and slave ======
        // ______ compute B on grid every 500 lines, 100 pixels
        // ______ in window for topo/defo ______
        int numpointsL = 20;                                  // grid for modelling
        int numpointsP = 10;                                  // grid for modelling
        double dlines = (toposlave.getCurrentWindow().linehi - toposlave.getCurrentWindow().linelo) /
                (numpointsL - 1);
        double dpixels = (toposlave.getCurrentWindow().pixhi - toposlave.getCurrentWindow().pixlo) /
                (numpointsP - 1);

        double[] LINENUMBER = new double[numpointsL * numpointsP];
        double[] PIXELNUMBER = new double[numpointsL * numpointsP];

        int k = 0;

        // define the grid
        for (int i = 0; i < numpointsL; ++i) {
            for (int j = 0; j < numpointsP; ++j) {
                LINENUMBER[k] = topounwrappedinterf.win.linelo + i * dlines;     // line coordinate
                PIXELNUMBER[k] = topounwrappedinterf.win.pixlo + j * dpixels;     // pixel coordinate
                ++k;
            }
        }

        // model baselines
        Baseline topo_baseline = new Baseline();
        topo_baseline.model(master, defoslave, masterorbit, defoorbit);

        Baseline defo_baseline = new Baseline();
        defo_baseline.model(master, toposlave, masterorbit, toposlaveorbit);


        double lastline = -1.0;

        double[] Bperptopo = new double[LINENUMBER.length];
        double[] Bperpdefo = new double[LINENUMBER.length];

        for (int i = 0; i < LINENUMBER.length; ++i) {
            Bperptopo[i] = topo_baseline.getBperp(LINENUMBER[i], PIXELNUMBER[i], 0);
            Bperpdefo[i] = defo_baseline.getBperp(LINENUMBER[i], PIXELNUMBER[i],0);
        }

        // ______ Now model ratio Bperpdefo/Bperptopo as linear ______
        // ______ r(l,p) = a00 + a10*l + a01*p ______
        // ______ give stats on max. error ______

        double[] Ratio = new double[Bperpdefo.length];

        for (int l = 0; l < Bperpdefo.length; l++) {
            Ratio[l] = Bperpdefo[l]/Bperptopo[l];
        }

        // ______ Set designmatrix, compute normalmatrix, righthandside ______
        DoubleMatrix A = new DoubleMatrix(Ratio.length, 3);
        for (int i = 0; i < A.rows; ++i) {
            A.put(i, 0, 1.0);
            A.put(i, 2, normalize2(PIXELNUMBER[i], minP, maxP));
            A.put(i, 1, normalize2(LINENUMBER[i], minL, maxL));
        }

        DoubleMatrix N = LinearAlgebraUtils.matTxmat(A, A);
        DoubleMatrix rhs = LinearAlgebraUtils.matTxmat(A, new DoubleMatrix(Ratio));

        rhs = Solve.solve(N, rhs);


//        // ______ Test inverse (thus stability cholesky) ______
//        for (uint i = 0; i < Qx_hat.lines(); i++)
//            for (uint j = 0; j < i; j++)
//                Qx_hat(j, i) = Qx_hat(i, j);// repiar only stored Lower tri
//        const real8 maxdev = max(abs(N * Qx_hat - eye(real8(Qx_hat.lines()))));
//        INFO << "dinsar: max(abs(N*inv(N)-I)) = " << maxdev;
//        INFO.print();
//        if (maxdev > 0.01) {
//            ERROR << ". Too large, normalization factors <-> crop?";
//            PRINT_ERROR(ERROR.get_str())
//            throw (some_error);
//        } else if (maxdev > 0.001) {
//            WARNING.print("Deviation quite large.  careful!");
//        }
//
//
//        // ______ Some other stuff for logfile ______
//        //  matrix<real8> Qy_hat        = A * (matxmatT(Qx_hat,A));
//        matrix<real8> y_hat = A * rhs;
//        matrix<real8> e_hat = Ratio - y_hat;
//
//        uint pos, dummy;
//        const real8 maxerrorratio = max(abs(e_hat), pos, dummy);
//        const real8 maxrelerror = 100.0 * maxerrorratio / Ratio(pos, 0);
//        INFO << "maximum error for l,p: " << LINENUMBER(pos, 0) << ","
//                << PIXELNUMBER(pos, 0) << "; Ratio=" << Ratio(pos, 0)
//                << " estimate=" << y_hat(pos, 0) << "; rel. err=" << maxrelerror << "%. ";
//        INFO.print();
//        if (maxrelerror < 5.0) {
//            INFO.print("max err OK");
//        } else {
//            WARNING.print("max err quite large");
//            WARNING.print("Error in deformation vector larger than 5% due to mismodeling baseline!");
//        }


        // ====== Per Tile
        boolean writescaledtopo = false;

        // scale TOPO and subtract from DEFO

//        int numlines = (int) (defointerferogram.win.lines() / defointerferogram.multilookL);
//        int numpixels = (int) (defointerferogram.win.pixels() / defointerferogram.multilookP);
//        float firstline = (float) (defointerferogram.win.linelo + (defointerferogram.multilookL - 1.) / 2.);
//        float firstpixel = (float) (defointerferogram.win.pixlo + (defointerferogram.multilookP - 1.) / 2.);
//
//        DoubleMatrix ratioline = new DoubleMatrix(1, numpixels);
//
//        for (int i = 0; i < numpixels; ++i) {
//            ratioline.put(0, i, firstpixel + i * defointerferogram.multilookP);
//        }
//
//        normalize(ratioline, real4(minP), real4(maxP));
//        ratioline *= real4(rhs(2, 0));         //     a01*p
//        ratioline += real4(rhs(0, 0));         // a00+a01*p
//
//        // ______ read in matrices line by line, correct phase ______
//        ComplexDoubleMatrix DEFO = new ComplexDoubleMatrix(1, numpixels);    // buffer
//
//        // test if reading phase was ok...
//        //cerr << "test: writing hgt phase and cint.\n";
//        //ofstream oftest1("TOPO.raw", ios::out | ios::binary | ios::trunc);
//        //ofstream oftest2("CINT.raw", ios::out | ios::binary | ios::trunc);
//
//        for (int i = 0; i < numlines; ++i) {
//
//            // ______ ratio=a00+a10*l+a01*p ______
//            double line = firstline + i * defointerferogram.multilookL;
//            DoubleMatrix ratio = ratioline + real4(rhs(1, 0)) * normalize(line, real4(minL), real4(maxL));
//
//            // ______ read from file, correct, write to file ______
//            ifdefocint >> DEFO;                         // read next full line
//            const window filewin(i + 1, i + 1, 1, numpixels);
//            matrix<real4> TOPO = topounwrappedinterf.readphase(filewin);
//
//            // seems faster, but how to check for unwrapping?
//            //TOPO *= ratio;    // scaled topo to defo baseline
//            //DEFO *= complr4(cos(TOPO),-sin(TOPO));
//            // better matrix <int32> index = TOPO.find(NaN); later reset??
//            // and topo=topo*ratio; and defo(index)=(0,0);
//            // but how to implement this best in matrixclass?
//            // BK 24-Oct-2000
//            for (int32 j = 0; j < numpixels; ++j) {
//                (TOPO(0, j) == NaN) ?                // if unwrapping ok, then subtract scaled phase
//                        DEFO(0, j) = complr4(0.0, 0.0) :
//                        DEFO(0, j) *= complr4(fast_cos(ratio(0, j) * TOPO(0, j)), fast_min_sin(ratio(0, j) * TOPO(0, j)));
//            }
//            ofcint << DEFO;             // now topo-corrected phase
//
//        }

    }
}

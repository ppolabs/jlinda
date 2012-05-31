package org.jdoris.core.geocode;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.Logger;
import org.jblas.DoubleMatrix;
import org.jblas.Solve;
import org.jdoris.core.*;
import org.jdoris.core.todo_classes.todo_classes;
import org.jdoris.core.utils.LinearAlgebraUtils;
import org.jdoris.core.utils.MathUtils;
import org.jdoris.core.utils.PolyUtils;
import org.slf4j.LoggerFactory;

import static org.jdoris.core.Constants.PI;
import static org.jdoris.core.Constants.SOL;
import static org.jdoris.core.utils.LinearAlgebraUtils.matTxmat;
import static org.jdoris.core.utils.PolyUtils.normalize2;

public class Slant2Height {

    private static final Logger logger = (Logger) LoggerFactory.getLogger(Slant2Height.class);

    /**
     * slant2h-eight (schwabisch)
     * <p/>
     * compute height in radar coded system (master):
     * <p/>
     * <p/>
     * 1.  compute reference phase for h=0,2000,4000 in Npoints
     * <p/>
     * 2.  solve system: h(phi) = a_0 + a_1*phi + a_2*phi*phi
     * (2nd degree 1D polynomial) for all Npoints
     * <p/>
     * 3.  compute a_i (l,p) = DEGREE2D 2D polynomial
     * <p/>
     * 4.0 set offset to the one of first pixel , add this to all
     * this step is skipped, phase is w.r.t. h=0, ref. is subtracted
     * <p/>
     * 4.1 evaluate polynomial of 3. for all points (l,p) of
     * (multilooked) unwrapped interferogram
     * <p/>
     * Note: solution to system for betas seems not be very stable!??
     */
    public void schwabisch(todo_classes.inputgeneral generalinput,
                           todo_classes.input_slant2h slant2hinput,
                           Ellipsoid ellips,
                           SLCImage master, SLCImage slave,
                           todo_classes.productinfo unwrappedinterf,
                           Orbit masterorbit, Orbit slaveorbit) throws Exception {

        logger.setLevel(Level.DEBUG);
        logger.trace("slant2h Schwabisch (PM 01-Apr-2011)");

        final int MAXITER = 10;
        final double CRITERPOS = 1e-6;
        final double CRITERTIM = 1e-10;
        final double m_minpi4cdivlambda = (-4. * PI * SOL) / master.getRadarWavelength();
        final double s_minpi4cdivlambda = (-4. * PI * SOL) / slave.getRadarWavelength();

        final int Npoints = slant2hinput.Npoints; // where ref.phase is evaluated
        final int DEGREE1D = slant2hinput.degree1d; // only possible now.
        final int DEGREE2D = slant2hinput.degree2d;
        final int Nheights = slant2hinput.Nheights;

        final int MAXHEIGHT = 5000; // max hei for ref.phase
        final int TEN = 10; // used in pointer
        if (DEGREE1D - 1 > TEN) {
            logger.error("panic, programmers problem: increase TEN.");
            throw new IllegalArgumentException();
        }

        final int HEIGHTSTEP = MAXHEIGHT / (Nheights - 1); // heights to eval ref.phase

        // Matrices for storing phase for all ref. ellipsoids
        //  PHASE(i,0)  phase for height 0
        //  PHASE(i,1)  phase for height Heigthsep * 1
        //  PHASE(i,Nh) phase for height 4000
        DoubleMatrix PHASE = new DoubleMatrix(Npoints, Nheights); // pseudo-observation

        // Distribute points in original master system (not multilooked)
        // (i,0): line, (i,1): pixel, (i,2) flagfromdisk (not used here)
        int[][] positionArray = MathUtils.distributePoints(Npoints, unwrappedinterf.win);

        DoubleMatrix Position = new DoubleMatrix(Npoints, 2);
        for (int i = 0; i < Npoints; i++) {
            Position.put(i, 1, positionArray[i][0]);
            Position.put(i, 2, positionArray[i][1]);
        }

        /** -- STEP 1 : compute reference phase in N points for nHeights ---------------*/
        // =============================================================================
        // Compute reference phase in N points for height (numheight)
        logger.debug("S2H: schwabisch: STEP1: compute reference phase for Nheights.");
        for (int numheight = 0; numheight < Nheights; numheight++) {
            int HEIGHT = numheight * HEIGHTSTEP;

            // Compute delta r for all points
            for (int i = 0; i < Npoints; i++) {

                int line = (int) Position.get(i, 0);
                int pixel = (int) Position.get(i, 1);
                double m_trange = master.pix2tr(pixel);

                // Compute xyz of point P on ELLIPS for this line,pixel
                Point xyzMaster = masterorbit.lph2xyz(line, pixel, HEIGHT, master);

                // Compute xyz of slave satelite in orbit_slave from P
                Point timeSlave = slaveorbit.xyz2t(xyzMaster, slave);

                PHASE.put(i, numheight, m_minpi4cdivlambda * m_trange - s_minpi4cdivlambda * timeSlave.x);
            }
        }

        //  Subtract ref. phase at h=0 for all point
        //  this is the same as adding reference phase for all in uint
        for (int i = 0; i < Npoints; ++i) {
            double offset = PHASE.get(i, 0);
            PHASE.put(i, 0, 0.d);
            for (int j = 1; j < Nheights; ++j) {
                PHASE.put(i, j, offset);
            }
        }

        /** -- STEP 2 : compute alpha coefficients of polynomials for these points ---*/

        logger.debug("S2H: schwabisch: STEP2: estimate coefficients 1d polynomial.");

        DoubleMatrix DESIGN = new DoubleMatrix(Nheights, DEGREE1D + 1); // design matrix
        DoubleMatrix ALPHAS = new DoubleMatrix(Npoints, DEGREE1D + 1); // pseudo-observation
        DoubleMatrix HEI = new DoubleMatrix(Nheights, 1);
        for (int i = 0; i < Nheights; i++)
            HEI.put(i, 0, i * HEIGHTSTEP); // 0, .., 5000

        // ______ normalize data to [0,1] ______
        double minphi = PHASE.min();
        double maxphi = PHASE.max();
        normalize(PHASE, minphi, maxphi); // regrid data

        for (int i = 0; i < Npoints; i++) // solve system for all points
        {
            // Set up design matrix
            for (int j = 0; j < Nheights; j++) {
                for (int k = 0; k <= DEGREE1D; k++) {
                    DESIGN.put(j, k, Math.pow(PHASE.get(i, j), k)); // PHASE is normalized
                }
            }

            // Solve by cholesky (even the exactly determined case)
            DoubleMatrix N = LinearAlgebraUtils.matTxmat(DESIGN, DESIGN);
            DoubleMatrix rhs = LinearAlgebraUtils.matTxmat(DESIGN, HEI);
            DoubleMatrix Qx_hat = N;

            rhs = Solve.solvePositive(Qx_hat, rhs);
            Qx_hat = Solve.solvePositive(Qx_hat, DoubleMatrix.eye(Qx_hat.getRows()));

            double maxdev = (N.mmul(Qx_hat).sub(DoubleMatrix.eye(Qx_hat.getRows()))).normmax();
            logger.info("s2h schwaebisch: max(abs(N*inv(N)-I)) = {}", maxdev);
            if (maxdev > .01)
                logger.warn("wrong solution for 1d polynomial? (decrease d1d or nhei)");

            // Scale back unknowns: alpha_i <= alpha_i * (scale)^i
            // Store solution in ALPHAS
            for (int alfa = 0; alfa <= DEGREE1D; alfa++) {
                ALPHAS.put(i, alfa, rhs.get(alfa, 0));
            }
        } // loop over all points

        /** -- STEP 3 : Compute alpha_i coefficients of polynomials as function of (l,p) --*/

        logger.debug("S2H: schwabisch: STEP3: estimate coefficients for 2d polynomial.");
        // Compute alpha_i coefficients of polynomials as function of (l,p)
        // ... alpha_i = sum(k,l) beta_kl l^k p^l;
        // ... Solve simultaneous for all betas
        // ... this does not seem to be possibly with my routine, so do per alfa_i
        final int Nunk = PolyUtils.numberOfCoefficients(DEGREE2D); // Number of unknowns

        // ______ Check redundancy is done before? ______
        if (Npoints < Nunk) {
            logger.error("slant2hschwabisch: N_observations<N_unknowns (increase S2H_NPOINTS or decrease S2H_DEGREE2D.");
            throw new IllegalArgumentException();
        }

        DoubleMatrix A = new DoubleMatrix(Npoints, Nunk); // designmatrix

        // Set up system of equations
        // .... Order unknowns: B00 B10 B01 B20 B11 B02 B30 B21 B12 B03 for degree=3
        double minL = Position.getColumn(1).min();
        double maxL = Position.getColumn(1).max();
        double minP = Position.getColumn(2).min();
        double maxP = Position.getColumn(2).max();
        for (int i = 0; i < Npoints; i++) {
            // ______ normalize coordinates ______
            double posL = normalize2(Position.get(i, 0), minL, maxL);
            double posP = normalize2(Position.get(i, 1), minP, maxP);

            int index = 0;
            for (int j = 0; j <= DEGREE2D; j++) {
                for (int k = 0; k <= j; k++) {
                    A.put(i, index, Math.pow(posL, j - k) * Math.pow(posP, k));
                    index++;
                }
            }
        }

        // Solve 2d polynomial system for alfas at these points
        DoubleMatrix N = matTxmat(A, A);
        DoubleMatrix rhs = matTxmat(A, ALPHAS);
        DoubleMatrix Qx_hat = N;
        // choles(Qx_hat); // Cholesky factorisation normalmatrix
        // solvechol(Qx_hat,rhs);        // Estimate of unknowns (betas) in rhs, NOT OK!

        // Solve the normal equations for all alpha_i
        // Simultaneous solution doesn't work somehow
        for (int i = 0; i < rhs.getColumns(); ++i) {
            DoubleMatrix rhs_alphai = rhs.getColumn(i);
            rhs.putColumn(i, Solve.solveSymmetric(Qx_hat, rhs_alphai));
//    		solvechol(Qx_hat, rhs_alphai); // Solution in rhs_alphai
//    		rhs.setcolumn(i, rhs_alphai); // place solution back
        }

        // Test solution by inverse
        Qx_hat = Solve.solveSymmetric(Qx_hat, DoubleMatrix.eye(Qx_hat.getRows()));
        double maxdev = (N.mmul(Qx_hat).sub(DoubleMatrix.eye(Qx_hat.getRows()))).normmax();
        logger.debug("s2h schwaebisch: max(abs(N*inv(N)-I)) = {}", maxdev);
        if (maxdev > 0.01) {
            logger.warn("slant2h: possibly wrong solution. deviation from unity AtA*inv(AtA) = {} > 0.01", maxdev);
        }

        /** -- STEP 4 : compute height for all pixels in N points for nHeights ---------*/

        logger.debug("S2H: schwabisch: STEP4: compute height for all pixels.");
        // Evaluate for all points interferogram h=f(l,p,phase)
        //  .....recon with multilook, degree1D, degree2D free
        //  .....Multilook factors
        double multiL = unwrappedinterf.multilookL;
        double multiP = unwrappedinterf.multilookP;

        // Number of lines/pixels of multilooked unwrapped interferogram
        int mllines = (int) (Math.floor((unwrappedinterf.win.linehi - unwrappedinterf.win.linelo + 1) / multiL));
        int mlpixels = (int) (Math.floor((unwrappedinterf.win.pixhi - unwrappedinterf.win.pixlo + 1) / multiP));

        // Line/pixel of first point in original master coordinates
        double veryfirstline = (double) (unwrappedinterf.win.linelo) + (multiL - 1.) / 2.;
        double firstpixel = (double) (unwrappedinterf.win.pixlo) + (multiP - 1.) / 2.;

        // ant axis of pixel coordinates ______
        DoubleMatrix p_axis = new DoubleMatrix(mlpixels, 1);
        for (int i = 0; i < mlpixels; i++) {
            p_axis.put(i, 0, firstpixel + i * multiP);
        }
        normalize(p_axis, minP, maxP);

//        int NUMMAT = 1 + (1 + DEGREE1D); // number of heavy matrices
        int bufferlines = mllines; //generalinput.memory / (NUMMAT * (mlpixels * sizeof(real4)));
        if (bufferlines > mllines) // whole image fits in BUFFER
            bufferlines = mllines;

        int FULLBUFFERS = mllines / bufferlines;
        int RESTLINES = mllines % bufferlines;
        int EXTRABUFFER = (RESTLINES > 0) ? 1 : 0;

        // ______ Window to be read into BUFFER from file in multilooked system ______
        int dummy = 999999; // large to force error if not ok
        Window bufferwin = new Window(1, bufferlines, 1, mlpixels); // initial
        Window offsetbuffer = new Window(1, dummy, 1, dummy); // dummy not used in readfile, no offset

        /** Process BUFFERS */
        for (int buffer = 1; buffer <= FULLBUFFERS + EXTRABUFFER; buffer++) {

            // In original master coordinate sytem
            double firstline = veryfirstline + (buffer - 1) * bufferlines * multiL;

            // ______ Set indices for loading / check last buffer ______
            bufferwin.linelo = 1 + (buffer - 1) * bufferlines; // Update window to be read from file
            if (buffer == FULLBUFFERS + 1) {
                bufferlines = RESTLINES;
                //BUFFER.resize(bufferlines,mlpixels);
            }
            bufferwin.linehi = bufferwin.linelo + bufferlines - 1; // window 2b read from file

            // ______ Read in phase buffer of unwrapped interferogram ______
            DoubleMatrix BUFFER = unwrappedinterf.readphase(bufferwin);

            // Evaluate polynomial coefficients for these points
            // Compute first line of current buffer in master coordinates
            DoubleMatrix l_axis = new DoubleMatrix(bufferlines, 1);
            for (int k = 0; k < bufferlines; k++) {
                l_axis.put(k, 0, firstline + k * multiL);
            }
            normalize(l_axis, minL, maxL);

            // ---> Lookup table because not known in advance what DEGREE1D is
            DoubleMatrix[] pntALPHA = new DoubleMatrix[TEN];
            for (int k = 0; k <= DEGREE1D; k++) {
                DoubleMatrix beta = new DoubleMatrix(PolyUtils.numberOfCoefficients(DEGREE2D), 1);
                for (int l = 0; l < PolyUtils.numberOfCoefficients(DEGREE2D); l++) {
                    beta.put(l, 0, rhs.get(l, k)); // solution stored in rhs
                }
                pntALPHA[k] = PolyUtils.polyval(l_axis, p_axis, beta, DEGREE2D);
            }

            // Evaluate h=f(l,p,phi) for all points in grid in BUFFER
            DoubleMatrix coeff_thispoint = new DoubleMatrix(DEGREE1D + 1, 1);

            for (int line = 0; line < BUFFER.getRows(); line++) {
                for (int pixel = 0; pixel < BUFFER.getColumns(); pixel++) {
                    // Check if unwrapped ok, else compute h
                    if (BUFFER.get(line, pixel) != Double.NaN) // else leave NaN
                    {
                        for (int k = 0; k < DEGREE1D + 1; k++) {
                            coeff_thispoint.put(k, 0, pntALPHA[k].get(line, pixel));
                        }
                        BUFFER.put(line, pixel, PolyUtils.polyVal1D(PolyUtils.normalize2(BUFFER.get(line, pixel), minphi, maxphi), coeff_thispoint.toArray()));
                    }
                }
            }

            // Write computed heights to file
//            ofile << BUFFER;
            // ______ new Matrix should be deleted ______
            // ______ if errors occur sigsegv maybe because of this ______
            logger.debug("deleting new matrix, memory errors could be caused by this");
//            for (int k = 0; k <= DEGREE1D; k++)
//                delete pntALPHA[ k];// correct?
        }// loop over BUFFERS

    } // END slant2h


    private void normalize(DoubleMatrix data, double min, double max) {
        data.mini(.5 * (min + max));
        data.divi(.25 * (max - min));
    }

}
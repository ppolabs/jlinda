package org.jdoris.core.filtering;

import org.apache.log4j.Logger;
import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jdoris.core.MathUtilities;
import org.jdoris.core.SLCImage;
import org.jdoris.core.todo_classes.todo_classes;

import static org.jblas.MatrixFunctions.pow;

/**
 * User: pmar@ppolabs.com
 * Date: 4/8/11
 * Time: 5:01 PM
 */
public class AzimuthFilter {

    static Logger logger = Logger.getLogger(AzimuthFilter.class.getName());

    /**
     * *************************************************************
     * azimuthfilter                                             *
     * Loop over whole master and slave image and filter out        *
     * part of the spectrum that is not common.                     *
     * Only do zero doppler freq. offset.                           *
     * do not use a polynomial from header for now.                 *
     * (next we will, but assume image are almost coreg. in range,  *
     * so f_dc polynomial can be eval. same)                       *
     * Per block in azimuth [1024] use a certain overlap with the   *
     * next block so that same data is partially used for spectrum  *
     * (not sure if this is requried).                              *
     * Filter is composed of: DE-hamming, RE-hamming (for correct   *
     * new size and center of the spectrum).                        *
     * Trick in processor.c: First call routine as:                 *
     * (generalinput,filtaziinput,master,slave)                    *
     * in order to process the master, and then as:                 *
     * (generalinput,filtaziinput,slave,master)                    *
     * to filter the slave slc image.
     */
    public static void azimuthfilter(final todo_classes.inputgeneral generalinput,
                                     final todo_classes.input_filtazi fitaziinput,
                                     SLCImage master, // not const, fdc possibly reset here?
                                     SLCImage slave) {

    }

    /**
     * *************************************************************
     * azimuth filter per block                                  *
     * Input is matrix of SIZE (e.g. 1024) lines, and N range pixs. *
     * Input is SLC of master. slave_info gives fDC polynomial      *
     * for slave + coarse offset. HAMMING is alpha for myhamming f. *
     * Filtered OUTPUT is same size as input block.                 *
     * Because of overlap (azimuth), only write to disk in calling  *
     * routine part (in matrix coord.) [OVERLAP:SIZE-OVERLAP-1]     *
     * = SIZE-(2*OVERLAP);  // number of output pixels              *
     * *
     * Filtering is performed in the spectral domain                *
     * (1DFFT over azimuth for all columns at once)                 *
     * Filter is different for each column due to shift in fd_c     *
     * doppler centroid frequency.                                  *
     * *
     * ! It should still be studied if real4 matrices are accurate  *
     * enough, but I guess it is (BK).                              *
     * *
     */
    public ComplexDoubleMatrix blockazifilt(
            final ComplexDoubleMatrix SLCIMAGE,
            final SLCImage master,          // PRF, BW, fd0
            final SLCImage slave,           // PRF, BW, fd0
            final double HAMMING) throws Exception {

//        return null;

        final long SIZE = SLCIMAGE.rows;  // fftlength
        final long NCOLS = SLCIMAGE.columns; // width
        if (NCOLS != master.getCurrentWindow().pixels())
            logger.warn("this will crash, size input matrix not ok...");

        // ______ Compute fDC_master, fDC_slave for all columns ______
        // ______ Create axis to evaluate fDC polynomial for master/slave ______
        // ______ fDC(column) = fdc_a0 + fDC_a1*(col/RSR) + fDC_a2*(col/RSR)^2 ______
        // ______ fDC = y = Ax ______
        // ______ Capitals indicate matrices (FDC_M <-> fDC_m) ______
        logger.debug("Filtering data by evaluated polynomial fDC for each column.");
        DoubleMatrix xaxis = new DoubleMatrix(1, (int) master.getCurrentWindow().pixels());         // lying

        // TODO: refactor to more efficient jblass call :: "linspace"
        for (long i = master.getCurrentWindow().pixlo; i <= master.getCurrentWindow().pixhi; ++i)
            xaxis.put(0, (int) (i - master.getCurrentWindow().pixlo), i - 1.0);

        xaxis.divi(master.getRsr2x() / 2.0);

        // TODO: better use SLCImage.pix2fdc()
        DoubleMatrix FDC_M = xaxis.mul(master.getF_DC_a1());
        FDC_M.addi(master.getF_DC_a0());
        FDC_M.addi(pow(xaxis, 2).mmul(master.getF_DC_a2()));

        // ______ fDC_slave for same(!) columns (coarse offset). ______
        // ______ offset defined as: cols=colm+offsetP ______
        for (long i = master.getCurrentWindow().pixlo; i <= master.getCurrentWindow().pixhi; ++i)
            xaxis.put(0, (int) master.getCurrentWindow().pixlo, i - 1.0 + slave.getCoarseOffsetP());

        xaxis.divi(slave.getRsr2x() / 2.0);
        DoubleMatrix FDC_S = xaxis.mul(slave.getF_DC_a1());
        FDC_S.addi(slave.getF_DC_a0());
        FDC_S.addi(pow(xaxis, 2).mmul(slave.getF_DC_a2()));

        logger.debug("Dumping matrices fDC_m, fDC_s (__DEBUG defined)");
        logger.debug("fDC_m: " + FDC_M.toString());
        logger.debug("fDC_s: " + FDC_S.toString());

        // ______ Axis for filter in frequencies ______
        // TODO check, rather shift, test matlab... or wshift,1D over dim1
        // use fft properties to shift...

        final boolean dohamming = (HAMMING < 0.9999) ? true : false;
        final double PRF = master.getPRF();               // pulse repetition freq. [Hz]
        final double ABW = master.getAzimuthBandwidth();  // azimuth band width [Hz]

        final float deltaf = (float) (PRF / SIZE);
        final float fr = (float) (-PRF / 2.0);
        DoubleMatrix freqaxis = new DoubleMatrix(1, (int) SIZE);
        for (int i = 0; i < SIZE; ++i)
            freqaxis.put(0, i, fr + (i * deltaf)); // [-fr:df:fr-df]

        DoubleMatrix FILTER;                         // i.e., the filter per column
        DoubleMatrix FILTERMAT = new DoubleMatrix((int) SIZE, (int) NCOLS);          // i.e., THE filter

        for (long i = 0; i < NCOLS; ++i) {
            final double fDC_m = FDC_M.get(0, (int) i);          // zero doppler freq. [Hz]
            final double fDC_s = FDC_S.get(0, (int) i);          // zero doppler freq. [Hz]
            final double fDC_mean = 0.5 * (fDC_m + fDC_s);   // mean doppler centroid freq.
            final double ABW_new = Math.max(1.0, 2.0 * (0.5 * ABW - Math.abs(fDC_m - fDC_mean)));       // new bandwidth > 1.0

            if (dohamming) {
                // ______ NOT a good implementation for per col., cause wshift AND fftshift.
                // ______ DE-weight spectrum at centered at fDC_m ______
                // ______ spectrum should be periodic! (use wshift) ______
                DoubleMatrix inversehamming = MathUtilities.myhamming(freqaxis, ABW, PRF, HAMMING);
                for (long ii = 0; ii < SIZE; ++ii)
                    inversehamming.put(0, (int) ii, (float) (1.0 / inversehamming.get(0, (int) ii)));

                // ______ Shift this circular by myshift pixels ______
                long myshift = (long) (Math.rint((SIZE * fDC_m / PRF)));// round
                MathUtilities.wshift(inversehamming, (int) -myshift);          // center at fDC_m

                // ______ Newhamming is scaled and centered around new mean ______
                myshift = (long) (Math.rint((SIZE * fDC_mean / PRF)));// round
                FILTER = MathUtilities.myhamming(freqaxis, ABW_new, PRF, HAMMING);         // fftshifted
                MathUtilities.wshift(FILTER, (int) -myshift);                  // center at fDC_mean
                FILTER.mmuli(inversehamming);
            } else {       // no weighting, but center at fDC_mean, size ABW_new
                long myshift = (long) (Math.rint((SIZE * fDC_mean / PRF)));// round
                FILTER = MathUtilities.myrect(freqaxis.divi((float) ABW_new)); // fftshifted
                MathUtilities.wshift(FILTER, (int) -myshift);                  // center at fDC_mean
            }

            MathUtilities.ifftshift(FILTER);                          // fftsh works on data!

            FILTERMAT.putColumn((int) i, FILTER);
        } // foreach column


        // ______ Filter slcdata ______
        ComplexDoubleMatrix FILTERED = SLCIMAGE.dup();
        MathUtilities.fft(FILTERED, 1);                              // fft foreach column
        FILTERED.mmuli(new ComplexDoubleMatrix(FILTERMAT));
        MathUtilities.ifft(FILTERED, 1);                             // ifft foreach column
        return FILTERED;

    }

}

package org.jdoris.core.filtering;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;
import org.apache.log4j.Logger;
import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
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
                DoubleMatrix inversehamming = myhamming(freqaxis, ABW, PRF, HAMMING);
                for (long ii = 0; ii < SIZE; ++ii)
                    inversehamming.put(0, (int) ii, (float) (1.0 / inversehamming.get(0, (int) ii)));

                // ______ Shift this circular by myshift pixels ______
                long myshift = (long) (Math.rint((SIZE * fDC_m / PRF)));// round
                wshift(inversehamming, (int) -myshift);          // center at fDC_m

                // ______ Newhamming is scaled and centered around new mean ______
                myshift = (long) (Math.rint((SIZE * fDC_mean / PRF)));// round
                FILTER = myhamming(freqaxis, ABW_new, PRF, HAMMING);         // fftshifted
                wshift(FILTER, (int) -myshift);                  // center at fDC_mean
                FILTER.mmuli(inversehamming);
            } else {       // no weighting, but center at fDC_mean, size ABW_new
                long myshift = (long) (Math.rint((SIZE * fDC_mean / PRF)));// round
                FILTER = myrect(freqaxis.divi((float) ABW_new)); // fftshifted
                wshift(FILTER, (int) -myshift);                  // center at fDC_mean
            }

            ifftshift(FILTER);                          // fftsh works on data!

            FILTERMAT.putColumn((int) i, FILTER);
        } // foreach column


        // ______ Filter slcdata ______
        ComplexDoubleMatrix FILTERED = SLCIMAGE.dup();
        fft(FILTERED, 1);                              // fft foreach column
        FILTERED.mmuli(new ComplexDoubleMatrix(FILTERMAT));
        ifft(FILTERED, 1);                             // ifft foreach column
        return FILTERED;

    }

    /**
     * *************************************************************
     * ifft(A,dim)                                                  *
     * inverse 1dfft over dim of A is returned in A by reference *
     * if dim=1 ifft is over all columns of A, if 2 over rows.   *
     * data is stored major row order in memory, so dim=2 is     *
     * probably much faster.                                     *
     * fftlength should be power of 2                            *
     * Bert Kampes, 22-Mar-2000                                  *
     * **************************************************************
     */
    private void ifft(ComplexDoubleMatrix A, int dimension) {
        int i;
        final int iopt = -1;                        // inverse FFT (scaled)

        switch (dimension) {
            case 1: {
                logger.debug("1d ifft over columns");
                int fftlength = A.rows;
                for (i = 0; i < A.columns; ++i) {
                    ComplexDoubleMatrix VECTOR = A.getColumn(i);
                    four1(VECTOR, fftlength, iopt);
                    A.putColumn(i, VECTOR);
                }
                break;
            }
            case 2: {
                logger.debug("1d ifft over rows");
                int fftlength = A.columns;

                for (i = 0; i < A.rows; ++i) {
                    ComplexDoubleMatrix VECTOR = A.getRow(i);
                    four1(VECTOR, fftlength, iopt);
                    A.putRow(i, VECTOR);
                }
                break;
            }
            default:
                logger.error("ifft: dimension != {1,2}");
        }
    } // END ifft


    private void fft(ComplexDoubleMatrix A, int dimension) {

        int i;
        final int iopt = 1;                         // forward FFT
        int fftlength;
        ComplexDoubleMatrix VECTOR = null;
        switch (dimension) {
            case 1: {
                fftlength = A.rows;
                for (i = 0; i < A.columns; ++i) {
                    VECTOR = A.getColumn(i);
                    four1(VECTOR, fftlength, iopt);// but generic.
                    A.putColumn(i, VECTOR);
                    // perhaps can be used directly to do this w/o data copying...
                    // four1(A.getColumn(i), fftlength, iopt);// but generic.
                }
                break;
            }
            case 2: {
                fftlength = A.columns;
                for (i = 0; i < A.rows; ++i) {
                    VECTOR = A.getRow(i);
                    four1(VECTOR, fftlength, iopt);// but generic.
                    A.putRow(i, VECTOR);
//                    four1(A.getRow(i), fftlength, iopt);
                }
                break;
            }
            default:
                logger.error("fft: dimension != {1,2}");
        }

    }

    /**
     * *************************************************************
     * four1(complr4 *, length, isign)                              *
     * four1(&A[0][0], 128, 1)                                      *
     * either based on numerical recipes or veclib                  *
     * helper function for other fft routines, if no veclib         *
     * cooley-turkey, power 2, replaces input,                      *
     * isign=1: fft , isign=-1: ifft                                *
     * handling vectors should be simpler (lying, standing)         *
     * note that this is not a good implementation, only to get     *
     * doris software working without veclib.                       *
     * *
     * define SAMEASVECIB if you want the order of the coefficients *
     * of the fft the same as veclib. it seems this is not required *
     * for a good version of Doris, but in case of problems this    *
     * may be the solution.                                         *
     * *
     * VECLIB defines the FT same as matlab:                        *
     * N-1                                                 *
     * X(k) = sum  x(n)*exp(-j*2*pi*k*n/N), 0 <= k <= N-1.        *
     * n=0                                                 *
     * *
     * FFTW defines the same as Matlab, but inv. not normalized.    *
     * I don't know if the matrix must be allocated somehow, so for *
     * now we try only 1d ffts to build 2d too.                     *
     * *
     * **************************************************************
     */
    private void four1(ComplexDoubleMatrix vector, int fftlength, int direction) {

        DoubleFFT_1D fft = new DoubleFFT_1D(fftlength);
        double[] tempDoubleArray;

        switch (direction) {
            case 1:
                tempDoubleArray = vector.toDoubleArray();
                fft.complexForward(tempDoubleArray);
                vector.copy(new ComplexDoubleMatrix(tempDoubleArray));

            case -1:
                tempDoubleArray = vector.toDoubleArray();
                fft.complexInverse(tempDoubleArray, false);
                vector.copy(new ComplexDoubleMatrix(tempDoubleArray));
        }

    }

    /**
     * *************************************************************
     * ifftshift(A)                                                 *
     * ifftshift of vector A is returned in A by reference       *
     * undo effect of fftshift. ?p=floor(m/2); A=A[p:m-1 0:p-1]; *
     * **************************************************************
     */
    private static void ifftshift(DoubleMatrix A) throws Exception {

        if (!A.isVector()) {
            logger.error("ifftshift: only vectors");
            throw new Exception();
        }

        DoubleMatrix Res = new DoubleMatrix(A.rows, A.columns);
        final int start = (int) (Math.floor((float) (A.length) / 2));

        System.arraycopy(A.data, start, Res.data, 0, A.length - start);
        System.arraycopy(A.data, 0, Res.data, A.length - start, start);

        A.copy(Res);

    } // END ifftshift

    /**
     * *************************************************************
     * wshift(A,n)                                                  *
     * circular shift of vector A by n pixels. positive n for    *
     * right to left shift.                                      *
     * implementation: WSHIFT(A,n) == WSHIFT(A,n-sizeA);         *
     * A is changed itself!                                      *
     * **************************************************************
     */

    public static void wshift(DoubleMatrix A, int n) throws Exception {

        if (n >= A.length) {
            System.err.println("wshift: shift larger than matrix not implemented.");
            throw new Exception();
        }

        if (!A.isVector()) {
            System.err.println("wshift: only vectors");
            throw new Exception();
        }

        // positive only, use rem!  n = n%A.nsize;
        if (n == 0) return;
        if (n < 0) n += A.length;

        DoubleMatrix Res = new DoubleMatrix(A.rows, A.columns);

        //  n always >0 here
        System.arraycopy(A.data, n, Res.data, 0, A.length - n);
        System.arraycopy(A.data, 0, Res.data, A.length - n, n);

        A.copy(Res);

    }


    private static DoubleMatrix myrect(DoubleMatrix X) throws Exception {

        if (X.rows != 1) {
            System.err.println("myrect: only lying vectors.");
            throw new Exception();
        }

        DoubleMatrix Res = new DoubleMatrix(1, X.rows);

        for (int i = 0; i < X.rows; ++i) {
            if (Math.abs(X.get(i)) <= 0.5) {
                Res.put(i, (float) 1.);
            }
        }

        return Res;
    }

    /**
     * *************************************************************
     * myhamming                                                    *
     * hamming window, lying vector                                 *
     * w = (a + (1.-a).*cos((2.*pi/fs).*fr)) .* myrect(fr./Br);     *
     * scale/shift filter by g(x)=f((x-xo)/s)                       *
     * alpha==1 yields a myrect window                              *
     * Bert Kampes, 31-Mar-2000                                  *
     * **************************************************************
     */
    private static DoubleMatrix myhamming(final DoubleMatrix fr, double RBW, double RSR, double alpha) throws Exception {

        if (fr.rows != 1) {
            System.err.println("myhamming: only lying vectors.");
            throw new Exception();
//            throw (argument_error);
        }

        if (alpha < 0.0 || alpha > 1.0) {
            System.err.println("myhamming: !alpha e{0..1}.");
            throw new Exception();
        }

        if (RBW > RSR) {
            System.err.println("myhamming: RBW>RSR.");
            throw new Exception();
        }

        DoubleMatrix Res = new DoubleMatrix(1, fr.columns);
        for (int i = 0; i < fr.columns; ++i) {
            if (Math.abs(fr.get(i)) < 0.5) {   // rect window
                Res.put(i, (float) (alpha + (1 - alpha) * Math.cos((2 * Math.PI / RSR) * fr.get(i))));
            }
        }
        return Res;
    }

}

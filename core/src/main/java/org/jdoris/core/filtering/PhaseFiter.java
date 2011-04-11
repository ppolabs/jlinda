package org.jdoris.core.filtering;

import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jdoris.core.todo_classes.todo_classes;

/**
 * User: pmar@ppolabs.com
 * Date: 4/8/11
 * Time: 5:01 PM
 */
public class PhaseFiter {

    //TODO: make template classes for generalInput, operatorInput, and ProductMetadata class

    /**
     * *************************************************************
     * phasefilter                                               *
     * goldsteins method, see routine goldstein and smooth.         *
     * After Goldstein and Werner, Radar interferogram filtering    *
     * for geophysical applications. GRL 25-21 pp 4035-4038, 1998.  *
     * and: ESA Florence 1997, vol2, pp969-972, Goldstein & Werner  *
     * "Radar ice motion interferometry".                           *
     * *************************************************************
     */
    public static void phasefilter(
            final todo_classes.inputgeneral generalinput,
            final todo_classes.productinfo interferogram,
            final todo_classes.input_filtphase filtphaseinput) {
    }

    /*
   *    phasefilter goldstein                                     *
   * Input is matrix of SIZE (e.g. 32) lines, and N range pixels. *
   * Filtered OUTPUT is same size as input block.                 *
   * Because of overlap, only write to disk in calling routine    *
   * part (in matrix coord.) [OVERLAP:SIZE-OVERLAP-1]             *
   *                                                              *
   * Smoothing of the amplitude of the spectrum is performed by   *
   * spatial convolution with a block kernel of size 2*SMOOTH+1.  *
   * (Which is done by FFT's). e.g. a spatial moving average with *
   * kernel (1d) k=[1 1 1 1 1]/5; kernel2d = transpose(k)*k.      *
   * Blocks in range direction.                                   *
   *                                                              *
   * After Goldstein and Werner, Radar interferogram filtering    *
   * for geophysical applications. GRL 25-21 pp 4035-4038, 1998.  *
   * and: ESA Florence 1997, vol2, pp969-972, Goldstein & Werner  *
   * "Radar ice motion interferometry".                           *
   * */
    public static ComplexDoubleMatrix goldstein(
            final ComplexDoubleMatrix CINT,
            final float ALPHA,
            final long OVERLAP,
            final DoubleMatrix smoothkernel) { // lying down
        return null;
    }

    /*
     *    phasefilter buffer by spatial conv. with kernel.          *
     * Input is matrix of SIZE (e.g. 256) lines, and N range pixels.*
     * Filtered OUTPUT is same size as input block.                 *
     * Because of overlap, only write to disk in calling routine    *
     * part (in matrix coord.) [OVERLAP:SIZE-OVERLAP-1]             *
     * (in line direction)                                          *
     * spatial convolution with a kernel function, such as a block  *
     * function 111 (1D) (By FFT's).                                *
     * Processing is done in blocks in range direction.             *
     * For the first block the part [0:OVERLAP-1] is set to 0.      *
     * For the last block the part [NPIX-1-OVERLAP:NPIX-1] is 0.    *
     *                                                              *
     * Input:                                                       *
     *  - matrix to be filtered of blocklines * numpixs             *
     *  - kernel2d: fft2 of 2d spatial kernel.                      *
     *  - overlap: half of the kernel size, e.g., 1 for 111.        *
     * Output:                                                      *
     *  - filtered matrix.                                          *
     *    ifft2d(BLOCK .* KERNEL2D) is returned, so if required for *
     *    non symmetrical kernel, offer the conj(KERNEL2D)!         *
     *                                                              *
    */
    public static ComplexDoubleMatrix convbuffer(
            final ComplexDoubleMatrix CINT,
            final ComplexDoubleMatrix KERNEL2D,
            final long OVERLAP) {         // overlap in column direction

        return null;

    }

    /**
     * *************************************************************
     * spatialphasefilt                                          *
     * For the first block the part [0:OVERLAP-1] is set to 0.      *
     * For the last block the part [NPIX-1-OVERLAP:NPIX-1] is 0.    *
     * *************************************************************
     */
    public static void spatialphasefilt(
            final todo_classes.inputgeneral generalinput,
            final todo_classes.productinfo interferogram,
            final todo_classes.input_filtphase filtphaseinput) {

    }

    /**
     * phasefilter spectral                                      *
     * Input is matrix of SIZE (e.g. 32) lines, and N range pixels. *
     * Filtered OUTPUT is same size as input block.                 *
     * Because of overlap, only write to disk in calling routine    *
     * part (in matrix coord.) [OVERLAP:SIZE-OVERLAP-1]             *
     * *
     * Filtering is performed by pointwise multiplication of the    *
     * spectrum per block by the KERNEL2D (input).                  *
     * Blocks in range direction,                                   *
     */

    public static ComplexDoubleMatrix spectralfilt(
            final ComplexDoubleMatrix CINT,
            final ComplexDoubleMatrix KERNEL2D,
            final long OVERLAP) {

        return null;
    }

    /**
     * B = smooth(A,KERNEL)                                         *
     * (circular) spatial moving average with a (2N+1,2N+1) block.  *
     * See also matlab script smooth.m for some tests.              *
     * implementation as convolution with FFT's                     *
     * input: KERNEL is the FFT of the kernel (block)
     */
    public static ComplexDoubleMatrix smooth(
            final DoubleMatrix A,
            final ComplexDoubleMatrix KERNEL2D) {
        return null;
    }

    /**
     * B = smooth(A,blocksize)                                      *
     * (circular) spatial moving average with a (2N+1,2N+1) block.  *
     * See also matlab script smooth.m for some tests.
     */
    public static DoubleMatrix smooth(
            final DoubleMatrix A,
            long N) {
        return null;
    }

    /**
     * phasefilterspectral                                       *
     * loop over whole file and multiply spectrum of interferogram  *
     * with kernel specified in input file.
     */
    public static void phasefilterspectral(
            final todo_classes.inputgeneral generalinput,
            final todo_classes.productinfo interferogram,
            final todo_classes.input_filtphase filtphaseinput) {


    }

}

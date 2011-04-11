package org.jdoris.core.filtering;

import org.jblas.ComplexDoubleMatrix;
import org.jdoris.core.SLCImage;
import org.jdoris.core.todo_classes.todo_classes;

/**
 * User: pmar@ppolabs.com
 * Date: 4/8/11
 * Time: 5:01 PM
 */
public class AzimuthFilter {

    /****************************************************************
     *    azimuthfilter                                             *
     * Loop over whole master and slave image and filter out        *
     * part of the spectrum that is not common.                     *
     * Only do zero doppler freq. offset.                           *
     * do not use a polynomial from header for now.                 *
     * (next we will, but assume image are almost coreg. in range,  *
     *  so f_dc polynomial can be eval. same)                       *
     * Per block in azimuth [1024] use a certain overlap with the   *
     * next block so that same data is partially used for spectrum  *
     * (not sure if this is requried).                              *
     * Filter is composed of: DE-hamming, RE-hamming (for correct   *
     * new size and center of the spectrum).                        *
     * Trick in processor.c: First call routine as:                 *
     *  (generalinput,filtaziinput,master,slave)                    *
     * in order to process the master, and then as:                 *
     *  (generalinput,filtaziinput,slave,master)                    *
     * to filter the slave slc image.                               */
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
    public static ComplexDoubleMatrix blockazifilt(
            final ComplexDoubleMatrix SLCIMAGE,
            final SLCImage master,          // PRF, BW, fd0
            final SLCImage slave,           // PRF, BW, fd0
            final double HAMMING) {

        return null;

    }

}

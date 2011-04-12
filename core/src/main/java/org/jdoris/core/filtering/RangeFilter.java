package org.jdoris.core.filtering;

import org.apache.log4j.Logger;
import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.ranges.IntervalRange;
import org.jdoris.core.MathUtilities;
import org.jdoris.core.Orbit;
import org.jdoris.core.SLCImage;
import org.jdoris.core.todo_classes.input_ell;
import org.jdoris.core.todo_classes.todo_classes;

import static org.jblas.MatrixFunctions.pow;

public class RangeFilter {

    static Logger logger = Logger.getLogger(RangeFilter.class.getName());


    //TODO: make template classes for generalInput, operatorInput, and ProductMetadata class
    public static void rangefilter(final todo_classes.inputgeneral input_gen,
                                   final SLCImage master,
                                   final SLCImage slave,
                                   final todo_classes.productinfo interferogram,
                                   final todo_classes.input_filtrange inputfiltrange) {
    }

    /**
     * rfilterblock
     * Computes powerspectrum of complex interferogram.
     * (product oversampled master. conj oversampled slave in range)
     * A peak in this spectrum corresponds to frequency shift.
     * The master and slave are LPF filtered for this shift.
     * *
     * Optionally the oversampling can be turned off, since no use
     * if only small baseline, and flat terrain.
     * The powerspectrum can be weighted to give more influence to
     * higher frequencies (conv. of 2 blocks, should be hamming).
     * The peak is detected by taking a mean of nlmean lines (odd).
     * Filtering is applied if the SNR (N*power peak / power rest)
     * is above a user supplied threshold.
     * At LPF filtering of the master/slave a hamming window may
     * be applied first to deweight, then to re-weight the spectrum
     * <p/>
     * Should filter based on zero terrain slope if below SNR,
     * but this requried knowledge of orbits, pixel,line coordinate
     * #%// BK 13-Nov-2000
     * *
     * Input:
     * - MASTER: block of master, that will be filtered
     * - SLAVE:  block of slave, that will be filtered
     * Output:
     * - MASTER (SLAVE): filtered from indeces[0:numl-1]
     * (nlmean-1)/2 to numlines-(nlmean-1)/2-1
     */
    public static void rfilterblock(ComplexDoubleMatrix master, // updated
                                    ComplexDoubleMatrix slave,  // updated
                                    long nlmean, float SNRthreshold,
                                    float RSR, // in MHz
                                    float RBW, // in MHz
                                    float alphahamming,
                                    long osfactor,
                                    boolean doweightcorrel,
                                    double meanSNR, // returned
                                    double percentnotfiltered) throws Exception { // returned


        long numlines = master.rows;
        long numpixs = master.columns;
        long outputlines = numlines - nlmean + 1;
        long firstline = (long) ((nlmean - 1) / 2);        // indices in matrix system
        long lastline = firstline + outputlines - 1;
        boolean dohamming = (alphahamming < 0.9999) ? true : false;
        // use oversampling before int. gen.
        boolean dooversample = (osfactor != 1) ? true : false;
        int notfiltered = 0;                                  // counter

        if (!MathUtilities.isodd(nlmean)) {
            logger.error("nlmean has to be odd.");
            throw new Exception();
        }
        if (!MathUtilities.ispower2(numpixs)) {

            logger.error("numpixels (FFT) has to be power of 2.");

            throw new Exception();
        }
        if (!MathUtilities.ispower2(osfactor)) {
            logger.error("oversample factor (FFT) has to be power of 2.");
            throw new Exception();
        }
        if (slave.rows != numlines) {
            logger.error("slave not same size as master.");
            throw new Exception();
        }
        if (slave.columns != numpixs) {
            logger.error("slave not same size as master.");
            throw new Exception();
        }
        if (outputlines < 1) {
            logger.warn("no outputlines, continuing.");
        }


        // SHIFT PARAMETERS
        int i, j;
        final double deltaf = RSR / numpixs;
        final double fr = -RSR / 2.;
        DoubleMatrix freqaxis = new DoubleMatrix(1, (int) numpixs);
        for (i = 0; i < numpixs; ++i) {
            freqaxis.put(0, i, fr + (i * deltaf));
        }

        DoubleMatrix inversehamming = null;
        if (dohamming) {
            inversehamming = MathUtilities.myhamming(freqaxis, RBW, RSR, alphahamming);
            for (i = 0; i < numpixs; ++i)
                if (inversehamming.get(0, i) != 0.)
                    inversehamming.put(0, i, 1. / inversehamming.get(0, i));
        }

        // COMPUTE CPLX IFG ON THE FLY -> power
        ComplexDoubleMatrix cint;
        if (dooversample) {
            cint = MathUtilities.dotmult(MathUtilities.oversample(slave.conj(), 1, (int) osfactor),
                    MathUtilities.oversample(master, 1, (int) osfactor));
        } else {
            cint = MathUtilities.dotmult(master, slave.conj());
        }


        long fftlength = cint.columns;

        logger.debug("is real4 accurate enough?");// seems so

        MathUtilities.fft(cint, 2);                 // cint=fft over rows
        DoubleMatrix power = MathUtilities.intensity(cint);      // power=cint.*conj(cint);

        // ______ Use weighted correlation due to bias in normal definition ______
        // ______ Actually better deweight with autoconvoluted hamming.
        // ______ No use a triangle for #points used for correlation estimation
        // ______ not in combination with dooversample...
        if (doweightcorrel) {

            // TODO: refactor this call to use arrays instead of loops

            //matrix<real4> weighting(1,fftlength);
            //for (i=0; i<fftlength; ++i)
            //  weighting(0,i) = abs(i)...;
            //power *= weightingfunction...== inv.triangle

            // weigth = numpoints in spectral convolution for fft squared for power...
            int indexnopeak = (int) ((1. - (RBW / RSR)) * (float) (numpixs));
            for (j = 0; j < fftlength; ++j) {

                long npnts = Math.abs(numpixs - j);
                // oversample: numpixs==fftlength/2
                double weight = (npnts < indexnopeak) ? Math.pow(numpixs, 2) : Math.pow(npnts, 2); // ==zero
                for (i = 0; i < numlines; ++i) {
//                  power(i, j) /= weight;
                    power.put(i, j, power.get(i, j) / weight);
                }
            }
        }

        // ______ Average power to reduce noise ______
        MathUtilities.fft(master, 2);                                // master=fft over rows
        MathUtilities.fft(slave, 2);                                 // slave=fft over rows
        logger.trace("Took FFT over rows of master, slave.");


        IntervalRange rangeRows = new IntervalRange(0, (int) (nlmean - 1));
        IntervalRange rangeColumns = new IntervalRange(0, (int) (fftlength - 1));


//        DoubleMatrix nlmeanpower = sum(power(0,nlmean-1, 0,fftlength-1),1);
        DoubleMatrix nlmeanpower = pow(power.get(rangeRows, rangeColumns), 2).rowSums();

        long shift = 0;                          // returned by max
        long dummy = 0;                          // returned by max
        meanSNR = 0.;
        double meanSHIFT = 0.;


        // ______ Start actual filtering ______
        // ______ oline is index in matrix system ______
        for (long oline = firstline; oline <= lastline; ++oline) {

            DoubleMatrix totalp = nlmeanpower.columnSums();        // 1x1 matrix ...
            double totalpower = totalp.get(0, 0);

            // TODO: check algorithmically this step
            // double maxvalue = max(nlmeanpower, dummy, shift);      // shift returned
            double maxvalue = nlmeanpower.max();

            long lastshift = shift;     // use this if current shift not ok.
            double SNR = fftlength * (maxvalue / (totalpower - maxvalue));
            meanSNR += SNR;

            // ______ Check for negative shift ______
            boolean negshift = false;
            if (shift > (int) (fftlength / 2)) {
                shift = (int) fftlength - shift;
                lastshift = shift; // use this if current shift not OK.
                negshift = true;
            }

            // ______ Do actual filtering ______
            if (SNR < SNRthreshold) {
                notfiltered++;                                    // update counter
                shift = lastshift;
                logger.warn("using last shift for filter");
            }
            meanSHIFT += shift;
            DoubleMatrix filter;

            if (dohamming) {
                // ______ Newhamming is scaled and centered around new mean ______
                filter = MathUtilities.myhamming(freqaxis.subi(0.5 * shift * deltaf),
                        RBW - (shift * deltaf),
                        RSR, alphahamming);          // fftshifted
                filter.mmul(inversehamming);
            } else { // no weighting of spectra
                filter = MathUtilities.myrect((freqaxis.subi(.5 * shift * deltaf)).divi((RBW - shift * deltaf)));   // fftshifted
            }

            // ______ Use freq. as returned by fft ______
            // ______ Note that filter_s = fliplr(filter_m) ______
            // ______ and that this is also valid after ifftshift ______
            MathUtilities.ifftshift(filter);                                // fftsh works on data!

            // ====== Actual spectral filtering ======
            // ______ Decide which side to filter, may be dependent on ______
            // ______ definition of FFT, this is ok for VECLIB ______
            // ______ if you have trouble with this step, either check your FFT
            // ______ (or use uinternal one, or add a card for changing false to true below
            if (negshift == false) {
//                dotmult(master[oline], filter, 1);
                MathUtilities.dotmult(master.getRow((int) oline), new ComplexDoubleMatrix(filter));
                MathUtilities.fliplr(filter);
                MathUtilities.dotmult(slave.getRow((int) oline), new ComplexDoubleMatrix(filter));
            } else {
                MathUtilities.dotmult(slave.getRow((int) oline), new ComplexDoubleMatrix(filter));
                MathUtilities.fliplr(filter);
                MathUtilities.dotmult(master.getRow((int) oline), new ComplexDoubleMatrix(filter));
            }
            // following is removed, we now always filter with last know spectral offset
            //      } // SNR>threshold
            //    else
            //      {
            //      notfiltered++;                                  // update counter
            //      }


            // ______ Update 'walking' mean ______
            if (oline != lastline) {                        // then breaks
                DoubleMatrix line1 = power.getRow((int) (oline - firstline));
                DoubleMatrix lineN = power.getRow((int) (oline - firstline + nlmean));
                nlmeanpower.add(lineN.sub(line1));
            }

        } // loop over outputlines

        // ______ IFFT of spectrally filtered data, and return these ______
        MathUtilities.ifft(master, 2);                               // master=ifft over rows
        MathUtilities.ifft(slave, 2);                                // slave=ifft over rows

        // ______ Return these to main ______
        meanSHIFT /= (outputlines - notfiltered);
        meanSNR /= outputlines;
        percentnotfiltered = 100. * (float) (notfiltered) / (float) outputlines;


        // ______ Some info for this block ______
        double meanfrfreq = meanSHIFT * deltaf;    // Hz?
        logger.debug("mean SHIFT for block"
                + ": " + meanSHIFT
                + " = " + meanfrfreq / 1e6 + " MHz (fringe freq.).");

        logger.debug("mean SNR for block: " + meanSNR);
        logger.debug("filtered for block"
                + ": " + (100.00 - percentnotfiltered) + "%");

        if (percentnotfiltered > 60.0) {
            logger.warn("more then 60% of singnal filtered?!?");
        }

    } // END rfilterblock

    // TODO: refactor InputEllips to "Ellipsoid" class of "org.esa.beam.framework.dataop.maptransf.Ellipsoid" and use GeoUtils of NEST;
    @Deprecated
    public static void rangefilterorbits(final todo_classes.inputgeneral generalinput,
                                         final todo_classes.input_filtrange inputfiltrange,
                                         final input_ell ellips,
                                         final SLCImage master,
                                         final SLCImage slave,
                                         Orbit masterorbit,
                                         Orbit slaveorbit) {
    }

}

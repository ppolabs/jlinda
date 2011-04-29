package org.jdoris.core.utils;

import org.apache.log4j.Logger;
import org.jblas.ComplexDouble;
import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;
import org.jdoris.core.Window;

public class SarUtils {

    static Logger logger = Logger.getLogger(SarUtils.class.getName());


    /**
     * B=oversample(A, factorrow, factorcol);
     * 2 factors possible, extrapolation at end.
     * no vectors possible.
     */
    public static ComplexDoubleMatrix oversample(ComplexDoubleMatrix AA, final int factorrow, final int factorcol) throws Exception {

        ComplexDoubleMatrix A = AA.dup(); // copy, AA is changed by in-place fft;

        final int l = A.rows;
        final int p = A.columns;
        final int halfl = l / 2;
        final int halfp = p / 2;
        final int L2 = factorrow * l;      // numrows of output matrix
        final int P2 = factorcol * p;      // columns of output matrix


        if (A.isVector()) {
            logger.error("OVERSAMPLE: only 2d matrices.");
            throw new Exception();
        }
        if (!MathUtils.ispower2(l) && factorrow != 1) {
            logger.error("OVERSAMPLE: numlines != 2^n.");
        }
        if (!MathUtils.ispower2(p) && factorcol != 1) {
            logger.error("OVERSAMPLE: numcols != 2^n.");
        }

        final ComplexDouble half = new ComplexDouble(0.5);

        ComplexDoubleMatrix Res = new ComplexDoubleMatrix(L2, P2);

//        int i, j;
        if (factorrow == 1) {

            // 1d fourier transform per row
            SpectralUtils.fft(A, 2);

            // divide by 2 because even fftlength
            A.putColumn(halfp, A.getColumn(halfp).mmuli(half));
//            for (i=0; i<l; ++i) {
//                A.put(i, halfp, A.get(i, halfp).mul(half));
//            }

            // zero padding windows
            Window winA1 = new Window(0, l - 1, 0, halfp);
            Window winA2 = new Window(0, l - 1, halfp, p - 1);
            Window winR2 = new Window(0, l - 1, P2 - halfp, P2 - 1);

            // prepare data
            LinearAlgebraUtils.setdata(Res, winA1, A, winA1);
            LinearAlgebraUtils.setdata(Res, winR2, A, winA2);

            // inverse fft per row
            SpectralUtils.ifft(Res, 2);

        } else if (factorcol == 1) {

            // 1d fourier transform per column
            SpectralUtils.fft(A, 1);

            // divide by 2 'cause even fftlength
            A.putRow(halfl, A.getRow(halfl).mmuli(half));
//            for (i=0; i<p; ++i){
//                A(halfl,i) *= half;
//            }

            // zero padding windows
            Window winA1 = new Window(0, halfl, 0, p - 1);
            Window winA2 = new Window(halfl, l - 1, 0, p - 1);
            Window winR2 = new Window(L2 - halfl, L2 - 1, 0, p - 1);

            // prepare data
            LinearAlgebraUtils.setdata(Res, winA1, A, winA1);
            LinearAlgebraUtils.setdata(Res, winR2, A, winA2);

            // inverse fft per row
            SpectralUtils.ifft(Res, 1);

        } else {

            // A=fft2d(A)
            SpectralUtils.fft2d(A);

            // divide by 2 'cause even fftlength
            A.putColumn(halfp, A.getColumn(halfp).mmuli(half));
            A.putRow(halfl, A.getRow(halfl).mmuli(half));
//            for (i=0; i<l; ++i) {
//                A(i,halfp) *= half;
//            }
//            for (i=0; i<p; ++i) {
//                A(halfl,i) *= half;
//            }

            // zero padding windows
            Window winA1 = new Window(0, halfl, 0, halfp);   // zero padding windows
            Window winA2 = new Window(0, halfl, halfp, p - 1);
            Window winA3 = new Window(halfl, l - 1, 0, halfp);
            Window winA4 = new Window(halfl, l - 1, halfp, p - 1);
            Window winR2 = new Window(0, halfl, P2 - halfp, P2 - 1);
            Window winR3 = new Window(L2 - halfl, L2 - 1, 0, halfp);
            Window winR4 = new Window(L2 - halfl, L2 - 1, P2 - halfp, P2 - 1);

            // prepare data
            LinearAlgebraUtils.setdata(Res, winA1, A, winA1);
            LinearAlgebraUtils.setdata(Res, winR2, A, winA2);
            LinearAlgebraUtils.setdata(Res, winR3, A, winA3);
            LinearAlgebraUtils.setdata(Res, winR4, A, winA4);

            // inverse back in 2d
            SpectralUtils.ifft2d(Res);
        }

        // scale
        Res.mmul((double) (factorrow * factorcol));
        return Res;

    }

    public static DoubleMatrix intensity(ComplexDoubleMatrix cint) {
        return MatrixFunctions.pow(cint.real(), 2).add(MatrixFunctions.pow(cint.imag(), 2));
    }

    public static DoubleMatrix magnitude(ComplexDoubleMatrix A) {
        return null;  //To change body of created methods use File | Settings | File Templates.
    }
}

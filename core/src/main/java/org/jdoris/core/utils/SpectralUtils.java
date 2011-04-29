package org.jdoris.core.utils;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_2D;
import org.apache.log4j.Logger;
import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;

public class SpectralUtils {

    static Logger logger = Logger.getLogger(MathUtils.class.getName());

    /**
     * ifft(A,dim)
     * inverse 1dfft over dim of A is returned in A by reference
     * if dim=1 ifft is over all columns of A, if 2 over rows.
     * data is stored major row order in memory, so dim=2 is
     * probably much faster.
     * fftlength should be power of 2
     */
    public static void ifft(ComplexDoubleMatrix A, int dimension) {
        int i;
        // inverse FFT (scaled)
        final int iOpt = -1;

        switch (dimension) {
            case 1: {
                logger.debug("1d ifft over columns");
                final int fftLength = A.rows;
                for (i = 0; i < A.columns; ++i) {
                    ComplexDoubleMatrix VECTOR = A.getColumn(i);
                    four1D(VECTOR, fftLength, iOpt);
                    A.putColumn(i, VECTOR);
                }
                break;
            }
            case 2: {
                logger.debug("1d ifft over rows");
                final int fftLength = A.columns;
                for (i = 0; i < A.rows; ++i) {
                    ComplexDoubleMatrix VECTOR = A.getRow(i);
                    four1D(VECTOR, fftLength, iOpt);
                    A.putRow(i, VECTOR);
                }
                break;
            }
            default:
                logger.error("ifft: dimension != {1,2}");
        }
    } // END ifft

    public static void fft(ComplexDoubleMatrix A, int dimension) {

        int i;
        // forward FFT
        final int intOption = 1;
        int fftLength;
        ComplexDoubleMatrix VECTOR;
        switch (dimension) {
            case 1: {
                fftLength = A.rows;
                for (i = 0; i < A.columns; ++i) {
                    VECTOR = A.getColumn(i);
                    four1D(VECTOR, fftLength, intOption);// but generic.
                    A.putColumn(i, VECTOR);
                    // perhaps can be used directly to do this w/o data copying...
                    // four1(A.getColumn(i), fftlength, iopt);// but generic.
                }
                break;
            }
            case 2: {
                fftLength = A.columns;
                for (i = 0; i < A.rows; ++i) {
                    VECTOR = A.getRow(i);
                    four1D(VECTOR, fftLength, intOption);// but generic.
                    A.putRow(i, VECTOR);
                }
                break;
            }
            default:
                logger.error("fft: dimension != {1,2}");
        }

    }

    /**
     * four1(complr4 *, length, isign)
     * four1(&A[0][0], 128, 1)
     * either based on numerical recipes or veclib
     * helper function for other fft routines, if no veclib
     * cooley-turkey, power 2, replaces input,
     * isign=1: fft , isign=-1: ifft
     * handling vectors should be simpler (lying, standing)
     * note that this is not a good implementation, only to get
     * doris software working without veclib.
     * *
     * define SAMEASVECIB if you want the order of the coefficients
     * of the fft the same as veclib. it seems this is not required
     * for a good version of Doris, but in case of problems this
     * may be the solution.
     * *
     * VECLIB defines the FT same as matlab:
     * N-1                                                 *
     * X(k) = sum  x(n)*exp(-j*2*pi*k*n/N), 0 <= k <= N-1.
     * n=0                                                 *
     * *
     * FFTW defines the same as Matlab, but inv. not normalized.
     * I don't know if the matrix must be allocated somehow, so for
     * now we try only 1d ffts to build 2d too.
     */
    public static void four1D(ComplexDoubleMatrix vector, int fftlength, int direction) {

        DoubleFFT_1D fft = new DoubleFFT_1D(fftlength);

        switch (direction) {
            case 1:
                fft.complexForward(vector.data);

            case -1:
                fft.complexInverse(vector.data, false);
        }

    }

    /**
     * ifftshift(A)                                                 *
     * ifftshift of vector A is returned in A by reference       *
     * undo effect of fftshift. ?p=floor(m/2); A=A[p:m-1 0:p-1]; *
     */
    public static void ifftshift(DoubleMatrix A) throws Exception {

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

    public static void fft2d(ComplexDoubleMatrix A) {
        DoubleFFT_2D fft2d = new DoubleFFT_2D(A.rows, A.columns);
        fft2d.complexForward(A.data);
    }

    // TODO: check declaration of FFT for Real 2D input and realForwardFull
    public static void fft2d(DoubleMatrix A) {
        DoubleFFT_2D fft2d = new DoubleFFT_2D(A.rows, A.columns);
        fft2d.realForwardFull(A.data);
    }

    public static void ifft2d(ComplexDoubleMatrix A) {
        DoubleFFT_2D fft2d = new DoubleFFT_2D(A.rows, A.columns);
        fft2d.complexInverse(A.data, true);
    }
}

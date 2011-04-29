package org.jdoris.core.utils;

import org.jblas.DoubleMatrix;

public class WeightWindows {

    public static DoubleMatrix myrect(DoubleMatrix X) throws Exception {

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
     * myhamming
     * hamming window, lying vector
     * w = (a + (1.-a).*cos((2.*pi/fs).*fr)) .* myrect(fr./Br);
     * scale/shift filter by g(x)=f((x-xo)/s)
     * alpha==1 yields a myrect window
     */
    public static DoubleMatrix myhamming(final DoubleMatrix fr, double RBW, double RSR, double alpha) throws Exception {

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

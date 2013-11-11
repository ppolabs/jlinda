package org.jlinda.core.unwrapping.mcf.utils;

import org.apache.commons.math.util.FastMath;
import org.jblas.DoubleMatrix;

public class JblasUtils {
    public static DoubleMatrix getMatrixFromIdx(DoubleMatrix in, DoubleMatrix idx, int offset) {

        DoubleMatrix out = new DoubleMatrix(idx.rows, idx.columns);
        for (int i = 0; i < idx.length; i++) {
            out.put(i, in.get((int) (FastMath.round(idx.get(i))) - offset));
        }
        return out;
    }

    public static DoubleMatrix getMatrixFromIdx(DoubleMatrix in, DoubleMatrix idx) {
        DoubleMatrix out = new DoubleMatrix(idx.rows, idx.columns);
        for (int i = 0; i < idx.length; i++) {
            out.put(i, in.get((int) FastMath.round(idx.get(i))));
        }
        return out;
    }

    public static DoubleMatrix getMatrixFromRange(int rowMin, int rowMax, int colMin, int colMax, DoubleMatrix in, int offset) {
        DoubleMatrix out = new DoubleMatrix(rowMax - rowMin + 1, colMax - colMin + 1);
        for (int i = rowMin - offset; i < rowMax; i++) {
            for (int j = colMin - offset; j < colMax; j++) {
                out.put(i, j, in.get(i, j));
            }
        }
        return out;
    }

    public static DoubleMatrix getMatrixFromRange(int rowMin, int rowMax, int colMin, int colMax, DoubleMatrix in) {
        DoubleMatrix out = new DoubleMatrix(rowMax - rowMin, colMax - colMin);
        for (int i = rowMin; i < rowMax; i++) {
            for (int j = colMin; j < colMax; j++) {
                out.put(i, j, in.get(i, j));
            }
        }
        return out;
    }

    public static DoubleMatrix setUpMatrixFromIdx(double nRows, double nCols, DoubleMatrix idxRow, DoubleMatrix idxCol) {
        return setUpMatrixFromIdx(nRows, nCols, idxRow, idxCol, 0, 1);
    }

    public static DoubleMatrix setUpMatrixFromIdx(double nRows, double nCols, DoubleMatrix idxRow, DoubleMatrix idxCol, int offset) {
        return setUpMatrixFromIdx(nRows, nCols, idxRow, idxCol, offset, 1);
    }

    public static DoubleMatrix setUpMatrixFromIdx(double nRows, double nCols, DoubleMatrix idxRow, DoubleMatrix idxCol, int offset, int value) {

        if (idxRow.length != idxCol.length) throw new IllegalArgumentException();

        DoubleMatrix out = new DoubleMatrix((int) nRows, (int) nCols);
        for (int i = 0; i < idxRow.length; i++) {
            out.put((int) (FastMath.round(idxRow.get(i))) - offset, (int) (FastMath.round(idxCol.get(i))) - offset, value);
        }
        return out;
    }
}

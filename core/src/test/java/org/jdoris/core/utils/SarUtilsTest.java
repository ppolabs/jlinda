package org.jdoris.core.utils;

import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jdoris.core.io.FlatBinary;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.awt.*;
import java.io.File;
import java.io.FileNotFoundException;

import static org.jblas.MatrixFunctions.cos;
import static org.jblas.MatrixFunctions.sin;

public class SarUtilsTest {

    static int nRows = 100;
    static int nCols = 200;

    static DoubleMatrix phaseMatrix;
    static DoubleMatrix magMatrix;
    static ComplexDoubleMatrix cplxData;
    static ComplexDoubleMatrix temp2;
    static ComplexDoubleMatrix temp1;
    private static final double DELTA = 1e-08;


    @Before
    public void setUpMultilookTestData() {
        // TODO: refactor to simulation package
        // multilook an array with a noisy phase trend (10 fringes),
        int numFringes = 10;
        phaseMatrix = MathUtils.ramp(nRows, nCols).muli(2 * Math.PI).muli(numFringes);//.add(DoubleMatrix.randn(nRows, nCols).mmul(0.25 * Math.PI));
        magMatrix = DoubleMatrix.ones(nRows, nCols);//.add(DoubleMatrix.randn(nRows, nCols).mmul(0.5));

        temp1 = new ComplexDoubleMatrix(magMatrix);//, DoubleMatrix.zeros(nRows, nCols));
        temp2 = new ComplexDoubleMatrix(cos(phaseMatrix), sin(phaseMatrix));
        cplxData = temp1.mul(temp2);

    }

    @Test
    public void testMultilook() throws Exception {
        ComplexDoubleMatrix cpxData_ml_ACTUAL = SarUtils.multilook(cplxData, 5, 5);
        DoubleMatrix cpxData_ml_real_EXPECTED = setData("/d2/test_real.out", cpxData_ml_ACTUAL.rows, cpxData_ml_ACTUAL.columns);
        DoubleMatrix cpxData_ml_imag_EXPECTED = setData("/d2/test_imag.out", cpxData_ml_ACTUAL.rows, cpxData_ml_ACTUAL.columns);
        ComplexDoubleMatrix cpxData_ml_EXPECTED = new ComplexDoubleMatrix(cpxData_ml_real_EXPECTED, cpxData_ml_imag_EXPECTED);

        Assert.assertArrayEquals(cpxData_ml_ACTUAL.real().toArray(), cpxData_ml_ACTUAL.real().toArray(), DELTA);
        Assert.assertArrayEquals(cpxData_ml_ACTUAL.imag().toArray(), cpxData_ml_ACTUAL.imag().toArray(), DELTA);


    }

    private void dumpData(double[][] data, String fileName, int rows, int columns) throws FileNotFoundException {
        FlatBinary outRealFile = new FlatBinary();
        outRealFile.setFile(new File(fileName));
        outRealFile.setDimensions(new Rectangle(columns, rows));
        outRealFile.setOutStream();
        outRealFile.setData(data);
        outRealFile.writeDoubleToStream();
    }

    private DoubleMatrix setData(String fileName, int rows, int columns) throws FileNotFoundException {
        FlatBinary outRealFile = new FlatBinary();
        outRealFile.setFile(new File(fileName));
        outRealFile.setDimensions(new Rectangle(columns, rows));
        outRealFile.setInStream();
        outRealFile.readDoubleFromStream();
        return new DoubleMatrix(outRealFile.data);
    }

//    @Test
//    public void testIntensity() throws Exception {
//
//    }
//
//    @Test
//    public void testMagnitude() throws Exception {
//
//    }
//
//    @Test
//    public void testCoherence() throws Exception {
//
//    }
//
//    @Test
//    public void testMultilook() throws Exception {
//
//    }

}

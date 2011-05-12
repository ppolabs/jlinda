package org.jdoris.core.io;

import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import java.awt.*;
import java.io.File;
import java.nio.ByteOrder;

public class FlatBinaryDoubleTest {

    private static FlatBinaryDouble flatBinaryDoubleRead;
    private static FlatBinaryDouble flatBinaryDoubleWrite;
    private static Rectangle testDataRectangle;

    private static double[][] testData;

    private static FlatBinaryDouble flatBinaryDoubleLittleRead;
    private static FlatBinaryDouble flatBinaryDoubleLittleWrite;

    @BeforeClass
    public static void setupTestData() throws Exception {

        testDataRectangle = new Rectangle(123, 321);
        testData = new double[testDataRectangle.width][testDataRectangle.height];

        for (int i = 0; i < testDataRectangle.width; i++) {
            for (int j = 0; j < testDataRectangle.height; j++) {
                testData[i][j] = Math.random()*100;
            }
        }

        flatBinaryDoubleRead = new FlatBinaryDouble();
        flatBinaryDoubleWrite = new FlatBinaryDouble();

        flatBinaryDoubleRead.setFile(new File("test/test.in"));
        flatBinaryDoubleWrite.setFile(new File("test/test.out"));

        flatBinaryDoubleWrite.setOutStream();

        flatBinaryDoubleLittleRead = new FlatBinaryDouble();
        flatBinaryDoubleLittleRead.setByteOrder(ByteOrder.LITTLE_ENDIAN);

        flatBinaryDoubleLittleWrite = new FlatBinaryDouble();
        flatBinaryDoubleLittleWrite.setByteOrder(ByteOrder.LITTLE_ENDIAN);

        flatBinaryDoubleLittleRead.setFile(new File("test/test.in.swapped"));
        flatBinaryDoubleLittleWrite.setFile(new File("test/test.out.swapped"));

        flatBinaryDoubleLittleWrite.setOutStream();
    }

    @AfterClass
    public static void cleanTestData() {

//        if (!testFile.exists())
//            throw new IllegalArgumentException("Delete: no such file or directory: " + testFile.getName());
//
//        boolean success = testFile.delete();
//
//        if (!success)
//            throw new IllegalArgumentException("Delete: deletion of file" + testFile.getName() + " failed");
//
//        System.gc();

    }

    @Test
    public void testCreateAndCheck() throws Exception {
        flatBinaryDoubleWrite.create();
        Assert.assertEquals("File creating: ", true, flatBinaryDoubleWrite.checkExists());
    }

    @Test
    public void testWritingReadingData() throws Exception {

        flatBinaryDoubleWrite.setDimensions(new Rectangle(testDataRectangle));
        flatBinaryDoubleWrite.setData(testData);
        flatBinaryDoubleWrite.writeToStream();

        flatBinaryDoubleRead.setFile(flatBinaryDoubleWrite.file);
        flatBinaryDoubleRead.setDimensions(new Rectangle(testDataRectangle));
        flatBinaryDoubleRead.setInStream();
        flatBinaryDoubleRead.readFromStream();

        Assert.assertArrayEquals(testData, flatBinaryDoubleRead.data);

    }

    @Test
    public void testWritingReadingLittleEndianData() throws Exception {

        flatBinaryDoubleLittleWrite.setDimensions(new Rectangle(testDataRectangle));
        flatBinaryDoubleLittleWrite.setData(testData);
        flatBinaryDoubleLittleWrite.writeToStream();

        flatBinaryDoubleLittleRead.setFile(flatBinaryDoubleLittleWrite.file);
        flatBinaryDoubleLittleRead.setDimensions(new Rectangle(testDataRectangle));
        flatBinaryDoubleLittleRead.setInStream();
        flatBinaryDoubleLittleRead.readFromStream();

        Assert.assertArrayEquals(testData, flatBinaryDoubleLittleRead.data);

    }

/*
    @Ignore
    @Test
    public void testCanRead() throws Exception {
    }

    @Ignore
    @Test
    public void testCanWrite() throws Exception {
    }
*/
}

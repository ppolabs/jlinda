package org.jdoris.core.io;

import org.junit.*;

import java.awt.*;
import java.io.File;

public class FlatBinaryTest {

    private static FlatBinary flatBinaryRead;
    private static FlatBinary flatBinaryWrite;
    private static Rectangle testDataRectangle;

    private static double[][] testData;

    @BeforeClass
    public static void setupTestData() throws Exception {

        testDataRectangle = new Rectangle(123, 321);
        testData = new double[testDataRectangle.height][testDataRectangle.width];

        for (int i = 0; i < testDataRectangle.height; i++) {
            for (int j = 0; j < testDataRectangle.width; j++) {
                testData[i][j] = Math.random();
            }
        }

        flatBinaryRead = new FlatBinary();
        flatBinaryWrite = new FlatBinary();

        flatBinaryRead.setFile(new File("test.in"));
        flatBinaryWrite.setFile(new File("test.out"));

        flatBinaryWrite.setOutStream();

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
        flatBinaryWrite.create();
        Assert.assertEquals("File creating: ", true, flatBinaryWrite.checkExists());
    }

    // TODO: brake up this test in smaller test blocks
    @Test
    public void testWritingReadingData() throws Exception {

        flatBinaryWrite.setDimensions(new Rectangle(testDataRectangle));
        flatBinaryWrite.setData(testData);
        flatBinaryWrite.writeDoubleToStream();

        flatBinaryRead.setFile(flatBinaryWrite.file);
        flatBinaryRead.setDimensions(new Rectangle(testDataRectangle));
        flatBinaryRead.setInStream();
        flatBinaryRead.readDoubleFromStream();

        Assert.assertArrayEquals(flatBinaryRead.data, testData);
    }

    @Ignore
    @Test
    public void testCanRead() throws Exception {
    }

    @Ignore
    @Test
    public void testCanWrite() throws Exception {
    }
}

package org.jdoris.core.io;

import org.jdoris.core.Window;
import org.junit.*;

import java.io.File;

public class FlatBinaryTest {

    private static FlatBinary flatBinaryRead;
    private static FlatBinary flatBinaryWrite;

    private static double[][] testData;
    private static Window testDataWindow;

    @BeforeClass
    public static void setupTestData() throws Exception {

        testDataWindow = new Window(0, 123, 0, 321);
        int lines = (int) testDataWindow.lines() - 1;
        int pixels = (int) testDataWindow.pixels() - 1;

        testData = new double[lines][pixels];

        flatBinaryRead = new FlatBinary();
        flatBinaryWrite = new FlatBinary();

//        flatBinaryRead.setFile(new File("test.in"));
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

        flatBinaryWrite.setDataWindow(new Window(testDataWindow));
        flatBinaryWrite.setData(testData);
        flatBinaryWrite.writeDoubleToStream();

        flatBinaryRead.setFile(flatBinaryWrite.file);
        flatBinaryRead.setDataWindow(new Window(testDataWindow));
        flatBinaryRead.setInStream();
        flatBinaryRead.readDoubleFromStream();

        Assert.assertArrayEquals(testData, flatBinaryRead.data);
    }

//    @Ignore
//    @Test
//    public void testCanRead() throws Exception {
//    }
//
//    @Ignore
//    @Test
//    public void testCanWrite() throws Exception {
//    }
}

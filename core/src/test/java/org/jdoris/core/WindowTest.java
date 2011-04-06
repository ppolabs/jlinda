package org.jdoris.core;

import org.junit.Test;
import static junit.framework.Assert.assertEquals;

/**
 * User: pmar@ppolabs.com
 * Date: 2/18/11
 * Time: 7:00 PM
 */
public class WindowTest {

    Window refWin = new Window(11, 21, 103, 114);

    @Test
    public void testSetWindow() throws Exception {

        Window testWin = new Window();
        testWin.setWindow(refWin);

        assertEquals(refWin.linelo, testWin.linelo);
        assertEquals(refWin.linehi, testWin.linehi);
        assertEquals(refWin.pixlo, testWin.pixlo);
        assertEquals(refWin.pixhi, testWin.pixhi);

    }

    @Test
    public void testCompareTo() throws Exception {
        Window testWin = new Window();
        testWin.setWindow(refWin);
        assertEquals(0, testWin.compareTo(refWin));
    }

    @Test
    public void testClone() throws Exception {
        Window testWin = (Window) refWin.clone();
        assertEquals(refWin,testWin);
    }

    @Test
    public void testLines() throws Exception {
        Window testWin = (Window) refWin.clone();
        assertEquals(refWin.lines(), testWin.lines());
    }

    @Test
    public void testPixels() throws Exception {
        Window testWin = (Window) refWin.clone();
        assertEquals(refWin.pixels(), testWin.pixels());
    }

}

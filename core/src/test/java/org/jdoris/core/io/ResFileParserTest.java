package org.jdoris.core.io;


import org.junit.After;
import org.junit.Before;
import org.junit.Test;

/**
 * User: pmar@ppolabs.com
 * Date: 4/1/11
 * Time: 5:38 PM
 */
public class ResFileParserTest {

    // template for testing
    String fileName = "/d2/test.res";

    @Before
    public void setUp() throws Exception {

    }

    @After
    public void tearDown() throws Exception {

    }

    @Test
    public void testBufferReading() throws Exception {

        String key = "key3";
//
//        ResFileParser.initializeLogger();
        StringBuffer bf = ResFileParser.createBuffer(fileName);

        System.out.println(bf.toString());




        System.out.println(ResFileParser.indexKey(bf, key, "start"));
        System.out.println(ResFileParser.indexKey(bf, key, "end"));

        StringBuffer tt;
        tt = new StringBuffer(bf.substring(34, 66));

        System.out.println(tt);

        System.out.println("---------");

        System.out.println(ResFileParser.createSubBuffer(bf, key, key));



//        ResFileParser.setLogging();

//        ArrayList values = ResFileParser.queryAsciiBuffer(bf, key, 3);
//
//        for (Object value : values) {
//            System.out.println("i = " + value);
//        }
//
////          DateTime dateAndTime = new DateTime();
//
//        ProductData.UTC dateAndTime = ProductData.UTC.parse(values.get(0).toString().trim());
//

//        ResFileParser.parseOrbit(bf);



    }

}
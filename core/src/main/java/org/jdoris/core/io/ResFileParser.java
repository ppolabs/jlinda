package org.jdoris.core.io;

import org.apache.log4j.*;
import org.esa.beam.framework.datamodel.ProductData;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.InputMismatchException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * User: pmar@ppolabs.com
 * Date: 4/1/11
 * Time: 5:01 PM
 */
public class ResFileParser {

    static Logger logger = Logger.getLogger(ResFileParser.class.getName());

    public enum IndexPositions {
        START, END
    }

//    public static void initializeLogger() {
//        BasicConfigurator.configure();
//        Logger.getRootLogger().setLevel(Level.ALL);
//        Layout layout = new PatternLayout("%d [%t] %-5p %c %x - %m%n");
//        Logger.getRootLogger().addAppender(new ConsoleAppender(layout));
//    }


    //        logger.setLevel(Level.INFO);
//    public static void setLogging() {
//
//        logger.info("test message!");
//
//    }

    // method to buffer doris res file
    public static StringBuffer createBuffer(String asciiFileName) {

        StringBuffer resFileContent = new StringBuffer();
        BufferedReader input = null;

        try {
            input = new BufferedReader(new FileReader(asciiFileName), 1);
            String line = null;
            while ((line = input.readLine()) != null) {
                resFileContent.append(line);
                resFileContent.append(System.getProperty("line.separator"));
            }
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        }

        return resFileContent;
    }

    public static StringBuffer createSubBuffer(StringBuffer inputBuffer, int idxStart, int idxEnd) {
        return new StringBuffer(inputBuffer.substring(idxStart, idxEnd));
    }

    public static StringBuffer createSubBuffer(StringBuffer inputBuffer, String key1, String key2) throws InputMismatchException {

        int idxStart = indexEndKey(inputBuffer, key1);
        int idxEnd = indexStartKey(inputBuffer, key2);

        logger.debug("Start key: " + idxStart);
        logger.debug("End key: " + idxEnd);
        return new StringBuffer(inputBuffer.substring(idxStart, idxEnd));

    }

    // define pattern: line starts with key, value(s) separated by ":"
    public static String createPattern(String key) {
        return "\\s*?(" + key + "):(\\s)(.*)";
    }

    private static Matcher createMatcher(StringBuffer buffer, String key) {
        String pattern = createPattern(key);
        Pattern dataPattern = Pattern.compile(pattern);
        return dataPattern.matcher(buffer);
    }

    // return index of first match of key
    public static int indexKey(StringBuffer inputBuffer, String key, String position) {

        int returnIndex = 0;
        Matcher match = createMatcher(inputBuffer, key);

        switch (IndexPositions.valueOf(position.toUpperCase())) {

            case START:
                if (match.find()) {
                    returnIndex = match.start();
                }
                break;

            case END:

                if (match.find()) {
                    returnIndex = match.end();
                }
                break;

            default:
                returnIndex = 0;
        }

        return returnIndex;
    }

    // return end index of first match of key
    public static int indexEndKey(StringBuffer inputBuffer, String key) {
        return indexKey(inputBuffer, key, "end");
    }

    // return end index of first match of key
    public static int indexStartKey(StringBuffer inputBuffer, String key) {
        return indexKey(inputBuffer, key, "start");
    }

    // method to query for keys in ascii file
    public static ArrayList queryKey(StringBuffer inputBuffer, String key, int groupToReturn) {

        ArrayList<String> valuesList = new ArrayList<String>();

        Matcher match = createMatcher(inputBuffer, key);
        while (match.find()) {
            valuesList.add(match.group(groupToReturn));

            // logging
            System.out.println("match.end() = " + match.end());
        }

        return valuesList;

    }

    // method to qyery for keys in acii file : returns group (3)
    public static ArrayList queryKey(StringBuffer inputBuffer, String key) {
        return queryKey(inputBuffer, key, 3);
    }

    public static ProductData.UTC parseDateTime(String dateTimeString) throws ParseException {
        return ProductData.UTC.parse(dateTimeString);
    }


//    public static void parseOrbit(StringBuffer inputBuffer) throws ParseException {
//
////        NUMBER_OF_DATAPOINTS:   5
////        22372.090719    2983988.510000     5512738.010000    3461111.080000
////        22376.158551    2997034.770000     5522094.860000    3434910.600000
////        22380.226384    3010032.990000     5531345.290000    3408648.260000
////        22384.294216    3022982.860000     5540489.160000    3382324.550000
////        22388.362049    3035884.100000     5549526.330000    3355939.930000
//
//
//        final String orbitStartKey = "NUMBER_OF_DATAPOINTS";
//
//        // get number of state vectors
//        // define pattNumStateVectors: line starts with key, value(s) separated by ":"
//        String numStateVectors = "\\s*?(" + orbitStartKey + "):(\\s)(.*)";
//        Pattern pattStateVectors = Pattern.compile(numStateVectors);
//        Matcher match = pattStateVectors.matcher(inputBuffer);
//
//
//
//        int numberOfStateVectors;
//        int idxOrbitBlockStart = 0;
//
//        if (match.find()) {
//            numberOfStateVectors = Integer.parseInt(match.group(3).trim());
//
//            // index of start of StateVectors block
//            match = dataPattern.matcher(inputBuffer);
//            idxOrbitBlockStart = match.end();
//        }
//
//        // get index of StateVectors block
////        pattNumStateVectors = "^(\\d+\\.\\d+)\\s*(\\d+\\.\\d+)\\s*(\\d+\\.\\d+)\\s*(\\d+\\.\\d+)";
////        System.out.println(inputBuffer.capacity());
//
//        String stateVectors = "\\s*?(\\d+\\.\\d+)\\s*?(\\d+\\.\\d+)\\s*?(\\d+\\.\\d+)\\s*?(\\d+\\.\\d+)";
//        pattStateVectors = Pattern.compile(stateVectors);
//
//        int stateVectorCounter = 1;
//        while (match.find() && stateVectorCounter <= 5) {
//
//            System.out.print(match.group(0));
//            stateVectorCounter++;
//        }

////        System.out.println("match.group(1) = " + match.group(1));

//        // numberOfElements 4*numberOfStateVectors = 20
////        // get
////        // define pattNumStateVectors: line starts with key, value(s) separated by ":"
////        String pattNumStateVectors = "\\s*?(" + orbitStartKey + "):(\\s)(.*)";
////        Pattern dataPattern = Pattern.compile(pattNumStateVectors);
////        Matcher match = dataPattern.matcher(inputBuffer);
//
//    }

}

package org.jdoris.core.io;

import org.apache.log4j.Logger;
import org.esa.beam.framework.datamodel.ProductData;

import java.io.*;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.InputMismatchException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ResFile {

    private static Logger logger = Logger.getLogger(ResFile.class.getName());

    /*
        public static void initializeLogger() {
            BasicConfigurator.configure();
            Logger.getRootLogger().setLevel(Level.ALL);
            Layout layout = new PatternLayout("%d [%t] %-5p %c %x - %m%n");
            Logger.getRootLogger().addAppender(new ConsoleAppender(layout));
        }
    */

    // fields
    private File resFile;

    // start/end index for subBuffers
    private int startIdx = 0;
    private int endIdx = 0;

    private StringBuffer buffer = new StringBuffer();

    public ResFile() {
    }

    private static enum IndexPositions {
        START, END;
    }

    public ResFile(File file) {

        resFile = file;

        try {
            BufferedReader input = new BufferedReader(new FileReader(file), 1);
            String line;
            while ((line = input.readLine()) != null) {
                buffer.append(line);
                buffer.append(System.getProperty("line.separator"));
            }
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        }

        startIdx = 0;
        endIdx = buffer.length();

    }

    // method to buffer doris res file
    public void streamBuffer(File fileName) {

        BufferedReader input;

        try {
            input = new BufferedReader(new FileReader(fileName), 1);
            String line;
            while ((line = input.readLine()) != null) {
                buffer.append(line);
                buffer.append(System.getProperty("line.separator"));
            }
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        }

    }

    public void dumpBuffer() throws IOException {

        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(resFile));
            out.write(buffer.toString());
            out.close();
        } catch (IOException e) {
            logger.error(" dumpBuffer() exception " + e.getMessage());
            e.printStackTrace();
        }

    }

    public void dumpBuffer(File newFile) throws IOException {

        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(newFile));
            out.write(buffer.toString());
            out.close();
        } catch (IOException e) {
            logger.error(" dumpBuffer() exception " + e.getMessage());
            e.printStackTrace();
        }

    }


    public void setSubBuffer(int start, int end) {
        if (start > end) {
            throw new InputMismatchException();
        } else {
            startIdx = start;
            endIdx = end;
        }
    }

    public void setSubBuffer(String key1, String key2) throws InputMismatchException {

        int idxStart = indexStartKey(key1);
        int idxEnd = indexEndKey(key2);

        if (idxStart > idxEnd) {
            throw new InputMismatchException();
        }

        setSubBuffer(idxStart, idxEnd);
    }

    public void resetSubBuffer() {
        startIdx = 0;
        endIdx = buffer.length();
    }

    // define pattern: line starts with key, value(s) separated by ":"
    private static String createPattern(String key) {
//        return "\\s*?(" + key + "):(\\s)(.*)";
        return "\\s*(" + key + "):?(.*)";
    }

    private static Matcher createMatcher(String buffer, String pattern) {
        Pattern dataPattern = Pattern.compile(pattern);
        return dataPattern.matcher(buffer);
    }

    // return index of first match of key
    private int indexKey(String key, String position) {

        int returnIndex = 0;
        Matcher match = createMatcher(buffer.substring(startIdx, endIdx), createPattern(key));

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
    private int indexEndKey(String key) {
        return indexKey(key, "end");
    }

    // return end index of first match of key
    private int indexStartKey(String key) {
        return indexKey(key, "start");
    }

    // method to query for keys in ascii file
    public ArrayList queryKey(String key, int groupToReturn) throws IndexOutOfBoundsException {

        ArrayList<String> valuesList = new ArrayList<String>();

        Matcher match = createMatcher(buffer.substring(startIdx, endIdx), createPattern(key));
        while (match.find()) {
            try {
                valuesList.add(match.group(groupToReturn));
            } catch (IndexOutOfBoundsException e) {
                logger.error("queryKey(key,group) : Exception handling regex : " + e.getLocalizedMessage());
            }
        }
        return valuesList;
    }

    // method to query for keys in acii file : returns group (2)
    public ArrayList queryKey(String key) throws IndexOutOfBoundsException {
        return queryKey(key, 2);
    }

    public String parseStringValue(String key) {
        return queryKey(key, 2).get(0).toString().trim();
    }

    public int parseIntegerValue(String key) {
        return Integer.parseInt(parseStringValue(key));
    }

    public double parseDoubleValue(String key) {
        return Double.parseDouble(parseStringValue(key));
    }

    public ProductData.UTC parseTimeValue(String key) throws ParseException {
        return ProductData.UTC.parse(parseStringValue(key));
    }

    public double[][] parseOrbit() throws Exception {

        // get number of state vectors
        final String numStateVectorsKey = "NUMBER_OF_DATAPOINTS";
        final int numberOfStateVectors = parseIntegerValue(numStateVectorsKey);

        final double[][] stateVectors = new double[numberOfStateVectors][4];

        String stateVectorsPattern = "\\s*?(\\d+\\.\\d+)\\s*?(\\d+\\.\\d+)\\s*?(\\d+\\.\\d+)\\s*?(\\d+\\.\\d+)";
        Matcher match = createMatcher(buffer.substring(startIdx, endIdx), stateVectorsPattern);

        int i = 0;
        while (match.find()) {
            try {

                stateVectors[i][0] = Double.parseDouble(match.group(1).trim());
                stateVectors[i][1] = Double.parseDouble(match.group(2).trim());
                stateVectors[i][2] = Double.parseDouble(match.group(3).trim());
                stateVectors[i][3] = Double.parseDouble(match.group(4).trim());

                i++;

            } catch (IndexOutOfBoundsException e) {
                logger.error("parseOrbit() : Exception handling regex : " + e.getLocalizedMessage());
            }
        }

        if (i != numberOfStateVectors) {
            logger.error("parseOrbit() : inconsistency in number of defined and parsed state vectors");
            throw new IllegalArgumentException("Cannot parse orbit : number of defined and parsed state vectors not the same");
        }

        return stateVectors;
    }

    // Getters & Setters
    // --------------------------------------------------------------------------
    public File getResFile() {
        return resFile;
    }

    public void setResFile(File resfile) {
        resFile = resfile;
    }

    public int getStartIdx() {
        return startIdx;
    }

    public void setStartIdx(int startidx) {
        startIdx = startidx;
    }

    public int getEndIdx() {
        return endIdx;
    }

    public void setEndIdx(int endidx) {
        endIdx = endidx;
    }

    public StringBuffer getBuffer() {
        return buffer;
    }

    public void setBuffer(StringBuffer inbuffer) {
        buffer = inbuffer;
    }

    @Override
    protected Object clone() throws CloneNotSupportedException {
        return super.clone();
    }


    /**
     * Computes a hash code for a double value, using the algorithm from
     * Joshua Bloch's book <i>Effective Java"</i>
     *
     * @param x
     * @return a hashcode for the double value
     */
    public static int hashCode(double x) {
        long f = Double.doubleToLongBits(x);
        return (int) (f ^ (f >>> 32));
    }

    /**
     * Gets a hashcode for this window.
     *
     * @return a hashcode for this window
     */
    public int hashCode() {
        //Algorithm from Effective Java by Joshua Bloch
        int result = 17;
        result = 37 * result + hashCode(startIdx);
        result = 37 * result + hashCode(endIdx);
        result = 37 * result + hashCode(buffer.length());
        return result;
    }

}

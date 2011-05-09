package org.jdoris.core.io;

import org.apache.log4j.Logger;

import java.awt.*;
import java.io.*;

public class FlatBinary {

    private static Logger logger = Logger.getLogger(FlatBinary.class.getName());

    File file;
    String format;
    long sizeBytes;
    Rectangle dimensions;
    DataOutputStream outStream;
    DataInputStream inStream;

    public double[][] data;

    public FlatBinary() {

    }

    public boolean checkExists() {
        return file.exists();
    }

    public boolean checkCanRead() throws FileNotFoundException {
        return file.canRead();
    }

    public boolean checkCanWrite() throws FileNotFoundException {
        return file.canWrite();
    }

    public void setInStream() throws FileNotFoundException {
        inStream = new DataInputStream(new BufferedInputStream(new FileInputStream(file.getAbsoluteFile())));
    }

    public void setOutStream() throws FileNotFoundException {
        outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file.getAbsoluteFile())));
    }

    public void readDoubleFromStream() throws FileNotFoundException {
        data = new double[dimensions.height][dimensions.width];
        for (int i = 0; i < dimensions.height; i++) {
            for (int j = 0; j < dimensions.width; j++) {
                try {
                    data[i][j] = inStream.readDouble();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

    }

    public void writeDoubleToStream() throws FileNotFoundException {
        for (int i = 0; i < dimensions.height; i++) {
            for (int j = 0; j < dimensions.width; j++) {
                try {
//                    outStream.writeLong(Double.doubleToLongBits(data[i][j]));
                    outStream.writeDouble(data[i][j]);
//                    System.out.println("data[" + i + "]" + "[" + j + "]:  " + data[i][j]);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

        try {
            this.outStream.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void setData(double[][] data) {
        this.data = data;
    }

    public void create() {
        try {
            file.createNewFile();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void create(File genericFile) {
        try {
            genericFile.createNewFile();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void create(String genericFileName) {
        this.create(new File(genericFileName));
    }

    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        sb.append("FlatBinary");
        sb.append("{file=").append(file.getAbsoluteFile());
        sb.append(", format='").append(format).append('\'');
        sb.append(", sizeBytes=").append(sizeBytes);
        sb.append(", dimensions=").append(dimensions);
        sb.append('}');
        return sb.toString();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        FlatBinary that = (FlatBinary) o;

        if (sizeBytes != that.sizeBytes) return false;
        if (dimensions != null ? !dimensions.equals(that.dimensions) : that.dimensions != null) return false;
        if (file != null ? !file.equals(that.file) : that.file != null) return false;
        if (format != null ? !format.equals(that.format) : that.format != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = file != null ? file.hashCode() : 0;
        result = 31 * result + (format != null ? format.hashCode() : 0);
        result = 31 * result + (int) (sizeBytes ^ (sizeBytes >>> 32));
        result = 31 * result + (dimensions != null ? dimensions.hashCode() : 0);
        return result;
    }

    public void setFile(File file) {
        this.file = file;
    }

    public void setFileName(String genericFileName) {
        this.file = new File(genericFileName);
    }

    public void setFormat(String format) {
        this.format = format;
    }

    public void setSizeBytes(long sizeBytes) {
        this.sizeBytes = sizeBytes;
    }

    public void setDimensions(Rectangle dimensions) {
        this.dimensions = dimensions;
    }

    public File getFile() {
        return file;
    }

    public String getFormat() {
        return format;
    }

    public long getSizeBytes() {
        return sizeBytes;
    }

    public Rectangle getDimensions() {
        return dimensions;
    }

}

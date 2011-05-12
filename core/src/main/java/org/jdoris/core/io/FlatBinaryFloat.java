package org.jdoris.core.io;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteOrder;

public final class FlatBinaryFloat extends FlatBinary {

    float [][] data;

    public FlatBinaryFloat() {
        this.byteOrder = ByteOrder.BIG_ENDIAN;
    }

    @Override
    public void readFromStream() throws FileNotFoundException {
        data = new float[dimensions.width][dimensions.height];
        for (int i = 0; i < dimensions.width; i++) {
            for (int j = 0; j < dimensions.height; j++) {
                try {
                    if (byteOrder == ByteOrder.LITTLE_ENDIAN) {
                        data[i][j] = ByteSwapper.swap(inStream.readFloat());
                    } else {
                        data[i][j] = inStream.readFloat();
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }

            }
        }
    }

    @Override
    public void writeToStream() throws FileNotFoundException {
        for (int i = 0; i < dimensions.width; i++) {
            for (int j = 0; j < dimensions.height; j++) {
                try {
                    if (byteOrder == ByteOrder.LITTLE_ENDIAN) {
                        outStream.writeFloat(ByteSwapper.swap(data[i][j]));
                    } else {
                        outStream.writeFloat(data[i][j]);
                    }
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
}

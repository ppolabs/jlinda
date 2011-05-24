package org.jdoris.core.filtering;

import org.jblas.ComplexDoubleMatrix;
import org.junit.Assert;
import org.junit.Test;

import static org.jdoris.core.io.DataReader.readCplxFloatData;

public class RangeFilterTest {

    private static final double DELTA_04 = 1e-04;

    // TODO: more robust tests -- for different options and params
    @Test
    public void testRangeFilterBlock() throws Exception {

        /// define parameters parameters
        final int nlMean = 15;
        final int snRthreshold = 5;
        final double RSR = 18962500.774137583;
        final int RBW = 15550000;
        final double alphaHamming = 0.75;
        final int ovsFactor = 2;
        final boolean doWeightCorrelFlag = true;

        /// load Input Data
        String fileMasterDataName = "test/testdata_cplx_MASTER_128_128.cr4.swap";
        ComplexDoubleMatrix masterCplx = readCplxFloatData(fileMasterDataName, 128, 128);

        String fileSlaveDataName = "test/testdata_cplx_SLAVE_128_128.cr4.swap";
        ComplexDoubleMatrix slaveCplx = readCplxFloatData(fileSlaveDataName, 128, 128);

        /// load Expected Data
        String fileMasterDataNameFiltered = "test/testdata_cplx_MASTER_RNGFILT_128_128.cr4.swap";
        ComplexDoubleMatrix masterCplx_rngFilter_EXPECTED = readCplxFloatData(fileMasterDataNameFiltered, 128, 128);

        String fileSlaveDataNameFiltered = "test/testdata_cplx_SLAVE_RNGFILT_128_128.cr4.swap";
        ComplexDoubleMatrix slaveCplx_rngFilter_EXPECTED = readCplxFloatData(fileSlaveDataNameFiltered, 128, 128);

        /// range filter data block
        RangeFilter.filterBlock(masterCplx, slaveCplx, nlMean, snRthreshold, RSR, RBW, alphaHamming, ovsFactor, doWeightCorrelFlag);

        /// assert
        Assert.assertArrayEquals(masterCplx_rngFilter_EXPECTED.toDoubleArray(), masterCplx.toDoubleArray(), DELTA_04);
        Assert.assertArrayEquals(slaveCplx_rngFilter_EXPECTED.toDoubleArray(), slaveCplx.toDoubleArray(), DELTA_04);

    }
}

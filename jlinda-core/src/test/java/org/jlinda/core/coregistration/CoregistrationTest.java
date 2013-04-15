package org.jlinda.core.coregistration;


import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.Logger;
import org.jblas.*;
import org.jlinda.core.SLCImage;
import org.jlinda.core.Window;
import org.jlinda.core.io.DataReader;
import org.jlinda.core.utils.MathUtils;
import org.jlinda.core.utils.PolyUtils;
import org.jlinda.core.utils.SarUtils;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;
import org.perf4j.StopWatch;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.nio.ByteOrder;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static org.jblas.MatrixFunctions.pow;
import static org.jblas.MatrixFunctions.sqrt;

/**
 * User: pmar@ppolabs.com
 * Date: 12/13/12
 * Time: 2:06 PM
 */
public class CoregistrationTest {

    private String dataPath;
    private String processingPath;
    private String masterFileName;
    private String slaveFileName;
    private String masterMagnitudeFileName;
    private String slaveMagnitudeFileName;
    private String correlFileName;

    private ComplexDoubleMatrix masterCplx;
    private ComplexDoubleMatrix slaveCplx;
    private FloatMatrix masterMagnitude;
    private FloatMatrix slaveMagnitude;
    private FloatMatrix correlMasterSlave;

    private ByteOrder littleEndian = ByteOrder.LITTLE_ENDIAN;
    private int rows;
    private int rowsAcc;
    private int cols;
    private int colsAcc;

    private static final Logger logger = (Logger) LoggerFactory.getLogger(CoregistrationTest.class);

    @Ignore
//    @Before
    public void setUpMagSpace() throws Exception {

        // CORRELATE
        // -------------------------------------------------------------------------------
        // declare file names
        processingPath = "/d2/test.processing/unit_tests/etna.volcano/process/crop/01486_21159.cpm/";
        dataPath = "coarseCorr/magSpace/";
//        path = "/d2/test.processing/unit_tests/etna.volcano/process/crop/01486_21159.cpm/";
        masterFileName = "master.cr4";
        slaveFileName = "slave.cr4";
        masterMagnitudeFileName = "masterMagnitude.r4";
        slaveMagnitudeFileName = "slaveMagnitude.r4";
        correlFileName = "correlMasterSlave.r4";

        // dimensions
        rows = 256 + 1;
        rowsAcc = 8;
        cols = 256 + 1;
        colsAcc = 8;

        // load all data
        int colsMaster = cols + 2 * rowsAcc;
        int rowsMaster = rows + 2 * colsAcc;
        int colsSlave = cols;
        int rowsSlave = rows;


        logger.setLevel(Level.TRACE);

        StopWatch clockReadingData = new StopWatch();
        clockReadingData.start();

        masterCplx = DataReader.readCplxFloatData(processingPath + dataPath + masterFileName, rowsMaster, colsMaster, littleEndian);
        slaveCplx = DataReader.readCplxFloatData(processingPath + dataPath + slaveFileName, rowsSlave, colsSlave, littleEndian);
        masterMagnitude = DataReader.readFloatData(processingPath + dataPath + masterMagnitudeFileName, rowsMaster, colsMaster, littleEndian);
        slaveMagnitude = DataReader.readFloatData(processingPath + dataPath + slaveMagnitudeFileName, rowsSlave, colsSlave, littleEndian);

        correlMasterSlave = DataReader.readFloatData(processingPath + dataPath + correlFileName, rowsMaster, colsMaster, littleEndian);

        clockReadingData.stop();
        logger.info("Time to read data  - new Way: {}", clockReadingData.getElapsedTime());


    }

    @Ignore
    @Test
    public void testMagSpace() throws Exception {

        // dump values to check whether everything is parsed right
        System.out.println(masterCplx.get(0, 0));
        System.out.println(slaveCplx.get(0, 0));
        System.out.println(masterMagnitude.get(0, 0));
        System.out.println(SarUtils.magnitude(masterCplx).get(0, 0));
        System.out.println(slaveMagnitude.get(0, 0));
        System.out.println(correlMasterSlave.get(0, 0));

        ComplexFloatMatrix masterCplx_FLOAT = castToComplexFloatMatrix(masterCplx);
        ComplexFloatMatrix slaveCplx_FLOAT = castToComplexFloatMatrix(slaveCplx);

        // casting FloatMatrix to DoubleMatrix
        DoubleMatrix masterMagnitudeDouble = castToDoubleMatrix(masterMagnitude);
        System.out.println(masterMagnitudeDouble.get(0, 0));
        DoubleMatrix slaveMagnitudeDouble = castToDoubleMatrix(slaveMagnitude);
        System.out.println(slaveMagnitudeDouble.get(0, 0));

        Coregistration coreg = new Coregistration();
        DoubleMatrix tempDouble = coreg.correlate(SarUtils.magnitude(masterCplx), SarUtils.magnitude(slaveCplx));

        FloatMatrix tempFloat = correlate(masterMagnitude, slaveMagnitude);
        FloatMatrix tempFloat_ESTIMATED = correlate(magnitude(masterCplx_FLOAT), magnitude(slaveCplx_FLOAT));

        System.out.println("Double: " + tempDouble.max());
        System.out.println("Float: " + tempFloat.max());
        System.out.println("Temp Float: " + tempFloat_ESTIMATED.max());
        System.out.println("Original: " + correlMasterSlave.max());

    }

    @Ignore
    @Test
    public void testMagFFT() throws Exception {

        // FFT
        // -------------------------------------------------------------------------------
        // declare file names
        processingPath = "/d2/test.processing/unit_tests/etna.volcano/process/crop/01486_21159.cpm/";
        dataPath = "coarseCorr/magFFT/";
        masterFileName = "master.cr4";
        slaveFileName = "slave.cr4";

        // dimensions
        rows = 256;
        cols = 256;

        masterCplx = DataReader.readCplxFloatData(processingPath + dataPath + masterFileName, rows, cols, littleEndian);
        slaveCplx = DataReader.readCplxFloatData(processingPath + dataPath + slaveFileName, rows, cols, littleEndian);


        // dump values to check whether everything is parsed right
        System.out.println(masterCplx.get(0, 0));
        System.out.println(slaveCplx.get(0, 0));

        // dump values to check whether everything is parsed right
        System.out.println("Master Dimensions: " + masterCplx.rows + ", " + masterCplx.columns);
        System.out.println("Slave Dimensions: " + slaveCplx.rows + ", " + slaveCplx.columns);

        Coregistration coreg = new Coregistration();
        int MasksizeP = 256;
        int MasksizeL = 256;
        double coherence = coreg.crosscorrelate(masterCplx, slaveCplx, 16, MasksizeL / 2, MasksizeP / 2, 0, 0);
        System.out.println("coherence = " + coherence);

    }


    @Ignore
    @Before
    public void setUpFineCoregData_magFFT() throws Exception {

        // CORRELATE
        // -------------------------------------------------------------------------------
        // declare file names
        processingPath = "/d2/test.processing/unit_tests/etna.volcano/process/crop/01486_21159.cpm/";
        dataPath = "fineCoreg/magFFT/";
        masterFileName = "master_patch.32x32.cr4";
        slaveFileName = "slave_patch.32x32.cr4";

        // dimensions
        rows = 32;
        cols = 32;

        // load all data
        masterCplx = DataReader.readCplxFloatData(processingPath + dataPath + masterFileName, rows, cols, littleEndian);
        slaveCplx = DataReader.readCplxFloatData(processingPath + dataPath + slaveFileName, rows, cols, littleEndian);


    }

    @Ignore
    @Test
    public void testFineMagFFT() throws Exception {

        // dump values to check whether everything is parsed right
        System.out.println(masterCplx.get(0, 0));
        System.out.println(slaveCplx.get(0, 0));

        // dump values to check whether everything is parsed right
        System.out.println("Master Dimensions: " + masterCplx.rows + ", " + masterCplx.columns);
        System.out.println("Slave Dimensions: " + slaveCplx.rows + ", " + slaveCplx.columns);

        Coregistration coreg = new Coregistration();
        int MasksizeL = rows;
        int MasksizeP = cols;
        int ovsfactor = 1;

        double coherence = coreg.crosscorrelate(masterCplx, slaveCplx, ovsfactor, MasksizeL / 2, MasksizeP / 2, 0, 0);
        System.out.println("coherence = " + coherence);

    }

    @Ignore
    @Test
    public void testShiftSpectrum() throws Exception {

        SLCImage minfo = new SLCImage();
        minfo.parseResFile(new File(processingPath + "01486.res"));

        SLCImage sinfo = new SLCImage();
        sinfo.parseResFile(new File(processingPath + "21159.res"));

        final int m_pixlo = 1212;
        final double s_pixlo = 1214;// neg.shift -> 0
        final double mPrf = minfo.getPRF();
        final double sPrf = sinfo.getPRF();
        final double mRsr2x = minfo.getRsr2x();
        final double sRsr2x = sinfo.getRsr2x();

        double[] mFdc = new double[3];
        mFdc[0] = minfo.doppler.getF_DC_a0();
        mFdc[1] = minfo.doppler.getF_DC_a1();
        mFdc[2] = minfo.doppler.getF_DC_a2();

        double[] sFdc = new double[3];
        sFdc[0] = sinfo.doppler.getF_DC_a0();
        sFdc[1] = sinfo.doppler.getF_DC_a1();
        sFdc[2] = sinfo.doppler.getF_DC_a2();

        Coregistration coreg = new Coregistration();

        coreg.shiftazispectrum(masterCplx, mPrf, mRsr2x, mFdc, -m_pixlo);// shift from fDC to zero
        coreg.shiftazispectrum(slaveCplx, sPrf, sRsr2x, sFdc, -s_pixlo);// shift from fDC to zero

    }

    @Ignore
    @Test
    public void testShiftSpectrumAndMagFFT() throws Exception {

        SLCImage minfo = new SLCImage();
        minfo.parseResFile(new File(processingPath + "01486.res"));

        SLCImage sinfo = new SLCImage();
        sinfo.parseResFile(new File(processingPath + "21159.res"));

        final int m_pixlo = 1212;
        final double s_pixlo = 1214;// neg.shift -> 0

        final int ovsFactorCorrelate = 8;
        final int ovsFactor = 2;

        double mPrf = minfo.getPRF();
        double sPrf = sinfo.getPRF();
        double mRsr2x = minfo.getRsr2x();
        double sRsr2x = sinfo.getRsr2x();

        double[] mFdc = new double[3];
        mFdc[0] = minfo.doppler.getF_DC_a0();
        mFdc[1] = minfo.doppler.getF_DC_a1();
        mFdc[2] = minfo.doppler.getF_DC_a2();

        double[] sFdc = new double[3];
        sFdc[0] = sinfo.doppler.getF_DC_a0();
        sFdc[1] = sinfo.doppler.getF_DC_a1();
        sFdc[2] = sinfo.doppler.getF_DC_a2();

        Coregistration coreg = new Coregistration();

        coreg.shiftazispectrum(masterCplx, mPrf, mRsr2x, mFdc, -m_pixlo);// shift from fDC to zero
        coreg.shiftazispectrum(slaveCplx, sPrf, sRsr2x, sFdc, -s_pixlo);// shift from fDC to zero

        final int AccL = masterCplx.rows / 2;
        final int AccP = masterCplx.columns / 2;

        masterCplx = SarUtils.oversample(masterCplx, ovsFactor, ovsFactor);
        slaveCplx = SarUtils.oversample(slaveCplx, ovsFactor, ovsFactor);

        int offsetL = 0;
        int offsetP = 0;

        double coherence = coreg.crosscorrelate(masterCplx, slaveCplx, ovsFactorCorrelate / ovsFactor, 2 * AccL, 2 * AccP, offsetL, offsetP);

        System.out.println("coherence = " + coherence);

    }

    @Ignore
    @Before
    public void setUpFineCoregData_CoherenceSpace() throws Exception {

        // CORRELATE
        // -------------------------------------------------------------------------------
        // declare file names
        processingPath = "/d2/test.processing/unit_tests/etna.volcano/process/crop/01486_21159.cpm/";
        dataPath = "fineCoreg/magSpace/";
        masterFileName = "master_patch.48x48.cr4";
        slaveFileName = "slave_patch.48x48.cr4";

        // dimensions
        rows = 48;
        cols = 48;

        // load all data
        masterCplx = DataReader.readCplxFloatData(processingPath + dataPath + masterFileName, rows, cols, littleEndian);
        slaveCplx = DataReader.readCplxFloatData(processingPath + dataPath + slaveFileName, rows, cols, littleEndian);


    }

    @Ignore
    @Test
    public void testCoherenceSpace() throws Exception {

        final int ovsFactor = 8;
        final int AccL = 8;
        final int AccP = 8;
        int offsetL = 0;
        int offsetP = 0;

        Coregistration coreg = new Coregistration();
        coreg.coherencespace(AccL, AccP, ovsFactor, masterCplx, slaveCplx, offsetL, offsetP);// shift from fDC to zero

    }

//    @Ignore
    @Test
    public void testResampling() throws Exception {

        logger.setLevel(Level.TRACE);
        logger.trace("Start Resampling [development code]");

        // PARAMETERS
        // ----------------------------------
        processingPath = "/d2/test.processing/unit_tests/etna.volcano/process/crop/01486_21159.cpm/";
        dataPath = "rsmp/";

        SLCImage master = new SLCImage();
        master.parseResFile(new File(processingPath + "01486.res"));
        master.setOriginalWindow(new Window(1, 26292, 1, 4900));

        SLCImage slave = new SLCImage();
        slave.parseResFile(new File(processingPath + "21159.res"));
/*
        // Estimated during CPM step: where inverse estimation is performed :
        // ...not really clear why estimation is performed?
        // ...I can just invert polynomials and apply them to slave?
        Deltaline_slave00_poly:                    	1.38941008e+02
        Deltapixel_slave00_poly:                   	-2.19844746e+00
        Deltaline_slave0N_poly:                    	1.38856333e+02
        Deltapixel_slave0N_poly:                   	-2.39968790e+00
        Deltaline_slaveN0_poly:                    	1.38911145e+02
        Deltapixel_slaveN0_poly:                   	-2.19893253e+00
        Deltaline_slaveNN_poly:                    	1.38936348e+02
        Deltapixel_slaveNN_poly:                   	-2.36860850e+00
*/
        slave.setSlaveMasterOffset(1.38941008e+02, -2.19844746e+00, 1.38856333e+02, -2.39968790e+00,
                1.38911145e+02, -2.19893253e+00, 1.38936348e+02, -2.36860850e+00);

        boolean shiftazi = true;
        int memory = 500;
        String method = "cc6p";

        int polyOrder = 2;
        int polyCoeffs = PolyUtils.numberOfCoefficients(polyOrder);
        // polynomial : 5th degree : AZIMUTH direction
        double[] cpmL = new double[polyCoeffs + 1];
        cpmL[0] = -138.60102699999999;
        cpmL[1] = -0.61342283399999997;
        cpmL[2] = 0.31677070099999999;
        cpmL[3] = 0.086886668400000006;
        cpmL[4] = -0.58635061099999997;
        cpmL[5] = -0.016647044699999999;

        // polynomial : 5th degree : RANGE direction
        double[] cpmP = new double[polyCoeffs + 1];
        cpmP[0] = 2.2027052999999999;
        cpmP[1] = 0.32394611499999998;
        cpmP[2] = -0.60784787200000001;
        cpmP[3] = -0.52443622899999998;
        cpmP[4] = -0.16831470000000001;
        cpmP[5] = -0.684955745;

        // PROCESSING
        // ----------------------------------
        if (shiftazi == true) {
            logger.info("Shifting kernel_L to data fDC");
        }

        final int BUFFERMEMSIZE = memory;

        final int Npoints = extractNumber(method); // #pnts interpolator
        logger.debug("Number of kernel points: {}", Npoints);

        if (MathUtils.isOdd(Npoints)) {
            logger.error("Resample only even point interpolators, defined number of points: {}", Npoints);
            throw new IllegalArgumentException();
        }

        final int Npointsd2 = Npoints / 2;
        final int Npointsd2m1 = Npointsd2 - 1;
//        final int sizeofci16 = sizeof(compli16);
//        final int sizeofcr4 = sizeof(complr4);

        // Normalize data for polynomial
        final double minL = master.getOriginalWindow().linelo;
        final double maxL = master.getOriginalWindow().linehi;
        final double minP = master.getOriginalWindow().pixlo;
        final double maxP = master.getOriginalWindow().pixhi;

        logger.info("resample: polynomial normalized by factors [AZIMUTH]: {} {} to [-2,2]", minL, maxL);
        logger.info("resample: polynomial normalized by factors [RANGE]: {} {} to [-2,2]", minP, maxP);


        // For KNAB/Raised Cosine kernel if requested
        // ...Because kernel is same in az. and rg. min. must be used.
        float CHI_az = (float) (slave.getPRF() / slave.getAzimuthBandwidth());// oversampling factor az
        float CHI_rg = (float) ((slave.getRsr2x() / 2.0) / slave.getRangeBandwidth());// oversampling factor rg
        float CHI = Math.min(CHI_az, CHI_rg);// min. oversampling factor of data
        logger.info("Oversampling ratio azimuth (PRF/ABW): {}", CHI_az);
        logger.info("Oversampling ratio range (RSR/RBW): {}", CHI_rg);
        logger.info("KNAB/RC kernel uses: oversampling ratio: {}", CHI);

        if (CHI < 1.1) {
            logger.warn("Oversampling ratio: {} not optimal for KNAB/RC", CHI);
        }

        /** Create lookup table */
        // ........ e.g. four point interpolator
        // ........ interpolating point: p=6.4925
        // ........ required points: 5, 6, 7, 8
        // ........ kernel number from lookup table: floor(.4925*interval+.5)
        // ........  table[0]= 0 1 0 0 ;table[interval]= 0 0 1 0
        // ........ intervals in lookup table: dx
        // ........ for high doppler 100 is OK (fdc=3prf; 6pi --> 10deg error?)
        final int INTERVAL = 127; // precision: 1./interval [pixel]
        final int Ninterval = INTERVAL + 1; // size of lookup table
        final double dx = 1.0d / INTERVAL; // interval look up table
        logger.info("resample: lookup table size: " + Ninterval);

        /** Notes:
         *  ...Lookup table complex because of multiplication with complex
         *  ...Loopkup table for azimuth and range and
         *  ...shift spectrum of azi kernel with doppler centroid
         *  ...kernel in azimuth should be sampled higher
         *  ...and may be different from range due to different oversampling ratio and spectral shift (const)
         */
        LUT lut = new LUT(LUT.CC6P, Npoints);
        lut.constructLUT();
//        lut.overviewOfLut();

        final ComplexDoubleMatrix pntKernelAz = new ComplexDoubleMatrix(lut.getKernel());
        final ComplexDoubleMatrix pntKernelRg = new ComplexDoubleMatrix(lut.getKernel());
        final DoubleMatrix pntAxis = lut.getAxis();

        // Degree of coregistration polynomial
        final int degree_cpmL = PolyUtils.degreeFromCoefficients(cpmL.length);
        final int degree_cpmP = PolyUtils.degreeFromCoefficients(cpmP.length);

        // Compute overlap between master and slave
        Window overlap = Coregistration.getOverlap(master, slave, (double) Npointsd2, 0d, 0d);
//        Window oldOverlap = Coregistration.getOverlap(master, slave, cpmL, cpmP);
        logger.debug("[New method] Overlap between master and slave, in master coordinates: {}", overlap.toString());

        // get the data
        // -------------------------------------------------------------------------------
        // declare file names
        String fileName = "resample_buffer.cr4";

        int nRows = 2750;
        int nCols = 550;

        ComplexDoubleMatrix buffer = DataReader.readCplxFloatData(processingPath + dataPath + fileName, nRows, nCols, littleEndian);
        logger.info("Buffer size: {} rows, {} cols", nRows, nCols);


        // ----------------------------------------------------------------------------

        // overlap between master and slave in MASTER COORDINATE SYSTEM

        Window masterWindow = new Window();
        Window slaveWindow = new Window();

//        Window fullMasterCrop = new Window(minL, maxL, minP, maxP);


        int firstline = 0;
        int lastline = 0;

        int line;
        int pixel;

        int percent = 0;
        int tenpercent = (int) (Math.rint(overlap.lines() / 10.0)); // round
        if (tenpercent == 0)
            tenpercent = 1000; // avoid error: x%0


        // loop that does the job!
        for (line = (int) overlap.linelo; line <= overlap.linehi; line++) {
            // Progress messages
            if (((line - overlap.linelo) % tenpercent) == 0) {
                logger.info("RESAMPLE: {}%", percent);
                percent += 10;
            }

            /*
            * if buffer full write the data, empty the buffer
            * */

            /*
            * if new buffer is required ask for more data
            * */

            /*
            * if all data fiddling is in place then perform the actual resampling line be line, pixel by pixel
            * */

            // FOR TEST I DO HAVE THE FIRST BUFFER AS PULLED FROM DORIS THAT I SHOULD RUN RESAMPLE ON AND CHECK WHETHER IT MATCHES THE OUTPUT!


        }


    }

    /*public static ComplexDoubleMatrix readCplxIntData(final String fileName,
                                                      final int rows, final int columns,
                                                      final ByteOrder byteOrder) throws FileNotFoundException {

        final FlatBinaryInt inRealFile = new FlatBinaryInt();
        inRealFile.setFile(new File(fileName));
        inRealFile.setByteOrder(byteOrder);
        inRealFile.setDataWindow(new Window(0, (rows - 1), 0, (2 * columns - 1)));
        inRealFile.setInStream();
        inRealFile.readFromStream();

        // parse data from :: assume it is stored in "major-row order"
        DoubleMatrix realData = new DoubleMatrix(rows, columns);
        DoubleMatrix imgData = new DoubleMatrix(rows, columns);
        final float[][] data = inRealFile.getData();
        int cnt;
        for (int i = 0; i < rows; i++) {
            cnt = 0;
            for (int j = 0; j < 2 * columns; j = j + 2) {
                realData.put(i, cnt, data[i][j]);
                imgData.put(i, cnt, data[i][j + 1]);
                cnt++;
            }
        }

        return new ComplexDoubleMatrix(realData, imgData);

    }*/


    private int extractNumber(String line) {
        String numbers = new String();

        Pattern p = Pattern.compile("\\d+");
        Matcher m = p.matcher(line);
        while (m.find()) {
            numbers = numbers + m.group();
        }

        return Integer.parseInt(numbers);
    }


    private DoubleMatrix castToDoubleMatrix(FloatMatrix matrixFloat) {
        int rows = matrixFloat.rows;
        int columns = matrixFloat.columns;
        DoubleMatrix matrixDouble = new DoubleMatrix(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                matrixDouble.put(i, j, matrixFloat.get(i, j));
            }
        }
        return matrixDouble;
    }

    private FloatMatrix castToFloatMatrix(DoubleMatrix matrixDouble) {
        int rows = matrixDouble.rows;
        int columns = matrixDouble.columns;
        FloatMatrix matrixFloat = new FloatMatrix(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                matrixFloat.put(i, j, (float) matrixDouble.get(i, j));
            }
        }
        return matrixFloat;
    }

    private ComplexFloatMatrix castToComplexFloatMatrix(ComplexDoubleMatrix matrixComplexDouble) {
        int rows = matrixComplexDouble.rows;
        int columns = matrixComplexDouble.columns;
        ComplexFloatMatrix matrixComplexFloat = new ComplexFloatMatrix(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                ComplexDouble complexDouble = matrixComplexDouble.get(i, j);
                matrixComplexFloat.put(i, j, new ComplexFloat((float) complexDouble.real(), (float) complexDouble.imag()));
            }
        }
        return matrixComplexFloat;
    }


    public FloatMatrix correlate(FloatMatrix A, FloatMatrix Mask) {

        double varM = 0.; // variance of Mask
        Mask.subi(Mask.mean());

        for (int ii = 0; ii < Mask.length; ii++) {
            varM += Math.pow(Mask.get(ii), 2); // 1/N later
        }

        // Compute correlation at these points
        int beginl = (Mask.rows - 1) / 2; // floor
        int beginp = (Mask.columns - 1) / 2; // floor

        FloatMatrix Result = FloatMatrix.zeros(A.rows, A.columns); // init to 0
        FloatMatrix Am = new FloatMatrix(Mask.rows, Mask.columns);

        // First Window of A, updated at end of loop______
        Window winA = new Window(0, Mask.rows - 1, 0, Mask.columns - 1);
        Window windef = new Window();// defaults to total Am

        // Correlate part of Result______
        for (int i = beginl; i < A.rows - beginl; i++) {
            for (int j = beginp; j < A.columns - beginp; j++) {

                // Am.setdata(windef, A, winA); // Am no allocs.
                setdata(Am, windef, A, winA);

                Am.subi(Am.mean()); // center around mean
                float covAM = (float) 0.; // covariance A,Mask
                float varA = (float) 0.; // variance of A(part)

                for (int l = 0; l < Mask.length; l++) {
                    covAM += (Mask.get(l) * Am.get(l));
                    varA += Math.pow(Am.get(l), 2);
                }

                Result.put(i, j, (float) (covAM / Math.sqrt(varM * varA)));
                winA.pixlo++;
                winA.pixhi++;
            }
            winA.linelo++;
            winA.linehi++;
            winA.pixlo = 0;
            winA.pixhi = winA.pixlo + Mask.columns - 1;
        }
        return Result;

    }


    /**
     * setdata(outMatrix, outWin, inMatrix, inWin):
     * set outWin of outMatrix to inWin of inMatrix
     * if outWin==0 defaults to totalB, inWin==0 defaults to totalA
     * first line matrix =0 (?)
     */
    public static void setdata(FloatMatrix outMatrix, Window outWin, FloatMatrix inMatrix, Window inWin) {

        if (outWin.linehi == 0 && outWin.pixhi == 0) {
            outWin.linehi = outMatrix.rows - 1;
            outWin.pixhi = outMatrix.columns - 1;
        }
        if (inWin.linehi == 0 && inWin.pixhi == 0) {
            inWin.linehi = inMatrix.rows - 1;
            inWin.pixhi = inMatrix.columns - 1;
        }

        if (((outWin.linehi - outWin.linelo) != (inWin.linehi - inWin.linelo)) ||
                ((outWin.pixhi - outWin.pixlo) != (inWin.pixhi - inWin.pixlo))) {
            throw new IllegalArgumentException("setdata: wrong input.");

        }
        if (outWin.linehi < outWin.linelo || outWin.pixhi < outWin.pixlo) {
            throw new IllegalArgumentException("setdata: wrong input.1");
        }

        if ((outWin.linehi > outMatrix.rows - 1) ||
                (outWin.pixhi > outMatrix.columns - 1)) {
            throw new IllegalArgumentException("setdata: wrong input.2");
        }

        if ((inWin.linehi > inMatrix.rows - 1) ||
                (inWin.pixhi > inMatrix.columns - 1)) {
            throw new IllegalArgumentException("setdata: wrong input.3");
        }

        //// Fill data ////
        int sizeLin = (int) inWin.lines();
        for (int i = (int) outWin.pixlo, j = (int) inWin.pixlo; i <= outWin.pixhi; i++, j++) {

            int startOut = (int) (i * outMatrix.rows + outWin.linelo);
            int startIn = (int) (j * inMatrix.rows + inWin.linelo);

            System.arraycopy(inMatrix.data, startIn, outMatrix.data, startOut, sizeLin);

        }
    }

    public FloatMatrix intensity(final ComplexFloatMatrix inputMatrix) {
        return pow(inputMatrix.real(), 2).add(pow(inputMatrix.imag(), 2));
    }

    public FloatMatrix magnitude(final ComplexFloatMatrix inputMatrix) {
        return sqrt(intensity(inputMatrix));
    }


}

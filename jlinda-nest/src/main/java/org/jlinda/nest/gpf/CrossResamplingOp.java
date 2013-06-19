package org.jlinda.nest.gpf;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.Logger;
import com.bc.ceres.core.ProgressMonitor;
import org.apache.commons.lang.ArrayUtils;
import org.esa.beam.framework.datamodel.*;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.Tile;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.Parameter;
import org.esa.beam.framework.gpf.annotations.SourceProduct;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.nest.datamodel.AbstractMetadata;
import org.esa.nest.datamodel.Unit;
import org.esa.nest.gpf.OperatorUtils;
import org.jlinda.core.SLCImage;
import org.jlinda.core.coregistration.LUT;
import org.jlinda.core.coregistration.cross.CrossGeometry;
import org.slf4j.LoggerFactory;

import javax.media.jai.*;
import java.awt.*;
import java.awt.image.DataBuffer;
import java.awt.image.RenderedImage;
import java.awt.image.renderable.ParameterBlock;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Image resampling for Cross Interferometry
 */
@OperatorMetadata(alias = "CrossResampling",
        category = "InSAR Tools\\Interpolation",
        authors = "Petar Marinkovic",
        copyright = "Copyright (C) 2013 by Array Systems Computing Inc.",
        description = "Estimate Resampling Polynomial using SAR Image Geometry, and Resample Input Images")
public class CrossResamplingOp extends Operator {

    private static final Logger logger = (Logger) LoggerFactory.getLogger(CrossResamplingOp.class);

    @SourceProduct
    private Product sourceProduct;
    @TargetProduct
    private Product targetProduct;

    @Parameter(description = "The order of polynomial function", valueSet = {"1", "2", "3"}, defaultValue = "2",
            label = "Interpolation Polynomial Order")
    private int warpPolynomialOrder = 2;

    // only complex data accepted
    @Parameter(valueSet = {LUT.CC4P, LUT.CC6P, LUT.TS6P, LUT.TS8P, LUT.TS16P}, defaultValue = LUT.CC6P, label = "Interpolation Method")
    private String interpolationMethod = LUT.CC6P;

    // only complex data accepted
    @Parameter(valueSet = {"ERS 1/2", "Envisat ASAR"}, defaultValue = "ERS 1/2", label = "Target Geometry")
    private String targetGeometry = "ERS 1/2";

    private Interpolation interp = null;
    private InterpolationTable interpTable = null;

    private Band masterBand = null;
    private Band masterBand2 = null;
    private boolean complexCoregistration = true;
    private boolean warpDataAvailable = false;

    // Processing Variables
    // target
    private double targetPRF;
    private double targetRSR;
    // source
    private double sourcePRF;
    private double sourceRSR;

    private SLCImage slcMetadata = null;
    private double[] coeffsAz;
    private double[] coeffsRg;

    // PARAMETERS FOR JAI INTERPOLATION KERNEL
    private final static int SUBSAMPLE_BITS = 7;
    private final static int PRECISION_BITS = 32;

    // ERS NOMINAL PRF and RSR
    private final static double ERS_PRF_NOMINAL = 1679.902; // [Hz]
    private final static double ERS_RSR_NOMINAL = 18.962468 * 1000; // [Hz]

    // ASAR NOMINAL PRF and RSR
    private final static double ASAR_PRF_NOMINAL = 1652.4156494140625;       // [Hz]
    private final static double ASAR_RSR_NOMINAL = 19.20768 * 1000;  // [Hz]

    private final Map<Band, Band> sourceRasterMap = new HashMap<Band, Band>(10);
    private final Map<Band, Band> complexSrcMap = new HashMap<Band, Band>(10);
    private final Map<Band, WarpData> warpDataMap = new HashMap<Band, WarpData>(10);

    private String processedSlaveBand;
    private String[] masterBandNames = null;

    /**
     * Default constructor. The graph processing framework
     * requires that an operator has a default constructor.
     */
    public CrossResamplingOp() {
    }

    /**
     * Initializes this operator and sets the one and only target product.
     * <p>The target product can be either defined by a field of type {@link org.esa.beam.framework.datamodel.Product} annotated with the
     * {@link org.esa.beam.framework.gpf.annotations.TargetProduct TargetProduct} annotation or
     * by calling {@link #setTargetProduct} method.</p>
     * <p>The framework calls this method after it has created this operator.
     * Any client code that must be performed before computation of tile data
     * should be placed here.</p>
     *
     * @throws org.esa.beam.framework.gpf.OperatorException
     *          If an error occurs during operator initialisation.
     * @see #getTargetProduct()
     */

    @Override
    public void initialize() throws OperatorException {

        logger.setLevel(Level.TRACE);

        try {

            final MetadataElement absRoot = AbstractMetadata.getAbstractedMetadata(sourceProduct);
            slcMetadata = new SLCImage(absRoot);

            final String mission = slcMetadata.getMission();

            // arrange bands
            masterBand = sourceProduct.getBandAt(0);
            if (masterBand.getUnit() != null && masterBand.getUnit().equals(Unit.REAL) && sourceProduct.getNumBands() > 1) {
                complexCoregistration = true;
                masterBand2 = sourceProduct.getBandAt(1);
            }

            if (!mission.equals("ERS") || !mission.equals("ASAR")) {
                throw new OperatorException("The Cross Interferometry operator is for ERS 1/2 and Envisat ASAR products only");
            }

            if (mission.equals(targetGeometry)) {
                throw new OperatorException("Some smart warning / exception message here");
            }

            // declare source
            sourcePRF = slcMetadata.getPRF();
            sourceRSR = slcMetadata.getRsr2x();

            // declare target conditionally
            if (mission.equals("ERS")) {
                targetPRF = ASAR_PRF_NOMINAL;
                targetRSR = ASAR_RSR_NOMINAL;
            } else if (mission.equals("ASAR")) {
                targetPRF = ERS_PRF_NOMINAL;
                targetRSR = ERS_RSR_NOMINAL;
            }


            constructPolynomial();
            constructInterpolationTable(interpolationMethod);

            createTargetProduct();


        } catch (Throwable e) {
            OperatorUtils.catchOperatorException(getId(), e);
        }
    }

    private int extractNumber(String line) {
        String numbers = new String();

        Pattern p = Pattern.compile("\\d+");
        Matcher m = p.matcher(line);
        while (m.find()) {
            numbers = numbers + m.group();
        }

        return Integer.parseInt(numbers);
    }


    private void constructPolynomial() {

        CrossGeometry crossGeometry = new CrossGeometry();

        crossGeometry.setPrfOriginal(sourcePRF);
        crossGeometry.setRsrOriginal(sourceRSR);

        crossGeometry.setPrfTarget(targetPRF);
        crossGeometry.setRsrTarget(targetRSR);

        // used for normalization in estimation
        crossGeometry.setDataWindow(slcMetadata.getCurrentWindow());
        crossGeometry.computeCoefficients();

        coeffsAz = crossGeometry.getCoeffsAz();
        coeffsRg = crossGeometry.getCoeffsRg();

        // show polynomials
        logger.debug("coeffsAZ : estimated with PolyUtils.polyFit2D : {}", ArrayUtils.toString(coeffsAz));
        logger.debug("coeffsRg : estimated with PolyUtils.polyFit2D : {}", ArrayUtils.toString(coeffsRg));

    }

    private void constructInterpolationTable(String interpolationMethod) {

        // construct interpolation LUT
        int kernelLength = extractNumber(interpolationMethod);
        LUT lut = new LUT(interpolationMethod, kernelLength);
        lut.constructLUT();

        // get lut and cast to float
        double[] lutArrayDoubles = lut.getKernel().toArray();
        float lutArrayFloats[] = new float[lutArrayDoubles.length];
        int i = 0;
        for (double lutElement : lutArrayDoubles) {
            lutArrayFloats[i++] = (float) lutElement;
        }

        // construct interpolation table for JAI resampling
        int padding = kernelLength / 2 - 1;
        interpTable = new InterpolationTable(padding, kernelLength, SUBSAMPLE_BITS, PRECISION_BITS, lutArrayFloats);


    }

    /**
     * Create target product.
     */
    private void createTargetProduct() {

        targetProduct = new Product(sourceProduct.getName(),
                sourceProduct.getProductType(),
                sourceProduct.getSceneRasterWidth(),
                sourceProduct.getSceneRasterHeight());

        // coregistrated image should have the same geo-coding as the master image
        OperatorUtils.copyProductNodes(sourceProduct, targetProduct);
        updateTargetProductMetadata();
    }

    /**
     * Update metadata in the target product.
     */
    private void updateTargetProductMetadata() {
        final MetadataElement absTgt = AbstractMetadata.getAbstractedMetadata(targetProduct);
        AbstractMetadata.setAttribute(absTgt, AbstractMetadata.coregistered_stack, 1);
    }

    /**
     * Called by the framework in order to compute a tile for the given target band.
     * <p>The default implementation throws a runtime exception with the message "not implemented".</p>
     *
     * @param targetBand The target band.
     * @param targetTile The current tile associated with the target band to be computed.
     * @param pm         A progress monitor which should be used to determine computation cancelation requests.
     * @throws org.esa.beam.framework.gpf.OperatorException
     *          If an error occurs during computation of the target raster.
     */
    @Override
    public void computeTile(Band targetBand, Tile targetTile, ProgressMonitor pm) throws OperatorException {

        final Rectangle targetRectangle = targetTile.getRectangle();
        final int x0 = targetRectangle.x;
        final int y0 = targetRectangle.y;
        final int w = targetRectangle.width;
        final int h = targetRectangle.height;
        //System.out.println("CrossResamplingOperator: x0 = " + x0 + ", y0 = " + y0 + ", w = " + w + ", h = " + h);

        try {

            final Band srcBand = sourceRasterMap.get(targetBand);
            if (srcBand == null)
                return;
            Band realSrcBand = complexSrcMap.get(srcBand);
            if (realSrcBand == null)
                realSrcBand = srcBand;

            // create source image
            final Tile sourceRaster = getSourceTile(srcBand, targetRectangle);

            if (pm.isCanceled())
                return;

            final WarpData warpData = warpDataMap.get(realSrcBand);
            if (warpData.notEnoughGCPs)
                return;

            final RenderedImage srcImage = sourceRaster.getRasterDataNode().getSourceImage();

            // get warped image
            final RenderedOp warpedImage = createWarpImage(warpData.jaiWarp, srcImage);

            // copy warped image data to target
            final float[] dataArray = warpedImage.getData(targetRectangle).getSamples(x0, y0, w, h, 0, (float[]) null);

            targetTile.setRawSamples(ProductData.createInstance(dataArray));

        } catch (Throwable e) {
            OperatorUtils.catchOperatorException(getId(), e);
        } finally {
            pm.done();
        }
    }

    /**
     * Compute WARP polynomial function using master and slave GCP pairs.
     *
     * @param warpData            Stores the warp information per band.
     * @param warpPolynomialOrder The WARP polynimal order.
     * @param masterGCPGroup      The master GCPs.
     */
    public static void computeWARPPolynomial(
            final WarpData warpData, final int warpPolynomialOrder, final ProductNodeGroup<Placemark> masterGCPGroup) {

        warpData.computeWARP(warpPolynomialOrder);

    }

    /**
     * Create warped image.
     *
     * @param warp     The WARP polynomial.
     * @param srcImage The source image.
     * @return The warped image.
     */
    private RenderedOp createWarpImage(WarpPolynomial warp, final RenderedImage srcImage) {

        // reformat source image by casting pixel values from ushort to float
        final ParameterBlock pb1 = new ParameterBlock();
        pb1.addSource(srcImage);
        pb1.add(DataBuffer.TYPE_FLOAT);
        final RenderedImage srcImageFloat = JAI.create("format", pb1);

        // get warped image
        final ParameterBlock pb2 = new ParameterBlock();

        pb2.addSource(srcImageFloat);
        pb2.add(warp);
        pb2.add(interpTable);

        RenderedOp warpOutput = JAI.create("warp", pb2);

        return warpOutput;
    }

    public static class WarpData {

        public final java.util.List<Placemark> slaveGCPList = new ArrayList<Placemark>();
        private WarpPolynomial jaiWarp = null;
        public double[] xCoef = null;
        public double[] yCoef = null;

        public int numValidGCPs = 0;
        public boolean notEnoughGCPs = false;
        public float[] rms = null;
        public float[] rowResiduals = null;
        public float[] colResiduals = null;
        public float[] masterGCPCoords = null;
        public float[] slaveGCPCoords = null;

        public double rmsStd = 0;
        public double rmsMean = 0;
        public double rowResidualStd = 0;
        public double rowResidualMean = 0;
        public double colResidualStd = 0;
        public double colResidualMean = 0;

        public WarpData(ProductNodeGroup<Placemark> slaveGCPGroup) {
            for (int i = 0; i < slaveGCPGroup.getNodeCount(); ++i) {
                slaveGCPList.add(slaveGCPGroup.get(i));
            }
        }

        /**
         * Compute WARP function using master and slave GCPs.
         *
         * @param warpPolynomialOrder The WARP polynimal order.
         */
        public void computeWARP(final int warpPolynomialOrder) {

            // check if master and slave GCP coordinates are identical, if yes set the warp polynomial coefficients
            // directly, no need to compute them using JAI function because JAI produces incorrect result due to ill
            // conditioned matrix.
            float sum = 0.0f;
            for (int i = 0; i < slaveGCPCoords.length; i++) {
                sum += Math.abs(slaveGCPCoords[i] - masterGCPCoords[i]);
            }
            if (sum < 0.01) {
                switch (warpPolynomialOrder) {
                    case 1: {
                        xCoef = new double[3];
                        yCoef = new double[3];
                        xCoef[0] = 0;
                        xCoef[1] = 1;
                        xCoef[2] = 0;
                        yCoef[0] = 0;
                        yCoef[1] = 0;
                        yCoef[2] = 1;
                        break;
                    }
                    case 2: {
                        xCoef = new double[6];
                        yCoef = new double[6];
                        xCoef[0] = 0;
                        xCoef[1] = 1;
                        xCoef[2] = 0;
                        xCoef[3] = 0;
                        xCoef[4] = 0;
                        xCoef[5] = 0;
                        yCoef[0] = 0;
                        yCoef[1] = 0;
                        yCoef[2] = 1;
                        yCoef[3] = 0;
                        yCoef[4] = 0;
                        yCoef[5] = 0;
                        break;
                    }
                    case 3: {
                        xCoef = new double[10];
                        yCoef = new double[10];
                        xCoef[0] = 0;
                        xCoef[1] = 1;
                        xCoef[2] = 0;
                        xCoef[3] = 0;
                        xCoef[4] = 0;
                        xCoef[5] = 0;
                        xCoef[6] = 0;
                        xCoef[7] = 0;
                        xCoef[8] = 0;
                        xCoef[9] = 0;
                        yCoef[0] = 0;
                        yCoef[1] = 0;
                        yCoef[2] = 1;
                        yCoef[3] = 0;
                        yCoef[4] = 0;
                        yCoef[5] = 0;
                        yCoef[6] = 0;
                        yCoef[7] = 0;
                        yCoef[8] = 0;
                        yCoef[9] = 0;
                        break;
                    }
                    default:
                        throw new OperatorException("Incorrect WARP degree");
                }
                return;
            }

            // ToDo - pre and post scale source and destination with zero mean for numerical stability
            jaiWarp = WarpPolynomial.createWarp(slaveGCPCoords, //source
                    0,
                    masterGCPCoords, // destination
                    0,
                    2 * numValidGCPs,
                    1.0F,
                    1.0F,
                    1.0F,
                    1.0F,
                    warpPolynomialOrder);

            final float[] jaiXCoefs = jaiWarp.getXCoeffs();
            final float[] jaiYCoefs = jaiWarp.getYCoeffs();
            final int size = jaiXCoefs.length;
            xCoef = new double[size];
            yCoef = new double[size];
            for (int i = 0; i < size; ++i) {
                xCoef[i] = jaiXCoefs[i];
                yCoef[i] = jaiYCoefs[i];
            }
        }
    }

    /**
     * The SPI is used to register this operator in the graph processing framework
     * via the SPI configuration file
     * {@code META-INF/services/org.esa.beam.framework.gpf.OperatorSpi}.
     * This class may also serve as a factory for new operator instances.
     *
     * @see org.esa.beam.framework.gpf.OperatorSpi#createOperator()
     * @see org.esa.beam.framework.gpf.OperatorSpi#createOperator(java.util.Map, java.util.Map)
     */
    public static class Spi extends OperatorSpi {
        public Spi() {
            super(CrossResamplingOp.class);
//            super.setOperatorUI(CrossResamplingOpUI.class);
        }
    }

}

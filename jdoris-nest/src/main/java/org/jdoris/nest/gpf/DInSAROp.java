package org.jdoris.nest.gpf;

import com.bc.ceres.core.ProgressMonitor;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.MetadataElement;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.Tile;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.Parameter;
import org.esa.beam.framework.gpf.annotations.SourceProduct;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.nest.dataio.ReaderUtils;
import org.esa.nest.datamodel.AbstractMetadata;
import org.esa.nest.datamodel.Unit;
import org.esa.nest.gpf.OperatorUtils;
import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;
import org.jblas.Solve;
import org.jdoris.core.Baseline;
import org.jdoris.core.Orbit;
import org.jdoris.core.SLCImage;
import org.jdoris.core.Window;
import org.jdoris.core.utils.LinearAlgebraUtils;
import org.jdoris.core.utils.MathUtils;
import org.jdoris.core.utils.PolyUtils;
import org.jdoris.core.utils.SarUtils;
import org.jdoris.nest.utils.BandUtilsDoris;
import org.jdoris.nest.utils.CplxContainer;
import org.jdoris.nest.utils.ProductContainer;
import org.jdoris.nest.utils.TileUtilsDoris;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;

import static org.jdoris.core.utils.PolyUtils.normalize2;


@OperatorMetadata(alias = "DInSAR",
        category = "InSAR\\Products",
        description = "Differential Interferometry", internal = false)
public class DInSAROp extends Operator {

    @SourceProduct
    private Product sourceProduct;

    @TargetProduct
    private Product targetProduct;

    @Parameter(interval = "(1, 10]",
            description = "Degree of orbit interpolation polynomial",
            defaultValue = "3",
            label = "Orbit Interpolation Degree")
    private int orbitDegree = 3;


    // source maps
    private HashMap<Integer, CplxContainer> masterMap = new HashMap<Integer, CplxContainer>();
    private HashMap<Integer, CplxContainer> slaveMap = new HashMap<Integer, CplxContainer>();

    // target maps
    private HashMap<String, ProductContainer> targetMap = new HashMap<String, ProductContainer>();

    // operator tags
    private static final boolean CREATE_VIRTUAL_BAND = true;
    private static final String PRODUCT_NAME = "dinsar";
    public static final String PRODUCT_TAG = "_dinsar";

    private int sourceImageWidth;
    private int sourceImageHeight;
    //    private static final int RATIO_DEGREE = 3;
    private DoubleMatrix baselineRatioPolynomial;
    private CplxContainer master;
    private CplxContainer slaveDefo;
    private CplxContainer slaveTopo;


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
        try {

            constructSourceMetadata();
            constructTargetMetadata();
            createTargetProduct();

            getSourceImageDimension();

            // ratio polynomial

            baselineRatio();

        } catch (Exception e) {
            throw new OperatorException(e);
        }
    }

    private void baselineRatio() throws Exception {
        int[] slaveKeys = new int[2];
        int masterKey = 0;

        for (Integer keyMaster : masterMap.keySet()) {

            masterKey = keyMaster;

            int i = 0;
            for (Integer keySlave : slaveMap.keySet()) {

                slaveKeys[i++] = keySlave;
            }


        }

        master = masterMap.get(masterKey);
        slaveDefo = slaveMap.get(slaveKeys[0]);
        slaveTopo = slaveMap.get(slaveKeys[1]);
        estimateRatioPolynomial(master.metaData, master.orbit, slaveDefo.metaData, slaveDefo.orbit, slaveTopo.metaData, slaveTopo.orbit);
    }

    private void getSourceImageDimension() {
        sourceImageWidth = sourceProduct.getSceneRasterWidth();
        sourceImageHeight = sourceProduct.getSceneRasterHeight();
    }

    private void createTargetProduct() {

        // construct target product
        targetProduct = new Product(PRODUCT_NAME,
                sourceProduct.getProductType(),
                sourceProduct.getSceneRasterWidth(),
                sourceProduct.getSceneRasterHeight());

        OperatorUtils.copyProductNodes(sourceProduct, targetProduct);

        for (final Band band : targetProduct.getBands()) {
            targetProduct.removeBand(band);
        }

        for (String key : targetMap.keySet()) {

            String targetBandName_I = targetMap.get(key).targetBandName_I;
            String targetBandName_Q = targetMap.get(key).targetBandName_Q;
            targetProduct.addBand(targetBandName_I, ProductData.TYPE_FLOAT64);
            targetProduct.addBand(targetBandName_Q, ProductData.TYPE_FLOAT64);

            final String tag0 = targetMap.get(key).sourceMaster.date;
            final String tag1 = targetMap.get(key).sourceSlave.date;
            if (CREATE_VIRTUAL_BAND) {
                String countStr = "_" + PRODUCT_TAG + "_" + tag0 + "_" + tag1;
                ReaderUtils.createVirtualIntensityBand(targetProduct, targetProduct.getBand(targetBandName_I), targetProduct.getBand(targetBandName_Q), countStr);
                ReaderUtils.createVirtualPhaseBand(targetProduct, targetProduct.getBand(targetBandName_I), targetProduct.getBand(targetBandName_Q), countStr);
            }

        }

        // For testing: the optimal results with 1024x1024 pixels tiles, not clear whether it's platform dependent?
        // targetProduct.setPreferredTileSize(512, 512);


    }

    private void constructTargetMetadata() {

        for (Integer keyMaster : masterMap.keySet()) {

            CplxContainer master = masterMap.get(keyMaster);

            int counter = 0;

            for (Integer keySlave : slaveMap.keySet()) {

                if (counter == 0) {
                    // generate name for product bands
                    final String productName = keyMaster.toString() + "_" + keySlave.toString();

                    final CplxContainer slave = slaveMap.get(keySlave);
                    final ProductContainer product = new ProductContainer(productName, master, slave, true);

                    product.targetBandName_I = "i_" + PRODUCT_TAG + "_" + master.date + "_" + slave.date;
                    product.targetBandName_Q = "q_" + PRODUCT_TAG + "_" + master.date + "_" + slave.date;

                    // put ifg-product bands into map
                    targetMap.put(productName, product);

                    counter++;

                }

            }
        }
    }


    private void constructSourceMetadata() throws Exception {

        // define sourceMaster/sourceSlave name tags
        final String masterTag = "ifg";
        final String slaveTag = "dummy";

        // get sourceMaster & sourceSlave MetadataElement
        final MetadataElement masterMeta = AbstractMetadata.getAbstractedMetadata(sourceProduct);
        final String slaveMetadataRoot = AbstractMetadata.SLAVE_METADATA_ROOT;

        /* organize metadata */

        // put sourceMaster metadata into the masterMap
        metaMapPut(masterTag, masterMeta, sourceProduct, masterMap);

        // pug sourceSlave metadata into slaveMap
        MetadataElement[] slaveRoot = sourceProduct.getMetadataRoot().getElement(slaveMetadataRoot).getElements();
        for (MetadataElement meta : slaveRoot) {
            metaMapPut(slaveTag, meta, sourceProduct, slaveMap);
        }
    }

    private void metaMapPut(final String tag,
                            final MetadataElement root,
                            final Product product,
                            final HashMap<Integer, CplxContainer> map) throws Exception {

        // TODO: include polarization flags/checks!
        // pull out band names for this product
        final String[] bandNames = product.getBandNames();
        final int numOfBands = bandNames.length;

        // map key: ORBIT NUMBER
        int mapKey = root.getAttributeInt(AbstractMetadata.ABS_ORBIT);

        // metadata: construct classes and define bands
        final String date = OperatorUtils.getAcquisitionDate(root);
        final SLCImage meta = new SLCImage(root);
        final Orbit orbit = new Orbit(root, orbitDegree);

        // TODO: resolve multilook factors
        meta.setMlAz(1);
        meta.setMlRg(1);

        Band bandReal = null;
        Band bandImag = null;

        for (int i = 0; i < numOfBands; i++) {
            String bandName = bandNames[i];
            if (bandName.contains(tag) && bandName.contains(date)) {
                final Band band = product.getBandAt(i);
                if (BandUtilsDoris.isBandReal(band)) {
                    bandReal = band;
                } else if (BandUtilsDoris.isBandImag(band)) {
                    bandImag = band;
                }
            }
        }
        try {
            map.put(mapKey, new CplxContainer(date, meta, orbit, bandReal, bandImag));
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    private void estimateRatioPolynomial(SLCImage masterMeta, Orbit masterOrbit, SLCImage slaveMeta, Orbit slaveOrbit, SLCImage topoMeta, Orbit topoOrbit) throws Exception {

        // ______ Normalization factors for polynomial ______
        final int minL = 0;
        final int maxL = sourceImageHeight;
        final int minP = 0;
        final int maxP = sourceImageWidth;

        // ====== Model perpendicular baseline for master and slave ======
        // ______ compute B on grid every 500 lines, 100 pixels
        // ______ in window for topo/defo ______
//        int numpointsL = 20;                                  // grid for modelling
//        int numpointsP = 10;                                  // grid for modelling

//        final int numberOfCoefficients = PolyUtils.numberOfCoefficients(RATIO_DEGREE);

        final int numberOfPoints = 200;

        final int[][] position = MathUtils.distributePoints(numberOfPoints, new Window(minL, maxL, minP, maxP));

        double[] LINENUMBER = new double[numberOfPoints];
        double[] PIXELNUMBER = new double[numberOfPoints];

        for (int i = 0; i < position.length; i++) {
            LINENUMBER[i] = position[i][0];
            PIXELNUMBER[i] = position[i][1];
        }

        // model baselines
        Baseline defoBaseline = new Baseline();
        defoBaseline.model(masterMeta, slaveMeta, masterOrbit, slaveOrbit);

        Baseline topoBaseline = new Baseline();
        topoBaseline.model(masterMeta, topoMeta, masterOrbit, topoOrbit);

        double lastline = -1.0;

        double[] bperpTopo = new double[LINENUMBER.length];
        double[] bperpDefo = new double[LINENUMBER.length];

        for (int i = 0; i < LINENUMBER.length; ++i) {
            bperpTopo[i] = topoBaseline.getBperp(LINENUMBER[i], PIXELNUMBER[i], 0);
            bperpDefo[i] = defoBaseline.getBperp(LINENUMBER[i], PIXELNUMBER[i], 0);
        }

        double[] baselineRatio = new double[bperpDefo.length];

        for (int l = 0; l < bperpDefo.length; l++) {
            baselineRatio[l] = bperpDefo[l] / bperpTopo[l];
        }

        // ______ Set designmatrix, compute normalmatrix, righthandside ______
        DoubleMatrix A = new DoubleMatrix(baselineRatio.length, 3);
        for (int i = 0; i < A.rows; ++i) {
            A.put(i, 0, 1.0);
            A.put(i, 2, normalize2(PIXELNUMBER[i], minP, maxP));
            A.put(i, 1, normalize2(LINENUMBER[i], minL, maxL));
        }

        DoubleMatrix N = LinearAlgebraUtils.matTxmat(A, A);
        DoubleMatrix rhs = LinearAlgebraUtils.matTxmat(A, new DoubleMatrix(baselineRatio));


        baselineRatioPolynomial = Solve.solve(N, rhs);

//        // setup observation and design matrix
//        DoubleMatrix y = new DoubleMatrix(srpNumberPoints);
//        DoubleMatrix A = new DoubleMatrix(srpNumberPoints, numberOfCoefficients);


    }


    private void checkUserInput() throws OperatorException {
        // check for the logic in input paramaters
        final MetadataElement masterMeta = AbstractMetadata.getAbstractedMetadata(sourceProduct);
        final int isCoregStack = masterMeta.getAttributeInt(AbstractMetadata.coregistered_stack);
        if (isCoregStack != 1) {
            throw new OperatorException("Input should be a coregistered SLC stack");
        }
    }


    /**
     * Called by the framework in order to compute a tile for the given target band.
     * <p>The default implementation throws a runtime exception with the message "not implemented".</p>
     *
     * @param targetTileMap   The target tiles associated with all target bands to be computed.
     * @param targetRectangle The rectangle of target tile.
     * @param pm              A progress monitor which should be used to determine computation cancelation requests.
     * @throws org.esa.beam.framework.gpf.OperatorException
     *          If an error occurs during computation of the target raster.
     */
    @Override
    public void computeTileStack(Map<Band, Tile> targetTileMap, Rectangle targetRectangle, ProgressMonitor pm)
            throws OperatorException {
        try {

            int y0 = targetRectangle.y;
            int yN = y0 + targetRectangle.height - 1;
            int x0 = targetRectangle.x;
            int xN = targetRectangle.x + targetRectangle.width - 1;
//            final Window tileWindow = new Window(y0, yN, x0, xN);

            Band targetBand_I;
            Band targetBand_Q;
            ComplexDoubleMatrix complexDefo = null;
            ComplexDoubleMatrix complexTopo = null;

            ProductContainer product = null;
            for (String ifgKey : targetMap.keySet()) {

                product = targetMap.get(ifgKey);

                /// check out results from source ///
                Tile tileReal = getSourceTile(product.sourceMaster.realBand, targetRectangle);
                Tile tileImag = getSourceTile(product.sourceMaster.imagBand, targetRectangle);
                complexDefo = TileUtilsDoris.pullComplexDoubleMatrix(tileReal, tileImag);

                for (int idxBand = 0; idxBand < sourceProduct.getBands().length; idxBand++) {

                    Band bandAt = sourceProduct.getBandAt(idxBand);
                    if (bandAt.getName().contains("unw") & bandAt.getUnit().contains(Unit.REAL)) {

                        /// check out results from source ///
                        tileReal = getSourceTile(bandAt, targetRectangle);
                        tileImag = getSourceTile(sourceProduct.getBandAt(idxBand + 1), targetRectangle);
                        complexTopo = TileUtilsDoris.pullComplexDoubleMatrix(tileReal, tileImag);

                    }
                }

            }


            if (baselineRatioPolynomial.length > 0) {

                // normalize range and azimuth axis
                DoubleMatrix rangeAxisNormalized = DoubleMatrix.linspace(x0, xN, complexDefo.columns);
                rangeAxisNormalized = normalizeDoubleMatrix(rangeAxisNormalized);

                DoubleMatrix azimuthAxisNormalized = DoubleMatrix.linspace(y0, yN, complexDefo.rows);
                azimuthAxisNormalized = normalizeDoubleMatrix(azimuthAxisNormalized);

                // estimate the phase on the grid
                DoubleMatrix ratio =
                        PolyUtils.polyval(azimuthAxisNormalized, rangeAxisNormalized,
                                baselineRatioPolynomial, PolyUtils.degreeFromCoefficients(baselineRatioPolynomial.length));

                // compute the reference phase
                ComplexDoubleMatrix complexReferencePhase =
                        new ComplexDoubleMatrix(MatrixFunctions.cos(ratio),
                                MatrixFunctions.sin(ratio));


                complexTopo.muli(complexReferencePhase); // no conjugate here!
            }

            SarUtils.computeIfg_inplace(complexDefo, complexTopo.conji());

            /// commit to target ///
            targetBand_I = targetProduct.getBand(product.targetBandName_I);
            Tile tileOutReal = targetTileMap.get(targetBand_I);
            TileUtilsDoris.pushDoubleMatrix(complexDefo.real(), tileOutReal, targetRectangle);

            targetBand_Q = targetProduct.getBand(product.targetBandName_Q);
            Tile tileOutImag = targetTileMap.get(targetBand_Q);
            TileUtilsDoris.pushDoubleMatrix(complexDefo.imag(), tileOutImag, targetRectangle);


        } catch (Throwable e) {
            OperatorUtils.catchOperatorException(getId(), e);

        }
    }

    private DoubleMatrix normalizeDoubleMatrix(DoubleMatrix matrix) {
        matrix.subi(0.5 * (1 + sourceImageWidth));
        matrix.divi(0.25 * (sourceImageWidth - 1));
        return matrix;
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
            super(DInSAROp.class);
        }
    }

}

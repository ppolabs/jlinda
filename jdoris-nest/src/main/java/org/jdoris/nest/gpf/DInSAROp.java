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
import org.esa.beam.framework.gpf.annotations.SourceProducts;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.nest.datamodel.AbstractMetadata;
import org.esa.nest.datamodel.Unit;
import org.esa.nest.gpf.OperatorUtils;
import org.esa.nest.gpf.ReaderUtils;
import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;
import org.jdoris.core.Baseline;
import org.jdoris.core.Orbit;
import org.jdoris.core.SLCImage;
import org.jdoris.core.utils.MathUtils;
import org.jdoris.core.utils.PolyUtils;
import org.jdoris.nest.utils.BandUtilsDoris;
import org.jdoris.nest.utils.CplxContainer;
import org.jdoris.nest.utils.ProductContainer;
import org.jdoris.nest.utils.TileUtilsDoris;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;

/*
 See DInSAR class in jdoris-core for implementation details
 ...core class is not being used because of threading issues
 ...core class will be integrated soon in the operator

 Known issues: polyFit will trigger warning that maxError is larger then allowed/expected,
    see DInSAR.java for more info and solution for this problem. It is on TODO list.
*/

// TODO: support for multiple DEFO interferometric pairs in sourceProduct

@OperatorMetadata(alias = "DInSAR",
        category = "InSAR\\Products",
        description = "Differential Interferometry", internal = false)
public class DInSAROp extends Operator {

    @SourceProducts(description = "Source products: InSAR DEFO pair product, and InSAR (unwrapped) TOPO pair product")
    private Product[] sourceProducts;

    @TargetProduct
    private Product targetProduct;

    @Parameter(interval = "(1, 10]",
            description = "Degree of orbit interpolation polynomial",
            defaultValue = "3",
            label = "Orbit Interpolation Degree")
    private int orbitDegree = 3;


    // source maps
    private HashMap<Integer, CplxContainer> masterDefoMap = new HashMap<Integer, CplxContainer>();
    private HashMap<Integer, CplxContainer> slaveDefoMap = new HashMap<Integer, CplxContainer>();

    private HashMap<Integer, CplxContainer> masterTopoMap = new HashMap<Integer, CplxContainer>();
    private HashMap<Integer, CplxContainer> slaveTopoMap = new HashMap<Integer, CplxContainer>();

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

    private Product topoProduct = null;
    private Product defoProduct = null;

    private final static int DEGREE_BASELINE_RATIO_MODEL = 1; // 1st degree polynomial is usually sufficient

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

            if (sourceProducts.length != 2) {
                throw new OperatorException("DInSAROp: requires 2 InSAR products, 'TOPO pair' and 'DEFO pairs' products.");
            }

            // work out which is which product: loop through source product and check which product has a 'unwrapped phase' band
            sortOutSourceProducts(); // -> declares topoProduct and defoProduct

            constructSourceMetadata();
            constructTargetMetadata();
            createTargetProduct();
            getSourceImageDimension();

            // compute topo/defo pair ratio polynomial needed for scaling of topoPair
            baselineRatio();

        } catch (Exception e) {
            throw new OperatorException(e);
        }
    }

    private void sortOutSourceProducts() {

        for (Product product : sourceProducts) {
            if (productHasABSPhase(product)) {
                topoProduct = product;
            } else {
                defoProduct = product;
            }
        }
        if (topoProduct == null || defoProduct == null) {
            throw new OperatorException("DInSAROp requires one product with Unwrapped Phase and one with Wrapped Phase");
        }
    }

    private boolean productHasABSPhase(Product product) {
        for (Band band : product.getBands()) {
            if (band.getUnit().contains(Unit.ABS_PHASE)) {
                return true;
            }
        }
        return false;
    }

    private void baselineRatio() throws Exception {
        int[] slaveKeys = new int[2];
        int masterKey = 0;

        for (Integer keyMaster : masterDefoMap.keySet()) {
            masterKey = keyMaster;
            int i = 0;
            for (Integer keySlave : slaveDefoMap.keySet()) {
                slaveKeys[i++] = keySlave;
            }
        }

        master = masterDefoMap.get(masterKey);
        slaveTopo = slaveTopoMap.get(slaveKeys[0]);
        slaveDefo = slaveDefoMap.get(slaveKeys[0]);

        estimateBperpRatioPolynomial(master.metaData, master.orbit,
                slaveDefo.metaData, slaveDefo.orbit,
                slaveTopo.metaData, slaveTopo.orbit);
    }

    private void getSourceImageDimension() {
        sourceImageWidth = defoProduct.getSceneRasterWidth();
        sourceImageHeight = defoProduct.getSceneRasterHeight();
    }

    private void createTargetProduct() {

        // construct target product
        targetProduct = new Product(PRODUCT_NAME,
                defoProduct.getProductType(),
                defoProduct.getSceneRasterWidth(),
                defoProduct.getSceneRasterHeight());

        OperatorUtils.copyProductNodes(defoProduct, targetProduct);

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

    }

    private void constructTargetMetadata() {

        for (Integer keyMaster : masterDefoMap.keySet()) {

            CplxContainer master = masterDefoMap.get(keyMaster);

            int counter = 0;

            for (Integer keySlave : slaveDefoMap.keySet()) {

                if (counter == 0) {
                    // generate name for product bands
                    final String productName = keyMaster.toString() + "_" + keySlave.toString();

                    final CplxContainer slave = slaveDefoMap.get(keySlave);
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

        /** --- DEFO PRODUCT -----*/
        // define sourceMaster/sourceSlave name tags
        String masterTag;
        String slaveTag;
        MetadataElement masterMeta;
        final String slaveMetadataRoot = AbstractMetadata.SLAVE_METADATA_ROOT;
        MetadataElement[] slaveRoot;

        masterTag = "ifg";
        slaveTag = "dummy";

        // get sourceMaster & sourceSlave MetadataElement
        masterMeta = AbstractMetadata.getAbstractedMetadata(defoProduct);

        /* organize metadata */

        // put sourceMaster metadata into the masterDefoMap
        metaMapPut(masterTag, masterMeta, defoProduct, masterDefoMap);

        // pug sourceSlave metadata into slaveDefoMap
        slaveRoot = defoProduct.getMetadataRoot().getElement(slaveMetadataRoot).getElements();
        for (MetadataElement meta : slaveRoot) {
            metaMapPut(slaveTag, meta, defoProduct, slaveDefoMap);
        }

        /** --- TOPO PRODUCT -----*/

        // get sourceMaster & sourceSlave MetadataElement
        masterMeta = AbstractMetadata.getAbstractedMetadata(topoProduct);

        /* organize metadata */

        // put sourceMaster metadata into the masterDefoMap
        metaMapPut(masterTag, masterMeta, topoProduct, masterTopoMap);

        // pug sourceSlave metadata into slaveDefoMap
        slaveRoot = topoProduct.getMetadataRoot().getElement(slaveMetadataRoot).getElements();
        for (MetadataElement meta : slaveRoot) {
            metaMapPut(slaveTag, meta, topoProduct, slaveTopoMap);
        }

    }

    private void metaMapPut(final String tag,
                            final MetadataElement root,
                            final Product product,
                            final HashMap<Integer, CplxContainer> map) throws Exception {

        // pull out band names for this product
        final String[] bandNames = product.getBandNames();
        final int numOfBands = bandNames.length;

        // map key: ORBIT NUMBER
        int mapKey = root.getAttributeInt(AbstractMetadata.ABS_ORBIT);

        // metadata: construct classes and define bands
        final String date = OperatorUtils.getAcquisitionDate(root);
        final SLCImage meta = new SLCImage(root);
        final Orbit orbit = new Orbit(root, orbitDegree);

        // TODO: mlook factores are hard-coded for now
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


    private void estimateBperpRatioPolynomial(SLCImage masterMeta, Orbit masterOrbit,
                                              SLCImage slaveDefoMeta, Orbit slaveDefoOrbit,
                                              SLCImage topoSlaveMeta, Orbit topoSlaveOrbit) throws Exception {

        // Normalization factors for polynomial
        double minL = 0;
        double maxL = sourceImageHeight;
        double minP = 0;
        double maxP = sourceImageWidth;

        // Model perpendicular baseline for master and defo
        // .....compute B on grid every 500 lines, 100 pixels
        // .....in window for topoData/defoData
        final int nPoints = 200;
        final int[][] positionArray = MathUtils.distributePoints(nPoints, masterMeta.getCurrentWindow());
        DoubleMatrix position = new DoubleMatrix(nPoints, 2);
        for (int i = 0; i < nPoints; i++) {
            position.put(i, 0, positionArray[i][0]);
            position.put(i, 1, positionArray[i][1]);
        }

        // model baselines
        Baseline topoBaseline = new Baseline();
        topoBaseline.model(masterMeta, topoSlaveMeta, masterOrbit, topoSlaveOrbit);

        Baseline defoBaseline = new Baseline();
        defoBaseline.model(masterMeta, slaveDefoMeta, masterOrbit, slaveDefoOrbit);

        DoubleMatrix bPerpTopo = new DoubleMatrix(nPoints);
        DoubleMatrix bPerpDefo = new DoubleMatrix(nPoints);

        for (int i = 0; i < nPoints; ++i) {
            bPerpTopo.put(i, topoBaseline.getBperp(position.get(i, 0), position.get(i, 1), 0));
            bPerpDefo.put(i, defoBaseline.getBperp(position.get(i, 0), position.get(i, 1), 0));
        }

        // model ratio bPerpDefo/bPerpTopo as linear polynomial
        //   ...r(l,p) = a00 + a10*l + a01*p
        //   ...give stats on max. error
        DoubleMatrix baselineRatio = bPerpDefo.div(bPerpTopo);

        DoubleMatrix normalizedAzimuthAxis = normalize(position.getColumn(0), minL, maxL);
        DoubleMatrix normalizedRangeAxis = normalize(position.getColumn(1), minP, maxP);

        baselineRatioPolynomial = new DoubleMatrix(PolyUtils.polyFit2D(normalizedAzimuthAxis, normalizedRangeAxis,
                baselineRatio, DEGREE_BASELINE_RATIO_MODEL));
    }

    private void checkUserInput() throws OperatorException {
        // check for the logic in input paramaters
        final MetadataElement masterMeta = AbstractMetadata.getAbstractedMetadata(defoProduct);
        final int isCoregStack = masterMeta.getAttributeInt(AbstractMetadata.coregistered_stack);
        if (isCoregStack != 1) {
            throw new OperatorException("Input should be a coregistered SLC stack");
        }
    }

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
            ComplexDoubleMatrix complexDefoPair = null;
            DoubleMatrix doubleTopoPair = null;

            ProductContainer product = null;
            for (String ifgKey : targetMap.keySet()) {

                product = targetMap.get(ifgKey);

                /// check out results from source ///
                Tile tileReal = getSourceTile(product.sourceMaster.realBand, targetRectangle);
                Tile tileImag = getSourceTile(product.sourceMaster.imagBand, targetRectangle);
                complexDefoPair = TileUtilsDoris.pullComplexDoubleMatrix(tileReal, tileImag);

            }

            // always pull from topoProduct
            for (Band band : topoProduct.getBands()) {
                if (band.getName().contains("unw") || band.getUnit().contains(Unit.ABS_PHASE)) {
                    /// check out results from source ///
                    Tile tileReal = getSourceTile(band, targetRectangle);
                    doubleTopoPair = TileUtilsDoris.pullDoubleMatrix(tileReal);
                }
            }

            ComplexDoubleMatrix scaledComplexTopoPair = null;
            if (baselineRatioPolynomial.length > 0) {

                // normalize azimuth and range axis
                DoubleMatrix azimuthAxisNormalized = DoubleMatrix.linspace(y0, yN, doubleTopoPair.rows);
                azimuthAxisNormalized = normalizeDoubleMatrix(azimuthAxisNormalized);

                DoubleMatrix rangeAxisNormalized = DoubleMatrix.linspace(x0, xN, doubleTopoPair.columns);
                rangeAxisNormalized = normalizeDoubleMatrix(rangeAxisNormalized);

                // evaluate scaling polynomial on tile grid
                DoubleMatrix baselineRatio = PolyUtils.polyval(azimuthAxisNormalized, rangeAxisNormalized,
                        baselineRatioPolynomial, PolyUtils.degreeFromCoefficients(baselineRatioPolynomial.length));

                // scale topo pair
                DoubleMatrix scaledDoubleTopoPair = doubleTopoPair.mul(baselineRatio);

                // convert to complex
                // NOTE(!!): sinus part in converting to complex numbers is .neg(), cheaper like this then with .conj()
                scaledComplexTopoPair = new ComplexDoubleMatrix(MatrixFunctions.cos(scaledDoubleTopoPair), MatrixFunctions.sin(scaledDoubleTopoPair).neg());
            }

            // subtract scaledComplexTopoPair from complexDefoPair
            // no .conj() since topoPair is converted to complex with .sin().neg()
            complexDefoPair.muli(scaledComplexTopoPair);

            /// commit complexDefoPair back to target ///
            targetBand_I = targetProduct.getBand(product.targetBandName_I);
            Tile tileOutReal = targetTileMap.get(targetBand_I);
            TileUtilsDoris.pushDoubleMatrix(complexDefoPair.real(), tileOutReal, targetRectangle);

            targetBand_Q = targetProduct.getBand(product.targetBandName_Q);
            Tile tileOutImag = targetTileMap.get(targetBand_Q);
            TileUtilsDoris.pushDoubleMatrix(complexDefoPair.imag(), tileOutImag, targetRectangle);

        } catch (Throwable e) {
            OperatorUtils.catchOperatorException(getId(), e);

        }
    }

    private DoubleMatrix normalizeDoubleMatrix(DoubleMatrix matrix) {
        matrix.subi(0.5 * (1 + sourceImageWidth));
        matrix.divi(0.25 * (sourceImageWidth - 1));
        return matrix;
    }

    private DoubleMatrix normalize(DoubleMatrix in, double min, double max) {
        DoubleMatrix out = in.dup();
        out.subi(.5 * (min + max));
        out.divi(.25 * (max - min));
        return out;
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
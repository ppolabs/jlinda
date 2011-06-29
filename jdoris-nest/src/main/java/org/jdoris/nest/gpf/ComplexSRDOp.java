package org.jdoris.nest.gpf;

import com.bc.ceres.core.ProgressMonitor;
import org.esa.beam.framework.datamodel.*;
import org.esa.beam.framework.dataop.dem.ElevationModel;
import org.esa.beam.framework.dataop.dem.ElevationModelDescriptor;
import org.esa.beam.framework.dataop.dem.ElevationModelRegistry;
import org.esa.beam.framework.dataop.resamp.Resampling;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.Tile;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.Parameter;
import org.esa.beam.framework.gpf.annotations.SourceProduct;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.nest.dataio.ReaderUtils;
import org.esa.nest.dataio.dem.FileElevationModel;
import org.esa.nest.datamodel.AbstractMetadata;
import org.esa.nest.datamodel.Unit;
import org.esa.nest.gpf.OperatorUtils;
import org.jblas.ComplexDoubleMatrix;
import org.jdoris.core.Orbit;
import org.jdoris.core.SLCImage;
import org.jdoris.nest.utils.BandUtilsDoris;
import org.jdoris.nest.utils.CplxContainer;
import org.jdoris.nest.utils.ProductContainer;
import org.jdoris.nest.utils.TileUtilsDoris;

import javax.media.jai.BorderExtender;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

@OperatorMetadata(alias = "ComplexSRD",
        category = "InSAR Products",
        description = "Compute and subtract TOPO phase", internal = false)
public final class ComplexSRDOp extends Operator {

    @SourceProduct
    private Product sourceProduct;

    @TargetProduct
    private Product targetProduct;

    @Parameter(interval = "[1, 10]",
            description = "Degree of orbit interpolation polynomial",
            defaultValue = "3",
            label = "Orbit poly degree")
    private int ORBIT_DEGREE = 3;

    @Parameter(valueSet = {"ACE", "GETASSE30", "SRTM 3Sec", "ASTER 1sec GDEM"},
            description = "The digital elevation model.", defaultValue = "SRTM 3Sec", label = "Digital Elevation Model")
    private String demName = "SRTM 3Sec";

    @Parameter(description = "The elevation band name.", defaultValue = "elevation", label = "Elevation Band Name")
    private String elevationBandName = "elevation";

    @Parameter(description = "The external DEM file.", defaultValue = " ", label = "External DEM")
    private String externalDEM = " ";

    @Parameter(valueSet = { NEAREST_NEIGHBOUR, BILINEAR, CUBIC }, defaultValue = BILINEAR,
                label="Resampling Method")
    private String resamplingMethod = BILINEAR;

    static final String NEAREST_NEIGHBOUR = "Nearest Neighbour";
    static final String BILINEAR = "Bilinear Interpolation";
    static final String CUBIC = "Cubic Convolution";

    private FileElevationModel fileElevationModel = null;
    private ElevationModel dem = null;
    private Band elevationBand = null;
    private float noDataValue = 0;

    // source
    private HashMap<Integer, CplxContainer> masterMap = new HashMap<Integer, CplxContainer>();
    private HashMap<Integer, CplxContainer> slaveMap = new HashMap<Integer, CplxContainer>();

    // target
    private HashMap<String, ProductContainer> targetMap = new HashMap<String, ProductContainer>();

    private static final boolean CREATE_VIRTUAL_BAND = true;

    private static final String PRODUCT_NAME = "srd_ifgs";
    public static final String PRODUCT_TAG = "srd";

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

            checkUserInput();
            constructSourceMetadata();
            constructTargetMetadata();
            createTargetProduct();

            // dem part
            defineDEM();

        } catch (Exception e) {
            throw new OperatorException(e);
        }
    }

    private void defineDEM() throws IOException {
        // dem part
        final ElevationModelRegistry elevationModelRegistry = ElevationModelRegistry.getInstance();
        final ElevationModelDescriptor demDescriptor = elevationModelRegistry.getDescriptor(demName);
        if (demDescriptor == null)
            throw new OperatorException("The DEM '" + demName + "' is not supported.");
        if (demDescriptor.isInstallingDem())
            throw new OperatorException("The DEM '" + demName + "' is currently being installed.");

        Resampling resampling = Resampling.NEAREST_NEIGHBOUR;
        if (externalDEM != null && !externalDEM.trim().isEmpty()) {
            fileElevationModel = new FileElevationModel(new File(externalDEM), resampling);
            noDataValue = fileElevationModel.getNoDataValue();
        } else {
            dem = demDescriptor.createDem(resampling);
            if (dem == null)
                throw new OperatorException("The DEM '" + demName + "' has not been installed.");
            noDataValue = dem.getDescriptor().getNoDataValue();
        }

        elevationBand = targetProduct.addBand(elevationBandName, ProductData.TYPE_FLOAT32);
        elevationBand.setSynthetic(true);
        elevationBand.setNoDataValue(noDataValue);
        elevationBand.setUnit(Unit.METERS);
        elevationBand.setDescription(demDescriptor.getName());
    }

    private void checkUserInput() {
        // TODO: use jdoris input.coherence class to check user input
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
        final Orbit orbit = new Orbit(root, ORBIT_DEGREE);
        Band bandReal = null;
        Band bandImag = null;

        // TODO: boy this is one ugly construction!?
        // loop through all band names(!) : and pull out only one that matches criteria
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

    private void constructTargetMetadata() {

        // this means there is only one slave! but still do it in the loop
        // loop through masters
        for (Integer keyMaster : masterMap.keySet()) {

            CplxContainer master = masterMap.get(keyMaster);

            for (Integer keySlave : slaveMap.keySet()) {

                // generate name for product bands
                String productName = keyMaster.toString() + "_" + keySlave.toString();

                final CplxContainer slave = slaveMap.get(keySlave);
                final ProductContainer product = new ProductContainer(productName, master, slave, false);

                product.targetBandName_I = "i_" + PRODUCT_TAG + "_" + master.date + "_" + slave.date;
                product.targetBandName_Q = "q_" + PRODUCT_TAG + "_" + master.date + "_" + slave.date;

                // put ifg-product bands into map
                targetMap.put(productName, product);

            }

        }
    }

    private void createTargetProduct() {

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

            if (CREATE_VIRTUAL_BAND) {
                final String tag0 = targetMap.get(key).sourceMaster.date;
                final String tag1 = targetMap.get(key).sourceSlave.date;
                String countStr = "_" + PRODUCT_TAG + "_" + tag0 + "_" + tag1;
                ReaderUtils.createVirtualIntensityBand(targetProduct, targetProduct.getBand(targetBandName_I), targetProduct.getBand(targetBandName_Q), countStr);
                ReaderUtils.createVirtualPhaseBand(targetProduct, targetProduct.getBand(targetBandName_I), targetProduct.getBand(targetBandName_Q), countStr);
            }

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
    public void computeTileStack(Map<Band, Tile> targetTileMap, Rectangle targetRectangle, ProgressMonitor pm) throws OperatorException {
        try {

            final GeoCoding geoCoding = targetProduct.getGeoCoding();
            final Rectangle rect = new Rectangle(targetRectangle);
            System.out.println("Original: x0 = " + rect.x + ", y = " + rect.y + ", w = " + rect.width + ", h = " + rect.height);
            int winRg = 1;
            int winAz = 1;
            final int x0 = rect.x - (winRg - 1) / 2;
            final int y0 = rect.y - (winAz - 1) / 2;
            final int w = rect.width + winRg - 1;
            final int h = rect.height + winAz - 1;
            rect.x = x0;
            rect.y = y0;
            rect.width = w;
            rect.height = h;

//            System.out.println("Shifted: x0 = " + x0 + ", y0 = " + y0 + ", w = " + w + ", h = " + h);
            final BorderExtender border = BorderExtender.createInstance(BorderExtender.BORDER_ZERO);
            Band targetBand_I;
            Band targetBand_Q;

            // anyhow do here for elevation
            Tile elevationTile = targetTileMap.get("elevation");
            final ProductData trgData = elevationTile.getDataBuffer();
            pm.beginTask("Computing elevations from " + demName + "...", h);
            try {
                final GeoPos geoPos = new GeoPos();
                final PixelPos pixelPos = new PixelPos();
                float elevation;

                for (int y = y0; y < y0 + h; ++y) {
                    for (int x = x0; x < x0 + w; ++x) {

                        pixelPos.setLocation(x + 0.5f, y + 0.5f);
                        geoCoding.getGeoPos(pixelPos, geoPos);
                        try {
                            if (fileElevationModel != null) {
                                elevation = fileElevationModel.getElevation(geoPos);
                            } else {
                                elevation = dem.getElevation(geoPos);
                            }
                        } catch (Exception e) {
                            elevation = noDataValue;
                        }

                        trgData.setElemDoubleAt(elevationTile.getDataBufferIndex(x, y), (short) Math.round(elevation));
                    }
                    pm.worked(1);
                }
            } finally {
                pm.done();
            }

            for (String ifgKey : targetMap.keySet()) {

                final ProductContainer product = targetMap.get(ifgKey);
                // check out from source
                Tile tileReal = getSourceTile(product.sourceMaster.realBand, rect, border);
                Tile tileImag = getSourceTile(product.sourceMaster.imagBand, rect, border);
                final ComplexDoubleMatrix dataMaster = TileUtilsDoris.pullComplexDoubleMatrix(tileReal, tileImag);// check out from source

                // check in to target
                targetBand_I = targetProduct.getBand(product.targetBandName_I);
                Tile tileOutReal = targetTileMap.get(targetBand_I);
                TileUtilsDoris.pushDoubleMatrix(dataMaster.real(), tileOutReal, rect);

                targetBand_Q = targetProduct.getBand(product.targetBandName_Q);
                Tile tileOutImag = targetTileMap.get(targetBand_Q);
                TileUtilsDoris.pushDoubleMatrix(dataMaster.imag(), tileOutImag, rect);

            }

        } catch (Exception e) {
            throw new OperatorException(e);
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
            super(ComplexSRDOp.class);
            setOperatorUI(ComplexSRDOpUI.class);
        }
    }

}

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
import org.esa.beam.util.ProductUtils;
import org.esa.nest.dataio.ReaderUtils;
import org.esa.nest.datamodel.AbstractMetadata;
import org.esa.nest.gpf.OperatorUtils;
import org.jblas.*;
import org.jdoris.core.Orbit;
import org.jdoris.core.SLCImage;
import org.jdoris.core.utils.SarUtils;
import org.jdoris.nest.utils.BandUtilsDoris;
import org.jdoris.nest.utils.CplxContainer;
import org.jdoris.nest.utils.IfgContainer;
import org.jdoris.nest.utils.TileUtilsDoris;

import javax.media.jai.BorderExtender;
import java.awt.*;
import java.util.HashMap;
import java.util.Map;


@OperatorMetadata(alias = "ComplexIfg",
        category = "InSAR Products",
        description = "Compute interferograms from stack of coregistered images : JBLAS implementation", internal = false)
public class ComplexIfgOp extends Operator {

    @SourceProduct
    private Product sourceProduct;

    @TargetProduct
    private Product targetProduct;

    @Parameter(valueSet = {"1", "2", "3", "4", "5", "6", "7", "8"},
            description = "Order of 'Flat earth phase' polynomial",
            defaultValue = "5",
            label = "Degree of \"Flat Earth\" polynomial")
    private int srpPolynomialDegree = 5;

    @Parameter(valueSet = {"301", "401", "501", "601", "701", "801", "901", "1001"},
            description = "Number of points for the 'flat earth phase' polynomial estimation",
            defaultValue = "501",
            label = "Number of 'Flat earth' estimation points")
    private int srpNumberPoints = 501;


    @Parameter(valueSet = {"1", "2", "3", "4", "5"},
            description = "Degree of orbit (polynomial) interpolator",
            defaultValue = "3",
            label = "Orbit interpolation degree")
    private int orbitPolynomialDegree = 3;

//        @Parameter(description = "Orbit interpolation method", valueSet = {"polynomial"},
//                defaultValue = "5",
//                label="SRP Polynomial Order")
//        private int srpPolynomialDegree = 5;

    private HashMap<Integer, CplxContainer> masterMap = new HashMap<Integer, CplxContainer>();
    private HashMap<Integer, CplxContainer> slaveMap = new HashMap<Integer, CplxContainer>();
    private HashMap<String, IfgContainer> ifgMap = new HashMap<String, IfgContainer>();

    private static final int ORBIT_DEGREE = 3; // hardcoded
    private static final boolean CREATE_VIRTUAL_BAND = true;

    private static final int TILE_OVERLAP_X = 0;
    private static final int TILE_OVERLAP_Y = 0;


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

            // getSourceImageGeocodings();
            // estimateFlatEarthPolynomial();
            // updateTargetProductMetadata();
            // updateTargetProductGeocoding();

            createTargetProduct();

        } catch (Exception e) {
            throw new OperatorException(e);
        }
    }

    private void constructTargetMetadata() {

        for (Integer keyMaster : masterMap.keySet()) {

            CplxContainer master = masterMap.get(keyMaster);

            String nameMaster_I = master.realBand.getName();
            String nameMaster_Q = master.imagBand.getName();

            int counter = 0;

            for (Integer keySlave : slaveMap.keySet()) {

                // generate name
                final CplxContainer slave = slaveMap.get(keySlave);

                String nameIfg = keyMaster.toString() + "_" + keySlave.toString();

                String nameSlave_I = slave.realBand.getName();
                String nameSlave_Q = slave.imagBand.getName();

                final String tagIfg = "ifg";

                String iBandName = "i_" + tagIfg + counter + "_" + nameMaster_I + "_" + nameSlave_I;
                String qBandName = "q_" + tagIfg + counter + "_" + nameMaster_Q + "_" + nameSlave_Q;

                // initialize container class
                final IfgContainer ifg = new IfgContainer(nameIfg, master, slave);

                ifg.iBandName = iBandName;
                ifg.qBandName = qBandName;

                ifgMap.put(nameIfg, ifg);

            }
        }
    }

    private void constructSourceMetadata() throws Exception {

//        masterMap = new HashMap<Integer, CplxContainer>();
//        slaveMap = new HashMap<Integer, CplxContainer>();

        final MetadataElement masterMeta = AbstractMetadata.getAbstractedMetadata(sourceProduct);

        metaMapPut("mst", masterMeta, sourceProduct, masterMap);

        final String slaveMetadataRoot = AbstractMetadata.SLAVE_METADATA_ROOT;
        MetadataElement[] slaveRoot = sourceProduct.getMetadataRoot().getElement(slaveMetadataRoot).getElements();
        for (MetadataElement meta : slaveRoot) {
            metaMapPut("slv", meta, sourceProduct, slaveMap);
        }

    }

    private void metaMapPut(final String tag,
                            final MetadataElement root, final Product product,
                            final HashMap<Integer, CplxContainer> map) throws Exception {

        final String[] bandNames = product.getBandNames();
        int numOfBands = bandNames.length;

        int key = root.getAttributeInt(AbstractMetadata.ABS_ORBIT);

        // metadata
        String date = OperatorUtils.getAcquisitionDate(root);
        SLCImage meta = new SLCImage(root);
        Orbit orbit = new Orbit(root, ORBIT_DEGREE);
        Band bandReal = null;
        Band bandImag = null;

        // TODO: include polarization flags/checks!
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
            map.put(key, new CplxContainer(date, meta, orbit, bandReal, bandImag));
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void createTargetProduct() throws Exception {

        // construct target product
        targetProduct = new Product("interferogram",
                sourceProduct.getProductType(),
                sourceProduct.getSceneRasterWidth(),
                sourceProduct.getSceneRasterHeight());


//        targetProduct.setPreferredTileSize(2,2);

        // copy product nodes
        OperatorUtils.copyProductNodes(sourceProduct, targetProduct);

        for (String key : ifgMap.keySet()) {

            // generate REAL band
            final IfgContainer ifg = ifgMap.get(key);
            final Band targetBandI = targetProduct.addBand(ifg.iBandName, ProductData.TYPE_FLOAT32);
            ProductUtils.copyRasterDataNodeProperties(ifg.master.realBand, targetBandI);

            // generate IMAGINARY band
            final Band targetBandQ = targetProduct.addBand(ifg.qBandName, ProductData.TYPE_FLOAT32);
            ProductUtils.copyRasterDataNodeProperties(ifg.master.imagBand, targetBandQ);

            // generate virtual bands
            if (CREATE_VIRTUAL_BAND) {
                ReaderUtils.createVirtualIntensityBand(targetProduct, targetBandI, targetBandQ, ("_" + key));
                ReaderUtils.createVirtualPhaseBand(targetProduct, targetBandI, targetBandQ, ("_" + key));
            }

        }
    }

    private void checkUserInput() {
        // check for the logic in input paramaters
    }

    private void updateTargetProductMetadata() {
        // update metadata of target product for the estimated polynomial
    }

    private void updateTargetProductGeocoding() {
        // update metadata of target product for the estimated polynomial
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

//            int x0 = targetRectangle.x;
//            int y0 = targetRectangle.y;
//            int w = targetRectangle.width;
//            int h = targetRectangle.height;
//            System.out.println("x0 = " + x0 + ", y0 = " + y0 + ", w = " + w + ", h = " + h);

            // target
            Band targetBand;
            Tile targetTile;

            // source
            Tile tileReal;
            Tile tileImag;

            final BorderExtender border = BorderExtender.createInstance(BorderExtender.BORDER_ZERO);

            final Rectangle rect = new Rectangle(targetRectangle);
            rect.width += TILE_OVERLAP_X;
            rect.height += TILE_OVERLAP_Y;

            for (String ifgTag : ifgMap.keySet()) {

                final IfgContainer ifg = ifgMap.get(ifgTag);

                // get master data
                tileReal = getSourceTile(ifg.master.realBand, rect, border);
                tileImag = getSourceTile(ifg.master.imagBand, rect, border);

                final ComplexDoubleMatrix masterMatrix = TileUtilsDoris.pullComplexDoubleMatrix(tileReal, tileImag);

                // get slave data
                tileReal = getSourceTile(ifg.slave.realBand, rect, border);
                tileImag = getSourceTile(ifg.slave.imagBand, rect, border);
                final ComplexDoubleMatrix slaveMatrix = TileUtilsDoris.pullComplexDoubleMatrix(tileReal, tileImag);

                // compute ifg
                final ComplexDoubleMatrix cplxIfg = SarUtils.computeIfg(masterMatrix, slaveMatrix);

                // save real band of computation
                targetBand = targetProduct.getBand(ifg.iBandName);
                tileReal = targetTileMap.get(targetBand);
                TileUtilsDoris.pushFloatMatrix(cplxIfg.real(), tileReal, targetRectangle);

                targetBand = targetProduct.getBand(ifg.qBandName);
                tileImag = targetTileMap.get(targetBand);
                TileUtilsDoris.pushFloatMatrix(cplxIfg.imag(), tileImag, targetRectangle);

//                // save imag band of computation : this is somehow too slow?
//                targetBand = targetProduct.getBand(ifg.qBandName);
//                tileReal = targetTileMap.get(targetBand);
//                targetBand = targetProduct.getBand(ifg.iBandName);
//                tileImag = targetTileMap.get(targetBand);
//                TileUtilsDoris.pushComplexFloatMatrix(cplxIfg, tileReal, tileImag, targetRectangle);


            }

        } catch (Throwable e) {
            OperatorUtils.catchOperatorException(getId(), e);
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
            super(ComplexIfgOp.class);
        }
    }


}

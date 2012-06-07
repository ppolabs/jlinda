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
import org.esa.nest.datamodel.AbstractMetadata;
import org.esa.nest.datamodel.Unit;
import org.esa.nest.gpf.OperatorUtils;
import org.jblas.DoubleMatrix;
import org.jblas.Solve;
import org.jdoris.core.Orbit;
import org.jdoris.core.Point;
import org.jdoris.core.SLCImage;
import org.jdoris.core.Window;
import org.jdoris.core.utils.MathUtils;
import org.jdoris.core.utils.PolyUtils;
import org.jdoris.nest.utils.BandUtilsDoris;
import org.jdoris.nest.utils.CplxContainer;
import org.jdoris.nest.utils.ProductContainer;
import org.jdoris.nest.utils.TileUtilsDoris;

import java.awt.*;
import java.util.HashMap;

import static org.jdoris.core.Constants.PI;
import static org.jdoris.core.Constants.SOL;
import static org.jdoris.core.utils.LinearAlgebraUtils.matTxmat;
import static org.jdoris.core.utils.PolyUtils.normalize2;
import static org.jdoris.core.utils.PolyUtils.polyFit;

@OperatorMetadata(alias = "Slant2Height",
        category = "InSAR\\Products",
        description = "Slant to Height conversion",
        internal = false)
public class Slant2HeightOp extends Operator {

    @SourceProduct(description = "Source product that contains unwrapped phase.")
    private Product sourceProduct;

    @TargetProduct
    private Product targetProduct;

    @Parameter(valueSet = {"100", "200", "300", "400", "500"},
            description = "Number of points for evaluation of flat earth phase at different altitudes",
            defaultValue = "200",
            label = "Number of estimation points")
    private int nPoints; // where ref.phase is evaluated

    @Parameter(valueSet = {"2", "3", "4", "5"},
            description = "Number of height samples in range [0,5000)",
            defaultValue = "3",
            label = "Number of height samples")
    private int nHeights;

    @Parameter(valueSet = {"1", "2", "3", "4", "5"},
            description = "Degree of the 1D polynomial to fit reference phase through.",
            defaultValue = "2",
            label = "Degree of 1D polynomial")
    private int degree1D; // only possible now.

    @Parameter(valueSet = {"1", "2", "3", "4", "5", "6", "7", "8"},
            description = "Degree of the 2D polynomial to fit reference phase through.",
            defaultValue = "5",
            label = "Degree of 2D polynomial")
    private int degree2D; // only possible now.

    @Parameter(valueSet = {"2", "3", "4", "5"},
            description = "Degree of orbit (polynomial) interpolator",
            defaultValue = "3",
            label = "Orbit interpolation degree")
    private int orbitDegree = 3;

    // source maps
    private HashMap<Integer, CplxContainer> masterMap = new HashMap<Integer, CplxContainer>();
    private HashMap<Integer, CplxContainer> slaveMap = new HashMap<Integer, CplxContainer>();

    // target maps
    private HashMap<String, ProductContainer> targetMap = new HashMap<String, ProductContainer>();

    // operator tags
    private static final boolean CREATE_VIRTUAL_BAND = true;
    private static final String PRODUCT_NAME = "slant2h";
    public static final String PRODUCT_TAG = "slant2h";

    private int sourceImageWidth;
    private int sourceImageHeight;

    private double minPhi;
    private double maxPhi;


    private static final int MAXHEIGHT = 5000;
    private static final int TEN = 10;

    private Band referenceBand = null;
    private DoubleMatrix rhs;

    @Override
    public void initialize() throws OperatorException {

        try {

            // work out which is which product: loop through source product and check which product has a 'unwrapped phase' band
            sortOutSourceProducts(); // -> declares topoProduct and defoProduct

            constructSourceMetadata();
            constructTargetMetadata();
            createTargetProduct();
            getSourceImageDimension();

            scwabisch();

        } catch (Exception e) {
            throw new OperatorException(e);
        }

    }

    private void scwabisch() throws Exception {

        final int heightStep = MAXHEIGHT / (nHeights - 1); // heights to eval ref.refPhase

        // Matrices for storing refPhase for all ref. ellipsoids
        //  refPhase(i,0)  refPhase for height 0
        //  refPhase(i,1)  refPhase for height Heigthsep * 1
        //  refPhase(i,Nh) refPhase for height 4000
//        Map<Integer, DoubleMatrix> refPhaseMap = Maps.newLinkedHashMap();

        DoubleMatrix refPhaseMatrix = new DoubleMatrix(nPoints, nHeights);

        // Distribute points in original master system (not multilooked)
        // (i,0): line, (i,1): pixel, (i,2) flagfromdisk (not used here)
        int[][] positionArray = MathUtils.distributePoints(nPoints, new Window(1, sourceImageHeight, 1, sourceImageWidth));

        DoubleMatrix Position = new DoubleMatrix(nPoints, 2);
        for (int i = 0; i < nPoints; i++) {
            Position.put(i, 0, positionArray[i][0]);
            Position.put(i, 1, positionArray[i][1]);
        }


        for (Integer keyMaster : masterMap.keySet()) {
            CplxContainer master = masterMap.get(keyMaster);
            for (Integer keySlave : slaveMap.keySet()) {
                CplxContainer slave = slaveMap.get(keySlave);


                /** ----------------------------------------------------------------------------*/
                /** -- STEP 1 : compute reference refPhase in N points for nHeights ------------*/
                /** ----------------------------------------------------------------------------*/

                // Compute reference refPhase in N points for height (numheight)
                //        logger.debug("S2H: schwabisch: STEP1: compute reference refPhase for nHeights.");
                DoubleMatrix refPhaseZero = new DoubleMatrix(nPoints);
                for (int heightIdx = 0; heightIdx < nHeights; heightIdx++) {

                    int height = heightIdx * heightStep;

                    DoubleMatrix refPhase = new DoubleMatrix(nPoints); // pseudo-observation

                    // Compute delta r for all points
                    for (int i = 0; i < nPoints; i++) {
                        double phase = computeReferencePhase((double) positionArray[i][0], (double) positionArray[i][1],
                                height,
                                master.metaData, slave.metaData, master.orbit, slave.orbit);
                        refPhase.put(i, phase);
                    }

                    // store refPhase at h = 0
                    if (height == 0) {
                        refPhaseZero = refPhase;
                    }

                    //  Subtract ref. refPhase at h=0 for all point
                    //  this is the same as adding reference refPhase for all in uint
                    refPhaseMatrix.putColumn(heightIdx, refPhase.sub(refPhaseZero));
                }

                /** ----------------------------------------------------------------------------*/
                /** -- STEP 2 : compute alpha coefficients of polynomials for these points -----*/
                /** ----------------------------------------------------------------------------*/

//                logger.debug("S2H: schwabisch: STEP2: estimate coefficients 1d polynomial.");
//
                //        DoubleMatrix design = new DoubleMatrix(nHeights, degree1D + 1); // design matrix
                DoubleMatrix alphas = new DoubleMatrix(nPoints, degree1D + 1); // pseudo-observation
                DoubleMatrix hei = new DoubleMatrix(nHeights, 1);
                for (int i = 0; i < nHeights; i++) {
                    hei.put(i, 0, i * heightStep); // 0, .., 5000
                }

                // normalize tile to [0,1]
                minPhi = refPhaseMatrix.min();
                maxPhi = refPhaseMatrix.max();
                normalize(refPhaseMatrix, minPhi, maxPhi);

                for (int i = 0; i < nPoints; i++) {// solve system for all points
                    alphas.putRow(i, new DoubleMatrix(polyFit(refPhaseMatrix.getRow(i), hei, degree1D)));
                } // loop over all points

                /** -------------------------------------------------------------------------------*/
                /** -- STEP 3 : Compute alpha_i coefficients of polynomials as function of (l,p) --*/
                /** -------------------------------------------------------------------------------*/

//                logger.debug("S2H: schwabisch: STEP3: estimate coefficients for 2d polynomial.");
                // Compute alpha_i coefficients of polynomials as function of (l,p)
                // ... alpha_i = sum(k,l) beta_kl l^k p^l;
                // ... Solve simultaneous for all betas
                // ... this does not seem to be possibly with my routine, so do per alfa_i
                final int Nunk = PolyUtils.numberOfCoefficients(degree2D); // Number of unknowns

                // ______ Check redundancy is done before? ______
                if (nPoints < Nunk) {
//                    logger.error("slant2hschwabisch: N_observations<N_unknowns (increase S2H_NPOINTS or decrease S2H_DEGREE2D.");
                    throw new IllegalArgumentException();
                }

                DoubleMatrix A = new DoubleMatrix(nPoints, Nunk); // designmatrix

                // Set up system of equations
                // .... Order unknowns: B00 B10 B01 B20 B11 B02 B30 B21 B12 B03 for degree=3
                double minL = 0; //Position.getColumn(0).min();
                double maxL = sourceImageHeight; //Position.getColumn(0).max();
                double minP = 0; //Position.getColumn(1).min();
                double maxP = sourceImageWidth; //Position.getColumn(1).max();

                for (int i = 0; i < nPoints; i++) {
                    // normalize coordinates
                    double posL = normalize2(Position.get(i, 0), minL, maxL);
                    double posP = normalize2(Position.get(i, 1), minP, maxP);

                    int index = 0;
                    for (int j = 0; j <= degree2D; j++) {
                        for (int k = 0; k <= j; k++) {
                            A.put(i, index, Math.pow(posL, j - k) * Math.pow(posP, k));
                            index++;
                        }
                    }
                }

                // Solve 2d polynomial system for alfas at these points
                DoubleMatrix N = matTxmat(A, A);
                rhs = matTxmat(A, alphas);
                DoubleMatrix Qx_hat = N;

                // Solve the normal equations for all alpha_i
                // Simultaneous solution doesn't work somehow
                for (int i = 0; i < rhs.getColumns(); ++i) {
                    DoubleMatrix rhs_alphai = rhs.getColumn(i);
                    rhs.putColumn(i, Solve.solveSymmetric(Qx_hat, rhs_alphai));
                }

//                // Test solution by inverse
//                Qx_hat = Solve.solveSymmetric(Qx_hat, DoubleMatrix.eye(Qx_hat.getRows()));
//                double maxdev = (N.mmul(Qx_hat).sub(DoubleMatrix.eye(Qx_hat.getRows()))).normmax();
////                logger.debug("s2h schwaebisch: max(abs(N*inv(N)-I)) = {}", maxdev);
//                if (maxdev > 0.01) {
////                    logger.warn("slant2h: possibly wrong solution. deviation from unity AtA*inv(AtA) = {} > 0.01", maxdev);
//                }

            }
        }


    }

    private void getSourceImageDimension() {
        sourceImageWidth = sourceProduct.getSceneRasterWidth();
        sourceImageHeight = sourceProduct.getSceneRasterHeight();
    }

    private void constructTargetMetadata() {

        for (Integer keyMaster : masterMap.keySet()) {

            CplxContainer master = masterMap.get(keyMaster);

            for (Integer keySlave : slaveMap.keySet()) {

                // generate name for product bands
                String productName = keyMaster.toString() + "_" + keySlave.toString();

                final CplxContainer slave = slaveMap.get(keySlave);
                final ProductContainer product = new ProductContainer(productName, master, slave, false);

                product.targetBandName_I = PRODUCT_TAG + "_" + master.date + "_" + slave.date;

                // put ifg-product bands into map
                targetMap.put(productName, product);

            }
        }
    }

    private void sortOutSourceProducts() {
        for (Band band : sourceProduct.getBands()) {
            if (band.getUnit().equals(Unit.ABS_PHASE)) {
                referenceBand = band;
            }
        }
        if (referenceBand == null) {
            throw new OperatorException("Slant2HeightOp requires minimum one 'unwrapped' phase band");
        }
    }

    private void constructSourceMetadata() throws Exception {

        // define sourceMaster/sourceSlave name tags
        final String masterTag = "ifg";
        final String slaveTag = "dummy";
        final MetadataElement masterMeta = AbstractMetadata.getAbstractedMetadata(sourceProduct);
        final String slaveMetadataRoot = AbstractMetadata.SLAVE_METADATA_ROOT;
        MetadataElement[] slaveRoot;

        /* organize metadata */
        // put sourceMaster metadata into the masterMap
        metaMapPut(masterTag, masterMeta, sourceProduct, masterMap);

        // pug sourceSlave metadata into slaveDefoMap
        slaveRoot = sourceProduct.getMetadataRoot().getElement(slaveMetadataRoot).getElements();
        for (MetadataElement meta : slaveRoot) {
            metaMapPut(slaveTag, meta, sourceProduct, slaveMap);
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
                } else if (product.getBandAt(i).getUnit().contains(Unit.ABS_PHASE)) {
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
            String bandName = targetMap.get(key).targetBandName_I;
            targetProduct.addBand(bandName, ProductData.TYPE_FLOAT32);
            targetProduct.getBand(bandName).setUnit(Unit.METERS);
        }

    }

    @Override
    public void computeTile(Band targetBand, Tile targetTile, ProgressMonitor pm) throws OperatorException {
        try {

            final Rectangle rect = targetTile.getRectangle();
//            System.out.println("Original: x0 = " + rect.x + ", y = " + rect.y + ", w = " + rect.width + ", h = " + rect.height);
            final int x0 = rect.x;
            final int y0 = rect.y;
            final int w = rect.width;
            final int h = rect.height;

            final int minL = 0;
            final int maxL = sourceImageHeight;
            final int minP = 0;
            final int maxP = sourceImageWidth;

            Window tileWindow = new Window(y0, y0 + h - 1, x0, x0 + w - 1);

            for (String absPhaseKey : targetMap.keySet()) {

                final ProductContainer product = targetMap.get(absPhaseKey);

                if (targetBand.getName().equals(product.targetBandName_I)) {

                    // check out from source
                    Tile tileRealMaster = getSourceTile(product.sourceMaster.realBand, rect);
                    final DoubleMatrix dataMaster = TileUtilsDoris.pullDoubleMatrix(tileRealMaster);// check out from source

//                    DoubleMatrix cohMatrix = SarUtils.coherence2(dataMaster, dataSlave, winAz, winRg);

                    double mlFacL = 1;//unwrappedinterf.multilookL;
                    double mlFacP = 1;//unwrappedinterf.multilookP;

                    // Number of lines/pixels of multilooked unwrapped interferogram
                    int mlLines = (int) (Math.floor((tileWindow.linehi - tileWindow.linelo + 1) / mlFacL));
                    int mlPixels = (int) (Math.floor((tileWindow.pixhi - tileWindow.pixlo + 1) / mlFacP));

                    // Line/pixel of first point in original master coordinates
                    double firstLine = (double) (tileWindow.linelo) + (mlFacL - 1.) / 2.;
                    double firstPixel = (double) (tileWindow.pixlo) + (mlFacP - 1.) / 2.;

                    // ant axis of pixel coordinates ______
                    DoubleMatrix p_axis = new DoubleMatrix(mlPixels, 1);
                    for (int i = 0; i < tileWindow.pixels(); i++) {
                        p_axis.put(i, 0, firstPixel + i * mlFacP);
                    }
                    normalize(p_axis, minP, maxP);

                    // ant axis for azimuth coordinates ______
                    DoubleMatrix l_axis = new DoubleMatrix(mlLines, 1);
                    for (int k = 0; k < tileWindow.lines(); k++) {
                        l_axis.put(k, 0, firstLine + k * mlFacL);
                    }
                    normalize(l_axis, minL, maxL);

                    // ---> Lookup table because not known in advance what degree1D is
                    DoubleMatrix[] pntALPHA = new DoubleMatrix[TEN];
                    for (int k = 0; k <= degree1D; k++) {
                        DoubleMatrix beta = new DoubleMatrix(PolyUtils.numberOfCoefficients(degree2D), 1);
                        for (int l = 0; l < PolyUtils.numberOfCoefficients(degree2D); l++) {
                            beta.put(l, 0, rhs.get(l, k)); // solution stored in rhs
                        }
                        pntALPHA[k] = PolyUtils.polyval(l_axis, p_axis, beta, degree2D);
                    }

                    // Evaluate h=f(l,p,phi) for all points in grid in BUFFER
                    double[] coeffThisPoint = new double[degree1D + 1]; //DoubleMatrix(degree1D + 1, 1);

                    DoubleMatrix BUFFER = dataMaster; //unwrappedinterf.readphase(bufferwin);

                    for (int line = 0; line < mlLines; line++) {
                        for (int pixel = 0; pixel < mlPixels; pixel++) {
                            // Check if unwrapped ok, else compute h
                            if (BUFFER.get(line, pixel) != Double.NaN) // else leave NaN
                            {
                                for (int k = 0; k < degree1D + 1; k++) {
                                    coeffThisPoint[k] = pntALPHA[k].get(line, pixel);
                                }
                                double data = BUFFER.get(line, pixel);
                                double x = PolyUtils.normalize2(data, minPhi, maxPhi);
                                double value = PolyUtils.polyVal1D(x, coeffThisPoint);
                                BUFFER.put(line, pixel, value);
                            }
                        }
                    }


                    TileUtilsDoris.pushDoubleMatrix(BUFFER, targetTile, targetTile.getRectangle());
//                    TileUtilsDoris.pushDoubleMatrix(dataMaster, targetTile, targetTile.getRectangle());

                }

            }

        } catch (Exception e) {
            throw new OperatorException(e);
        }
    }

    public static class Spi extends OperatorSpi {
        public Spi() {
            super(Slant2HeightOp.class);
        }
    }


    private double computeReferencePhase(final double line, final double pixel, final double height,
                                         final SLCImage master, final SLCImage slave,
                                         final Orbit masterOrbit, final Orbit slaveOrbit) throws Exception {

        double mTimeRange = master.pix2tr(pixel);

        // Compute xyz of point P on ELLIPS for this line,pixel
        Point xyzMaster = masterOrbit.lph2xyz(line, pixel, height, master);

        // Compute xyz of slave satelite in orbit_slave from P
        Point timeSlave = slaveOrbit.xyz2t(xyzMaster, slave);

        return mTimeRange * ((-4. * PI * SOL) / master.getRadarWavelength()) - timeSlave.x * ((-4. * PI * SOL) / slave.getRadarWavelength());
    }

    private void normalize(DoubleMatrix data, double min, double max) {
        data.subi(.5 * (min + max));
        data.divi(.25 * (max - min));
    }

}

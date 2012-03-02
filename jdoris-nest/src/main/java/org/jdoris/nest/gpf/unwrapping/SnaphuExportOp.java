package org.jdoris.nest.gpf.unwrapping;

import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.Parameter;
import org.esa.beam.framework.gpf.annotations.SourceProduct;

@OperatorMetadata(alias = "SnaphuExport",
        category = "Unwrapping",
        description = "Export data and prepare conf file for SNAPHU processing", internal = false)
public class SnaphuExportOp extends Operator {

    @SourceProduct
    private Product sourceProduct;

    @Parameter(interval = "TOPO, DEFO, SMOOTH, NOSTATCOSTS",
            description = "Size of coherence estimation window in Azimuth direction",
            defaultValue = "DEFO",
            label = "Statistical-cost mode")
    private String statCostMode = "DEFO";

    @Parameter(interval = "MST, MCF",
            description = "Algorithm used for initialization of wrapped phase values",
            defaultValue = "DEFO",
            label = "Initial method")
    private String initMethod = "MST";

    // no target : I am lazy to setup a writer!??!
    // @TargetProduct
    // private Product targetProduct;

    @Override
    public void initialize() throws OperatorException {

        // as a source use "SNAPHU product" : only float information for phase and coherence!

        // prepare snaphu-parameters class (part of jDoris core)

        // prepare file locations and dirs: see how it is done in envi writer!

        // create snaphu.conf file in corresponding directory


    }


    public static class Spi extends OperatorSpi {
        public Spi() {
            super(SnaphuExportOp.class);
        }
    }

}

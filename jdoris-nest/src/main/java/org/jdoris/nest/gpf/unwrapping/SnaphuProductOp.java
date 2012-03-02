package org.jdoris.nest.gpf.unwrapping;

import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.SourceProduct;

@OperatorMetadata(alias = "SnaphuProduct",
        category = "Unwrapping",
        description = "Create Snpahu product from coherence and interferometric products", internal = false)
public class SnaphuProductOp extends Operator {

    @SourceProduct
    private Product sourceProduct;

    // Perhaps as input only complex vs float data

    @Override
    public void initialize() throws OperatorException {
    }

    public static class Spi extends OperatorSpi {
        public Spi() {
            super(SnaphuExportOp.class);
        }
    }

}

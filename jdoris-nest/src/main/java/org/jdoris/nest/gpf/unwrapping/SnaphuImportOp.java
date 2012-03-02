package org.jdoris.nest.gpf.unwrapping;

import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;

@OperatorMetadata(alias = "SnaphuImport",
        category = "Unwrapping",
        description = "Import result of Snaphu to InSAR product", internal = false)
public class SnaphuImportOp extends Operator {

    @Override
    public void initialize() throws OperatorException {
    }


    public static class Spi extends OperatorSpi {
        public Spi() {
            super(SnaphuImportOp.class);
        }
    }

}

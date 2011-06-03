package org.jdoris.core.filtering;

import org.jdoris.core.SLCImage;

/**
 * User: pmar@ppolabs.com
 * Date: 6/3/11
 * Time: 12:36 PM
 */
public class ProductDataFilter extends SlcDataFilter {

    SLCImage metadata1;

    public SLCImage getMetadata1() {
        return metadata1;
    }

    public void setMetadata1(SLCImage metadata1) {
        this.metadata1 = metadata1;
    }
}

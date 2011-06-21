package org.jdoris.nest.utils;

import org.esa.beam.framework.datamodel.Band;
import org.jdoris.core.Orbit;
import org.jdoris.core.SLCImage;

/**
* User: pmar@ppolabs.com
* Date: 6/20/11
* Time: 11:16 PM
*/
public class CplxContainer {

    String date;
    SLCImage metaData;
    Orbit orbit;
    public Band realBand;
    public Band imagBand;
    int type;

    public CplxContainer(String date, SLCImage metaData, Orbit orbit, Band realBand, Band imagBand) {
        this.date = date;
        this.metaData = metaData;
        this.orbit = orbit;
        this.realBand = realBand;
        this.imagBand = imagBand;
        this.type = realBand.getDataType();
    }

}

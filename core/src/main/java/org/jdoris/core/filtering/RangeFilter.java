package org.jdoris.core.filtering;

import org.jblas.ComplexDoubleMatrix;
import org.jdoris.core.Orbit;
import org.jdoris.core.SLCImage;
import org.jdoris.core.todo_classes.todo_classes;

/**
 * User: pmar@ppolabs.com
 * Date: 4/8/11
 * Time: 5:01 PM
 */
public class RangeFilter {

    //TODO: make template classes for generalInput, operatorInput, and ProductMetadata class
    public static void rangefilter(final todo_classes.inputgeneral input_gen,
                                   final SLCImage master,
                                   final SLCImage slave,
                                   final todo_classes.productinfo interferogram,
                                   final todo_classes.input_filtrange inputfiltrange) {
    }

    public static void rfilterblock(ComplexDoubleMatrix MASTER, // updated
                                    ComplexDoubleMatrix SLAVE,  // updated
                                    long nlmean, float SNRthreshold,
                                    float RSR, // in MHz
                                    float RBW, // in MHz
                                    float hammingalpha, long oversamplefactor, boolean docorrectcorrel,
                                    double meanSNR, // returned
                                    double percentnotfiltered)  { // returned

    }

    // TODO: refactor InputEllips to "Ellipsoid" class of "org.esa.beam.framework.dataop.maptransf.Ellipsoid" and use GeoUtils of NEST;
    public static void rangefilterorbits(final todo_classes.inputgeneral generalinput,
                                         final todo_classes.input_filtrange inputfiltrange,
                                         final todo_classes.input_ell ellips,
                                         final SLCImage master,
                                         final SLCImage slave,
                                         Orbit masterorbit,
                                         Orbit slaveorbit) {
    }
}

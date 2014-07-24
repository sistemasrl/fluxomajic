/*
 *
 *   Copyright (C) 2013 Geobeyond Srl
 *
 *   This library is free software; you can redistribute it and/or modify it under
 *   the terms of the GNU Lesser General Public License as published by the Free
 *   Software Foundation; either version 2.1 of the License, or (at your option)
 *   any later version.
 *
 *   This library is distributed in the hope that it will be useful, but WITHOUT
 *   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *   FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 *   details.
 *
 *   You should have received a copy of the GNU Lesser General Public License along
 *   with this library; if not, write to the Free Software Foundation, Inc., 59
 *   Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *   
 *
 */

package org.geotools.filter.function;

import com.vividsolutions.jts.algorithm.distance.DistanceToPoint;
import com.vividsolutions.jts.algorithm.distance.PointPairDistance;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;

import static org.geotools.filter.capability.FunctionNameImpl.*;

import org.geotools.filter.FunctionExpressionImpl;
import org.geotools.filter.capability.FunctionNameImpl;
import org.opengis.filter.capability.FunctionName;

import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geomgraph.Position;
import com.vividsolutions.jts.operation.buffer.BufferOp;
import com.vividsolutions.jts.operation.buffer.BufferParameters;
import com.vividsolutions.jts.operation.buffer.OffsetCurveBuilder;
import com.vividsolutions.jts.operation.linemerge.LineMerger;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.geotools.geometry.DirectPosition2D;
import org.geotools.geometry.jts.JTS;
import org.geotools.geometry.jts.ReferencedEnvelope;
import org.geotools.referencing.CRS;
import org.geotools.referencing.GeodeticCalculator;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.opengis.filter.expression.Expression;
import org.opengis.filter.expression.Literal;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.TransformException;

/**
 * @see <a
 *      href="http://docs.geotools.org/stable/userguide/tutorial/factory.html">http://docs.geotools.org/stable/userguide/tutorial/factory.html</a>
 */
public class FluxoFilterFunction extends FunctionExpressionImpl implements GeometryTransformation {

    public final static double radiantsOfOneDegree = 0.01745329251994;// PiGreek/180
    public final static double earthRadiusMeters = 6371005.076123;// earth
                                                                  // radius (m)

    private static double mitreLimit = 10.0;

    public static FunctionName NAME = new FunctionNameImpl("fluxo", Geometry.class, parameter("geometry", Geometry.class), parameter("offset", Double.class), parameter("width", Double.class), parameter("driveMode", Integer.class), parameter("quadseg", Integer.class), parameter("endcap", Integer.class), parameter("join", Integer.class), parameter("scalingWidth", Integer.class), parameter("outputCRS", CoordinateReferenceSystem.class), parameter("outputWidth", Integer.class), parameter("outputHeight", Integer.class), parameter("outputBBOX", ReferencedEnvelope.class));

    public FluxoFilterFunction() {

        super(NAME);
    }

    public FluxoFilterFunction(List<Expression> parameters, Literal fallback) {

        super(NAME);
        setParameters(parameters);
        setFallbackValue(fallback);

    }

    @Override
    public Object evaluate(Object feature) {
       
        Geometry geom = null;

        CoordinateReferenceSystem outCRS = null;
        ReferencedEnvelope outBBox = null;
        try {
            geom = getExpression(0).evaluate(feature, Geometry.class);

            Double offsetPx = getExpression(1).evaluate(feature, Double.class);

            Double widthPx = getExpression(2).evaluate(feature, Double.class);

            Integer dMode = getExpression(3).evaluate(feature, Integer.class);

            Integer quadseg = getExpression(4).evaluate(feature, Integer.class);

            Integer endcap = getExpression(5).evaluate(feature, Integer.class);

            Integer join = getExpression(6).evaluate(feature, Integer.class);

            Integer scalingWidth = getExpression(7).evaluate(feature, Integer.class);

            outCRS = getExpression(8).evaluate(feature, CoordinateReferenceSystem.class);

            Integer wmsWidth = getExpression(9).evaluate(feature, Integer.class);

            Integer wmsHeight = getExpression(10).evaluate(feature, Integer.class);

            outBBox = getExpression(11).evaluate(feature, ReferencedEnvelope.class);

            //debug_sdc
//            Geometry tmp_ret=TestGeodetic.getGeodeticLineBuf(geom, widthPx);
//            return tmp_ret;
            return buildGeometryToReturn(outBBox, geom, offsetPx, widthPx, endcap, join, quadseg, dMode, scalingWidth, wmsWidth, wmsHeight);
        }
        catch (Exception ex) {
            //debug_sdc
            if(true){
                Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.INFO, "Output CRS:{0}", "" + outCRS);
                if (outCRS != null && outCRS.getCoordinateSystem() != null && outCRS.getCoordinateSystem().getIdentifiers() != null) {
                    Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.INFO, "SRS identifier of the output coordinate system:{0}", outCRS.getCoordinateSystem().getIdentifiers().toString());
                }
                if (geom != null && geom.getUserData() != null) {
                    Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.INFO, "CRS value read by the geometry of the db:{0}", ((CoordinateReferenceSystem) geom.getUserData()).toString());
                }
                Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.INFO, "Output Boundind Box:{0}", "" + outBBox);
                Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.SEVERE, null, ex);
                
            }
        }
        return geom;

    }

    /**
     * Build the geometry to return from the evaluate filter function
     * 
     * @param outBBox
     * @param geom
     *            Input geometry
     * @param offsetPx
     * @param widthPx
     * @param endcap
     * @param join
     * @param dMode
     * @param scalingWidth
     * @param bbox_width_pixel
     * @param bbox_heigth_pixel
     * @param quadseg
     * @return
     * @throws TransformException
     * @throws FactoryException
     * @throws ExecutionException
     */
    Geometry buildGeometryToReturn(ReferencedEnvelope outBBox, Geometry geom, double offsetPx, double widthPx, int endcap, int join, int quadseg, int dMode, int scalingWidth, int bbox_width_pixel, int bbox_heigth_pixel) throws TransformException, FactoryException, ExecutionException {

//      ReferencedEnvelope tre = null;
//        // Run a CRS transformation if necessary
//        if (outBBox.getCoordinateReferenceSystem().equals(geom.getUserData()))
//            tre = outBBox;
//        else
//            tre = transfEnvelope(outBBox, (CoordinateReferenceSystem) geom.getUserData());

        MathTransform transform = null;
        transform = CRS.findMathTransform((CoordinateReferenceSystem) geom.getUserData(),outBBox.getCoordinateReferenceSystem());
        Geometry newGeom=JTS.transform(geom, transform);
        newGeom.setUserData(outBBox.getCoordinateReferenceSystem());

        
        Double crsUnitPerPixel = null;

        double bbox_width_crs = outBBox.getWidth();
        double bbox_heigth_crs = outBBox.getHeight();
        if (FluxoFunctionFactory.useCache)
            crsUnitPerPixel = FluxoFunctionFactory.cache_width.getIfPresent(bbox_width_crs + ";" + bbox_heigth_crs);
        if (crsUnitPerPixel == null) {
            crsUnitPerPixel = getCRSUnitInOnePixel(outBBox, bbox_width_pixel, bbox_heigth_pixel);
            FluxoFunctionFactory.cache_width.put(bbox_width_crs + ";" + bbox_heigth_crs, crsUnitPerPixel);
        }

        double offsetCrs = offsetPx * crsUnitPerPixel;

        //double widthMt;
        double widthCrs= widthPx * crsUnitPerPixel;
        
//        double offsetCrs=0;
//        double widthCrs =0;
//        if(offsetPx>0){
//            offsetCrs = distanceInCrs(offsetMt, geom, (CoordinateReferenceSystem) geom.getUserData());
//            if(widthPx>0)
//                widthCrs = (offsetCrs / offsetMt) * widthMt;
//            double CRSPerPixel = offsetCrs / offsetPx;
//        }else {
//            if(widthPx>0)
//                widthCrs = distanceInCrs(widthMt, geom, (CoordinateReferenceSystem) geom.getUserData());         
//        }
//        
        
        BufferParameters bufferparams = new BufferParameters();
        bufferparams.setSingleSided(true);

        bufferparams.setEndCapStyle(endcap);
        bufferparams.setJoinStyle(join);
        bufferparams.setQuadrantSegments(quadseg);
        bufferparams.setMitreLimit(mitreLimit);

        Geometry ret = null;
        if (dMode == 1) {
            ret = bufferWithParams(offsetCurve(newGeom, -offsetCrs, bufferparams, false, bufferparams.getQuadrantSegments()), widthCrs, false, bufferparams.getQuadrantSegments(), endcap, join, mitreLimit, scalingWidth, crsUnitPerPixel);
        }
        else {
            ret = bufferWithParams(offsetCurve(newGeom, offsetCrs, bufferparams, false, bufferparams.getQuadrantSegments()), widthCrs, false, bufferparams.getQuadrantSegments(), endcap, join, mitreLimit, scalingWidth, crsUnitPerPixel);
        }

        transform = CRS.findMathTransform((CoordinateReferenceSystem) newGeom.getUserData(),(CoordinateReferenceSystem) geom.getUserData());
        Geometry reverseRet=JTS.transform(ret, transform);
        reverseRet.setUserData( geom.getUserData());

        
        return reverseRet;
    }

    private Geometry offsetCurve(Geometry geometry, double d, BufferParameters parameters, Boolean roughOffsetCurve, Integer qS) throws ExecutionException {

        if(d==0) return (Geometry) geometry.clone();
        
        GeometryFactory gf = geometry.getFactory();
        // If "geometry" is a surface, process its boundary
        if (geometry.getDimension() == 2) {
            geometry = geometry.getBoundary();
        }

        Collection<LineString> offsetCurves = new ArrayList<LineString>();
        if (roughOffsetCurve) {
            addRoughOffsetCurves(offsetCurves, geometry, parameters, d);
        }
        else {
            addCleanOffsetCurves(offsetCurves, geometry, parameters, d, qS);
        }
        return gf.buildGeometry(offsetCurves);
    }

    private void addCleanOffsetCurves(Collection<LineString> offsetCurves, Geometry sourceCurve, BufferParameters parameters, Double offsetDistance, Integer qS) throws ExecutionException {

        parameters.setSingleSided(true);
        parameters.setQuadrantSegments(qS);

        Geometry sidedBuffer = BufferOp.bufferOp(sourceCurve, offsetDistance, parameters).getBoundary();

        Collection<LineString> offsetSegments = new ArrayList<LineString>();
        double offsetPositiveDistance = Math.abs(offsetDistance);
        // Segments located entirely under this distance are excluded
        double lowerBound = offsetPositiveDistance * Math.sin(Math.PI / (4 * qS));
        // Segments located entirely over this distance are included
        // note that the theoretical approximation made with quadrantSegments
        // is offset*cos(PI/(4*quadrantSegments) but
        // offset*cos(PI/(2*quadrantSegments)
        // is used to make sure to include segments located on the boundary
        double upperBound = offsetPositiveDistance * Math.cos(Math.PI / (2 * qS));
        for (int i = 0; i < sidedBuffer.getNumGeometries(); i++) {
            Coordinate[] cc = sidedBuffer.getGeometryN(i).getCoordinates();
            PointPairDistance ppd = new PointPairDistance();
            DistanceToPoint.computeDistance(sourceCurve, cc[0], ppd);

            double dj = ppd.getDistance();
            for (int j = 1; j < cc.length; j++) {
                double di = dj;
                ppd = new PointPairDistance();
                DistanceToPoint.computeDistance(sourceCurve, cc[j], ppd);
                dj = ppd.getDistance();
                // segment along or touching the source geometry : exclude it
                if (Math.max(di, dj) < lowerBound || di == 0 || dj == 0) {
                    continue;
                } // segment along the buffer boundary : include it
                else if (Math.min(di, dj) > upperBound) {
                    LineString segment = sourceCurve.getFactory().createLineString(new Coordinate[] { cc[j - 1], cc[j] });
                    offsetSegments.add(segment);
                } // segment entirely located inside the buffer : exclude it
                else if (Math.min(di, dj) > lowerBound && Math.max(di, dj) < upperBound) {
                    continue;
                } // segment with a end at the offset distance and the other
                  // located within the buffer : divide it
                else {
                    // One of the coordinates is closed to but not on the source
                    // curve and the other is more or less closed to offset
                    // distance
                    divide(offsetSegments, sourceCurve, cc[j - 1], cc[j], di, dj, lowerBound, upperBound);
                }
            }
        }

        offsetCurves.addAll(merge(offsetSegments));
    }

    /**
     * Recursive function to split segments located on the single-side buffer
     * boundary, but having a part of them inside the full buffer.
     */
    private void divide(Collection<LineString> offsetSegments, Geometry sourceCurve, Coordinate c1, Coordinate c2, double d1, double d2, double lb, double ub) {

        // I stop recursion for segment < 2*lb to exclude small segments
        // perpendicular but very close to the boundary
        if (c1.distance(c2) < 2 * lb) {
            return;
        }

        Coordinate c = new Coordinate((c1.x + c2.x) / 2.0, (c1.y + c2.y) / 2.0);
        PointPairDistance ppd = new PointPairDistance();
        DistanceToPoint.computeDistance(sourceCurve, c, ppd);
        double d = ppd.getDistance();
        if (Math.max(d1, d) < lb) {}
        else if (Math.min(d1, d) > lb && Math.max(d1, d) < ub) {}
        else if (Math.min(d1, d) > ub) {
            LineString segment = sourceCurve.getFactory().createLineString(new Coordinate[] { c1, c });
            offsetSegments.add(segment);
        }
        else {
            divide(offsetSegments, sourceCurve, c1, c, d1, d, lb, ub);
        }
        if (Math.max(d, d2) < lb) {}
        else if (Math.min(d, d2) > lb && Math.max(d, d2) < ub) {}
        else if (Math.min(d, d2) > ub) {
            LineString segment = sourceCurve.getFactory().createLineString(new Coordinate[] { c, c2 });
            offsetSegments.add(segment);
        }
        else {
            divide(offsetSegments, sourceCurve, c, c2, d, d2, lb, ub);
        }
    }

    private void addRoughOffsetCurves(Collection<LineString> offsetCurves, Geometry sourceCurve, BufferParameters parameters, Double offsetDistance) {

        OffsetCurveBuilder builder = new OffsetCurveBuilder(sourceCurve.getFactory().getPrecisionModel(), parameters);

        for (int i = 0; i < sourceCurve.getNumGeometries(); i++) {
            if (sourceCurve.getGeometryN(i) instanceof LineString) {
                LineString lineString = (LineString) sourceCurve.getGeometryN(i);
                Coordinate[] cc = lineString.getCoordinates();
                if (lineString.isClosed()) {
                    offsetCurves.add(lineString.getFactory().createLineString(builder.getRingCurve(cc, offsetDistance > 0 ? Position.LEFT : Position.RIGHT, Math.abs(offsetDistance))));
                }
                else {
                    offsetCurves.add(lineString.getFactory().createLineString(builder.getOffsetCurve(cc, offsetDistance)));
                }
            }
        }
    }

    @SuppressWarnings("unchecked")
    private Collection<LineString> merge(Collection<LineString> linestrings) {

        LineMerger merger = new LineMerger();
        merger.add(linestrings);
        return merger.getMergedLineStrings();
    }

    /**
     * Returns a buffered geometry with old shapes in the center of new ones. If
     * the buffer is issued at single side then a negative offset renders the
     * shape on the left while a positive offset on the right
     */
    public static Geometry bufferWithParams(Geometry geometry, Double offset, boolean singleSided, Integer quadrantSegments, Integer capStyle, Integer joinStyle, Double mitreLimit, Integer scalingWidth, double meterPerPixel) {

        double widthCRS = 0;
        if (offset != null) {
            if (scalingWidth == 1) {
                // TODO

            }
            else {
                widthCRS = offset;
            }
        }

        BufferParameters bufferparameters = new BufferParameters();

        // Custom code to be able to draw only on the side of the offset curve
        bufferparameters.setSingleSided(singleSided);

        if (quadrantSegments != null) {
            bufferparameters.setQuadrantSegments(quadrantSegments.intValue());
        }
        if (capStyle != null) {
            bufferparameters.setEndCapStyle(capStyle.intValue());
        }
        if (joinStyle != null) {
            bufferparameters.setJoinStyle(joinStyle.intValue());
        }
        if (mitreLimit != null) {
            bufferparameters.setMitreLimit(mitreLimit.doubleValue());
        }
        
        return BufferOp.bufferOp(geometry, widthCRS, bufferparameters);
    }

    /**
     * Returns an translated rendering envelope if the offsets are not using
     * feature attributes. If the offsets are feature dependent the user will
     * have to expand the rendering area via the renderer buffer parameter
     */
    @Override
    public ReferencedEnvelope invert(ReferencedEnvelope renderingEnvelope) {

        Envelope bufferedEnvelope = JTS.toGeometry((Envelope) renderingEnvelope).getEnvelopeInternal();
        return new ReferencedEnvelope(bufferedEnvelope, renderingEnvelope.getCoordinateReferenceSystem());
    }

    private double distanceInCrs(double inMeters, Geometry geom, CoordinateReferenceSystem crs) {

        try {
            double[] refXY = { geom.getCentroid().getCoordinates()[0].x, geom.getCentroid().getCoordinates()[0].y };
            return distanceInCrs(inMeters, refXY, crs);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return 0d;
    }

    private double distanceInCrs(double distanceInMeters, double[] refXY, CoordinateReferenceSystem crs) throws TransformException {

        double dist = 0;

        // calculate the distance in meters of 0.01 * refY in the ref CRS
        double[] sp = { refXY[0], refXY[1] };
        double[] dp = { refXY[0], refXY[1] * 1.01 };

        GeodeticCalculator gc = new GeodeticCalculator(crs);

        gc.setStartingPosition(new DirectPosition2D(crs, sp[0], sp[1]));
        gc.setDestinationPosition(new DirectPosition2D(crs, dp[0], dp[1]));

        double refY01InMeters = gc.getOrthodromicDistance();

        // now, calculate the CRS distance as a proportional of 0.01 * refY
        dist = distanceInMeters * (refXY[1] * 0.01) / refY01InMeters;

        return dist;
    }

    /**
     * 
     * @param bbox_crs
     * @param bbox_width_pixel
     * @param bbox_heigth_pixel
     * @return
     */
    private double getCRSUnitInOnePixel(ReferencedEnvelope bbox_crs, int bbox_width_pixel, int bbox_heigth_pixel) {

        double diag_distance_pixel = Math.sqrt((bbox_width_pixel * bbox_width_pixel) + (bbox_heigth_pixel * bbox_heigth_pixel));
        //double diag_distance_m = getGeodeticSegmentLength(bbox_crs.getMinX(), bbox_crs.getMinY(), bbox_crs.getMaxX(), bbox_crs.getMaxY());
        double diag_distance_CRS = Math.sqrt((bbox_crs.getWidth() * bbox_crs.getWidth()) + (bbox_crs.getHeight() * bbox_crs.getHeight()));

        return diag_distance_CRS / diag_distance_pixel;
    }

    /**
     * Compute the distance in meters between the two specified points
     * 
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @return
     */
    private static double getDistance(double x1, double y1, double x2, double y2) {

        double xFactor = Math.cos(y1) * radiantsOfOneDegree;
        double leng = Math.sqrt((x1 - x2) * xFactor * (x1 - x2) * xFactor + (y1 - y2) * (y1 - y2));
        leng = leng * radiantsOfOneDegree * earthRadiusMeters;
        return leng;
    }

    private static double getGeodeticSegmentLength(double minx, double miny, double maxx, double maxy) {

        final GeodeticCalculator calculator = new GeodeticCalculator(DefaultGeographicCRS.WGS84);
        double rminx = rollLongitude(minx);
        double rminy = rollLatitude(miny);
        double rmaxx = rollLongitude(maxx);
        double rmaxy = rollLatitude(maxy);
        calculator.setStartingGeographicPoint(rminx, rminy);
        calculator.setDestinationGeographicPoint(rmaxx, rmaxy);
        return calculator.getOrthodromicDistance();
    }

    protected static double rollLongitude(final double x) {

        double rolled = x - (((int) (x + Math.signum(x) * 180)) / 360) * 360.0;
        return rolled;
    }

    protected static double rollLatitude(final double x) {

        double rolled = x - (((int) (x + Math.signum(x) * 90)) / 180) * 180.0;
        return rolled;
    }

    private ReferencedEnvelope transfEnvelope(ReferencedEnvelope re, CoordinateReferenceSystem targetCRS) throws TransformException, FactoryException {

        ReferencedEnvelope result = null;
        result = re.transform(targetCRS, true, 10);
        return result;
    }

}

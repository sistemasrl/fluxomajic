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

import org.geotools.geometry.jts.JTS;
import org.geotools.geometry.jts.ReferencedEnvelope;
import org.geotools.referencing.CRS;
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

    public static FunctionName NAME = new FunctionNameImpl("fluxo", Geometry.class, parameter("geometry", Geometry.class), parameter("offset", Double.class), parameter("width", Double.class), parameter("driveMode", Integer.class), parameter("quadseg", Integer.class), parameter("endcap", Integer.class), parameter("join", Integer.class), parameter("scalingWidth", Integer.class), parameter("maxPixelLengthToDraw", Integer.class), parameter("outputCRS", CoordinateReferenceSystem.class), parameter("outputWidth", Integer.class), parameter("outputHeight", Integer.class), parameter("outputBBOX", ReferencedEnvelope.class));

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

        @SuppressWarnings("unused")
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

            int maxPixelLengthToDraw = getExpression(8).evaluate(feature,Integer.class);

            outCRS = getExpression(9).evaluate(feature, CoordinateReferenceSystem.class);

            Integer wmsWidth = getExpression(10).evaluate(feature, Integer.class);

            Integer wmsHeight = getExpression(11).evaluate(feature, Integer.class);

            outBBox = getExpression(12).evaluate(feature, ReferencedEnvelope.class);

            BufferParameters bufferParameters = new BufferParameters();
            bufferParameters.setSingleSided(true);
            bufferParameters.setEndCapStyle(endcap);
            bufferParameters.setJoinStyle(join);
            bufferParameters.setQuadrantSegments(quadseg);
            bufferParameters.setMitreLimit(mitreLimit);

            // Update offset and width in order to respect the input values.
            // Indeed to drawn a correct feature is necessary to run a
            // both-sided buffer
            // that doubled the current width
            widthPx = widthPx / 2;
            offsetPx += widthPx;

            return buildGeometryToReturn(outBBox, geom, offsetPx, widthPx, bufferParameters, maxPixelLengthToDraw, dMode, scalingWidth, wmsWidth, wmsHeight);
        }
        catch (Exception ex) {
            Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.SEVERE, null, ex);
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
    Geometry buildGeometryToReturn(ReferencedEnvelope outBBox, Geometry geom, double offsetPx, double widthPx, BufferParameters bufferParameters, int maxPixelLengthToDraw, int dMode, int scalingWidth, int bbox_width_pixel, int bbox_heigth_pixel) throws TransformException, FactoryException, ExecutionException {

        // Compute CRS units in one pixel
        Double crsUnitPerPixel = null;

        double bbox_width_crs = outBBox.getWidth();
        double bbox_heigth_crs = outBBox.getHeight();
        if (FluxoFunctionFactory.useCache)
            crsUnitPerPixel = FluxoFunctionFactory.cache_width.getIfPresent(bbox_width_crs + ";" + bbox_heigth_crs);
        if (crsUnitPerPixel == null) {
            crsUnitPerPixel = getCRSUnitInOnePixel(outBBox, bbox_width_pixel, bbox_heigth_pixel);
            FluxoFunctionFactory.cache_width.put(bbox_width_crs + ";" + bbox_heigth_crs, crsUnitPerPixel);
        }

        // Transorf the geometry according to the output CRS
        Geometry geom_transformed = geom;
        boolean runTransformation = false;
        MathTransform transform = null;
        if (!((CoordinateReferenceSystem) geom.getUserData()).equals(outBBox.getCoordinateReferenceSystem())) {
            runTransformation = true;
            transform = CRS.findMathTransform((CoordinateReferenceSystem) geom.getUserData(), outBBox.getCoordinateReferenceSystem());
            geom_transformed = JTS.transform(geom, transform);
            geom_transformed.setUserData(outBBox.getCoordinateReferenceSystem());
        }

        // Simplify geometry
        Geometry geom_simplified = geom_transformed;
        if(maxPixelLengthToDraw>0){
            double xcoo1, ycoo1, xcoo2, ycoo2;
            double geom_bbox_width, geom_bbox_height;
            if (geom_transformed.getCoordinates().length > 2) {
                Envelope envGeom = geom_transformed.getEnvelopeInternal();
                xcoo1 = envGeom.getMinX();
                ycoo1 = envGeom.getMinY();
                xcoo2 = envGeom.getMaxX();
                ycoo2 = envGeom.getMaxY();
                geom_bbox_width = envGeom.getWidth();
                geom_bbox_height = envGeom.getHeight();
            }
            else {
                Coordinate[] coord = geom_transformed.getCoordinates();
                xcoo1 = coord[0].x;
                ycoo1 = coord[0].y;
                xcoo2 = coord[1].x;
                ycoo2 = coord[1].y;
                geom_bbox_width = Math.abs(xcoo1 - xcoo2);
                geom_bbox_height = Math.abs(ycoo1 - ycoo2);
            }
            double diag_distance_crs = Math.sqrt(geom_bbox_width * geom_bbox_width + geom_bbox_height * geom_bbox_height);
            double diag_distance_pixel = diag_distance_crs / crsUnitPerPixel;
            // check if diagonal length of bounding box of the geometry is minor
            // than 3 pixel
            if (diag_distance_pixel < maxPixelLengthToDraw) {// simplify geometry...
                return null;
                // geom_simplified=geom.getFactory().createPoint(new
                // Coordinate(xcoo1, ycoo1));
                // geom_simplified=geom.getFactory().createLineString(new
                // Coordinate[]{new Coordinate(xcoo1, ycoo1), new Coordinate(xcoo2,
                // ycoo2)});
            }
            
        }

        if (scalingWidth == 1) {
            double unitPerPixel=crsUnitPerPixel;
            if(outBBox.getCoordinateReferenceSystem().getCoordinateSystem().getAxis(0).getUnit().getStandardUnit().toString().equals("rad")){
                unitPerPixel=crsUnitPerPixel*earthRadiusMeters;
            }
            double multiplier=Math.pow(10/unitPerPixel,1.0/12);
            widthPx=widthPx*multiplier;
        }
        
        double offsetCrs = offsetPx * crsUnitPerPixel;
        double widthCrs = widthPx * crsUnitPerPixel;

        bufferParameters.setSingleSided(true);

        Geometry geom_offseted = null;
        // Offset
        if (dMode == 0) {
            geom_offseted = offsetCurve(geom_simplified, -offsetCrs, bufferParameters, false);//right side
        }
        else {
            geom_offseted = offsetCurve(geom_simplified, offsetCrs, bufferParameters, false);//left side
        }

        // Buffer
        bufferParameters.setSingleSided(false);
        Geometry geom_buffered = BufferOp.bufferOp(geom_offseted, widthCrs, bufferParameters);

        // Reverse the geometry trasformation returnig to the original CRS
        Geometry geom_toReturn = geom_buffered;
        if (runTransformation) {
            transform = CRS.findMathTransform(outBBox.getCoordinateReferenceSystem(), (CoordinateReferenceSystem) geom.getUserData());
            geom_toReturn = JTS.transform(geom_buffered, transform);
            geom_toReturn.setUserData(geom.getUserData());
        }

        return geom_toReturn;
    }

    private Geometry offsetCurve(Geometry geometry, double d, BufferParameters parameters, Boolean roughOffsetCurve) throws ExecutionException {

        if (d == 0)
            return (Geometry) geometry.clone();

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
            addCleanOffsetCurves(offsetCurves, geometry, parameters, d);
        }
        return gf.buildGeometry(offsetCurves);
    }

    private void addCleanOffsetCurves(Collection<LineString> offsetCurves, Geometry sourceCurve, BufferParameters parameters, Double offsetDistance) throws ExecutionException {

        Geometry sidedBuffer = BufferOp.bufferOp(sourceCurve, offsetDistance, parameters).getBoundary();

        Collection<LineString> offsetSegments = new ArrayList<LineString>();
        double offsetPositiveDistance = Math.abs(offsetDistance);
        // Segments located entirely under this distance are excluded
        double lowerBound = offsetPositiveDistance * Math.sin(Math.PI / (4 * parameters.getQuadrantSegments()));
        // Segments located entirely over this distance are included
        // note that the theoretical approximation made with quadrantSegments
        // is offset*cos(PI/(4*quadrantSegments) but
        // offset*cos(PI/(2*quadrantSegments)
        // is used to make sure to include segments located on the boundary
        double upperBound = offsetPositiveDistance * Math.cos(Math.PI / (2 * parameters.getQuadrantSegments()));
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
     * Returns an translated rendering envelope if the offsets are not using
     * feature attributes. If the offsets are feature dependent the user will
     * have to expand the rendering area via the renderer buffer parameter
     */
    @Override
    public ReferencedEnvelope invert(ReferencedEnvelope renderingEnvelope) {

        Envelope bufferedEnvelope = JTS.toGeometry((Envelope) renderingEnvelope).getEnvelopeInternal();
        return new ReferencedEnvelope(bufferedEnvelope, renderingEnvelope.getCoordinateReferenceSystem());
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
        double diag_distance_CRS = Math.sqrt((bbox_crs.getWidth() * bbox_crs.getWidth()) + (bbox_crs.getHeight() * bbox_crs.getHeight()));
        return diag_distance_CRS / diag_distance_pixel;
    }
}

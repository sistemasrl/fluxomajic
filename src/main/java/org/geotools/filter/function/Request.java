/**
 * CacheFeature.java
 */
package org.geotools.filter.function;

import org.geotools.geometry.jts.ReferencedEnvelope;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

import com.vividsolutions.jts.geom.Geometry;

/**
 * @author simone.decristofaro
 * 18/lug/2014
 */
public class Request{
	
	//ATTIBUTES
    private FluxoFilterFunction fluxoFilterFunction;
	private Geometry geom;
	private double offsetPx;
	private double widthPx;
	private int dMode;
	private int quadseg;
	private int endcap;
	private int join;
	private int scalingWidth;
	private CoordinateReferenceSystem outCRS;
	private int wmsWidth;
	private int wmsHeight;
	private ReferencedEnvelope outBBox;
	
	
	//CONSTRUCTORS
	/**
	 * @param geom
	 * @param offsetPx
	 * @param widthPx
	 * @param dMode
	 * @param quadseg
	 * @param endcap
	 * @param join
	 * @param scalingWidth
	 * @param outCRS
	 * @param wmsWidth
	 * @param wmsHeight
	 * @param outBBox
	 */
	public Request(FluxoFilterFunction fluxoFilterFunction, Geometry geom, double offsetPx, double widthPx, int dMode,
			int quadseg, int endcap, int join, int scalingWidth,
			CoordinateReferenceSystem outCRS, int wmsWidth, int wmsHeight,
			ReferencedEnvelope outBBox) {
	    this.fluxoFilterFunction=fluxoFilterFunction;
		this.geom = geom;
		this.offsetPx = offsetPx;
		this.widthPx = widthPx;
		this.dMode = dMode;
		this.quadseg = quadseg;
		this.endcap = endcap;
		this.join = join;
		this.scalingWidth = scalingWidth;
		this.outCRS = outCRS;
		this.wmsWidth = wmsWidth;
		this.wmsHeight = wmsHeight;
		this.outBBox = outBBox;
	}
	
	//METHODS
	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		if(obj==null || !obj.getClass().equals(Request.class)) return false;
		Request request=(Request) obj;
		
		if(!geom.equals(request.geom)) return false;
		if(!outCRS.equals(request.outCRS)) return false;
		if(!outBBox.equals(request.outBBox)) return false;
		
		return offsetPx==request.offsetPx && widthPx==request.widthPx
				&& dMode==request.dMode && quadseg==request.quadseg
				&& endcap==request.endcap && join==request.join
				&& scalingWidth==request.scalingWidth
				&& wmsWidth==request.wmsWidth && wmsHeight==request.wmsHeight;
	}

	/* (non-Javadoc)
     * @see java.lang.Object#clone()
     */
    @Override
    protected Object clone() throws CloneNotSupportedException {

        Request r=new Request(fluxoFilterFunction,(Geometry) geom.clone(), offsetPx, widthPx, dMode, quadseg, endcap, join, scalingWidth, outCRS, wmsWidth, wmsHeight, outBBox);
        return r;
    }

    @Override
	public int hashCode() {
        return 0;
	};
        
	/**
     * @return the fluxoFilterFunction
     */
    public FluxoFilterFunction getFluxoFilterFunction() {
    
        return fluxoFilterFunction;
    }

    
    /**
     * @return the geom
     */
    public Geometry getGeom() {
    
        return geom;
    }

    
    /**
     * @return the offsetPx
     */
    public double getOffsetPx() {
    
        return offsetPx;
    }

    
    /**
     * @return the widthPx
     */
    public double getWidthPx() {
    
        return widthPx;
    }

    
    /**
     * @return the dMode
     */
    public int getdMode() {
    
        return dMode;
    }

    
    /**
     * @return the quadseg
     */
    public int getQuadseg() {
    
        return quadseg;
    }

    
    /**
     * @return the endcap
     */
    public int getEndcap() {
    
        return endcap;
    }

    
    /**
     * @return the join
     */
    public int getJoin() {
    
        return join;
    }

    
    /**
     * @return the scalingWidth
     */
    public int getScalingWidth() {
    
        return scalingWidth;
    }

    
    /**
     * @return the outCRS
     */
    public CoordinateReferenceSystem getOutCRS() {
    
        return outCRS;
    }

    
    /**
     * @return the wmsWidth
     */
    public int getWmsWidth() {
    
        return wmsWidth;
    }

    
    /**
     * @return the wmsHeight
     */
    public int getWmsHeight() {
    
        return wmsHeight;
    }

    
    /**
     * @return the outBBox
     */
    public ReferencedEnvelope getOutBBox() {
    
        return outBBox;
    }

}

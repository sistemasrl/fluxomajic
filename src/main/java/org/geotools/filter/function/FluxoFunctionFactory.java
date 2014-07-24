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

import org.geotools.filter.function.FluxoFilterFunction;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;

import org.geotools.feature.NameImpl;
import org.geotools.filter.FunctionFactory;
import org.opengis.feature.type.Name;
import org.opengis.filter.capability.FunctionName;
import org.opengis.filter.expression.Expression;
import org.opengis.filter.expression.Function;
import org.opengis.filter.expression.Literal;

import com.google.common.cache.Cache;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.io.WKBReader;
import com.vividsolutions.jts.io.WKBWriter;
import com.vividsolutions.jts.io.WKTWriter;
import com.vividsolutions.jts.operation.buffer.BufferOp;

public class FluxoFunctionFactory implements FunctionFactory {

    static Cache<String, Double> cache_width;
    //debug_sdc
    static boolean useCache=false;

    
    public List<FunctionName> getFunctionNames() {

        List<FunctionName> functionList = new ArrayList<FunctionName>();
        functionList.add(FluxoFilterFunction.NAME);
        // initialize cache
        cache_width = CacheBuilder.newBuilder().maximumSize(1000000).build();

        functionList.add(FluxoFilterFunctionOld.NAME);
        functionList.add(FluxoFilterFunction_noGeoc.NAME);
        return Collections.unmodifiableList(functionList);
    }

    public Function function(String name, List<Expression> args, Literal fallback) {

        return function(new NameImpl(name), args, fallback);
    }

    public Function function(Name name, List<Expression> args, Literal fallback) {

        if (FluxoFilterFunction.NAME.getFunctionName().equals(name)) {
            return new FluxoFilterFunction(args, fallback);
        }
        else if (FluxoFilterFunctionOld.NAME.getFunctionName().equals(name)) {
            return new FluxoFilterFunctionOld(args, fallback);
        }
        else if (FluxoFilterFunction_noGeoc.NAME.getFunctionName().equals(name)) {
            return new FluxoFilterFunction_noGeoc(args, fallback);
        }
        return null; // we do not implement that function
    }
}

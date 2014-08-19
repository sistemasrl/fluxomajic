Fluxomajic
=======

Fluxomajic is a GeoServer plugin consisting in a new pluggable SLD geometry transformation. It renders the state of flows related to linestrings and their associated information (i.e.: traffic state).

It takes the following input parameters:

- Linestring collection
- Offset
- Width
- Drive side
- Number of quadrants
- ENDCAP style
- JOIN style
- Scaling Widht (YES or NO)

The output is a new polygonal geometry that follows the shape of the source linestrings; directly in your WMS request:

![fluxomajic shaped traffic layer](https://raw.github.com/geobeyond/fluxomajic/master/img/fluxomajic.jpg "fluxomajic behavior")


## Versions

Fluxomajic version number is 3.x-geoserver.a.b.c

- 3 is the main version (Fluxomajic comes from a long path...)
- x is the minor version (bugs and minor improvement)
- geoserver.a.b.c is the supported version of Geoserver


## How to build

+ Clone Fluxomajic repository.

+ Make sure each dependency in the POM is aligned with the Geotools version used in your GeoServer. As an example version 9.0 for GeoServer 2.3:

```xml
    ...
	<properties>
		<project.build.sourceEncoding>UTF-8</project.build.sourceEncoding>
		<project.reporting.outputEncoding>UTF-8</project.reporting.outputEncoding>
		<geoserver-version>9.0</geoserver-version>
	</properties>
    ...
	<dependency>
		<groupId>org.geotools</groupId>
		<artifactId>gt-xxx</artifactId>
		<version>${geoserver-version}</version>
	</dependency>
```

See also [GeoTools and GeoServer release schedule](http://geoserver.org/display/GEOS/GeoTools+and+GeoServer+release+schedule).

+ Use Maven to build

```bash
	mvn clean install
```


## How to install

+ Copy the generated JAR file into the WEB-INF/lib folder of Geoserver.


## How to use

Fluxomajic adds a function "fluxo" that can be used in SLD. Here is an example of how to use it.

```xml
<PolygonSymbolizer uom="http://www.opengeospatial.org/se/units/pixel">
    <Geometry>
        <ogc:Function name="fluxo">
            <ogc:PropertyName>shap</ogc:PropertyName>  <!-- Layer's attribute name -->
            <ogc:Literal>5</ogc:Literal>              <!-- Offset; expressed in pixel -->
            <ogc:Literal>5</ogc:Literal>              <!-- Width; expressed in pixel -->
            <ogc:Literal>0</ogc:Literal>              <!-- Drive side;   0 = RIGHT(default),
                                                                         1 = LEFT -->
            <ogc:Literal>4</ogc:Literal>              <!-- Number of quandrants (number of facets
                                                           into which to divide a fillet of
                                                           90 degrees);  16(default) -->
            <ogc:Literal>3</ogc:Literal>              <!-- ENDCAP style; 1 = ROUND(default),
                                                                         2 = FLAT,
                                                                         3 = SQUARE -->
            <ogc:Literal>1</ogc:Literal>              <!-- JOIN style;   1 = ROUND(default),
                                                                         2 = MITRE,
                                                                         3 = BEVEL -->
            <ogc:Literal>2</ogc:Literal>              <!-- Scaling Width (if the width has to scale
                                                           according to zoom level);
                                                                         1 = YES(default),
                                                                         2 = NO -->
            <ogc:Literal>3</ogc:Literal>            <!-- Minimum length (in pixel) of the diagonal 
                                                      of the bounding box of the single geometry
                                                      to run the drawn process;  3(default) -->
            <ogc:Function name="env">                 <!-- envs (fixed and mandatory) -->
                <ogc:Literal>wms_crs</ogc:Literal>
            </ogc:Function>
            <ogc:Function name="env">
                <ogc:Literal>wms_width</ogc:Literal>
            </ogc:Function>
            <ogc:Function name="env">
                <ogc:Literal>wms_height</ogc:Literal>
            </ogc:Function>
            <ogc:Function name="env">
                <ogc:Literal>wms_bbox</ogc:Literal>
            </ogc:Function>
        </ogc:Function>
    </Geometry>
    <Fill>
        <CssParameter name="fill">
            <ogc:Literal>#33BA2E</ogc:Literal>
        </CssParameter>
    </Fill>
    </PolygonSymbolizer>
```


## How to setup development environment

+ Clone Geoserver repository and setup following the [GeoServer Developer Manual](http://docs.geoserver.org/stable/en/developer/).

+ Clone Fluxomajic repository into the community folder of your Geoserver installation.

+ Add the following profile to the community/pom.xml.

```xml
    <profile>
        <id>fluxo</id>
        <modules>
            <module>"fluxo-project-folder-name"</module>
        </modules>
    </profile>
```

+ Add the following dependency to the web/app/pom.xml.

```xml
    <dependency>
        <groupId>org.geotools</groupId>
        <artifactId>fluxo</artifactId>
        <version>"current-version"</version>
    </dependency>
```

+ Build using Maven both fluxo and gs-web-app modules.

+ As a confirmation, if you are using Eclipse, check in the gs-web-app properties whether fluxo is among the dependencies.

## Geotools buffer behavior

Since the geotools buffer function not consider end cap in case of single side buffer, is necessary to run a buffer on both shape size.
But in this way the width is double (respect to the input parameter), so is necessary to internally correct the two parameters offset and width in the following way:
+ Or = Oi + 0.5Wi
+ Wr = Wi/2
where:
+ Oi = input offset
+ Wi = input width
+ Or = offset to use for the geotools buffer function
+ Wr = width to use for the geotools buffer function  
![fluxomajic workflow](/img/workflow.png?raw=true "getools buffer behavior")

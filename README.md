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
    <sld:PolygonSymbolizer uom="http://www.opengeospatial.org/se/units/pixel">
    <sld:Geometry>
      <ogc:Function name="fluxo">
        <ogc:PropertyName>xxx</ogc:PropertyName>     <!-- Layer's attribute name -->
        <ogc:Literal>5</ogc:Literal>                 <!-- Offset expressed in pixel -->                  
        <ogc:Literal>5</ogc:Literal>                 <!-- Width expressed in pixel -->                 
        <ogc:Literal>1</ogc:Literal>                 <!-- Drive side; 0 = RIGHT(default), 1 = LEFT -->
        <ogc:Literal>4</ogc:Literal>                 <!-- Number of quandrants; default = 16 -->
        <ogc:Literal>3</ogc:Literal>                 <!-- ENDCAP style; 1 = ROUND(default), 2 = FLAT, 3 = SQUARE -->
        <ogc:Literal>1</ogc:Literal>                 <!-- JOIN style; 1 = ROUND(default), 2 = MITRE, 3 = BEVEL -->
        <ogc:Function name="env">                    <!-- envs (fixed and mandatory) -->
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
    </sld:Geometry>
    <sld:Fill>
      <sld:CssParameter name="fill">
          <ogc:Literal>#33BA2E</ogc:Literal>
      </sld:CssParameter>
    </sld:Fill>
    </sld:PolygonSymbolizer>
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
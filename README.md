# Many-to-Many Polygon Matching à la Jaccard

The code in this repository implements the algorithm proposed in the paper 'Many-to-Many Polygon Matching à la Jaccard'.

## Dependencies

The code makes use of the following libraries:

* [boost](https://www.boost.org/)
* [CGAL](https://www.cgal.org/)
* [shapelib](https://github.com/OSGeo/shapelib)
* [gurobi](https://www.gurobi.com/)

## Installation

1. Clone https://github.com/anaumannigg/m-n-Polygon-Matching-Jaccard.git
2. Clone https://github.com/OSGeo/shapelib.git into the subdirectory 'lib'.
3. Put the FindGUROBI.cmake file into the project folder
4. make project:

 ```bash
mkdir build
cd build
cmake ..
make -j 8
```

5. Run via

```bash
 ./build/PolygonMatching
```

## Running the Code
The code can be parameterized with the following command line arguments:

* -d dataset_name   | Provide the name of the dataset.
* -l lambda         | Specify the parameter for lambda.
* -t threads        | (optional, default 1) Specify the number of threads that will be used for computations on connected components. (We recommend 8)
* -k time_limit     | (optional, default is 1000) Specify the time limit in [s] for the computation per connected component.
* -s size_limit     | (optional, deactivated on default) Set a limit on the set size on cumulative vertex computation.
* -h, --help        | Display help message.

The program expects two Shapefiles containing single Polygons without holes, which are located in 'input/data_name'. The files should be named in the way 'data_name_osm' and 'data_name_atkis'.
Should a polygon have holes, only its outer boundary will be considered.
The code will store its results in 'export/data_name/lambdavalue*10' in the form of Shapefiles containing the original polygon sets per input file, where each polygon has the attribute 'match_id', which encodes the match,
that the polygon has been assigned to. Negative match IDs are used for polygons, which have not been matched. Additionally, two csv-Files are exported into the same directory. 

* 'analysis' - provides the sizes of the components after the precomputation steps as well as information about timing, and if the optimal solution could be computed within the given time limit. 
* 'data' provides the information about the optimal solution, where for each polygon, its match ID as well as the quality of the match is listed.

Within this table 'map' = 0 refers to the osm data and 'map' = 1 to atkis.

This repository contains one exemplary dataset of Auerberg, a district of Bonn. It consists of the [OpenStreeMap](https://www.openstreetmap.org/)-Data, fetched from the Overpass-API and the [official cadastral data](https://www.opengeodata.nrw.de/produkte/geobasis/lk/akt/gru_xml/), publicly available as a free download. (both data sets were fetched on November 9th, 2023)

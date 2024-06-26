# Many-to-Many Polygon Matching à la Jaccard

The code in this repository implements the algorithm proposed in the paper 'Many-to-Many Polygon Matching à la Jaccard'.

## Installation

1. Clone https://github.com/anaumannigg/m-n-Polygon-Matching-Jaccard.git
2. Clone https://github.com/OSGeo/shapelib.git into the subdirectory 'extern'.
3. make project:

 ```bash
mkdir build
cd build
cmake ..
make
```

3. Run via

```bash
 ./build/main
```

## Running the Code

The program expects two Shapefiles containing Polygons, which are located in 'input/data_name'. The files should be named in the way 'data_name_osm' and 'data_name_alkis'.
The code will store its results in 'export/data_name/lambdavalue*10' in the form of Shapefiles containing the original polygon sets per input file, where each polygon has the attribute 'match_id', which encodes the match,
that the polygon has been assigned to. Negative match IDs are used for polygons, which have not been matched. Additionally, two csv-Files are exported into the same directory. 'analysis' provides the sizes of the components
after the precomputation steps as well as information about timing, and if the optimal solution could be computed within the given time limit. 'data' provides the information about the optimal solution, where for each polygon, its match ID
as well as the quality of the match is listed. Within this table 'map' = 0 refers to the osm data and 'map' = 1 to alkis.

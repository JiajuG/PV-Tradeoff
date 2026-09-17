# PV-Tradeoff

This repository provides supplementary information related to Sentinel-2 processing on Google Earth Engine, the MaxEnt suitability model, and the Pareto-constrained ecological-development co-optimization framework used in this study.

## Sentinel Data Processing Code on Google Earth Engine (GEE)

To facilitate reproducibility, we provide the code used to process and download Sentinel data on the GEE platform.
<img width="1430" height="744" alt="image" src="https://github.com/user-attachments/assets/f6c58d96-2f47-4c67-bac2-1b3658d1aef3" />


```javascript

var ASSET_PV = 'projects/global-phenology/assets/RCR-PV-PV2020';
var ASSET_BUFFER = 'projects/global-phenology/assets/RCR-BUFFER-PV2020';
var ID_FIELD = '编号';

var START_YEAR = 2017;
var END_YEAR = 2023;
var BUILD_YEAR = 2020;
var SCALE = 10;
var MIN_VALID_OBS = 5;
var MAX_SCENE_CLOUD = 80;
var PARALLEL_SCALE = 8;
var TILE_SCALE = 16;
var EXPORT_FOLDER = 'RCR_Sentinel2_PV2020_Yearly';

var pv = ee.FeatureCollection(ASSET_PV);
var buffer = ee.FeatureCollection(ASSET_BUFFER);
var analysisRegion = buffer.geometry().bounds();

print('PV feature count:', pv.size());
print('Buffer feature count:', buffer.size());
print('PV unique station IDs:', pv.aggregate_count_distinct(ID_FIELD));
print('Buffer unique station IDs:', buffer.aggregate_count_distinct(ID_FIELD));

Map.centerObject(pv, 4);
Map.addLayer(pv, {color: 'red'}, 'PV', false);
Map.addLayer(buffer, {color: 'blue'}, 'Buffer', false);

// ============================================================
// 1. Mask invalid observations and calculate NDVI/EVI
// ============================================================

function maskAndAddIndices(image) {
  image = image.select(['B2', 'B4', 'B8', 'SCL']);

  var scl = image.select('SCL');
  var validMask = scl.neq(0)
    .and(scl.neq(1))
    .and(scl.neq(3))
    .and(scl.neq(6))
    .and(scl.neq(8))
    .and(scl.neq(9))
    .and(scl.neq(10))
    .and(scl.neq(11));

  var blue = image.select('B2').multiply(0.0001);
  var red = image.select('B4').multiply(0.0001);
  var nir = image.select('B8').multiply(0.0001);
  var reflectanceMask = blue.gt(0).and(red.gt(0)).and(nir.gt(0));

  var ndvi = nir.subtract(red).divide(nir.add(red))
    .rename('NDVI')
    .updateMask(nir.add(red).abs().gt(0.0001))
    .updateMask(validMask)
    .updateMask(reflectanceMask);

  var eviDenominator = nir.add(red.multiply(6))
    .subtract(blue.multiply(7.5)).add(1);

  var evi = nir.subtract(red).multiply(2.5)
    .divide(eviDenominator)
    .rename('EVI')
    .updateMask(eviDenominator.abs().gt(0.0001))
    .updateMask(validMask)
    .updateMask(reflectanceMask);

  ndvi = ndvi.updateMask(ndvi.gte(-1).and(ndvi.lte(1)));
  evi = evi.updateMask(evi.gte(-1).and(evi.lte(1.5)));

  return ee.Image.cat([ndvi, evi])
    .copyProperties(image, ['system:time_start', 'system:index']);
}

// ============================================================
// 2. Pixel-wise annual maximum
// ============================================================

function makeAnnualTmax(year) {
  var start = ee.Date.fromYMD(year, 1, 1);
  var end = start.advance(1, 'year');

  var collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(analysisRegion)
    .filterDate(start, end)
    .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', MAX_SCENE_CLOUD))
    .select(['B2', 'B4', 'B8', 'SCL'])
    .map(maskAndAddIndices);

  var temporalReducer = ee.Reducer.max().combine({
    reducer2: ee.Reducer.count(),
    sharedInputs: true
  });

  var annualStats = collection.select(['NDVI', 'EVI'])
    .reduce(temporalReducer, PARALLEL_SCALE);

  var ndviCount = annualStats.select('NDVI_count');
  var eviCount = annualStats.select('EVI_count');

  var ndviTmax = annualStats.select('NDVI_max')
    .rename('NDVI_Tmax')
    .updateMask(ndviCount.gte(MIN_VALID_OBS));

  var eviTmax = annualStats.select('EVI_max')
    .rename('EVI_Tmax')
    .updateMask(eviCount.gte(MIN_VALID_OBS));

  return ee.Image.cat([ndviTmax, eviTmax]).set({
    year: year,
    relative_year: year - BUILD_YEAR,
    source_image_count: collection.size()
  });
}

// ============================================================
// 3. Polygon spatial mean, median and valid pixel count
// ============================================================

var spatialReducer = ee.Reducer.mean()
  .combine({
    reducer2: ee.Reducer.median(),
    sharedInputs: true
  })
  .combine({
    reducer2: ee.Reducer.count(),
    sharedInputs: true
  });

function reduceZone(annualImage, features, zoneType, year) {
  return annualImage.reduceRegions({
    collection: features,
    reducer: spatialReducer,
    scale: SCALE,
    tileScale: TILE_SCALE,
    maxPixelsPerRegion: 1e8
  }).map(function(feature) {
    return feature.set({
      construction_year: BUILD_YEAR,
      year: year,
      relative_year: year - BUILD_YEAR,
      zone_type: zoneType,
      source_image_count: annualImage.get('source_image_count'),
      min_valid_observations: MIN_VALID_OBS,
      temporal_method: 'pixelwise_annual_maximum',
      spatial_method: 'polygon_mean_median_count',
      spatial_scale_m: SCALE
    }).setGeometry(null);
  });
}

var selectors = [
  ID_FIELD,
  'construction_year',
  'year',
  'relative_year',
  'zone_type',
  'NDVI_Tmax_mean',
  'NDVI_Tmax_median',
  'NDVI_Tmax_count',
  'EVI_Tmax_mean',
  'EVI_Tmax_median',
  'EVI_Tmax_count',
  'source_image_count',
  'min_valid_observations',
  'temporal_method',
  'spatial_method',
  'spatial_scale_m'
];

// ============================================================
// 4. Create seven independent yearly export tasks
// ============================================================

for (var year = START_YEAR; year <= END_YEAR; year++) {
  var annualImage = makeAnnualTmax(year);
  var pvTable = reduceZone(annualImage, pv, 'PV', year);
  var bufferTable = reduceZone(annualImage, buffer, 'BUFFER', year);
  var yearlyTable = pvTable.merge(bufferTable);

  Export.table.toDrive({
    collection: yearlyTable,
    description: 'RCR_S2_Tmax_PV_BUFFER_' + year,
    folder: EXPORT_FOLDER,
    fileNamePrefix: 'RCR_S2_Tmax_PV_BUFFER_' + year,
    fileFormat: 'CSV',
    selectors: selectors
  });
}
```

## Pareto-Constrained Optimization of PV Layouts

The framework combines vegetation-greening probabilities predicted by XGBoost with PV development suitability predicted by MaxEnt. Candidate pathways are generated through a systematic search over the full weight domain. Pareto non-dominated solutions and their distances to the ideal point are then used to derive a fixed compromise weight for each vegetation-response scenario.

Uncertainty in vegetation responses is represented by the 10th, 50th, and 90th percentile scenarios.

```text
INPUTS:
    PV development suitability for pixel i: S_maxent[i]
    Annual energy yield for pixel i: Y[i]
    Area of pixel i: A[i]
    Vegetation-greening probability under scenario q: P_green[i, q]
    Vegetation-response class under scenario q: L[i, q]
    Vegetation-response scenarios Q = {10th, 50th, 90th}
    Candidate weights Alpha = {0.00, 0.01, ..., 1.00}
    Cumulative annual energy-yield evaluation targets T

FOR each vegetation-response scenario q in Q:

    Retain pixels with valid suitability, energy-yield,
    vegetation-response, and area data.

    FOR each valid pixel i:
        EcologicalCost[i, q] = 1 - P_green[i, q]
        DevelopmentCost[i]  = 1 - S_maxent[i]

    EcoRank[:, q] = percentile_rank(EcologicalCost[:, q])
    DevRank[:]    = percentile_rank(DevelopmentCost[:])

    Assign average ranks to tied values.
    CandidateSolutions = empty table

    FOR each alpha in Alpha:

        FOR each valid pixel i:
            PriorityScore[i] =
                alpha * DevRank[i]
                + (1 - alpha) * EcoRank[i, q]

        Sort all valid pixels in ascending order of PriorityScore
        to obtain one complete development sequence.

        Along this sequence, calculate at every position k:

            CumulativeEnergy[k] = SUM(Y[i]) / 10^9

            CumulativeBrowningArea[k] =
                SUM(A[i] * Indicator(L[i, q] = browning))

            CumulativeDevelopmentCost[k] =
                SUM(DevelopmentCost[i])

        FOR each cumulative annual energy-yield target t in T:
            Find the sequence position k nearest to t.

            Add the following record to CandidateSolutions:
                target = t
                weight = alpha
                browning_area = CumulativeBrowningArea[k]
                development_cost = CumulativeDevelopmentCost[k]

    TargetSpecificWeights = empty list

    FOR each cumulative annual energy-yield target t in T:

        Extract all candidate solutions evaluated at target t.

        Identify the Pareto non-dominated solutions by jointly
        minimizing cumulative browning area and cumulative
        development cost.

        Normalize both objectives within the Pareto set using
        min-max normalization.

        FOR each solution j in the normalized Pareto set:
            IdealDistance[j] = SQRT(
                NormalizedBrowningArea[j]^2
                + NormalizedDevelopmentCost[j]^2
            )

        Select the Pareto solution with the minimum IdealDistance.
        Record its weight as the target-specific compromise weight.

    Average the target-specific compromise weights over the focal
    planning evaluation nodes.

    Map the average to the nearest value in Alpha to obtain the
    scenario-specific fixed compromise weight alpha_fixed[q].

    Generate the final pathways:

        Scenario A, ecological priority:
            alpha = 0

        Scenario B, development priority:
            alpha = 1

        Scenario C, Pareto compromise:
            alpha = alpha_fixed[q]

    FOR each of Scenarios A, B, and C:
        Calculate the priority score using its fixed weight.
        Sort the valid pixels once in ascending score order.
        Use successive prefixes of this single sequence to represent
        increasing cumulative annual energy-yield levels.

        The resulting layouts are spatially nested: each layout at a
        lower energy-yield level is retained within layouts at higher
        energy-yield levels.

OUTPUTS:
    Scenario-specific fixed compromise weights
    Continuous development sequences
    Cumulative annual energy-yield trajectories
    Cumulative browning-area trajectories
    Cumulative development-cost trajectories
    Spatial PV-layout rasters
```

This repository documents the optimization procedure and does not contain the original spatial datasets or model-training data.

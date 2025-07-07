
# UAV Data Processing for Cassava Field Trials (2022–2024)

This repository contains three Jupyter notebooks designed to process UAV-derived raster data from cassava field trials conducted in Taiwan over three growing seasons: 2022, 2023, and 2024. The notebooks perform geospatial processing to extract quantitative plant traits from DEMs (Digital Elevation Models).

---

## Goal

The objective of these scripts is to extract per-plot values of key plant traits (canopy height and volume) from UAV-based raster outputs using georeferenced plant polygons. These outputs are formatted into tabular `.csv` files, ready for downstream phenotyping analysis.

---

## General Logic

Each notebook follows a similar multi-step processing pipeline:

1. **Load Inputs**:
   - Load plant ROI shapefile (with attributes like `plot_name`, `plot_number`).
   - Load raster files like DEMs, Crop Surface Models (CSMs), from specified folders.

2. **Per-plant Analysis**:
   - For each plant polygon, mask the raster and extract only the pixel values within that polygon.

3. **Trait Extraction**:
   - **Height**: Compute the 95th percentile of valid height pixels.
   - **Volume**: Compute the total sum of valid volume pixels.

4. **Filtering**:
   - Exclude plots labeled as "IH" (in-harvest) for later dates.
   - In 2024, entries starting with 2024_EF are excluded, as they were not part of the results presented.

5. **Export**:
   - The processed data is saved into `.csv` files for height and volume, named accordingly.

---

## Required Data Inputs

For each year, the following two inputs are required:

| Input | Description |
|-------|-------------|
| "add the path to the shapefile containing plant-level plot geometries" | GeoJSON or SHP file with individual plant polygons. Must include `plot_name`, `harvest`, and optionally `no_plant`. |
| "add the path to your DEM raster files" | Folder containing Digital Elevation Model (DEM) rasters per date. |

From these two inputs:
- Canopy Surface Models (CSMs) and volume rasters are automatically computed using the DEMs.
- Then, plant height (95th percentile) and volume (summed values) are extracted per plant geometry and saved as yearly CSV outputs.

---

## Notes per Year

- **2022 & 2023**: Each plot is processed independently.
- **2024**: Entries starting with `2024_EF` excluded from the outputs, as they were not included in the plots used in the paper.

---

## Output Files

Each notebook creates two `.csv` files named according to the year:

- `height_YYYY.csv`: Contains plot-level height metrics.
- `vol_YYYY.csv`: Contains plot-level volume estimates.

Each row includes:
- `plot_name`
- `date` (formatted as `YYYYMMDD`)
- `height_quantile_95` or `volume`

---

## How to Run

1. Open each notebook (e.g., `uav_data_processing_2022.ipynb`) in JupyterLab or VSCode.
2. Modify the placeholder path strings to match your data folders.
3. Choose the data to process by setting the `dates_to_process` variable in the pop up window (options are `"all"` or a list of dates in `DDMM` format).
4. Run all cells sequentially.
5. Outputs will be saved to the specified destinations.

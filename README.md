# DamTools

The DamTools.pyt is a Python Toolbox for ArcGis Pro (ESRI Inc.).
It contains several scripts to estimate and analyse dams and their reservoirs.

## Automated reservoir estimation (2022.03)

Tool to estimate a hydropower dam's reservoir size and shape using the dam position planned metrics, like dam wall height and installed hydropower capacity, in combination with the surrounding elevation model and river network.

### Geoprocessing Parameters:

- **Dam Positions:** The dam positions used for the estimation process as point feature class. **Notice:** Attribute table fields for the dam positions need to be defined *in* the Python Code of the tool itself to be able to use the given data (see below).
- **SRTMGL3 DEM:** The elevation model used for the estimation process. It should cover the surrounding area of the dam, but using a continent or global DEM is fine when a BasinATLAS is defined for the tool. The SRTMGL3 DEM (https://lpdaac.usgs.gov/products/srtmgl3v003/) with a resolujtion of ~90 Meters is recommended and tested. When a different DEM is used the resolution in the Python code may needs to be adjusted. The "Combine DEM rasters" tool can be used to combine the downloaded DEM rasters. The SRTMGL3 raster for Africa is available here: https://myshare.uni-osnabrueck.de/f/84c2d95f3875430d9560/
- **RiverATLAS river segments:** The river segments with information about discharges to be used for the water drop height estimation required to reach the planned hydropower capacity. The RiverATLAS (https://www.hydrosheds.org/page/hydroatlas) is recommended and a continental or global river atlas is fine, when used with a BasinATLAS. When a different river network is used the discharge field names need to be adjusted in the Python code.
- **Height Step [Meters]:** The vertical resolution used for the estimation. E.g. a Height step of 2 meters results in reservoir calculations for every other meter to be then checked for accordance with the given dam parameters like dam wall height or hydropower capacity.
- **Number of upriver river segments defining the initial search area:** The DEM is cut to an initial search area using X upriver river segments of the dam position. The standard value of 700 is set to encompass an area way larger than the expected reservoir will have. If other river networks are used the value may needs to be adjusted, if the resulting watershed is noticably to big or small.
- **Dam ID List: (optional)** A comma seperated list of dam IDs can be given to the tool to calculate only a subset of the given dam postions.
- **Project Name:** The name of the current tool run, it will create a new folder with the given name in the **Output Folder**.
- **Output Folder:** The path to the folder, where all new tool run folders will be created.
- **BasinATLAS: (optional, but recommended)** A database containing river basins to highly increase reservoir caluculation when multiple dams are calculated in one run or continental or global datasets are used for the elevation model or river network. The basins are used to partition dam positions and environmental data according to the basins to decrease dataset sizes for the calcualtion. The recommended BasinATLAS lev03 (https://www.hydrosheds.org/page/hydroatlas) yielded the best time savings in the testing phase.
- **Keep all temporary files for debugging?:** The tool deletes all temporary files and folders during or after the estimation process to highly decrease disk space of the results. If a reservoir yields unrealistic results they can be rerun (using the **Dam ID List**) with this parameter set to True, to analyse the files created in the process. 

### Python Parameters:

The Python Parameters need to be defined *in* the Python code itself.
They are set in the *class Automated_reservoir_estimation():* in the function *def \__init\__()* and they all start with *self.xxx*.
Soem of the parameters need to be adjusted, when a new dam position layer is used, while the other can be left unchanged if the recommended data sources are used.
The once that need to be changed are:
- **self.height_field:** Needs to be set to the name of the attribute table field containing the projected dam wall heights in meters (None or Null values in the attribute table are fine here).
- **self.cap_field:** Needs to be set to the field name of the hyopower capacity in MW.
- **self.id_field:** Needs to be set to the field name of *unique* dam IDs (when multiple dams share the same ID one of them will be overwritten).
- **self.lat_field:** Needs to be set to the field name of the lattitude value of the dam positions.

GIS Hydro - Geoprocessing Tasks CHANGELOG
===========================================================
- Created by: Javier Mardones
- Created on: 06/2021
- Purpose: GIS Hydro Geoprocessing Tasks
- Description: Provides the back end change log to the GIS Hydro application.

---------------------------------------------------------

## (02/27/2024)

- Developer: Javier Mardones
- Version 0.96
- Added MARISA option to precipitation and TR20 scripts

## (01/13/2024)

- Developer: Javier Mardones
- Version 0.954
- Fixed watersheds not being delineated outside of maryland
- Fixed snapping issues

## (08/01/2024)

- Developer: Javier Mardones
- Version 0.952
- Fixed FRRE Imperviousness value reference (stuck on 27.91)

## (07/12/2024)

- Developer: Javier Mardones
- Version 0.951
- Fixed TOTTIME output from the Velocity Method

## (07/10/2024)

- Developer: Javier Mardones
- Version 0.95
- Fixed AVGAREA output from the Velocity Method

## (06/26/2024)

- Developer: Javier Mardones
- Version 0.94
- Formatting changes to Python scripts (it will probable make the GP services run faster)
- Fix to flowpaths not running in server

## (05/29/2024)

- Developer: Javier Mardones
- Version#: 0.93
- Made numerous corrections to ensure a seamless interface and more effective error handling.
- Improved some GP Task scripts for better readibility and efficiency, with plans to improve all others in the future.
- Fixed issue where the FRRE would vary based on the user’s layer selection. Now, the FRRE consistently utilizes the same datasets (MDP2010, SSURGO, 30-meter). If a watershed extends into another state, the NLCD 2019 land use dataset will be employed.
- Removed LFP correction from TC options
- Added Database version to interface

## (05/28/2024)

- Developer: Javier Mardones
- Version#: 0.92
- Fixed BasinOuputs GP issues with local Land Use CLASS_NAME id not found
- Improved DataSelection and BasinOutputs code, also changed the GetParameter input index from "set_i" to the respective number/index

## (05/15/2024)

- Developer: Javier Mardones
- Version #: 0.9
- Work: FY24 task order changes before testing

## (06/??/2021)

- Developer: Javier Mardones
- Version #: 0.2
- Work: Initial Development


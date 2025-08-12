GISHydroWEB  - Geoprocessing Tasks
===========================================

## GISHydroWEB Version

0.954

## GISHydroWEB Overview
The geoprocessing tasks in this repo support the front-end GISHydroWEB application.
This repo contains multiple Python scripts collected into a GP Toolbox. The
scripts should be deployed to an ArcGIS Server instance and published as Geoprocessing

[Changelog](./changelog.md)

## Table of Contents
- [Requirements](#requirements)
- [Developer Documentation](#developer-documentation)
- [Deployment](#deployment)
- [Getting Started](#getting-started)
- [Acknowledgements](#acknowledgements)
- [License](#license)

## Requirements

- ArcGIS v10.5 or higher
- Machine with ArcGIS installed (for running and publishing the scripts)
- The Microsoft **azcopy** command v10+ must exist on the target machine's path


## Developer Documentation

It is recommended to follow [Semantic Commit Messages](https://gist.github.com/joshbuchea/6f47e86d2510bce28f8e7f42ae84c716)
to make updates easier  to categorize and understand.

- [Branching Strategy](./branching-strategy.md)

## Deployment
GISHydroWEB consists of a set of Geoprocessing Services, which are deployed
using an ArcGIS Server account.

A File Geodatabase (FGDB) of various supporting data need to be available for these
scripts to run properly. Please contact program developer (Javier Mardones, javier.mardones@mbakerintl.com)
for more information regarding the data sets.

## Needed environment setup
1. Pull scripts from repo to target machine( automated by build pipeline at SHA )
1. Configure environment variables
    - **GISHydro_Dir** -- to point to the parent folder where the data will live,
        e.g. _D:\GISHydroWEB_. The scripts have a hardcoded system variable path
        for the "data" folder beneath this directory.
    - **GIS_HYDRO_APP_ID** -- the Azure ID secret for the Azure ID which has
        permission to copy the files from Azure Blob Storage
    - **AZURE_TENANT_ID** -- is the tenant ID for Azure where the blob storage
        exists.
1. From a machine with ArcGIS installed (perhaps the server where these are published),
   navigate to the deployment_scripts folder using the Python
1. From with the deployment_scripts folder run the __publish_hydro_gp.py__ script with
	```
	>python ./publish_hydro_gp.py
	```
 This script automates the deployment of the GP Tasks to the ArcGIS server. Detailed manual instructions are provided below should there be a problem with the script.


**Important: The toolbox scripts should live under the same parent as the data and look similar to this:**

```
   _GISHydro_Dir
   |__data
   |__deployment_scripts
   |__projects
   |__scripts
```

#### How to create a system variable in Windows 10:

- In Search, search for and then select: System (Control Panel)
- Click the Advanced system settings link.
- Click Environment Variables. In the section System Variables create a new variable "GISHydro_DIR"
- In the Edit System Variable (or New System Variable) window, specify the path
  of your GISHydroWEB data to the newly created "GISHydro_Dir" environment variable.
  Click OK.
- Close all remaining windows by clicking OK.

### Manually Publishing the Geoprocessing Services:
The GISHydro tools consist of two subsets of tools:

1. GISHydro
2. Tools

The GISHydro tool has all the required scripts to make a GISHydroWEB project (e.g., Project creation, Delineation, etc.), while the Tools tool has optional scripts that will help the user in the process (e.g., loading layers to the framework)

The GISHydro tool scripts must be executed in a specific order:

1. DataSelection
2. BasinOutputs
3. Flowpaths
4. Outlets
5. Subsheds
6. SetTOC
7. Transect
8. Reservoir
9. PrecipitationDepths
10. SetTOC

### The Tools scripts can be executed in any order. However, they must be executed after the GISHydro tool is done.

- DesignStorm
- LoadLayer
- Tasker
- DelCheck

* Once the "DataSelection" has run, a project folder ("project") will be created
(if it was not created previously), and a new project will be created inside
the **project** folder.
The tools already contain default values included, however, most of the tools require an input folder (Project Name) that must be typed manually (or copied) from the newly created project in the "project" folder (e.g., "20210602_160047_My_Project").

To publish these services please [refer to ESRI's tutorial](https://enterprise.arcgis.com/en/server/10.8/publish-services/windows/a-quick-tour-of-publishing-a-geoprocessing-service.htm)

Once the geoprocessing services are published, a link will be created (REST API) that will allow to execute the tool from the Server URL. Use the link for each service in the GISHydro Web Application's JavaScript code.

Please refer to the GISHydroWEB documentation for the GISHydroWEB architecture.

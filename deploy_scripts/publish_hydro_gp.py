## Run and publish the gisHydro services to ArcGIS Server
## Template for approach: https://pro.arcgis.com/en/pro-app/latest/arcpy/functions/creategpsddraft.htm
import inspect
import os
import shutil
import tempfile
from azcopyfiles import azcopyfiles

from pprint import pprint

from arcpy import (
    ImportToolbox,
    mp,
    CreateGPSDDraft,
    StageService_server,
    UploadServiceDefinition_server,
)

# makes a complete fresh path and returns the name
def makePath(path):
    print(path)
    shutil.rmtree(path, ignore_errors=True)
    os.makedirs(path)
    return path

## Rough outline of tasks:
## loop through list of services
## create an sd for each
## create a draft for each
## combine to a single result
## stage service(need admin creds for this part)
## publish to server

## list of services we are running and publishing
hydro_services = (
    'DataSelection',
    'BasinOutputs',
    'Flowpaths',
    'Outlets',
    'Subsheds',
    'SetTOC',
    'Transect',
    'Reservoir',
    'PrecipitationDepths',
    'TR20Control',
    'DesignStorm',
    'LoadLayer',
    'FRREBasin',
    'FRRECalculator',
    'DelCheck',
    'UserStreams',
    'DownloadLayer'
)

project_not_required_services = [
    'DataSelection',
    'DesignStorm',
    'FRRECalculator',
    'DelCheck',
    'UserStreams'
]

## first copy the data files locally
## This ensure the new scripts have the latest data to work on

print('Beginning Data File Copy.')
azcopyfiles()
print('Completed Copy of Files.')

# get the machines temp directory
# on windows this is generally: C:\Users\<user_name>\AppData\Local\Temp\
work_path = r"" + str(os.environ['GISHydro_DIR'])

# The connection_path is picked up via env variable
# we have a centralized location for storing ags connections
connection_path = r"" + str(os.environ['LOCAL_AGS_ADMIN_CONN'])
temp_path = tempfile.gettempdir()
toolbox_path = os.path.join(work_path, 'GISHydroWeb.tbx')
# create a service_definitions sub-directory
sd_path = makePath(os.path.join(temp_path, 'gishydro', 'service_definitions'))
# create a service_draft sub-directory
sd_draft_path = makePath(os.path.join(temp_path, 'gishydro', 'service_drafts'))

toolbox = ImportToolbox(toolbox_path, 'GISHydro')
print('This is our toolbox:', toolbox)

# this list will hold all the results generated from running each service
results = []
PROJECT_NAME_IDX = 1
project_name = ''

for service_name in hydro_services:
    print(service_name)
    sd_file = os.path.join(sd_path, "" + service_name + ".sd")
    draft_file = os.path.join(sd_draft_path, "" + service_name + ".sddraft")
    hydro_draft = os.path.join(sd_draft_path, "gisHydro.sddraft")
    hydro_stage= os.path.join(sd_path, "gisHydro.sd")
    print('\tSD AND DRAFT:\n\t', sd_file, '\n\t', draft_file)

    # import the service and run
    service = getattr(toolbox, service_name)
    print('\tTOOL IS: ', service)
    params = inspect.getargspec(service).args
    for param in params:
        print('\t\t', param)

    # execute the tool (project name is set below)
    if service_name in project_not_required_services:
        result = service()
    else:
        print('CALLING WITH PROJECT NAME:', project_name)
        result = service(projectname=project_name)

    # this will just print a string representation of the object
    print('Result is:', result)
    # this prints the actual members of the object
    pprint(inspect.getmembers(result))

    # add result to list of all results
    results.append(result)

    # Once the First items runs we have a proper reference which the dependent
    # tools will use
    if service_name == 'DataSelection':
        project_name = result.getOutput(PROJECT_NAME_IDX)
        print('PROJECT NAME SET TO:', project_name)


# Not necessary, but lets see what happens
print('Results are:', '\n'.join(map(str,results)) )

# we create a combined service draft with the results of all the tools
# This will publish them to a single GP service
print ('Creating a combined draft for GIS Hydro')
## CreateGPSDDraft (result, out_sddraft, service_name, {server_type}, {connection_file_path},
## {copy_data_to_server}, {folder_name}, {summary}, {tags}, {executionType}, {resultMapServer},
## {showMessages}, {maximumRecords}, {minInstances}, {maxInstances}, {maxUsageTime}, {maxWaitTime},
## {maxIdleTime}, {capabilities}, {constantValues})
draft_results = CreateGPSDDraft(results, hydro_draft, 'GISHydro',
    server_type="ARCGIS_SERVER",
    connection_file_path=connection_path,
    executionType="Synchronous",
    copy_data_to_server=False,
    folder_name='GISHydro',
    ##summary='GIS Hydro Services', tags='GIS Hydro'
    )

# Enable this to allow the review of the settings
# xml_doc = xml.dom.minidom.parse(hydro_draft)
# settings_xml = xml_doc.getElementsByTagName('StagingSettings')[0].toprettyxml()
# resources_xml = xml_doc.getElementsByTagName('Resources')[0].toprettyxml()
# print(settings_xml, resources_xml)

# Next analyze for errors
# Script written to utilize Python 3, however "AnalyzeForSD" is not available
#   see: https://support.esri.com/en/bugs/nimbus/QlVHLTAwMDEyMTEyNg==
# so we poll the result object for errors
if draft_results['errors'] == {}:

    # Next we need to stage the service.
    print('Staging the GIS Hydro services')
    StageService_server(hydro_draft, hydro_stage)

    # This is where the service is pushed to the server
    # If odd 404 errors arise here, check that the PublishingToolsEx service
    # is running on the server (under the Systems folder)
    # see: https://community.esri.com/t5/arcgis-pro-questions/map-service-fails-to-publish-in-pro-2-4/td-p/164469/page/2
    print ('Beginning the publish of the services');
    # used for testing: in_folder_type='NEW', in_folder='hydrotest'
    upStatus = UploadServiceDefinition_server(hydro_stage, connection_path)
    print ('Up Status: ', upStatus)
    print ("Completed upload")

else:
    print('Error Occcured creating draft:\n')
    print(*draft_results['errors'], sep='\n')

print ('Exiting deployment script.')

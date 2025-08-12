# Run this only if the main deploy script is failing during the staging and
# publishing state. It will try to upload the existing draft without the need
# to re-run all of the service tools
import xml.dom.minidom
import os
import shutil
import tempfile

from pprint import pprint

from arcpy import (
    ImportToolbox,
    mp,
    CreateGPSDDraft,
    StageService_server,
    UploadServiceDefinition_server,
)

connection_path = r"" + str(os.environ['LOCAL_AGS_ADMIN_CONN'])
temp_path = tempfile.gettempdir()
sd_path = os.path.join(temp_path, 'gishydro', 'service_definitions')
sd_draft_path = os.path.join(temp_path, 'gishydro', 'service_drafts')
hydro_draft = os.path.join(sd_draft_path, "gisHydro.sddraft")
hydro_stage= os.path.join(sd_path, "gisHydro.sd")

xml_doc = xml.dom.minidom.parse(hydro_draft)
staging_settings = xml_doc.getElementsByTagName('StagingSettings')[0].toprettyxml()
staging_resources = xml_doc.getElementsByTagName('Resources')[0].toprettyxml()
print(staging_settings, staging_resources)

print ('Beginnng the publish of the services');

upStatus = UploadServiceDefinition_server(hydro_stage, connection_path)
print ('Up Status: ', upStatus)
print ("Completed upload")

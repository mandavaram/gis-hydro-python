"""Contains the azcopy command to move the GIS Hydro data files"""
## Assumes the two environment variables to be in place
## 1. GIS_HYDRO_APP_ID = A service ID created in azure with permission to the blob container
## 2. AZURE_TENANT_ID = the tenant id where the blob storage lives
import os

def azcopyfiles():
    """Uses azcopy to copy data files from Azure blob to local directory"""
    work_dir = os.environ["GISHydro_DIR"]
    data_dir = os.path.join(work_dir, 'data')

    print('Beginning copy of GIS Hydro Data Files.')
    # utilizing auto login on target server
    copy_process = (
        f"azcopy cp "
        f"https://storagisbackup.blob.core.windows.net/gis-hydro/data/* "
        f"{data_dir} --recursive=true --overwrite=ifSourceNewer"
    )
    copy_result = os.system(copy_process)
    if copy_result != 0:
        raise Exception('The azcopy failed. Check Blob container is correct.')
    print('Completed copy.')

if __name__ == "__main__":
    azcopyfiles()

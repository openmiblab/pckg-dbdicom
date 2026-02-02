import os
import shutil
import zipfile
from pathlib import Path
from typing import Union
from tqdm import tqdm
import numpy as np
import vreg

from dbdicom.dbd import DataBaseDicom




def open(path:str) -> DataBaseDicom:
    """Open a DICOM database

    Args:
        path (str): path to the DICOM folder

    Returns:
        DataBaseDicom: database instance.
    """
    return DataBaseDicom(path)

def to_json(path):
    """Summarise the contents of the DICOM folder in a json file

    Args:
        path (str): path to the DICOM folder
    """
    dbd = open(path)
    dbd.close()   

def print(path):
    """Print the contents of the DICOM folder

    Args:
        path (str): path to the DICOM folder
    """
    dbd = open(path)
    dbd.print()
    dbd.close()


def summary(path) -> dict:
    """Return a summary of the contents of the database.

    Args:
        path (str): path to the DICOM folder

    Returns:
        dict: Nested dictionary with summary information on the database.
    """
    dbd = open(path)
    s = dbd.summary()
    dbd.close()
    return s


def tree(path) -> dict:
    """Return the structure of the database as a dictionary tree.

    Args:
        path (str): path to the DICOM folder

    Returns:
        dict: Nested dictionary with summary information on the database.
    """
    dbd = open(path)
    s = dbd.register
    dbd.close()
    return s


def to_csv(path, csv_file) -> dict:
    """Write a summary of the contents of the database to csv.

    Args:
        path (str): path to the DICOM folder
        csv_file (str): path to the csv file
    """
    dbd = open(path)
    dbd.to_csv(csv_file)
    dbd.close()


def patients(path, name:str=None, contains:str=None, isin:list=None)->list:
    """Return a list of patients in the DICOM folder.

    Args:
        path (str): path to the DICOM folder
        name (str, optional): value of PatientName, to search for 
            individuals with a given name. Defaults to None.
        contains (str, optional): substring of PatientName, to 
            search for individuals based on part of their name. 
            Defaults to None.
        isin (list, optional): List of PatientName values, to search 
            for patients whose name is in the list. Defaults to None.

    Returns:
        list: list of patients fulfilling the criteria.
    """
    dbd = open(path)
    p = dbd.patients(name, contains, isin)
    dbd.close()
    return p


def studies(entity:str | list, desc:str=None, contains:str=None, isin:list=None)->list:
    """Return a list of studies in the DICOM folder.

    Args:
        entity (str or list): path to a DICOM folder (to search in 
            the whole folder), or a two-element list identifying a 
            patient (to search studies of a given patient).
        desc (str, optional): value of StudyDescription, to search for 
            studies with a given description. Defaults to None.
        contains (str, optional): substring of StudyDescription, to 
            search for studies based on part of their description. 
            Defaults to None.
        isin (list, optional): List of StudyDescription values, to search 
            for studies whose description is in a list. Defaults to None.

    Returns:
        list: list of studies fulfilling the criteria.
    """
    if isinstance(entity, str): # path = folder
        dbd = open(entity)
        s = dbd.studies(entity, desc, contains, isin)
        dbd.close()
        return s
    elif len(entity)==2: # path = patient
        dbd = open(entity[0])
        s = dbd.studies(entity, desc, contains, isin)
        dbd.close()
        return s
    else:
        raise ValueError(
            "The path must be a folder or a 2-element list "
            "with a folder and a patient name."
        )

def series(entity:str | list, desc:str=None, contains:str=None, isin:list=None)->list:
    """Return a list of series in the DICOM folder.

    Args:
        entity (str or list): path to a DICOM folder (to search in 
            the whole folder), or a list identifying a 
            patient or a study (to search series of a given patient 
            or study).
        desc (str, optional): value of SeriesDescription, to search for 
            series with a given description. Defaults to None.
        contains (str, optional): substring of SeriesDescription, to 
            search for series based on part of their description. 
            Defaults to None.
        isin (list, optional): List of SeriesDescription values, to search 
            for series whose description is in a list. Defaults to None.

    Returns:
        list: list of series fulfilling the criteria.
    """
    if isinstance(entity, str): # path = folder
        dbd = open(entity)
        s = dbd.series(entity, desc, contains, isin)
        dbd.close()
        return s
    elif len(entity) in [2,3]:
        dbd = open(entity[0])
        s = dbd.series(entity, desc, contains, isin)
        dbd.close()
        return s
    else:
        raise ValueError(
            "To retrieve a series, the entity must be a database, patient or study."
        )
    
def copy(from_entity:list, to_entity=None):
    """Copy a DICOM  entity (patient, study or series)

    Args:
        from_entity (list): entity to copy
        to_entity (list, optional): entity after copying. If this is not 
            provided, a copy will be made in the same study and returned.

    Returns:
        entity: the copied entity. If th to_entity is provided, this is 
        returned.
    """
    dbd = open(from_entity[0])
    from_entity_copy = dbd.copy(from_entity, to_entity)
    dbd.close()
    return from_entity_copy


def delete(entity:list, not_exists_ok=False):
    """Delete a DICOM entity

    Args:
        entity (list): entity to delete
        not_exists_ok (bool): By default, an exception is raised when attempting 
            to delete an entity that does not exist. Set this to True to pass over this silently.
    """
    dbd = open(entity[0])
    dbd.delete(entity, not_exists_ok)
    dbd.close()


def remove_duplicate_frames(entity:list, dims:list=None, verbose=1, filter:dict=None, dry_run=False):
    """Remove duplicate frames from the entity

    Args:
        entity (list): remove duplicates from this entity
        dims (list, optional): Dimensions of grid in which to search for duplicates. 
            If None are provided, the duplicate InstanceNumbers are removed.
        verbose (int, optional): Provide feedback of the computation. Defaults to 1.
        filter (dict, optional): keywords to filter the series. Duplicates will only 
            be removed from those that are not filtered out.
        dry_run (bool, optional): if True, files to be deleted are printed but not actually deleted.
            Default is False (delete for real!).
    Returns:
        list: deleted files
    """
    dbd = open(entity[0])
    files = dbd.remove_duplicate_frames(entity, dims, verbose, filter, dry_run)
    dbd.close()  
    return files 


def move(from_entity:list, to_entity:list):
    """Move a DICOM entity

    Args:
        entity (list): entity to move
    """
    dbd = open(from_entity[0])
    dbd.copy(from_entity, to_entity)
    dbd.delete(from_entity)
    dbd.close()

def split_series(series:list, attr:Union[str, tuple], key=None, verbose=1)->list:
    """
    Split a series into multiple series
    
    Args:
        series (list): series to split.
        attr (str or tuple): dicom attribute to split the series by. 
        key (function): split by by key(attr) 
        verbose (bool, optional): If set to 1, shows progress bar. Defaults to 1.
    Returns:
        list: list of two-element tuples, where the first element is
        is the value and the second element is the series corresponding to that value.      
    """
    dbd = open(series[0])
    split_series = dbd.split_series(series, attr, key, verbose)
    dbd.close()
    return split_series


def exists(entity:list):
    """Check if the entity exists in the database

    Args:
        entity (list): List of 1, 2, 3 or 4 elements.

    Returns:
        bool: True if the entity exists (i.e. has data), false otherwise.
    """
    dbd = open(entity[0])
    result = dbd.exists(entity)
    dbd.close()
    return result


def volume(series:list, dims:list=None, verbose=1, **kwargs) -> vreg.Volume3D:
    """Read volume from a series.

    Args:
        series (list, str): DICOM entity to read
        dims (list, optional): Non-spatial dimensions of the volume. Defaults to None.
        verbose (int, optional): If set to 1, shows progress bar. Defaults to 1.
        kwargs (dict, optional): keywords to filter the series.

    Returns:
        vreg.Volume3D.
    """
    dbd = open(series[0])
    vol = dbd.volume(series, dims, verbose, **kwargs)
    dbd.close()
    return vol


def slices(series:list, dims:list=None, verbose=1) -> vreg.Volume3D:
    """Read 2D volumes from the series

    Args:
        entity (list, str): DICOM series to read
        dims (list, optional): Non-spatial dimensions of the volume. Defaults to None.
        verbose (int, optional): If set to 1, shows progress bar. Defaults to 1.
        
    Returns:
        list of vreg.Volume3D
    """
    dbd = open(series[0])
    vol = dbd.slices(series, dims, verbose)
    dbd.close()
    return vol


# Obsolete API - phase out
def volumes_2d(*args, **kwargs):
    return slices(*args, **kwargs)


def values(series:list, *attr, dims:list=None, verbose=1, filter:dict=None) -> Union[np.ndarray, list]:
    """Read the values of some attributes from a DICOM series

    Args:
        series (list): DICOM series to read. 
        attr (tuple, optional): DICOM attributes to read.
        dims (list, optional): Dimensions to sort the values. 
            If dims is not provided, values are sorted by 
            InstanceNumber.
        filter (dict, optional): keywords to filter the series.

    Returns:
        tuple: arrays with values for the attributes.
    """
    dbd = open(series[0])
    values = dbd.values(series, *attr, dims=dims, verbose=verbose, filter=filter)
    dbd.close()
    return values



def write_volume(vol:Union[vreg.Volume3D, tuple], series:list, 
                 ref:list=None, append=False, verbose=1, **kwargs):
    """Write a vreg.Volume3D to a DICOM series

    Args:
        vol (vreg.Volume3D or tuple): Volume to write to the series.
        series (list): DICOM series to read
        dims (list, optional): Non-spatial dimensions of the volume. Defaults to None.
        append (bool): by default write_volume will only write to a new series, 
            and raise an error when attempting to write to an existing series. 
            To overrule this behaviour and add the volume to an existing series, set append to True. 
            Default is False.
        verbose (bool): if set to 1, a progress bar is shown. verbose=0 does not show updates.
        kwargs: Keyword-value pairs to be set on the fly
    """
    dbd = open(series[0])
    dbd.write_volume(vol, series, ref=ref, append=append, verbose=verbose, **kwargs)
    dbd.close()


def edit(series:list, new_values:dict, dims:list=None, verbose=1):
    """Edit attribute values in a DICOM series

    Warning: this function edits all values as requested. Please take care 
    when editing attributes that affect the DICOM file organisation, such as 
    UIDs, as this could corrupt the database.

    Args:
        series (list): DICOM series to edit
        new_values (dict): dictionary with attribute: value pairs to write to the series
        dims (list, optional): Non-spatial dimensions of the volume. Defaults to None.
        verbose (bool, optional): If set to 1, shows progress bar. Defaults to 1.
        
    """
    dbd = open(series[0])
    dbd.edit(series, new_values, dims=dims, verbose=verbose)
    dbd.close()


def to_npz(database, destination, dims:list=None, skip=False, overwrite=False):
    """Save a dicom database in a lightweight .npz (numpy) format

    This will save each series in a single .npz file following the 
    same folder structure and naming conventions as dicom 
    exports by dbdicom.

    Args:
        database (str): path to DICOM database
        destination (str): path to new npz database
        dims (list, optional): additional dimensions for the 
            volumes to be saved.
        skip (bool, optional): if skip=True, the function will 
            silently skip over and any series that 
            are not volumetric or cannot be arranged in the dimensions
            specified. If skip=False (default) an error will be 
            raised if such series are found.
        overwrite (bool, optional): if True, and existing 
            .npz files will be overwritten. If False (default) 
            they will be silently ignored and not overwritten.

    Note:
        At this point this feature is only available for 
        volumetric series. If multislice series are found, 
        or volumes that cannot be arranged in the dimensions 
        specified, an error will be raised unless skip is set to True.

        Most header information contained in DICOM data will NOT 
        be saved in the npz database. The only information 
        that is preserved is geometry (affine), dimensions and 
        coordinates of the volumes saved.
    """
    dbd = open(database)
    dbd.to_npz(destination, dims, skip, overwrite)
    dbd.close() 


def to_nifti(series:list, file:str, dims:list=None, verbose=1):
    """Save a DICOM series in nifti format.

    Args:
        series (list): DICOM series to read
        file (str): file path of the nifti file.
        dims (list, optional): Non-spatial dimensions of the volume. 
            Defaults to None.
        verbose (bool, optional): If set to 1, shows progress bar. Defaults to 1.
    """
    dbd = open(series[0])
    dbd.to_nifti(series, file, dims, verbose)
    dbd.close()

def from_nifti(file:str, series:list, ref:list=None):
    """Create a DICOM series from a nifti file.

    Args:
        file (str): file path of the nifti file.
        series (list): DICOM series to create
        ref (list): DICOM series to use as template.
    """
    dbd = open(series[0])
    dbd.from_nifti(file, series, ref)
    dbd.close()


def dir(entity:list) -> str:
    """Return directory of the entity

    Args:
        entity (list): DICOM entity as list (database, patient, study or series)

    Returns:
        list: directory
    """
    dbd = open(entity[0])
    dir = dbd.dir(entity)
    dbd.close()
    return dir


def files(entity:list) -> list:
    """Read the files in a DICOM entity

    Args:
        entity (list or str): DICOM entity to read. This can 
            be a path to a folder containing DICOM files, or a 
            patient or study to read all series in that patient or 
            study. 

    Returns:
        list: list of valid dicom files.
    """
    if isinstance(entity, str):
        entity = [entity]
    dbd = open(entity[0])
    files = dbd.files(entity)
    dbd.close()
    return files


def pixel_data(series:list, dims:list=None, verbose=1) -> tuple:
    """Read the pixel data from a DICOM series

    Args:
        series (list or str): DICOM series to read. This can also 
            be a path to a folder containing DICOM files, or a 
            patient or study to read all series in that patient or 
            study. In those cases a list is returned.
        dims (list, optional): Dimensions of the array.

    Returns:
        numpy.ndarray or tuple: numpy array with pixel values, with 
            at least 3 dimensions (x,y,z). 
    """
    if isinstance(series, str):
        series = [series]
    dbd = open(series[0])
    array = dbd.pixel_data(series, dims, verbose)
    dbd.close()
    return array


def unique(pars:list, entity:list) -> dict:
    """Return a list of unique values for a DICOM entity

    Args:
        pars (list, str/tuple): attribute or attributes to return.
        entity (list): DICOM entity to search (Patient, Study or Series)

    Returns:
        dict: if a pars is a list, this returns a dictionary with 
        unique values for each attribute. If pars is a scalar this returnes a list of values
    """
    dbd = open(entity[0])
    u = dbd.unique(pars, entity)
    dbd.close()
    return u


def archive(path, archive_path):
    dbd = open(path)
    dbd.archive(archive_path)
    dbd.close()


def restore(archive_path, path):
    _copy_and_extract_zips(archive_path, path)
    dbd = open(path)
    dbd.close()


def _copy_and_extract_zips(src_folder, dest_folder):
    if not os.path.exists(dest_folder):
        os.makedirs(dest_folder)

    # First pass: count total files
    total_files = sum(len(files) for _, _, files in os.walk(src_folder))

    with tqdm(total=total_files, desc="Copying and extracting") as pbar:
        for root, dirs, files in os.walk(src_folder):
            rel_path = os.path.relpath(root, src_folder)
            dest_path = os.path.join(dest_folder, rel_path)
            os.makedirs(dest_path, exist_ok=True)

            for file in files:
                src_file_path = os.path.join(root, file)
                dest_file_path = os.path.join(dest_path, file)

                if file.lower().endswith('.zip'):
                    try:
                        zip_dest_folder = dest_file_path[:-4]
                        if os.path.exists(zip_dest_folder):
                            continue
                        with zipfile.ZipFile(src_file_path, 'r') as zip_ref:
                            zip_ref.extractall(zip_dest_folder)
                        #tqdm.write(f"Extracted ZIP: {src_file_path}")
                        #_flatten_folder(zip_dest_folder) # still needed?
                    except zipfile.BadZipFile:
                        tqdm.write(f"Bad ZIP file skipped: {src_file_path}")
                else:
                    if os.path.exists(dest_file_path):
                        continue
                    shutil.copy2(src_file_path, dest_file_path)

                pbar.update(1)


def _flatten_folder(root_folder):
    for dirpath, dirnames, filenames in os.walk(root_folder, topdown=False):
        for filename in filenames:
            src_path = os.path.join(dirpath, filename)
            dst_path = os.path.join(root_folder, filename)
            
            # If file with same name exists, optionally rename or skip
            if os.path.exists(dst_path):
                base, ext = os.path.splitext(filename)
                counter = 1
                while os.path.exists(dst_path):
                    dst_path = os.path.join(root_folder, f"{base}_{counter}{ext}")
                    counter += 1

            shutil.move(src_path, dst_path)

        # Remove empty subdirectories (but skip the root folder)
        if dirpath != root_folder:
            try:
                os.rmdir(dirpath)
            except OSError:
                print(f"Could not remove {dirpath} — not empty or in use.")



if __name__=='__main__':
    pass
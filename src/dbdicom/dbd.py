import os
import shutil
import json
from typing import Union
import zipfile
import re
from copy import deepcopy

from tqdm import tqdm
import numpy as np
import vreg
from pydicom.dataset import Dataset
import pydicom

import dbdicom.utils.arrays
import dbdicom.dataset as dbdataset
import dbdicom.database as dbdatabase
import dbdicom.register as register
import dbdicom.const as const
from dbdicom.utils.pydicom_dataset import (
    get_values, 
    set_values,
    set_value,
    )



class DataBaseDicom():
    """Class to read and write a DICOM folder.

    Args:
        path (str): path to the DICOM folder.
    """

    def __init__(self, path):

        if not os.path.exists(path):
            os.makedirs(path)
        self.path = path

        file = self._register_file()
        if os.path.exists(file):
            try:
                with open(file, 'r') as f:
                    self.register = json.load(f)
                # remove the json file after reading it. If the database
                # is not properly closed this will prevent that changes
                # have been made which are not reflected in the json 
                # file on disk
                # os.remove(file)
            except Exception as e:
                # raise ValueError(
                #     f'Cannot open {file}. Please close any programs that are '
                #     f'using it and try again. Alternatively you can delete the file '
                #     f'manually and try again.'
                # )
                # If the file can't be read, delete it and load again
                os.remove(file)
                self.read()
        else:
            self.read()


    def read(self):
        """Read the DICOM folder again
        """
        self.register = dbdatabase.read(self.path)
        # For now ensure all series have just a single CIOD
        # Leaving this out for now until the issue occurs again.
        # self._split_series()
        return self

    def delete(self, entity, not_exists_ok=False):
        """Delete a DICOM entity from the database

        Args:
            entity (list): entity to delete
            not_exists_ok (bool): By default, an exception is raised when attempting 
                to delete an entity that does not exist. Set this to True to pass over this silently.
        """
        # delete datasets on disk
        try:
            removed = register.index(self.register, entity)
        except ValueError:
            if not_exists_ok:
                return self
            else:
                raise ValueError(
                    f"The entity you are trying to delete does not exist. \n"
                    f"You can set not_exists_ok=True in dbdicom.delete() to avoid this error."
                )
        for index in removed:
            file = os.path.join(self.path, index)
            if os.path.exists(file): 
                os.remove(file)
        # drop the entity from the register
        register.remove(self.register, entity)
        # cleanup empty folders
        remove_empty_folders(entity[0])
        return self
    

    def close(self): 
        """Close the DICOM folder
        
        This also saves changes in the header file to disk.
        """
        file = self._register_file()
        with open(file, 'w') as f:
            json.dump(self.register, f, indent=4)
        return self

    def _register_file(self):
        return os.path.join(self.path, 'dbtree.json') 
    

    def summary(self):
        """Return a summary of the contents of the database.

        Returns:
            dict: Nested dictionary with summary information on the database.
        """
        return register.summary(self.register)
    
    def to_csv(self, csv_file) -> dict:
        """Write a summary of the contents of the database to csv.

        Args:
            path (str): path to the DICOM folder
            csv_file (str): path to the csv file
        """
        register.to_csv(self.register, csv_file)
    

    def print(self):
        """Print the contents of the DICOM folder
        """
        register.print_tree(self.register)
        return self
    
    def patients(self, name=None, contains=None, isin=None):
        """Return a list of patients in the DICOM folder.

        Args:
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
        return register.patients(self.register, self.path, name, contains, isin)
    
    def studies(self, entity=None, desc=None, contains=None, isin=None):
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
        if entity == None:
            entity = self.path
        if isinstance(entity, str):
            studies = []
            for patient in self.patients():
                studies += self.studies(patient, desc, contains, isin)
            return studies
        elif len(entity)==1:
            studies = []
            for patient in self.patients():
                studies += self.studies(patient, desc, contains, isin)
            return studies
        else:
            return register.studies(self.register, entity, desc, contains, isin)
    
    def series(self, entity=None, desc=None, contains=None, isin=None):
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
        if entity == None:
            entity = self.path
        if isinstance(entity, str):
            series = []
            for study in self.studies(entity):
                series += self.series(study, desc, contains, isin)
            return series
        elif len(entity)==1:
            series = []
            for study in self.studies(entity):
                series += self.series(study, desc, contains, isin)
            return series            
        elif len(entity)==2:
            series = []
            for study in self.studies(entity):
                series += self.series(study, desc, contains, isin)
            return series
        else: # path = None (all series) or path = patient (all series in patient)
            return register.series(self.register, entity, desc, contains, isin)


    def volume(self, entity:Union[list, str], dims:list=None, verbose=1) -> vreg.Volume3D:
        """Read volume.

        Args:
            entity (list, str): DICOM series to read
            dims (list, optional): Non-spatial dimensions of the volume. Defaults to None.
            verbose (bool, optional): If set to 1, shows progress bar. Defaults to 1.

        Returns:
            vreg.Volume3D:
        """
        # if isinstance(entity, str): # path to folder
        #     return [self.volume(s, dims) for s in self.series(entity)]
        # if len(entity) < 4: # folder, patient or study
        #     return [self.volume(s, dims) for s in self.series(entity)]
        
        if dims is None:
            dims = []
        elif isinstance(dims, str):
            dims = [dims]
        else:
            dims = list(dims)
        dims = ['SliceLocation'] + dims

        # Read dicom files
        values = [[] for _ in dims]
        volumes = []

        files = register.files(self.register, entity)
        for f in tqdm(files, desc='Reading volume..', disable=(verbose==0)):
            ds = pydicom.dcmread(f)
            values_f = get_values(ds, dims)
            for d in range(len(dims)):
                values[d].append(values_f[d])
            volumes.append(dbdataset.volume(ds))

        # Format coordinates as mesh
        coords = [np.array(v) for v in values]
        coords, inds = dbdicom.utils.arrays.meshvals(coords)

        # Check that all slices have the same coordinates
        if len(dims) > 1:
            # Loop over all coordinates after slice location
            for c in coords[1:]:
                # Loop over all slice locations
                for k in range(1, c.shape[0]):
                    # Coordinate c of slice k
                    if not np.array_equal(c[k,...], c[0,...]):
                        raise ValueError(
                            "Cannot build a single volume. Not all slices "
                            "have the same coordinates."     
                        )

        # Build volumes
        vols = np.array(volumes)
        vols = vols[inds].reshape(coords[0].shape)

        # Infer spacing between slices from slice locations
        # Technically only necessary if SpacingBetweenSlices not set or incorrect
        vols = infer_slice_spacing(vols)

        # Join 2D volumes into 3D volumes
        try:
            vol = vreg.join(vols)
        except ValueError:
            # some vendors define the slice vector as -cross product 
            # of row and column vector. Check if that solves the issue.
            for v in vols.reshape(-1):
                v.affine[:3,2] = -v.affine[:3,2]
                # Then try again
            vol = vreg.join(vols)

        # For multi-dimensional volumes, set dimensions and coordinates
        if vol.ndim > 3:
            # Coordinates of slice 0
            c0 = [c[0,...] for c in coords[1:]]
            vol.set_coords(c0)
            vol.set_dims(dims[1:])
        return vol

    # Obsolete API - phase out
    def volumes_2d(self, *args, **kwargs):
        return self.slices(*args, **kwargs)
    

    def slices(self, entity:Union[list, str], dims:list=None, verbose=1) -> list:
        """Read 2D volumes from the series

        Args:
            entity (list, str): DICOM series to read
            dims (list, optional): Non-spatial dimensions of the volume. Defaults to None.
            verbose (bool, optional): If set to 1, shows progress bar. Defaults to 1.

        Returns:
            list of vreg.Volume3D
        """
        # if isinstance(entity, str): # path to folder
        #     return [self.volume(s, dims) for s in self.series(entity)]
        # if len(entity) < 4: # folder, patient or study
        #     return [self.volume(s, dims) for s in self.series(entity)]
        
        if dims is None:
            dims = []
        elif isinstance(dims, str):
            dims = [dims]
        else:
            dims = list(dims)
        dims = ['SliceLocation'] + dims

        # Read dicom files
        values = {}
        volumes = {}

        files = register.files(self.register, entity)
        for f in tqdm(files, desc='Reading volume..', disable=(verbose==0)):
            ds = pydicom.dcmread(f)
            values_f = get_values(ds, dims)
            vol = dbdataset.volume(ds, multislice=True)
            slice_loc = values_f[0]
            if slice_loc in volumes:
                volumes[slice_loc].append(vol)
                for d in range(len(dims)):
                    values[slice_loc][d].append(values_f[d])
            else:
                volumes[slice_loc] = [vol]
                values[slice_loc] = [[values_f[d]] for d in range(len(dims))]

        # Build a volume for each slice location
        volumes_2d = []
        for slice_loc in volumes.keys():
            vols_list = volumes[slice_loc]

            if values == {}:
                if len(vols_list) > 1:
                    raise ValueError(
                        "Cannot return a 2D volume - multiple slices at the same "
                        "location. \n Use InstanceNumber or another suitable DICOM "
                        "attribute as dimension to sort them.")
                volumes_2d.append(vols_list[0])
                continue

            # Sort by coordinata values
            vals_list = values[slice_loc]

            # Format coordinates as mesh
            coords = [np.array(v) for v in vals_list]
            coords, inds = dbdicom.utils.arrays.meshvals(coords)

            # Check that all slices have the same coordinates
            if len(dims) > 1:
                # Loop over all coordinates after slice location
                for c in coords[1:]:
                    # Loop over all slice locations
                    for k in range(1, c.shape[0]):
                        # Coordinate c of slice k
                        if not np.array_equal(c[k,...], c[0,...]):
                            raise ValueError(
                                "Cannot build a single volume. Not all slices "
                                "have the same coordinates."     
                            )

            # Build volumes, sort and reshape along the coordinates
            vols = np.array(vols_list)
            vols = vols[inds].reshape(coords[0].shape)
        
            # Join 2D volumes along the extra dimensions
            vol = vreg.join(vols[0,...].reshape((1,) + vols.shape[1:]))

            # For multi-dimensional volumes, set dimensions and coordinates
            if vol.ndim > 3:
                # Coordinates of slice 0
                c0 = [c[0,...] for c in coords[1:]]
                vol.set_coords(c0)
                vol.set_dims(dims[1:])

            volumes_2d.append(vol)
        
        # sort volumes by slice location
        volumes_2d.sort(key=lambda v: dbdataset.affine_to_slice_loc(v.affine))

        return volumes_2d


    def pixel_data(self, series:list, dims:list=None, verbose=1) -> np.ndarray:
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
        vols = self.volumes_2d(series, dims, verbose)
        for v in vols[1:]:
            if v.shape != vols[0].shape:
                raise ValueError(
                    "Cannot return a pixel array because slices have different shapes." \
                    "Instead try using volumes_2d to return a list of 2D volumes."
                )
        slices = [v.values for v in vols]
        pixel_array = np.concatenate(slices, axis=2)
        return pixel_array
        
    

    def values(self, series:list, *attr, dims:list=None, verbose=1) -> Union[dict, tuple]:
        """Read the values of some attributes from a DICOM series

        Args:
            series (list): DICOM series to read. 
            attr (tuple, optional): DICOM attributes to read.
            dims (list, optional): Non-spatial dimensions of the volume. Defaults to None.
            verbose (bool, optional): If set to 1, shows progress bar. Defaults to 1.

        Returns:
            tuple: arrays with values for the attributes.
        """
        # if isinstance(series, str): # path to folder
        #     return [self.values(s, attr, dims) for s in self.series(series)]
        # if len(series) < 4: # folder, patient or study
        #     return [self.values(s, attr, dims) for s in self.series(series)]

        if dims is None:
            dims = ['InstanceNumber']
        elif np.isscalar(dims):
            dims = [dims]
        else:
            dims = list(dims)
            
        # Read dicom files
        coord_values = [[] for _ in dims]
        attr_values = [[] for _ in attr]

        files = register.files(self.register, series)
        for f in tqdm(files, desc='Reading values..', disable=(verbose==0)):
            ds = pydicom.dcmread(f)
            coord_values_f = get_values(ds, dims)
            for d in range(len(dims)):
                coord_values[d].append(coord_values_f[d])
            attr_values_f = get_values(ds, attr)
            for a in range(len(attr)):
                attr_values[a].append(attr_values_f[a])

        # Format coordinates as mesh
        coords = [np.array(v) for v in coord_values]
        coords, inds = dbdicom.utils.arrays.meshvals(coords)

        # Sort values accordingly
        values = [np.array(v) for v in attr_values]
        values = [v[inds].reshape(coords[0].shape) for v in values]

        # Return values
        if len(values) == 1:
            return values[0]
        else:
            return tuple(values)


    def write_volume(
            self, vol:Union[vreg.Volume3D, tuple], series:list, 
            ref:list=None, append=False, verbose=1, **kwargs,
        ):
        """Write a vreg.Volume3D to a DICOM series

        Args:
            vol (vreg.Volume3D): Volume to write to the series.
            series (list): DICOM series to read
            ref (list): Reference series
            append (bool): by default write_volume will only write to a new series, 
               and raise an error when attempting to write to an existing series. 
               To overrule this behaviour and add the volume to an existing series, set append to True. 
               Default is False.
            verbose (bool): if set to 1, a progress bar is shown
            kwargs: Keyword-value pairs to be set on the fly
        """
        series_full_name = full_name(series)
        if series_full_name in self.series():
            if not append:
                raise ValueError(f"Series {series_full_name[-1]} already exists in study {series_full_name[-2]}.")

        if isinstance(vol, tuple):
            vol = vreg.volume(vol[0], vol[1])
        if ref is None:
            ds = dbdataset.new_dataset('MRImage')
            #ds = dbdataset.new_dataset('ParametricMap')
        else:
            if ref[0] == series[0]:
                ref_mgr = self
            else:
                ref_mgr = DataBaseDicom(ref[0])
            files = register.files(ref_mgr.register, ref)
            ref_mgr.close()
            ds = pydicom.dcmread(files[0]) 

        # Get the attributes of the destination series
        attr = self._series_attributes(series)
        n = self._max_instance_number(attr['SeriesInstanceUID'])

        if vol.ndim==3:
            slices = vol.split()
            for i, sl in tqdm(enumerate(slices), desc='Writing volume..', disable=verbose==0):
                dbdataset.set_volume(ds, sl)
                if kwargs != {}:
                    set_values(ds, list(kwargs.keys()), list(kwargs.values()))
                self._write_dataset(ds, attr, n + 1 + i)
        else:
            i=0
            vols = vol.separate().reshape(-1)
            for vt in tqdm(vols, desc='Writing volume..', disable=verbose==0):
                slices = vt.split()
                for sl in slices:
                    dbdataset.set_volume(ds, sl)
                    sl_coords = [c.ravel()[0] for c in sl.coords]
                    set_value(ds, sl.dims, sl_coords)
                    if kwargs != {}:
                        set_values(list(kwargs.keys()), list(kwargs.values()))
                    self._write_dataset(ds, attr, n + 1 + i)
                    i+=1
        return self
    

    def edit(
            self, series:list, new_values:dict, dims:list=None, verbose=1,
        ):
        """Edit attribute values in a new DICOM series

        Args:
            series (list): DICOM series to edit
            new_values (dict): dictionary with attribute: value pairs to write to the series
            dims (list, optional): Non-spatial dimensions of the volume. Defaults to None.
            verbose (bool, optional): If set to 1, shows progress bar. Defaults to 1.
        """

        if dims is None:
            dims = ['InstanceNumber']
        elif np.isscalar(dims):
            dims = [dims]
        else:
            dims = list(dims)
            
        # Check that all values have the correct nr of elements
        files = register.files(self.register, series)
        for a in new_values.values():
            if np.isscalar(a):
                pass
            elif np.array(a).size != len(files):
                raise ValueError(
                    f"Incorrect value lengths. All values need to have {len(files)} elements"
                )

        # Read dicom files to sort them
        coord_values = [[] for _ in dims]
        for f in tqdm(files, desc='Sorting series..', disable=(verbose==0)):
            ds = pydicom.dcmread(f)
            coord_values_f = get_values(ds, dims)
            for d in range(len(dims)):
                coord_values[d].append(coord_values_f[d])

        # Format coordinates as mesh
        coords = [np.array(v) for v in coord_values]
        coords, inds = dbdicom.utils.arrays.meshvals(coords)

        # Sort files accordingly
        files = np.array(files)[inds]

        # Now edit and write the files
        attr = self._series_attributes(series)
        n = self._max_instance_number(attr['SeriesInstanceUID'])

        # Drop existing attributes if they are edited
        attr = {a:attr[a] for a in attr if a not in new_values}

        # List instances to be edited
        to_drop = register.index(self.register, series)
         
        # Write the instances
        tags = list(new_values.keys())
        for i, f in tqdm(enumerate(files), desc='Writing values..', disable=(verbose==0)):
            ds = pydicom.dcmread(f)
            values = []
            for a in new_values.values():
                if np.isscalar(a):
                    values.append(a)
                else:
                    values.append(np.array(a).reshape(-1)[i])
            set_values(ds, tags, values)
            self._write_dataset(ds, attr, n + 1 + i)

        # Delete the originals files
        register.drop(self.register, to_drop)
        [os.remove(os.path.join(self.path, idx)) for idx in to_drop]

        return self


    def to_nifti(self, series:list, file:str, dims=None, verbose=1):
        """Save a DICOM series in nifti format.

        Args:
            series (list): DICOM series to read
            file (str): file path of the nifti file.
            dims (list, optional): Non-spatial dimensions of the volume. 
                Defaults to None.
            verbose (bool, optional): If set to 1, shows progress bar. Defaults to 1.
            
        """
        vol = self.volume(series, dims, verbose)
        vreg.write_nifti(vol, file)
        return self

    def from_nifti(self, file:str, series:list, ref:list=None):
        """Create a DICOM series from a nifti file.

        Args:
            file (str): file path of the nifti file.
            series (list): DICOM series to create
            ref (list): DICOM series to use as template.
        """
        vol = vreg.read_nifti(file)
        self.write_volume(vol, series, ref)
        return self
    

        


        

    def files(self, entity:list) -> list:
        """Read the files in a DICOM entity

        Args:
            entity (list or str): DICOM entity to read. This can 
                be a path to a folder containing DICOM files, or a 
                patient or study to read all series in that patient or 
                study. 

        Returns:
            list: list of valid dicom files.
        """
        if isinstance(entity, str): # path to folder
            files = []
            for s in self.series(entity):
                files += self.files(s)
            return files
        if len(entity) < 4: # folder, patient or study
            files = []
            for s in self.series(entity):
                files += self.files(s)
            return files

        return register.files(self.register, entity)

    
    
    def unique(self, pars:list, entity:list) -> dict:
        """Return a list of unique values for a DICOM entity

        Args:
            pars (list, str/tuple): attribute or attributes to return.
            entity (list): DICOM entity to search (Patient, Study or Series)

        Returns:
            dict: if a pars is a list, this returns a dictionary with 
            unique values for each attribute. If pars is a scalar 
            this returnes a list of values.
        """
        if not isinstance(pars, list):
            single=True
            pars = [pars]
        else:
            single=False

        v = self._values(pars, entity)

        # Return a list with unique values for each attribute
        values = []
        for a in range(v.shape[1]):
            va = v[:,a]
            # Remove None values
            va = va[[x is not None for x in va]]
            va = list(va)
            # Get unique values and sort
            va = [x for i, x in enumerate(va) if i==va.index(x)]
            try: 
                va.sort()
            except:
                pass
            values.append(va)

        if single:
            return values[0]
        else:
            return {p: values[i] for i, p in enumerate(pars)} 
        
    
    def copy(self, from_entity, to_entity=None):
        """Copy a DICOM  entity (patient, study or series)

        Args:
            from_entity (list): entity to copy
            to_entity (list, optional): entity after copying. If this is not 
                provided, a copy will be made in the same study and returned

        Returns:
            entity: the copied entity. If th to_entity is provided, this is 
            returned.
        """
        if len(from_entity) == 4:
            if to_entity is None:
                to_entity = deepcopy(from_entity)
                if isinstance(to_entity[-1], tuple):
                    to_entity[-1] = (to_entity[-1][0] + '_copy', 0)
                else:
                    to_entity[-1] = (to_entity[-1] + '_copy', 0)
                while to_entity in self.series():
                    to_entity[-1] = (to_entity[-1][0], to_entity[-1][1] + 1)
            if len(to_entity) != 4:
                raise ValueError(
                    f"Cannot copy series {from_entity} to series {to_entity}. "
                    f"{to_entity} is not a series (needs 4 elements)."
                )
            self._copy_series(from_entity, to_entity)
            return to_entity
        
        if len(from_entity) == 3:
            if to_entity is None:
                to_entity = deepcopy(from_entity)
                if isinstance(to_entity[-1], tuple):
                    to_entity[-1] = (to_entity[-1][0] + '_copy', 0)
                else:
                    to_entity[-1] = (to_entity[-1] + '_copy', 0)
                while to_entity in self.studies():
                    to_entity[-1][1] += 1
            if len(to_entity) != 3:
                raise ValueError(
                    f"Cannot copy study {from_entity} to study {to_entity}. "
                    f"{to_entity} is not a study (needs 3 elements)."
                )
            self._copy_study(from_entity, to_entity)
            return to_entity
        
        if len(from_entity) == 2:
            if to_entity is None:
                to_entity = deepcopy(from_entity)
                to_entity[-1] += '_copy'
                while to_entity in self.patients():
                    to_entity[-1] += '_copy'
            if len(to_entity) != 2:
                raise ValueError(
                    f"Cannot copy patient {from_entity} to patient {to_entity}. "
                    f"{to_entity} is not a patient (needs 2 elements)."
                )                
            self._copy_patient(from_entity, to_entity)
            return to_entity
        
        raise ValueError(
            f"Cannot copy {from_entity} to {to_entity}. "
        )
    
    def move(self, from_entity, to_entity):
        """Move a DICOM entity

        Args:
            entity (list): entity to move
        """
        self.copy(from_entity, to_entity)
        self.delete(from_entity)
        return self
    
    def split_series(self, series:list, attr:Union[str, tuple], key=None) -> list:
        """
        Split a series into multiple series
        
        Args:
            series (list): series to split.
            attr (str or tuple): dicom attribute to split the series by. 
            key (function): split by by key(attr)
        Returns:
            list: list of two-element tuples, where the first element is
            is the value and the second element is the series corresponding to that value.         
        """

        # Find all values of the attr and list files per value
        all_files = register.files(self.register, series)
        files = []
        values = []
        for f in tqdm(all_files, desc=f'Reading {attr}'):
            ds = pydicom.dcmread(f)
            v = get_values(ds, attr)
            if key is not None:
                v = key(v)
            if v in values:
                index = values.index(v)
                files[index].append(f)
            else:
                values.append(v)
                files.append([f])

        # Copy the files for each value (sorted) to new series
        split_series = []
        for index, v in tqdm(enumerate(values), desc='Writing new series'):
            series_desc = series[-1] if isinstance(series, str) else series[-1][0]
            series_desc = clean_folder_name(f'{series_desc}_{attr}_{v}')
            series_v = series[:3] + [(series_desc, 0)]
            self._files_to_series(files[index], series_v)
            split_series.append((v, series_v))
        return split_series


    def _values(self, attributes:list, entity:list):
        # Create a np array v with values for each instance and attribute
        # if set(attributes) <= set(dbdatabase.COLUMNS):
        #     index = register.index(self.register, entity)
        #     v = self.register.loc[index, attributes].values
        # else:
        files = register.files(self.register, entity)
        v = np.empty((len(files), len(attributes)), dtype=object)
        for i, f in enumerate(files):
            ds = pydicom.dcmread(f)
            v[i,:] = get_values(ds, attributes)
        return v

    def _copy_patient(self, from_patient, to_patient):
        from_patient_studies = register.studies(self.register, from_patient)
        for from_study in tqdm(from_patient_studies, desc=f'Copying patient {from_patient[1:]}'):
            # Count the studies with the same description in the target patient
            study_desc = from_study[-1][0]
            if to_patient[0]==from_patient[0]:
                cnt = len(self.studies(to_patient, desc=study_desc))
            else:
                mgr = DataBaseDicom(to_patient[0])
                cnt = len(mgr.studies(to_patient, desc=study_desc))
                mgr.close()    
            # Ensure the copied studies end up in a separate study with the same description
            to_study = to_patient + [(study_desc, cnt)]         
            self._copy_study(from_study, to_study)

    def _copy_study(self, from_study, to_study):
        from_study_series = register.series(self.register, from_study)
        for from_series in tqdm(from_study_series, desc=f'Copying study {from_study[1:]}'):
            # Count the series with the same description in the target study
            series_desc = from_series[-1][0]
            if to_study[0]==from_study[0]:
                cnt = len(self.series(to_study, desc=series_desc))
            else:
                mgr = DataBaseDicom(to_study[0])
                cnt = len(mgr.series(to_study, desc=series_desc))
                mgr.close()
            # Ensure the copied series end up in a separate series with the same description
            to_series = to_study + [(series_desc, cnt)]
            self._copy_series(from_series, to_series)

    def _copy_series(self, from_series, to_series):
        # Get the files to be exported
        from_series_files = register.files(self.register, from_series)
        if to_series[0] == from_series[0]:
            # Copy in the same database
            self._files_to_series(from_series_files, to_series)
        else:
            # Copy to another database
            mgr = DataBaseDicom(to_series[0])
            mgr._files_to_series(from_series_files, to_series)
            mgr.close()


    def _files_to_series(self, files, to_series):

        # Get the attributes of the destination series
        attr = self._series_attributes(to_series)
        n = self._max_instance_number(attr['SeriesInstanceUID'])
        
        # Copy the files to the new series 
        for i, f in tqdm(enumerate(files), total=len(files), desc=f'Copying series {to_series[1:]}'):
            # Read dataset and assign new properties
            ds = pydicom.dcmread(f)
            self._write_dataset(ds, attr, n + 1 + i)

    def _max_study_id(self, patient_id):
        for pt in self.register:
            if pt['PatientID'] == patient_id:
                # Find the largest integer StudyID
                n = []
                for st in pt['studies']:
                    try:
                        n.append(int(st['StudyID']))
                    except:
                        pass
                if n == []:
                    return 0
                else:
                    return int(np.amax(n))
        return 0
    
    def _max_series_number(self, study_uid):
        for pt in self.register:
            for st in pt['studies']:
                if st['StudyInstanceUID'] == study_uid:
                    n = [sr['SeriesNumber'] for sr in st['series']]
                    return int(np.amax(n))
        return 0

    def _max_instance_number(self, series_uid):
        for pt in self.register:
            for st in pt['studies']:
                for sr in st['series']:
                    if sr['SeriesInstanceUID'] == series_uid:
                        n = list(sr['instances'].keys())
                        return int(np.amax([int(i) for i in n]))
        return 0

    # def _attributes(self, entity):
    #     if len(entity)==4:
    #         return self._series_attributes(entity)
    #     if len(entity)==3:
    #         return self._study_attributes(entity)
    #     if len(entity)==2:
    #         return self._patient_attributes(entity)       


    def _patient_attributes(self, patient):
        try:
            # If the patient exists and has files, read from file
            files = register.files(self.register, patient)
            attr = const.PATIENT_MODULE
            ds = pydicom.dcmread(files[0])
            vals = get_values(ds, attr)
        except:
            # If the patient does not exist, generate values
            if patient in self.patients():
                raise ValueError(
                    f"Cannot create patient with id {patient[1]}."
                    f"The ID is already taken. Please provide a unique ID."
                )
            attr = ['PatientID', 'PatientName']
            vals = [patient[1], 'Anonymous']
        return {attr[i]:vals[i] for i in range(len(attr)) if vals[i] is not None}


    def _study_attributes(self, study):
        patient_attr = self._patient_attributes(study[:2])
        try:
            # If the study exists and has files, read from file
            files = register.files(self.register, study)
            attr = const.STUDY_MODULE
            ds = pydicom.dcmread(files[0])
            vals = get_values(ds, attr)
        except register.AmbiguousError as e:
            raise register.AmbiguousError(e)
        except:
            # If the study does not exist or is empty, generate values
            if study[:-1] not in self.patients():
                study_id = 1
            else:
                study_id = 1 + self._max_study_id(study[1])
            attr = ['StudyInstanceUID', 'StudyDescription', 'StudyID']
            study_uid = pydicom.uid.generate_uid()
            study_desc = study[-1] if isinstance(study[-1], str) else study[-1][0]
            #study_date = datetime.today().strftime('%Y%m%d')
            vals = [study_uid, study_desc, str(study_id)]
        return patient_attr | {attr[i]:vals[i] for i in range(len(attr)) if vals[i] is not None}


    def _series_attributes(self, series):
        study_attr = self._study_attributes(series[:3])
        try:
            # If the series exists and has files, read from file
            files = register.files(self.register, series)
            attr = const.SERIES_MODULE
            ds = pydicom.dcmread(files[0])
            vals = get_values(ds, attr)
        except register.AmbiguousError as e:
            raise register.AmbiguousError(e)
        except:
            # If the series does not exist or is empty, generate values
            try:
                study_uid = register.study_uid(self.register, series[:-1])
            except:
                series_number = 1
            else:
                series_number = 1 + self._max_series_number(study_uid)
            attr = ['SeriesInstanceUID', 'SeriesDescription', 'SeriesNumber']
            series_uid = pydicom.uid.generate_uid()
            series_desc = series[-1] if isinstance(series[-1], str) else series[-1][0]
            vals = [series_uid, series_desc, int(series_number)]
        return study_attr | {attr[i]:vals[i] for i in range(len(attr)) if vals[i] is not None}

        
    def _write_dataset(self, ds:Dataset, attr:dict, instance_nr:int):
        # Set new attributes 
        attr['SOPInstanceUID'] = pydicom.uid.generate_uid()
        attr['InstanceNumber'] = str(instance_nr)
        set_values(ds, list(attr.keys()), list(attr.values()))
        # Save results in a new file
        rel_dir = os.path.join(
            f"Patient__{attr['PatientID']}", 
            f"Study__{attr['StudyID']}__{attr['StudyDescription']}", 
            f"Series__{attr['SeriesNumber']}__{attr['SeriesDescription']}",
        )
        os.makedirs(os.path.join(self.path, rel_dir), exist_ok=True)
        rel_path = os.path.join(rel_dir, pydicom.uid.generate_uid() + '.dcm')
        dbdataset.write(ds, os.path.join(self.path, rel_path))
        # Add an entry in the register
        register.add_instance(self.register, attr, rel_path)


    def archive(self, archive_path):
        # TODO add flat=True option for zipping at patient level
        for pt in tqdm(self.register, desc='Archiving '):
            for st in pt['studies']:
                zip_dir = os.path.join(
                    archive_path,
                    f"Patient__{pt['PatientID']}", 
                    f"Study__{st['StudyID']}__{st['StudyDescription']}", 
                )
                os.makedirs(zip_dir, exist_ok=True)
                for sr in st['series']:
                    zip_file = os.path.join(
                        zip_dir, 
                        f"Series__{sr['SeriesNumber']}__{sr['SeriesDescription']}.zip",
                    )
                    if os.path.exists(zip_file):
                        continue
                    try:
                        with zipfile.ZipFile(zip_file, 'w') as zipf:
                            for rel_path in sr['instances'].values():
                                file = os.path.join(self.path, rel_path)
                                zipf.write(file, arcname=os.path.basename(file))
                    except Exception as e:
                        raise RuntimeError(
                            f"Error extracting series {sr['SeriesDescription']} "
                            f"in study {st['StudyDescription']} of patient {pt['PatientID']}."
                        )



def full_name(entity):

    if len(entity)==3: # study
        if isinstance(entity[-1], tuple):
            return entity
        else:
            full_name_study = deepcopy(entity)
            full_name_study[-1] = (full_name_study[-1], 0)
            return full_name_study
        
    elif len(entity)==4: # series
        full_name_study = full_name(entity[:3])
        series = full_name_study + [entity[-1]]
        if isinstance(series[-1], tuple):
            return series
        else:
            full_name_series = deepcopy(series)
            full_name_series[-1] = (full_name_series[-1], 0)
            return full_name_series
    else:
        return entity


def clean_folder_name(name, replacement="", max_length=255):
    # Strip leading/trailing whitespace
    name = name.strip()

    # Replace invalid characters (Windows, macOS, Linux-safe)
    illegal_chars = r'[<>:"/\\|?*\[\]\x00-\x1F\x7F]'
    name = re.sub(illegal_chars, replacement, name)

    # Replace reserved Windows names
    reserved = {
        "CON", "PRN", "AUX", "NUL",
        *(f"COM{i}" for i in range(1, 10)),
        *(f"LPT{i}" for i in range(1, 10))
    }
    name_upper = name.upper().split(".")[0]  # Just base name
    if name_upper in reserved:
        name = f"{name}_folder"

    # Truncate to max length (common max: 255 bytes)
    return name[:max_length] or "folder"



def remove_empty_folders(path):
    """
    Removes all empty subfolders from a given directory.

    This function walks through the directory tree from the bottom up.
    This is crucial because it allows child directories to be removed before
    their parents, potentially making the parent directory empty and
    eligible for removal in the same pass.

    Args:
        path (str): The absolute or relative path to the directory to scan.
    """
    # Walk the directory tree in a bottom-up manner (topdown=False)
    for dirpath, dirnames, filenames in os.walk(path, topdown=False):
        # A directory is considered empty if it has no subdirectories and no files
        if not dirnames and not filenames:
            try:
                shutil.rmtree(dirpath)
            except OSError as e:
                # This might happen due to permissions issues
                print(f"Error removing {dirpath}: {e}")



def infer_slice_spacing(vols):
    # In case spacing between slices is not (correctly) encoded in 
    # DICOM it can be inferred from the slice locations.

    shape = vols.shape
    vols = vols.reshape((shape[0], -1))
    slice_spacing = np.zeros(vols.shape[-1])

    for d in range(vols.shape[-1]):

        # For single slice volumes there is nothing to do
        if vols[:,d].shape[0]==1:
            continue

        # Get a normal slice vector from the first volume.
        mat = vols[0,d].affine[:3,:3]
        normal = mat[:,2]/np.linalg.norm(mat[:,2])

        # Get slice locations by projection on the normal.
        pos = [v.affine[:3,3] for v in vols[:,d]]
        slice_loc = [np.dot(p, normal) for p in pos]

        # Sort slice locations and take consecutive differences.
        slice_loc = np.sort(slice_loc)
        distances = slice_loc[1:] - slice_loc[:-1]

        # Round to 10 micrometer and check if unique
        distances = np.around(distances, 2)
        slice_spacing_d = np.unique(distances)

        # Check if slice spacings are unique - otherwise this is not a volume
        if len(slice_spacing_d) > 1:
            raise ValueError(
                'Cannot build a volume - spacings between slices are not unique.'
            )
        else:
            slice_spacing_d= slice_spacing_d[0]
        
        # Set correct slice spacing in all volumes
        for v in vols[:,d]:
            v.affine[:3,2] = normal * abs(slice_spacing_d)

        slice_spacing[d] = slice_spacing_d

    # Check slice_spacing is the same across dimensions
    # Not sure if this is possible as volumes are sorted by slice location
    slice_spacing = np.unique(slice_spacing)
    if len(slice_spacing) > 1:
        raise ValueError(
            'Cannot build a volume - spacings between slices are not unique.'
        )    

    return vols.reshape(shape)


import os
from pathlib import Path
import platform

from tqdm import tqdm
import numpy as np
import pandas as pd
import pydicom

import dbdicom.utils.dcm4che as dcm4che
from dbdicom.utils.pydicom_dataset import get_values


COLUMNS = [   
    # Identifiers (unique)
    'PatientID', 
    'StudyInstanceUID', 
    'SeriesInstanceUID', 
    'SOPInstanceUID', 
    # Human-readable identifiers (not unique)
    'PatientName', 
    'StudyDescription', 
    'StudyDate', 
    'StudyID',
    'SeriesDescription', 
    'SeriesNumber', 
    'InstanceNumber', 
]

def read(path):
    files = _all_files(path)
    tags = COLUMNS + ['NumberOfFrames'] # + ['SOPClassUID']
    array = []
    for i, file in tqdm(enumerate(files), total=len(files), desc='Reading DICOM folder'):
        try:
            ds = pydicom.dcmread(file, force=True, specific_tags=tags+['Rows'])
        except:
            pass
        else:
            if isinstance(ds, pydicom.dataset.FileDataset):
                if 'TransferSyntaxUID' in ds.file_meta:
                    if not 'Rows' in ds: # Image only
                        continue
                    row = get_values(ds, tags)
                    index = os.path.relpath(file, path)
                    p = Path(index)
                    parts = list(p.parts)
                    row = [parts] + row
                    array.append(row)
    df = pd.DataFrame(array, columns = ['rel_path'] + tags)
    df = _multiframe_to_singleframe(path, df) # needs updating and testing
    dbtree = _tree(df)
    return dbtree


def _all_files(path):
    files = [item.path for item in _scan_tree(path) if item.is_file()]
    # Windows has maximum path length of 260 - ignore any files that are longer
    if platform.system() == 'Windows':
        files = [f for f in files if len(f) <= 260]
    return files


def _scan_tree(directory):
    """Helper function: yield DirEntry objects for the directory."""

    for entry in os.scandir(directory):
        if entry.is_dir(follow_symlinks=False):
            yield from _scan_tree(entry.path)
        else:
            yield entry


def _multiframe_to_singleframe(path, df):
    """Converts all multiframe files in the folder into single-frame files.
    
    Reads all the multi-frame files in the folder,
    converts them to singleframe files, and delete the original multiframe file.
    """
    singleframe = df.NumberOfFrames.isnull() 
    multiframe = singleframe == False
    nr_multiframe = multiframe.sum()
    if nr_multiframe != 0: 
        raise ValueError(
            "dbdicom currently does not support multiframe data."
            "Please remove them from the database and try again."
            )
        for relpath in tqdm(df[multiframe].index.values, desc="Converting multiframe file " + relpath):
            filepath = [path] + [relpath]
            filepath = Path(*filepath)
            singleframe_files = dcm4che.split_multiframe(str(filepath)) 
            if singleframe_files != []:            
                # add the single frame files to the dataframe
                dfnew = read(singleframe_files, df.columns, path) # This needs fixing
                df = pd.concat([df, dfnew])
                # delete the original multiframe 
                os.remove(filepath)
            # drop the file also if the conversion has failed
            df.drop(index=relpath, inplace=True)
    df.drop('NumberOfFrames', axis=1, inplace=True)
    return df


def _tree(df):
    # A human-readable summary tree
    # TODO: Add version number

    df.sort_values(['PatientID','StudyInstanceUID','SeriesNumber'], inplace=True)
    df = df.fillna('None')
    summary = []

    for uid_patient in df.PatientID.unique():
        df_patient = df[df.PatientID == uid_patient]
        patient_name = df_patient.PatientName.values[0]
        patient = {
            'PatientName': patient_name,
            'PatientID': uid_patient,
            'studies': [],
        }
        summary.append(patient)
        for uid_study in df_patient.StudyInstanceUID.unique():
            df_study = df_patient[df_patient.StudyInstanceUID == uid_study]
            study_desc = df_study.StudyDescription.values[0]
            study_id = df_study.StudyID.values[0]
            study_date = df_study.StudyDate.values[0]
            study = {
                'StudyDescription': study_desc,
                'StudyDate': study_date,
                'StudyID': study_id,
                'StudyInstanceUID': uid_study,
                'series': [],
            }
            patient['studies'].append(study)
            for uid_sery in df_study.SeriesInstanceUID.unique():
                df_series = df_study[df_study.SeriesInstanceUID == uid_sery]
                series_desc = df_series.SeriesDescription.values[0]
                series_nr = int(df_series.SeriesNumber.values[0])
                series = {
                    'SeriesNumber': series_nr,
                    'SeriesDescription': series_desc,
                    'SeriesInstanceUID': uid_sery,
                    'instances': {},
                }
                study['series'].append(series)
                for uid_instance in df_series.SOPInstanceUID.unique():
                    df_instance = df_series[df_series.SOPInstanceUID == uid_instance]
                    instance_nr = int(df_instance.InstanceNumber.values[0])
                    relpath = df_instance.rel_path.values[0]
                    series['instances'][instance_nr]=relpath

    return summary
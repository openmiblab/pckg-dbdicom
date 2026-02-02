import os
from pathlib import Path
import vreg
import numpy as np


def _parse_patient_dir(dir):
    # Patient__2128_008
    patient_id = dir.name[9:]
    return patient_id

def _parse_study_dir(dir):
    # Study__1__Baseline
    nr_desc = dir.name[7:].split('__')
    study_id = nr_desc[0] # 1 
    study_desc = "".join(nr_desc[1:]) # Baseline
    return int(study_id), str(study_desc)

def _parse_series_file(file:Path):
    # Series__1__normalized_kidney_right.npz
    filename = file.name.split('.')
    filename = "".join(filename[:-1]) # Series__1__normalized_kidney_right
    nr_desc = filename[8:].split('__')  
    series_nr = nr_desc[0] # 1 
    series_desc = "".join(nr_desc[1:]) # normalized_kidney_right
    return int(series_nr), str(series_desc)

def _studies(database, patient_id, study_desc): # all studies with given desc
    patient_dir = os.path.join(
        database, 
        f"Patient__{patient_id}", 
    )
    if not os.path.exists(patient_dir):
        return []
    study_dirs = [s for s in Path(patient_dir).iterdir() if s.is_dir()]
    study_dirs = [s for s in study_dirs if study_desc == _parse_study_dir(s)[1]]
    return study_dirs

def _series(database, patient_id, study_desc, study_id, series_desc):
    study_dir = os.path.join(
        database, 
        f"Patient__{patient_id}", 
        f"Study__{study_id}__{study_desc}",
    )
    if not os.path.exists(study_dir):
        return []
    series_files = [s for s in Path(study_dir).iterdir() if s.is_file()]
    series_files = [s for s in series_files if series_desc == _parse_series_file(s)[1]]
    return series_files

def _file_for_reading(series):
    database = series[0]
    patient = series[1]
    study_desc, study_idx = series[2][0], series[2][1]
    series_desc, series_idx = series[3][0], series[3][1]

    dir = Path(database)

    # Get all subfolders that include the patient_id
    patient_dirs = [p for p in dir.iterdir() if p.is_dir()]
    patient_dirs = [p for p in patient_dirs if patient == _parse_patient_dir(p)] 
    if len(patient_dirs) == 0:
        raise ValueError(f"No patients {patient}")
    if len(patient_dirs) > 1:
        raise ValueError(f"Multiple patients {patient}")
    patient_dir = patient_dirs[0]

    # Get all subfolders that include study_desc
    study_dirs = [s for s in patient_dir.iterdir() if s.is_dir()]
    study_dirs = [s for s in study_dirs if study_desc == _parse_study_dir(s)[1]]
    if len(study_dirs) == 0:
        raise ValueError(f"No studies {study_desc} in patient{str(patient_dir)}")
    if len(study_dirs) > study_idx + 1:
        raise ValueError(f"Only {len(study_dirs)} studies {study_desc} in patient{str(patient_dir)}")
    study_dir = study_dirs[study_idx]

    # Get all files that include series_desc
    series_files = [s for s in study_dir.iterdir() if s.is_file()]
    series_files = [s for s in series_files if series_desc == _parse_series_file(s)[1]]
    if len(series_files) == 0:
        raise ValueError(f"No series {series_desc} in {study_desc} of patient{str(patient_dir)}")
    if len(series_files) > series_idx + 1:
        raise ValueError(f"Only {len(series_files)} series {series_desc} in study {study_desc} in patient{str(patient_dir)}")
    series_file = series_files[series_idx]

    return str(series_file)

def _max_study_id(database, patient_id):
    patient_dir = os.path.join(
        database, 
        f"Patient__{patient_id}", 
    )  
    if not os.path.exists(patient_dir):
        return 0
    study_dirs = [s for s in Path(patient_dir).iterdir() if s.is_dir()]  
    study_ids = [_parse_study_dir(s)[0] for s in study_dirs] 
    return np.max(study_ids)

def _max_series_number(study_dir):
    if not os.path.exists(study_dir):
        return 0
    series_files = [s for s in Path(study_dir).iterdir() if s.is_file()]
    series_nrs = [_parse_series_file(s)[0] for s in series_files]
    return np.max(series_nrs)

def _file_for_writing(series):
    studies_list = _studies(series[0], series[1], series[2][0])
    study_idx = series[2][1]
    if studies_list == []:
        study_id = 1 + _max_study_id(series[0], series[1])
        study_dir = os.path.join(
            series[0], 
            f"Patient__{series[1]}", 
            f"Study__{study_id}__{series[2][0]}",
        )        
    elif study_idx > len(studies_list) - 1:
        max_study_id, _ = _parse_study_dir(studies_list[-1])
        study_id = 1 + max_study_id
        study_dir = os.path.join(
            series[0], 
            f"Patient__{series[1]}", 
            f"Study__{study_id}__{series[2][0]}",
        )
    else:
        study_dir = studies_list[study_idx]
        study_id, _ = _parse_study_dir(study_dir)

    series_list = _series(series[0], series[1], series[2][0], study_id, series[3][0])
    series_idx = series[3][1]
    if series_list == []:
        series_nr = 1 + _max_series_number(study_dir)
        npz_file = os.path.join(
            series[0], 
            f"Patient__{series[1]}", 
            f"Study__{study_id}__{series[2][0]}",
            f"Series__{series_nr}__{series[3][0]}.npz"
        )
    elif series_idx > len(series_list) - 1:
        max_series_nr, _ = _parse_series_file(series_list[-1])
        series_nr = 1 + max_series_nr
        npz_file = os.path.join(
            series[0], 
            f"Patient__{series[1]}", 
            f"Study__{study_id}__{series[2][0]}",
            f"Series__{series_nr}__{series[3][0]}.npz"
        )
    else:
        npz_file = series_list[series_idx]
    return npz_file


def patients(database:str):
    dir = Path(database)
    patient_dirs = [p for p in dir.iterdir() if p.is_dir()]   
    return [[database, _parse_patient_dir(d)] for d in patient_dirs]

def studies(database:str):
    studies_list = []
    dir = Path(database)
    # Loop over patient directories
    patient_dirs = [p for p in dir.iterdir() if p.is_dir()]
    for patient_dir in patient_dirs:
        patient_id = _parse_patient_dir(patient_dir)
        # Study directores of the patient
        study_dirs = [p for p in patient_dir.iterdir() if p.is_dir()]
        # Sort by study id
        study_dirs = sorted(study_dirs, key=lambda d: _parse_study_dir(d)[0])
        # Read unique descriptions
        study_descs = np.unique([_parse_study_dir(d)[1] for d in study_dirs])
        # Build a list of studies for each description
        for study_desc in study_descs:
            study_dirs_desc = [d for d in study_dirs if _parse_study_dir(d)[1]==study_desc]
            for study_idx, study_dir in enumerate(study_dirs_desc):
                study = [database, patient_id, (str(study_desc), study_idx)]
                studies_list.append(study)
    return studies_list

def series(database:str):
    series_list = []
    dir = Path(database)
    # Loop over patient directories
    patient_dirs = [p for p in dir.iterdir() if p.is_dir()]
    for patient_dir in patient_dirs:
        patient_id = _parse_patient_dir(patient_dir)
        # Study directores of the patient
        study_dirs = [p for p in patient_dir.iterdir() if p.is_dir()]
        # Sort by study id
        study_dirs = sorted(study_dirs, key=lambda d: _parse_study_dir(d)[0])
        # Read unique descriptions
        study_descs = np.unique([_parse_study_dir(d)[1] for d in study_dirs])
        # Build a list of studies for each description
        for study_desc in study_descs:
            study_dirs_desc = [d for d in study_dirs if _parse_study_dir(d)[1]==study_desc]
            # Loop over studies with the same description
            for study_idx, study_dir in enumerate(study_dirs_desc):
                # Get the series files for the study
                series_files = [s for s in study_dir.iterdir() if s.is_file()]
                # Sort by series number
                series_files = sorted(series_files, key=lambda f: _parse_series_file(f)[0])
                # Build a list of unique series descriptions
                series_descs = np.unique([_parse_series_file(f)[1] for f in series_files])
                # Build a list of series for each description
                for series_desc in series_descs:
                    series_files_desc = [f for f in series_files if _parse_series_file(f)[1]==series_desc]
                    for series_idx, series_file in enumerate(series_files_desc):
                        series = [database, patient_id, (str(study_desc), study_idx), (str(series_desc), series_idx)]
                        series_list.append(series)
    return series_list


def identifiers(file):
    filepath = Path(file)
    series_nr, series_desc = _parse_series_file(filepath)
    study_id, study_desc = _parse_study_dir(filepath.parent)
    patient_id = _parse_patient_dir(filepath.parent.parent)
    database = filepath.parent.parent.parent
    return {
        'database': str(database),
        'PatientID': patient_id, 
        'StudyDescription': study_desc, 
        'SeriesDescription': series_desc,
        'StudyID': study_id,
        'SeriesNumber': series_nr,
    }

def files(database, PatientID=None, StudyDescription=None, SeriesDescription=None):
    # Get all files
    files = [f for f in Path(database).rglob("*") if f.is_file()]

    # Filter by patient id
    if PatientID is not None:
        if not isinstance(PatientID, list):
            PatientID = [PatientID]
        files = [f for f in files if _parse_patient_dir(f.parent.parent) in PatientID]

    # Filter by study desc
    if StudyDescription is not None:
        if not isinstance(StudyDescription, list):
            StudyDescription = [StudyDescription]
        files = [f for f in files if _parse_study_dir(f.parent)[1] in StudyDescription]

    # Filter by series desc
    if SeriesDescription is not None:
        if not isinstance(SeriesDescription, list):
            SeriesDescription = [SeriesDescription]
        files = [f for f in files if _parse_series_file(f)[1] in SeriesDescription]

    return [str(f) for f in files]


def file(series):
    return _file_for_reading(series)
 
def write_volume(vol, series):
    npz_file = _file_for_writing(series)
    vol.write_npz(npz_file)

def volume(series) -> vreg.Volume3D:
    npz_file = _file_for_reading(series)
    return vreg.read_npz(npz_file)

def exists(series):
    try:
        _file_for_reading(series)
    except:
        return False
    else:
        return True

def filepath(
        dir, 
        patient_id, 
        study_desc, 
        series_desc, 
        study_id=0, 
        series_nr=0,
) -> str:
    return os.path.join(
        dir, 
        f"Patient__{patient_id}", 
        f"Study__{study_id}__{study_desc}",
        f"Series__{series_nr}__{series_desc}.npz"
    )
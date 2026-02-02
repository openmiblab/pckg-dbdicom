import os
from pathlib import Path
import csv
import numpy as np


class AmbiguousError(Exception):
    pass


def _sort_study(st):
    if st['StudyID'] is not None:
        return st['StudyID']
    # Need a backup as StudyID is not required in DICOM
    return st['StudyInstanceUID']


def add_instance(dbtree:list, attr, rel_path):
    
    # Get patient and create if needed
    pts = [pt for pt in sorted(dbtree, key=lambda pt: pt['PatientID']) if pt['PatientID']==attr['PatientID']]
    if pts==[]:
        pt = {
            'PatientName': attr['PatientName'],
            'PatientID': attr['PatientID'],
            'studies': [],
        }
        dbtree.append(pt)
    else:
        pt = pts[0]
    
    # Get study and create if needed
    sts = [st for st in sorted(pt['studies'], key=_sort_study) if st['StudyInstanceUID']==attr['StudyInstanceUID']]
    if sts==[]:
        st = {
            'StudyDescription': attr['StudyDescription'],
            'StudyID': attr['StudyID'],
            'StudyInstanceUID': attr['StudyInstanceUID'],
            'series': [],
        }
        pt['studies'].append(st)
    else:
        st = sts[0]

    # Get series and create if needed
    srs = [sr for sr in sorted(st['series'], key=lambda sr: sr['SeriesNumber']) if sr['SeriesInstanceUID']==attr['SeriesInstanceUID']]
    if srs==[]:
        sr = {
            'SeriesNumber': attr['SeriesNumber'],
            'SeriesDescription': attr['SeriesDescription'],
            'SeriesInstanceUID': attr['SeriesInstanceUID'],
            'instances': {},
        }
        st['series'].append(sr)
    else:
        sr = srs[0]

    # Add instance
    sr['instances'][attr['InstanceNumber']] = rel_path

    return dbtree


def max_study_id(dbtree, patient_id):
    for pt in dbtree:
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

def max_series_number(dbtree, study_uid):
    for pt in dbtree:
        for st in pt['studies']:
            if st['StudyInstanceUID'] == study_uid:
                n = [sr['SeriesNumber'] for sr in st['series']]
                return int(np.amax(n))
    return 0

def max_instance_number(dbtree, series_uid):
    for pt in dbtree:
        for st in pt['studies']:
            for sr in st['series']:
                if sr['SeriesInstanceUID'] == series_uid:
                    n = list(sr['instances'].keys())
                    return int(np.amax([int(i) for i in n]))
    return 0

def study_id(dbtree, study): # absolute path to series
    uid = study_uid(dbtree, study) # TODO need a function study_list() to avoid unnecessary lookups
    for pt in dbtree:
        for st in pt['studies']:
            if st['StudyInstanceUID'] == uid:
                nr = st['StudyID']
                return int(nr)
    return 0

def series_number(dbtree, series): # absolute path to series
    uid = series_uid(dbtree, series)
    for pt in dbtree:
        for st in pt['studies']:
            for sr in st['series']:
                if sr['SeriesInstanceUID'] == uid:
                    nr = sr['SeriesNumber']
                    return int(nr)
                
def dir(dbtree, entity:list) -> list:
    dir = [entity[0]]
    if len(entity) == 1:
        return dir
    patient_id = entity[1]
    dir += [f"Patient__{patient_id}"]
    if len(entity) == 2:
        return dir
    study_desc = entity[2][0]
    studyid = study_id(dbtree, entity[:3])
    dir += [f"Study__{studyid}__{study_desc}"]
    if len(entity) == 3:
        return dir
    series_desc = entity[3][0]
    series_nr = series_number(dbtree, entity)
    dir += [f"Series__{series_nr}__{series_desc}"]
    return dir


def files(dbtree, entity):
    # Raises an error if the entity does not exist or has no files
    relpath = index(dbtree, entity)
    if relpath==[]:
        raise ValueError(f'No files in entity {entity}')
    if isinstance(entity, str):
        return [os.path.join(entity, str(Path(*f))) for f in relpath]
    else:
        return [os.path.join(entity[0], str(Path(*f))) for f in relpath]
    

def index(dbtree, entity):
    if isinstance(entity, str):
        idx = []
        for pt in sorted(dbtree, key=lambda pt: pt['PatientID']):
            for st in sorted(pt['studies'], key=_sort_study):
                for sr in sorted(st['series'], key=lambda sr: sr['SeriesNumber']):
                    idx += list(sr['instances'].values())
        return idx
    elif len(entity)==2:
        patient_id = entity[1]
        idx = []
        for pt in sorted(dbtree, key=lambda pt: pt['PatientID']):
            if pt['PatientID'] == patient_id:
                for st in sorted(pt['studies'], key=_sort_study):
                    for sr in sorted(st['series'], key=lambda sr: sr['SeriesNumber']):
                        idx += list(sr['instances'].values())
                return idx
        raise ValueError(f'Patient {patient_id} not found')
    elif len(entity)==3:
        study_uid = uid(dbtree, entity)
        idx = []
        for pt in sorted(dbtree, key=lambda pt: pt['PatientID']):
            for st in sorted(pt['studies'], key=_sort_study):
                if st['StudyInstanceUID'] == study_uid:
                    for sr in sorted(st['series'], key=lambda sr: sr['SeriesNumber']):
                        idx += list(sr['instances'].values())
                    return idx
        raise ValueError(f'Study {study_uid} not found')
    elif len(entity)==4:
        series_uid = uid(dbtree, entity)
        for pt in sorted(dbtree, key=lambda pt: pt['PatientID']):
            for st in sorted(pt['studies'], key=_sort_study):
                for sr in sorted(st['series'], key=lambda sr: sr['SeriesNumber']):
                    if sr['SeriesInstanceUID'] == series_uid:
                        return list(sr['instances'].values())
        raise ValueError(f'Series {series_uid} not found')
    
                    
def remove(dbtree, entity):
    if len(entity)==2:
        patient_id = entity[1]
        for pt in sorted(dbtree, key=lambda pt: pt['PatientID']):
            if pt['PatientID'] == patient_id:
                dbtree.remove(pt)
    elif len(entity)==3:
        study_uid = uid(dbtree, entity)
        for pt in sorted(dbtree, key=lambda pt: pt['PatientID']):
            for st in sorted(pt['studies'], key=_sort_study):
                if st['StudyInstanceUID'] == study_uid:
                    pt['studies'].remove(st)
    elif len(entity)==4:
        series_uid = uid(dbtree, entity)
        for pt in sorted(dbtree, key=lambda pt: pt['PatientID']):
            for st in sorted(pt['studies'], key=_sort_study):
                for sr in sorted(st['series'], key=lambda sr: sr['SeriesNumber']):
                    if sr['SeriesInstanceUID'] == series_uid:
                        st['series'].remove(sr)
    return dbtree
                    

def drop(dbtree, relpaths):
    for pt in dbtree[:]:
        for st in pt['studies'][:]:
            for sr in st['series'][:]:
                for nr, relpath in list(sr['instances'].items()):
                    if relpath in relpaths:
                        del sr['instances'][nr]
                if sr['instances'] == []:
                    st['series'].remove(sr)
            if st['series'] == []:
                pt['studies'].remove(st)
        if pt['studies'] == []:
            dbtree.remove(pt)
    return dbtree



def uid(dbtree, entity): # uid from entity
    if len(entity)==2:
        return entity[1]
    if len(entity)==3:
        return study_uid(dbtree, entity)
    if len(entity)==4:
        return series_uid(dbtree, entity)
    

def study_uid(dbtree, study):
    patient_id, study = study[1], study[2]
    for pt in sorted(dbtree, key=lambda pt: pt['PatientID']):
        if pt['PatientID'] == patient_id:

            studies = {}
            for st in sorted(pt['studies'], key=_sort_study):
                study_desc = st['StudyDescription']
                uid_study = st['StudyInstanceUID']
                if study_desc not in studies:
                    studies[study_desc] = [uid_study]
                else:
                    studies[study_desc].append(uid_study)

            if isinstance(study, str):
                if study not in studies:
                    raise ValueError(f"Series {study} not found in patient {patient_id}.")
                if len(studies[study]) == 1:
                    return studies[study][0]
                raise ValueError(
                    f"Multiple studies with name {study} in patient {patient_id}. "
                    f"Please specify the index along with the description. "
                    f"For instance ({study}, {len(studies)-1})'. "
                )
            else:
                if study[0] not in studies:
                    raise ValueError(f"Series {study[0]} not found in patient {patient_id}.")
                try:
                    return studies[study[0]][study[1]]
                except IndexError:
                    raise ValueError(
                        f"Index {study[1]} is out of bounds for series {study[0]}. "
                        f"The are {len(studies[study[0]])} series with description {study[0]} in study {study}."
                    )       


def series_uid(dbtree, series): # absolute path to series
    uid_study = study_uid(dbtree, series[:-1])
    study, sery = series[2], series[3]
    for pt in sorted(dbtree, key=lambda pt: pt['PatientID']):
        for st in sorted(pt['studies'], key=_sort_study):
            if st['StudyInstanceUID'] == uid_study:

                series = {}
                for sr in sorted(st['series'], key=lambda sr: sr['SeriesNumber']):
                    series_desc = sr['SeriesDescription']
                    uid_series = sr['SeriesInstanceUID']
                    if series_desc not in series:
                        series[series_desc] = [uid_series]
                    else:
                        series[series_desc].append(uid_series)

                if isinstance(sery, str):
                    if sery not in series:
                        raise ValueError(f"Series {sery} not found in study {study}.")
                    if len(series[sery]) == 1:
                        return series[sery][0]
                    raise ValueError(
                        f"Multiple series with name {sery}. "
                        f"Please specify the index along with the description. "
                        f"For instance ({sery}, {len(series)-1})'. "
                    )
                else:
                    if sery[0] not in series:
                        raise ValueError(f"Series {sery[0]} not found in study {study}.")
                    try:
                        return series[sery[0]][sery[1]]
                    except IndexError:
                        raise ValueError(
                            f"Index {sery[1]} is out of bounds for series {sery[0]}. "
                            f"The are {len(series[sery[0]])} series with description {sery[0]} in study {study}."
                        )                    


def patients(dbtree, database, name=None, contains=None, isin=None):

    patients = []
    for pt in sorted(dbtree, key=lambda pt: pt['PatientID']):
        patient_name = pt['PatientName']
        append = True
        if name is not None:
            append = append and (patient_name==name)
        if contains is not None:
            append = append and (contains in patient_name)
        if isin is not None:
            append = append and (patient_name in isin)
        if append:
            patients.append(pt['PatientID'])

    return [[database, p] for p in patients]


def studies(dbtree, pat, desc=None, contains=None, isin=None):
    database, patient_id = pat[0], pat[1]
    studies = []
    for pt in sorted(dbtree, key=lambda pt: pt['PatientID']):
        if pt['PatientID'] == patient_id:
            study_idx = {}
            for st in sorted(pt['studies'], key=_sort_study):
                study_desc = st['StudyDescription']
                if study_desc in study_idx:
                    study_idx[study_desc] += 1
                else:
                    study_idx[study_desc] = 0
                studies.append((study_desc, study_idx[study_desc]))

    # Apply filters
    if desc is not None:   
        studies = [s for s in studies if s[0]==desc] 
    if contains is not None:
        studies = [s for s in studies if contains in s[0]]
    if isin is not None:
        studies = [s for s in studies if s[0] in isin]

    # Return result
    return [[database, patient_id, study] for study in studies]



def series(dbtree, stdy, desc=None, contains=None, isin=None):
    database, patient_id, study = stdy[0], stdy[1], stdy[2]
    study_as_str = isinstance(study, str)
    if study_as_str:
        study = (study, 0)
    series = []
    for pt in sorted(dbtree, key=lambda pt: pt['PatientID']):
        if pt['PatientID'] == patient_id:
            study_idx = {}
            for st in sorted(pt['studies'], key=_sort_study):
                study_desc = st['StudyDescription']
                if study_desc in study_idx:
                    study_idx[study_desc] += 1
                else:
                    study_idx[study_desc] = 0
                if study[0] == study_desc:
                    if study_as_str:
                        if study_idx[study_desc] > 0:
                            raise AmbiguousError(
                                f"Multiple studies named {study_desc} in patient {patient_id}. Please provide an index along with the study description."
                            )
                if study == (study_desc, study_idx[study_desc]):
                    series_idx = {}
                    for sr in sorted(st['series'], key=lambda sr: sr['SeriesNumber']):
                        series_desc = sr['SeriesDescription']
                        if series_desc in series_idx:
                            series_idx[series_desc] += 1
                        else:
                            series_idx[series_desc] = 0
                        series.append((series_desc, series_idx[series_desc]))
                    if not study_as_str:
                        break

    # Apply filters (if any)
    if desc is not None:    
        series = [s for s in series if s[0]==desc]
    if contains is not None:    
        series = [s for s in series if contains in s[0]]
    if isin is not None:   
        series = [s for s in series if s[0] in isin] 

    # Return result
    return [[database, patient_id, study, s] for s in series] 
    



def print_tree(dbtree):
    tree = summary(dbtree)
    for patient, studies in tree.items():
        print(f"Patient: ({patient[0]}, {patient[1]})")
        for study, series in studies.items():
            print(f"  Study: ({study[0]}, {study[1]})")
            for s in series:
                print(f"    Series: ({s[0]}, {s[1]})")


def summary(dbtree):
    # A human-readable summary tree

    summary = {}

    for patient in sorted(dbtree, key=lambda pt: pt['PatientID']):
        pat_id, pat_name = patient['PatientID'], patient['PatientName']
        summary[pat_id, pat_name] = {}

        study_idx = {}
        for study in sorted(patient['studies'], key=_sort_study):
            study_desc = study['StudyDescription']
            if study_desc in study_idx:
                study_idx[study_desc] += 1
            else:
                study_idx[study_desc] = 0
            summary[pat_id, pat_name][study_desc, study_idx[study_desc]] = []

            series_idx = {}
            for series in sorted(study['series'], key=lambda sr: sr['SeriesNumber']):
                series_desc = series['SeriesDescription']
                if series_desc in series_idx:
                    series_idx[series_desc] += 1
                else:
                    series_idx[series_desc] = 0
                summary[pat_id, pat_name][study_desc, study_idx[study_desc]].append((series_desc, series_idx[series_desc]))
    
    return summary


def to_csv(dbtree, csv_file):
    
    columns = ['PatientID', 'PatientName', 'StudyDescription', 'study index', 'SeriesDescription', 'SeriesNumber', 'nr of files']

    with open(csv_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(columns)  

        for patient in sorted(dbtree, key=lambda pt: pt['PatientID']):
            pat_id, pat_name = patient['PatientID'], patient['PatientName']

            study_idx = {}
            for study in sorted(patient['studies'], key=_sort_study):
                study_desc = study['StudyDescription']
                if study_desc in study_idx:
                    study_idx[study_desc] += 1
                else:
                    study_idx[study_desc] = 0

                for series in sorted(study['series'], key=lambda sr: sr['SeriesNumber']):
                    n_instances = len(series['instances'])
                    row = [pat_id, pat_name, study_desc, study_idx[study_desc], series['SeriesDescription'], series['SeriesNumber'], n_instances]
                    writer.writerow(row) 




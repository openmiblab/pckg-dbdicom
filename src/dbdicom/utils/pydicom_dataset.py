import math
import datetime

import numpy as np
import pydicom



def get_values(ds, tags:list):
    """Return a list of values for a dataset"""

    # https://pydicom.github.io/pydicom/stable/guides/element_value_types.html
    
    row = []  
    for tag in tags:
        value = None

        # If the tag is provided as string
        if isinstance(tag, str):
            if hasattr(ds, tag):
                pydcm_value = ds[tag].value
                try:
                    VR = pydicom.datadict.dictionary_VR(tag)
                except:
                    VR = None
                value = to_set_type(pydcm_value, VR) # ELIMINATE THIS STEP - return pydicom datatypes

        # If the tag is a tuple of hexadecimal values
        else: 
            if tag in ds:
                try:
                    VR = pydicom.datadict.dictionary_VR(tag)
                except:
                    VR = None
                value = to_set_type(ds[tag].value, VR)

        # If a tag is not present in the dataset, check if it can be derived
        if value is None:
            value = derive_data_element(ds, tag)

        row.append(value)
    return row


def set_values(ds, tags, values, VR=None, coords=None):

    if np.isscalar(tags): 
        tags = [tags]
        values = [values]
        VR = [VR]
    elif VR is None:
        VR = [None] * len(tags)

    if coords is not None: # obsolete
        tags += list(coords.keys())
        values += list(coords.values())

    for i, tag in enumerate(tags):
                
        if values[i] is None:
            if isinstance(tag, str):
                if hasattr(ds, tag):
                    del ds[tag]
            else: # hexadecimal tuple
                if tag in ds:
                    del ds[tag]

        elif isinstance(tag, str):
            if hasattr(ds, tag):
                ds[tag].value = format_value(values[i], tag=tag)
            else:
                _add_new(ds, tag, values[i], VR=VR[i])

        else: # hexadecimal tuple
            if tag in ds:
                ds[tag].value = format_value(values[i], tag=tag)
            else:
                _add_new(ds, tag, values[i], VR=VR[i])

        #_set_derived_data_element(ds, tag, values[i])
                
    return ds


def value(ds, tags):
    # Same as get_values but without VR lookup

    # https://pydicom.github.io/pydicom/stable/guides/element_value_types.html
    if np.isscalar(tags): 
        return get_values(ds, [tags])[0]
            
    row = []  
    for tag in tags:
        value = None

        # If the tag is provided as string
        if isinstance(tag, str):

            if hasattr(ds, tag):
                value = to_set_type(ds[tag].value)

        # If the tag is a tuple of hexadecimal values
        else: 
            if tag in ds:
                value = to_set_type(ds[tag].value)

        # If a tag is not present in the dataset, check if it can be derived
        if value is None:
            value = derive_data_element(ds, tag)

        row.append(value)
    return row


def set_value(ds, tags, values):
    # Same as set_values but without VR lookup
    # This excludes new private tags - set those using add_private()
    if np.isscalar(tags): 
        tags = [tags]
        values = [values]

    for i, tag in enumerate(tags):
                
        if values[i] is None:
            if isinstance(tag, str):
                if hasattr(ds, tag):
                    del ds[tag]
            else: # hexadecimal tuple
                if tag in ds:
                    del ds[tag]

        elif isinstance(tag, str):
            if hasattr(ds, tag):
                ds[tag].value = check_value(values[i], tag)
            else:
                add_new(ds, tag, values[i])

        else: # hexadecimal tuple
            if tag in ds:
                ds[tag].value = check_value(values[i], tag)
            else:
                add_new(ds, tag, values[i])
                
    return ds


def check_value(value, tag):

    # If the change below is made (TM, DA, DT) then this needs to 
    # convert those to string before setting

    LO = [
        'SeriesDescription',
        'StudyDescription',
    ]
    TM = [
        'AcquisitionTime',
    ]

    if tag in LO:
        if len(value) > 64:
            return value[-64:]
    if tag in TM:
        return seconds_to_str(value)
    
    return value


def add_new(ds, tag, value):
    if not isinstance(tag, pydicom.tag.BaseTag):
        tag = pydicom.tag.Tag(tag)
    if tag.is_private:
        raise ValueError("if you want to add a private data element, use "
                         "dataset.add_private()")
   # Add a new data element
    value_repr = pydicom.datadict.dictionary_VR(tag)
    if value_repr == 'US or SS':
        if value >= 0:
            value_repr = 'US'
        else:
            value_repr = 'SS'
    elif value_repr == 'OB or OW':
        value_repr = 'OW'
    ds.add_new(tag, value_repr, format_value(value, value_repr))

def _add_new(ds, tag, value, VR='OW'):
    if not isinstance(tag, pydicom.tag.BaseTag):
        tag = pydicom.tag.Tag(tag)
    if not tag.is_private: # Add a new data element
        value_repr = pydicom.datadict.dictionary_VR(tag)
        if value_repr == 'US or SS':
            if value >= 0:
                value_repr = 'US'
            else:
                value_repr = 'SS'
        elif value_repr == 'OB or OW':
            value_repr = 'OW'
        ds.add_new(tag, value_repr, format_value(value, value_repr))
    else:
        if (tag.group, 0x0010) not in ds:
            ds.private_block(tag.group, 'dbdicom ' + str(tag.group), create=True)
        ds.add_new(tag, VR, format_value(value, VR))

def format_value(value, VR=None, tag=None):

    # If the change below is made (TM, DA, DT) then this needs to 
    # convert those to string before setting

    # Slow - dictionary lookup for every value write

    if VR is None:
        VR = pydicom.datadict.dictionary_VR(tag)

    if VR == 'LO':
        if len(value) > 64:
            return value[-64:]
            #return value[:64]
    if VR == 'TM':
        return seconds_to_str(value)
    if VR == 'DA':
        if not is_valid_dicom_date(value):
            return '99991231'
    
    return value


def is_valid_dicom_date(da_str: str) -> bool:
    if not isinstance(da_str, str) or len(da_str) != 8 or not da_str.isdigit():
        return False
    try:
        datetime.datetime.strptime(da_str, "%Y%m%d")
        return True
    except ValueError:
        return False


def to_set_type(value, VR=None):
    """
    Convert pydicom datatypes to the python datatypes used to set the parameter.
    """
    # Not a good idea to modify pydicom set/get values. confusing and requires extra VR lookups

    if VR == 'TM':
        # pydicom sometimes returns string values for TM data types
        if isinstance(value, str):
            return str_to_seconds(value)

    if value.__class__.__name__ == 'MultiValue':
        return [to_set_type(v, VR) for v in value]
    if value.__class__.__name__ == 'PersonName':
        return str(value)
    if value.__class__.__name__ == 'Sequence':
        return [ds for ds in value]
    if value.__class__.__name__ == 'TM': 
        return time_to_seconds(value) # return datetime.time
    if value.__class__.__name__ == 'UID': 
        return str(value) 
    if value.__class__.__name__ == 'IS': 
        return int(value)
    if value.__class__.__name__ == 'DT': 
        return datetime_to_str(value) # return datetime.datetime
    if value.__class__.__name__ == 'DA':  # return datetime.date
        return date_to_str(value)
    if value.__class__.__name__ == 'DSfloat': 
        return float(value)
    if value.__class__.__name__ == 'DSdecimal': 
        return int(value)
    
    return value


def derive_data_element(ds, tag):
    """Tags that are not required but can be derived from other required tags"""

    if tag == 'SliceLocation' or tag == (0x0020, 0x1041):
        if 'ImageOrientationPatient' in ds and 'ImagePositionPatient' in ds:
            image_orientation = ds['ImageOrientationPatient'].value
            image_position = ds['ImagePositionPatient'].value
            row_cosine = np.array(image_orientation[:3])    
            column_cosine = np.array(image_orientation[3:]) 
            slice_cosine = np.cross(row_cosine, column_cosine)
            return np.dot(np.array(image_position), slice_cosine)
    # To be extended ad hoc with other tags that can be derived


def add_private(ds, tag, value, VR):
    if not isinstance(tag, pydicom.tag.BaseTag):
        tag = pydicom.tag.Tag(tag)
    if (tag.group, 0x0010) not in ds:
        ds.private_block(tag.group, 'dbdicom ' + str(tag.group), create=True)
    ds.add_new(tag, VR, format_value(value, VR))




# UTILITIES
    

def str_to_seconds(dicom_tm):
    if dicom_tm is None:
        return None
    if dicom_tm == '':
        return None
    # dicom_tm is of the form 'HHMMSS.FFFFFF'
    # Split the seconds into seconds and fractional seconds
    seconds, fractional_seconds = dicom_tm.split('.')
    # Convert the hours, minutes, and seconds to integers
    hours = int(seconds[:2])
    minutes = int(seconds[2:4])
    seconds = int(seconds[4:])
    # Convert the fractional seconds to a decimal
    fractional_seconds = float('.' + fractional_seconds)
    # Create a datetime object representing the time
    seconds_since_midnight = (hours * 3600) + (minutes * 60) + seconds + fractional_seconds
    return seconds_since_midnight

def seconds_to_str(seconds_since_midnight):
    # if not isinstance(seconds_since_midnight, float): 
    #     return None
    if seconds_since_midnight is None:
        return None
    seconds_since_midnight = float(seconds_since_midnight)
    hours = math.floor(seconds_since_midnight/3600)
    minutes = math.floor((seconds_since_midnight-hours*3600)/60)
    seconds = math.floor(seconds_since_midnight-hours*3600-minutes*60)
    fractional_seconds = round(seconds_since_midnight-hours*3600-minutes*60-seconds, 6)
    hours = str(hours).zfill(2)
    minutes = str(minutes).zfill(2)
    seconds = str(seconds).zfill(2)
    fractional_seconds = str(fractional_seconds)
    fractional_seconds = fractional_seconds.split('.')
    if len(fractional_seconds) == 2:
        fractional_seconds = fractional_seconds[1].ljust(6,'0')
    else:
        fractional_seconds = '0'.ljust(6,'0')
    return hours + minutes + seconds + '.' + fractional_seconds

def time_to_seconds(tm):
    if tm is None:
        return None
    hours = tm.hour
    minutes = tm.minute
    seconds = tm.second
    fractional_seconds = tm.microsecond / 1000000.0
    seconds_since_midnight = (hours * 3600) + (minutes * 60) + seconds + fractional_seconds
    return seconds_since_midnight

def seconds_to_time(seconds_since_midnight):
    # if not isinstance(seconds_since_midnight, float): 
    #     return None
    if seconds_since_midnight is None:
        return None
    seconds_since_midnight = float(seconds_since_midnight)
    hours = math.floor(seconds_since_midnight/3600)
    minutes = math.floor((seconds_since_midnight-hours*3600)/60)
    seconds = math.floor(seconds_since_midnight-hours*3600-minutes*60)
    fractional_seconds = round(seconds_since_midnight-hours*3600-minutes*60-seconds, 6)
    microseconds = fractional_seconds*1000000.0
    return datetime.time(int(hours), int(minutes), int(seconds), int(microseconds))

def time_to_str(tm):
    if tm is None:
        return None
    hours = tm.hour
    minutes = tm.minute
    seconds = tm.second
    fractional_seconds = tm.microsecond / 1000000.0   
    hours = str(hours).zfill(2)
    minutes = str(minutes).zfill(2)
    seconds = str(seconds).zfill(2)
    fractional_seconds = str(fractional_seconds)
    _, fractional_seconds = fractional_seconds.split('.')
    fractional_seconds = fractional_seconds.ljust(6,'0')
    return hours + minutes + seconds + '.' + fractional_seconds 

def date_to_str(tm):
    if tm is None:
        return None
    year = str(tm.year).rjust(4, '0')
    month = str(tm.month).rjust(2, '0')
    day = str(tm.day).rjust(2, '0')
    return year + month + day

def datetime_to_str(dt):
    if dt is None:
        return None
    date = date_to_str(dt.date())
    time = time_to_str(dt.time())
    return date + time
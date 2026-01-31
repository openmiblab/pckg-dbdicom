import os
import shutil
from itertools import product

import numpy as np
import vreg

import dbdicom.dbd
import dbdicom as db
from dbdicom.utils.points import duplicates


# pytest test_duplicates.py




# ------------------------
# Basic functionality tests
# ------------------------
def test_duplicates_no_duplicates():
    x = [0, 1, 0, 1]
    y = [0, 0, 1, 1]
    assert duplicates([x, y]) == []


def test_duplicates_simple_duplicate():
    x = [0, 1, 0, 1, 0]
    y = [0, 0, 1, 1, 0]  # last point duplicates the first
    assert duplicates([x, y]) == [4]


def test_duplicates_multiple_duplicates():
    x = [0, 1, 0, 1, 0, 1]
    y = [0, 0, 1, 1, 0, 0]
    assert duplicates([x, y]) == [4, 5]


def test_duplicates_single_point():
    coords = [[42], ['A']]
    assert duplicates(coords) == []


# ------------------------
# Unhashable types
# ------------------------
def test_duplicates_lists():
    x = [0, 1, 0, 1]
    y = [[1,2],[1,2],[3,4],[3,4]]
    assert duplicates([x,y]) == []


def test_duplicates_lists_with_duplicate():
    x = [0, 1, 0, 1, 0]
    y = [[1,2],[1,2],[3,4],[3,4],[1,2]]
    assert duplicates([x,y]) == [4]


def test_duplicates_numpy_arrays():
    x = [0, 1, 0, 1]
    y = [np.array([1,2]), np.array([1,2]), np.array([3,4]), np.array([3,4])]
    assert duplicates([x,y]) == []


def test_duplicates_numpy_arrays_with_duplicate():
    x = [0, 1, 0, 1, 0]
    y = [np.array([1,2]), np.array([1,2]), np.array([3,4]), np.array([3,4]), np.array([1,2])]
    assert duplicates([x,y]) == [4]


# ------------------------
# Mixed types
# ------------------------
def test_duplicates_mixed_types():
    x_vals = [0, 1]
    y_vals = ['A', 'B']
    z_vals = [(1,2), (3,4)]
    w_vals = [np.array([1,2]), np.array([3,4])]
    q_vals = [['X','Y'], ['Z','W']]

    points = list(product(x_vals, y_vals, z_vals, w_vals, q_vals))
    x, y, z, w, q = zip(*points)
    # No duplicates in the full Cartesian product
    assert duplicates([x,y,z,w,q]) == []

    # Introduce a duplicate
    x = list(x) + [0]
    y = list(y) + ['A']
    z = list(z) + [(1,2)]
    w = list(w) + [np.array([1,2])]
    q = list(q) + [['X','Y']]
    assert duplicates([x,y,z,w,q])[-1] == len(x)-1  # last point is detected as duplicate


# ------------------------
# Edge case: empty input
# ------------------------
def test_duplicates_empty_input():
    assert duplicates([[], []]) == []



def test_full_name():

    tmp = os.path.join(os.getcwd(), 'tests', 'tmp')
    os.makedirs(tmp, exist_ok=True)
    shutil.rmtree(tmp)
    os.makedirs(tmp, exist_ok=True)

    values = 100*np.random.rand(128, 192, 20).astype(np.float32)
    vol = vreg.volume(values)
    series = [tmp, '007', 'dbdicom_test', 'ax']
    db.write_volume(vol, series)

    series_fn = [tmp, '007', ('dbdicom_test', 0), ('ax', 0)]
    series_fn_test = dbdicom.dbd.full_name(series)
    for i, a in enumerate(series_fn):
        a == series_fn_test[i]
    series_fn_test = dbdicom.dbd.full_name(series_fn_test)
    for i, a in enumerate(series_fn):
        a == series_fn_test[i]   

    study_fn = [tmp, '007', ('dbdicom_test', 0)]
    study_fn_test = dbdicom.dbd.full_name(series[:3])
    for i, a in enumerate(study_fn):
        a == study_fn_test[i]
    study_fn_test = dbdicom.dbd.full_name(study_fn_test)
    for i, a in enumerate(study_fn):
        a == study_fn_test[i]     

    shutil.rmtree(tmp)



if __name__=='__main__':

    test_full_name()

    print('All utils tests have passed!!!')
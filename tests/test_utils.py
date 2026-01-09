import os
import shutil
import numpy as np
import vreg

import dbdicom.utils.arrays
import dbdicom.dbd
import dbdicom as db


def test_meshvals():

    z = [0,1,2,0,1,2]
    p = ['A','A','A','B','B','B']
    coords, inds = dbdicom.utils.arrays.meshvals([z,p])
    assert coords[0].shape == (3,2)
    assert np.array_equal(coords[1][0,:], ['A','B'])
    assert np.array_equal(coords[0][0,:], [0,0])
    assert np.array_equal(coords[0][:,0], [0,1,2])
    assert np.array_equal(coords[1][:,0], ['A','A','A'])
    assert np.array_equal(inds, [0,3,1,4,2,5])

    z = [0,1,2,0,1,2]
    p = [['A','B'],['A','B'],['A','B'],['C','D'],['C','D'],['C','D']]
    coords, inds = dbdicom.utils.arrays.meshvals([z,p])
    assert coords[0].shape == (3,2)
    assert np.array_equal(coords[1][0,0], ['A','B'])
    assert np.array_equal(coords[1][0,1], ['C','D'])
    assert np.array_equal(coords[0][0,:], [0,0])
    assert np.array_equal(coords[0][:,0], [0,1,2])
    assert np.array_equal(coords[1][0,0], ['A','B'])
    assert np.array_equal(coords[1][1,0], ['A','B'])
    assert np.array_equal(coords[1][2,0], ['A','B'])
    assert np.array_equal(inds, [0,3,1,4,2,5])

    z = [0,1,2,0,1,2,3]
    p = ['A','A','A','B','B','B']
    try:
        dbdicom.utils.arrays.meshvals([z,p])
    except:
        assert True
    else:
        assert False

    z = [0,1,2,0,1,3]
    p = ['A','A','A','B','B','B']
    try:
        dbdicom.utils.arrays.meshvals([z,p])
    except:
        assert True
    else:
        assert False

    z = [0,1,2,0,1,2]
    p = ['A','A','A','B','B','A']
    try:
        dbdicom.utils.arrays.meshvals([z,p])
    except:
        assert True
    else:
        assert False


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

    test_meshvals()
    test_full_name()

    print('All utils tests have passed!!!')
import os
import shutil
import numpy as np
import dbdicom as db
import vreg


tmp = os.path.join(os.getcwd(), 'tests', 'tmp')
os.makedirs(tmp, exist_ok=True)
shutil.rmtree(tmp)
os.makedirs(tmp, exist_ok=True)




def test_write_volume():

    values = 100*np.random.rand(128, 192, 20).astype(np.float32)
    vol = vreg.volume(values)
    series = [tmp, '007', 'dbdicom_test', 'ax']
    db.write_volume(vol, series)

    values = np.zeros((256, 256, 16, 2))
    affine = np.eye(4)
    vol = vreg.volume(values, affine, coords=(['INPHASE', 'OUTPHASE'], ), dims=['ImageType'])
    series = [tmp, '007', 'dbdicom_test', 'dixon']
    db.write_volume(vol, series)

    # Writing to an existing series returns an error by default
    try:
        db.write_volume(vol, series)
    except:
        assert True
    else:
        assert False

    # Translate the volume in the z-direction over 10mm and append to the series
    # This creates a series with two volumes separated by a gap of 5 mm
    vol2 = vol.translate([0,0,20], coords='volume')
    db.write_volume(vol2, series, append=True)

    # Reading now throws an error as there are multiple volumes in the series
    try:
        db.volume(series, dims=['ImageType'])
    except:
        assert True
    else:
        assert False


    shutil.rmtree(tmp)


def test_volumes_2d():

    # Write one volume
    values = 100*np.random.rand(128, 192, 5).astype(np.float32)
    vol = vreg.volume(values)
    series = [tmp, '007', 'dbdicom_test', 'ax']
    db.write_volume(vol, series)

    # Shift it up to leave a gap and write to the same series
    vol2 = vol.translate([0,0,10], coords='volume')
    db.write_volume(vol2, series, append=True)

    # Trying to read as a single volume throws an error because of the gap
    try:
        db.volume(series)
    except:
        assert True
    else:
        assert False

    # But we can read them as 2D volumes, returning 10 2D volumes
    vols = db.volumes_2d(series)
    assert len(vols) == 10

    # Now 4D
    values = np.zeros((256, 256, 5, 2))
    affine = np.eye(4)
    vol = vreg.volume(values, affine, coords=(['INPHASE', 'OUTPHASE'], ), dims=['ImageType'])
    series = [tmp, '007', 'dbdicom_test', 'dixon']
    db.write_volume(vol, series)

    vol2 = vol.translate([0,0,10], coords='volume')
    db.write_volume(vol2, series, append=True)

    vols = db.volumes_2d(series, dims=['ImageType'])
    assert len(vols) == 10
    assert vols[-1].shape == (256, 256, 1, 2)

    shutil.rmtree(tmp)


def test_volume():

    # One slice
    values = 100*np.random.rand(128, 192, 1).astype(np.float32)
    vol = vreg.volume(values)
    series = [tmp, '007', 'test', 'slice']
    db.write_volume(vol, series)
    vol2 = db.volume(series)
    assert np.linalg.norm(vol2.values-vol.values) < 0.0001*np.linalg.norm(vol.values)
    assert np.linalg.norm(vol2.affine-vol.affine) == 0

    # 3D volume
    values = 100*np.random.rand(2, 192, 20).astype(np.float32)
    vol = vreg.volume(values)
    series = [tmp, '007', 'test', 'ax']
    db.write_volume(vol, series)
    vol2 = db.volume(series)
    assert np.linalg.norm(vol2.values-vol.values) < 0.0001*np.linalg.norm(vol.values)
    assert np.linalg.norm(vol2.affine-vol.affine) == 0

    # 4D volume
    # values = 100*np.random.rand(256, 256, 3, 2).astype(np.float32)
    values = 100*np.random.rand(2, 2, 2, 2).astype(np.float32)
    image_type = [['ORIGINAL', 'INPHASE'], ['ORIGINAL', 'OUTPHASE']]
    vol = vreg.volume(values, dims=['ImageType'], coords=(image_type, ), orient='coronal')
    series = [tmp, '007', 'dbdicom_test', 'dixon']
    db.write_volume(vol, series)
    vol2 = db.volume(series, dims=['ImageType'])
    assert np.linalg.norm(vol2.values-vol.values) < 0.0001 * np.linalg.norm(vol.values)
    assert np.linalg.norm(vol2.affine-vol.affine) == 0
    assert vol2.dims == vol.dims
    assert np.array_equal(vol2.coords[0], vol.coords[0])

    values = 100*np.random.rand(256, 256, 3, 2, 2).astype(np.float32)
    dims = ['FlipAngle','ImageType']
    coords = ([10, 20], image_type)
    vol = vreg.volume(values, dims=dims, coords=coords, orient='coronal')
    series = [tmp, '007', 'dbdicom_test', 'vfa_dixon']
    db.write_volume(vol, series)
    vol2 = db.volume(series, dims=dims)
    assert np.linalg.norm(vol2.values-vol.values) < 0.0001*np.linalg.norm(vol.values)
    assert np.linalg.norm(vol2.affine-vol.affine) == 0
    assert vol2.dims == vol.dims
    assert np.array_equal(vol2.coords[0], vol.coords[0])

    # Test filtering feature
    vol3 = db.volume(series, dims=['FlipAngle'], ImageType=['ORIGINAL', 'INPHASE'])
    assert np.linalg.norm(vol3.values-vol.values[...,0]) < 0.0001*np.linalg.norm(vol.values[...,0])

    shutil.rmtree(tmp)


def test_values():

    values = 100*np.random.rand(256, 256, 3, 2, 2).astype(np.float32)
    dims = ['FlipAngle','ImageType']
    coords = ([10, 20], ['INPHASE', 'OUTPHASE'])
    vol = vreg.volume(values, dims=dims, coords=coords, orient='coronal')
    series = [tmp, '007', 'dbdicom_test', 'vfa_dixon']
    db.write_volume(vol, series)

    # Read all slice locations as 1D array
    locs = db.values(series, 'SliceLocation')
    assert locs.shape == (12,)
    assert np.array_equal(locs[-3:], [0,1,2])

    locs, fa = db.values(series, 'SliceLocation', 'FlipAngle')
    assert np.array_equal(np.unique(fa), [10,20])

    locs = db.values(series, 'SliceLocation', dims=['SliceLocation', 'FlipAngle', 'ImageType'])
    assert locs.shape == (3,2,2)
    assert np.unique(locs[0,...]) == [0]

    locs, it = db.values(series, 'SliceLocation', 'ImageType', dims=['SliceLocation', 'FlipAngle', 'ImageType'])
    assert it.shape == (3,2,2)

    pn = db.values(series, 'PatientName', dims=['SliceLocation', 'FlipAngle', 'ImageType'])
    assert pn.shape == (3,2,2)
    assert np.unique(pn) == ['Anonymous']

    # Improper dimensions
    try:
        db.values(series, 'SliceLocation', 'ImageType', dims=['SliceLocation'])
    except:
        assert True
    else:
        assert False

    shutil.rmtree(tmp)


def test_edit():

    values = 100*np.random.rand(256, 256, 3, 2, 2).astype(np.float32)
    dims = ['FlipAngle','ImageType']
    coords = ([10, 20], ['INPHASE', 'OUTPHASE'])
    vol = vreg.volume(values, dims=dims, coords=coords, orient='coronal')
    series = [tmp, '007', 'dbdicom_test', 'vfa_dixon']
    db.write_volume(vol, series)

    shape = (3,2,2)
    dims = ('SliceLocation', 'FlipAngle', 'ImageType')
    new_tr = np.arange(np.prod(shape)).reshape(shape)
    new_pn = np.full(shape, 'James Bond').reshape(shape)
    new_values = {'RepetitionTime': new_tr, 'PatientName': new_pn}
    db.edit(series, new_values, dims=dims)
    tr, pn = db.values(series, 'RepetitionTime', 'PatientName', dims=dims)
    assert np.array_equal(tr, new_tr)
    assert np.array_equal(pn, new_pn)

    shutil.rmtree(tmp)


def test_write_database():
    values = 100*np.random.rand(16, 16, 4).astype(np.float32)
    vol = vreg.volume(values)
    db.write_volume(vol, [tmp, '007', 'test', 'ax'])    # create series ax
    try:
        db.write_volume(vol, [tmp, '007', 'test', 'ax'])    # add to it
    except:
        assert True
    else:
        assert False
    try:
        db.write_volume(vol, [tmp, '007', 'test', ('ax', 0)])   # add to it
    except:
        assert True
    else:
        assert False
    db.write_volume(vol, [tmp, '007', 'test', ('ax', 1)])   # create a new series ax
    db.write_volume(vol, [tmp, '007', 'test', ('ax', 3)])   # create a new series ax
    try:
        db.write_volume(vol, [tmp, '007', 'test', 'ax'])   # Ambiguous
    except:
        assert True
    else:
        assert False
    db.write_volume(vol, [tmp, '008', 'test', 'ax'])            # Create a new patient
    db.write_volume(vol, [tmp, '008', 'test', 'ax-2'])          # Add a new series
    try:
        db.write_volume(vol, [tmp, '008', ('test', 0), 'ax'])       # Add to the series ax 
    except:
        assert True
    else:
        assert False
    db.write_volume(vol, [tmp, '008', ('test', 1), 'ax'])       # Add to a new study
    try:
        db.write_volume(vol, [tmp, '008', 'test', 'ax'])       # Ambiguous
    except:
        assert True
    else:
        assert False

    series = db.series(tmp)
    [print(s) for s in series]

    assert ('ax', 2) in [s[-1] for s in series]
    assert [] == db.series(tmp, contains='b')
    assert 2 == len(db.patients(tmp))
    assert 2 == len(db.patients(tmp, name='Anonymous'))

    shutil.rmtree(tmp)

def test_copy():

    # Build some data
    tmp1 = os.path.join(tmp, 'dir1')
    tmp2 = os.path.join(tmp, 'dir2')
    os.makedirs(tmp1, exist_ok=True)
    os.makedirs(tmp2, exist_ok=True)
    values = 100*np.random.rand(16, 16, 4).astype(np.float32)
    vol = vreg.volume(values)
    db.write_volume(vol, [tmp1, '007', 'test', 'ax'])    # create series ax
    db.write_volume(vol, [tmp1, '007', 'test2', 'ax2'])    # create series ax

    # Copy to named entity
    db.copy([tmp1, '007', 'test2', 'ax2'], [tmp2, '007', 'test2', 'ax2'])
    db.copy([tmp1, '007', 'test2', 'ax2'], [tmp2, '007', 'test2', 'ax'])
    db.copy([tmp1, '007', 'test2', 'ax2'], [tmp2, '007', 'test2', 'ax'])
    copy_ax2 = db.copy([tmp1, '007', 'test2', 'ax2'])
    print('0')
    [print(s) for s in db.series(tmp2)]

    db.copy([tmp1, '007', 'test2'], [tmp2, '008', 'test2'])
    copy_test2 = db.copy([tmp1, '007', 'test2'])
    assert len(db.series(copy_test2)) == len(db.series([tmp1, '007', 'test2']))
    print('1')
    [print(s) for s in db.series(tmp2)]
    assert 2==len(db.patients(tmp2))
    assert 4==len(db.series(tmp2))
    db.copy([tmp1, '007', 'test2'], [tmp2, '008', 'test2']) 
    print('2')
    [print(s) for s in db.series(tmp2)]
    assert 6==len(db.series(tmp2))
    db.copy([tmp1, '007'], [tmp2, '008'])
    copy_007 = db.copy([tmp1, '007'])
    print('3')
    [print(s) for s in db.series(tmp2)]
    assert 11==len(db.series(tmp2))
    assert 5==len(db.studies(tmp2))
    assert len(db.series(copy_007)) == len(db.series([tmp1, '007']))

    shutil.rmtree(tmp)


def test_db_read():

    # Build some data
    tmp1 = os.path.join(tmp, 'dir')
    os.makedirs(tmp1, exist_ok=True)
    values = np.arange(16 * 16 * 4).reshape((16, 16, 4))
    vol = vreg.volume(values)
    series = [tmp1, '007', 'test', 'ax']
    db.write_volume(vol, series)

    # Delete the index file and read again
    idx = os.path.join(tmp1, 'index.json')
    os.remove(idx)
    vol_rec = db.volume(series)
    assert np.linalg.norm(vol_rec.values - vol.values) == 0



if __name__ == '__main__':

    test_write_volume()
    test_volumes_2d()
    test_values()
    test_edit()
    test_write_database()
    test_copy()
    test_volume()
    test_db_read()

    print('All api tests have passed!!!')
import os
import shutil
import numpy as np
import dbdicom as db
import vreg


tmp = os.path.join(os.getcwd(), 'tests', 'tmp')
os.makedirs(tmp, exist_ok=True)
shutil.rmtree(tmp)
os.makedirs(tmp, exist_ok=True)



def test_dti_volume():

    dti_series = [tmp, '007', 'diff_test', 'DTI']
    dims = ['DiffusionBValue', 'DiffusionGradientOrientation']

    # Define pixel values
    arr = np.ones(16 * 16 * 4).reshape((16, 16, 4, 1, 1))
    b0 = 0 * arr
    b700_1 = 1 * arr
    b700_2 = 2 * arr
    b700_3 = 3 * arr

    # Define coordinates
    b0_coords = ([0], [[0,0,0]])
    b700_1_coords = ([700], [[0,0,1]])
    b700_2_coords = ([700], [[0,1,0]])
    b700_3_coords = ([700], [[1,0,0]])

    # Build volumes
    b0_vol = vreg.volume(b0, dims=dims, coords=b0_coords)
    b700_1_vol = vreg.volume(b700_1, dims=dims, coords=b700_1_coords)
    b700_2_vol = vreg.volume(b700_2, dims=dims, coords=b700_2_coords)
    b700_3_vol = vreg.volume(b700_3, dims=dims, coords=b700_3_coords)

    # Save in a single dicom series
    db.write_volume(b0_vol, dti_series)
    db.write_volume(b700_1_vol, dti_series, append=True)
    db.write_volume(b700_2_vol, dti_series, append=True)
    db.write_volume(b700_3_vol, dti_series, append=True)

    # The series is not a single volume, so this fails:
    try:
        db.volume(dti_series)
    except:
        assert True
    else:
        assert False

    # Sorting by dimensions does not work either 
    # because b0 and b700 have different shape
    try:
        db.volume(dti_series, dims=dims)
    except:
        assert True
    else:
        assert False

    # What does work is use a filter to separate out b0 and b700:
    b0_vol_rec = db.volume(dti_series, dims=dims, DiffusionBValue=0)
    b700_vol_rec = db.volume(dti_series, dims=dims, DiffusionBValue=700)

    # Check the shapes and values
    assert b0_vol_rec.shape == (16, 16, 4, 1, 1)
    assert b700_vol_rec.shape == (16, 16, 4, 1, 3)

    # Check coordinates
    assert b0_vol_rec.coords[0][0,0] == 0
    assert b700_vol_rec.coords[0][0,2] == 700
    assert np.array_equal(b0_vol_rec.coords[1][0,0], [0,0,0])
    assert np.array_equal(b700_vol_rec.coords[1][0,2], [1,0,0])

    # Check values
    assert np.array_equal(b0_vol_rec.values, b0_vol.values)
    assert np.array_equal(b700_vol_rec.values[...,0], b700_1[...,0])
    assert np.array_equal(b700_vol_rec.values[...,1], b700_2[...,0])
    assert np.array_equal(b700_vol_rec.values[...,2], b700_3[...,0])

    # If we are not interested in reading all coordinates, 
    # we can also read differently as below:

    # b0 is just a single 3D volume so we don't have to specify dimensions
    b0_vol_rec = db.volume(dti_series, DiffusionBValue=0)

    # Now this has a 3D shape
    assert b0_vol_rec.shape == (16, 16, 4)
    assert b0_vol_rec.dims is None
    assert b0_vol_rec.coords is None
    assert np.array_equal(b0_vol_rec.values, b0_vol.values[...,0,0])

    # b700 has only one b-value so we just need to specify orientation as dimensions
    b700_vol_rec = db.volume(dti_series, dims='DiffusionGradientOrientation', DiffusionBValue=700)

    # This has 4D shape
    assert b700_vol_rec.shape == (16, 16, 4, 3)
    assert np.array_equal(b700_vol_rec.dims, ['DiffusionGradientOrientation'])
    assert np.array_equal(b700_vol_rec.coords[0][2], [1,0,0])
    assert np.array_equal(b700_vol_rec.values[...,2], b700_3[...,0,0])

    # We can also read just the b-values and b-vectors from the data
    bvals, bvecs, zloc = db.values(dti_series, 'DiffusionBValue', 'DiffusionGradientOrientation', 'SliceLocation')

    shutil.rmtree(tmp)


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

    shutil.rmtree(tmp)


def test_rw_series():

    # Write three series with the same desc in the same study
    tmp1 = os.path.join(tmp, 'dir1')
    os.makedirs(tmp1, exist_ok=True)
    values_1 = 1 * np.arange(16 * 16 * 4).reshape((16, 16, 4))
    values_2 = 2 * np.arange(16 * 16 * 4).reshape((16, 16, 4))
    values_3 = 3 * np.arange(16 * 16 * 4).reshape((16, 16, 4))
    vol_1 = vreg.volume(values_1)
    vol_2 = vreg.volume(values_2)
    vol_3 = vreg.volume(values_3)
    series_1 = [tmp1, '007', 'test', ('ax', 0)]
    series_2 = [tmp1, '007', 'test', ('ax', 1)]
    series_3 = [tmp1, '007', 'test', ('ax', 2)]

    db.write_volume(vol_1, series_1)
    db.write_volume(vol_2, series_2)
    db.write_volume(vol_3, series_3)

    try:
        db.write_volume(vol_2, [tmp1, '007', 'test', 'ax'])
    except:
        assert True
    else:
        assert False

    v_3 = db.volume(series_3)
    assert np.array_equal(v_3.values, values_3)

    v_3 = db.volume([tmp1, '007', 'test', ('ax', -1)])
    assert np.array_equal(v_3.values, values_3)

    shutil.rmtree(tmp)


def test_rw_studies():

    # Write three series with the same desc in the same study
    tmp1 = os.path.join(tmp, 'dir1')
    os.makedirs(tmp1, exist_ok=True)
    values_1 = 1 * np.arange(16 * 16 * 4).reshape((16, 16, 4))
    values_2 = 2 * np.arange(16 * 16 * 4).reshape((16, 16, 4))
    values_3 = 3 * np.arange(16 * 16 * 4).reshape((16, 16, 4))
    vol_1 = vreg.volume(values_1)
    vol_2 = vreg.volume(values_2)
    vol_3 = vreg.volume(values_3)
    series_1 = [tmp1, '007', ('test', 0), 'ax']
    series_2 = [tmp1, '007', ('test', 1), 'ax']
    series_3 = [tmp1, '007', ('test', 2), 'ax']

    db.write_volume(vol_1, series_1)
    db.write_volume(vol_2, series_2)
    db.write_volume(vol_3, series_3)

    try:
        db.write_volume(vol_2, [tmp1, '007', 'test', 'ax'])
    except:
        assert True
    else:
        assert False

    v_2 = db.volume([tmp1, '007', ('test', 1), 'ax'])
    assert np.array_equal(v_2.values, values_2)

    v_3 = db.volume([tmp1, '007', ('test', 2), 'ax'])
    assert np.array_equal(v_3.values, values_3)

    v_3 = db.volume([tmp1, '007', ('test', -1), 'ax'])
    assert np.array_equal(v_3.values, values_3)

    shutil.rmtree(tmp)



if __name__ == '__main__':

    test_dti_volume()
    test_write_volume()
    test_volumes_2d()
    test_values()
    test_edit()
    test_write_database()
    test_copy()
    test_volume()
    test_db_read()
    test_rw_series()
    test_rw_studies()

    print('All api tests have passed!!!')
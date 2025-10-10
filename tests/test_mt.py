import os
import shutil
import numpy as np
import dbdicom as db
import vreg


tmp = os.path.join(os.getcwd(), 'tests', 'tmp')
os.makedirs(tmp, exist_ok=True)
shutil.rmtree(tmp)
os.makedirs(tmp, exist_ok=True)


def test_write():

    # Currently works as expected
    values = 100*np.random.rand(256, 256, 16, 2).astype(np.float32)
    dims = ['ScanOptions']
    coords = (['PFP', 'PFP'], )
    vol = vreg.volume(values, dims=dims, coords=coords, orient='coronal')
    series = [tmp, '007', 'dbdicom_test', 'mt1']
    db.write_volume(vol, series)
    
    # Currently does not work
    values = 100*np.random.rand(256, 256, 16, 2).astype(np.float32)
    dims = ['ScanOptions']
    coords = (['PFP', ['PFP', 'MT']], )
    vol = vreg.volume(values, dims=dims, coords=coords, orient='coronal')
    series = [tmp, '007', 'dbdicom_test', 'mt2']
    db.write_volume(vol, series)

    shutil.rmtree(tmp)


if __name__ == '__main__':

    test_write()

    print('All api tests have passed!!!')

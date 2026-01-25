# Test with real data - not included in the distribution

import dbdicom as db

def dfoot_dti():

    foot = r'C:\Users\md1spsx\Downloads\diabetes_foot'

    dti = db.series(foot)[0]
    manuf = db.unique('Manufacturer', dti)[0]

    print(f'Manufacturer: {manuf}')

    bvals = 'DiffusionBValue'
    bvecs = 'DiffusionGradientOrientation'

    dims = ['SliceLocation', bvals, bvecs, 'InstanceNumber']
    vals = db.unique(dims, dti)
    
    print(f"Number of slices: {len(vals['SliceLocation'])}")
    print(f"Number of unique b-values: {len(vals[bvals])}")
    print(f"Number of unique b-vectors: {len(vals[bvecs])}")
    print(f"Number of images: {len(vals['InstanceNumber'])}")

    z0 = vals['SliceLocation'][0]
    b1 = vals[bvals][0]
    v1 = db.values(dti, bvecs, filter={bvals:b1, 'SliceLocation': z0})
    print(f"Number of b-vectors per slice for the first b-value: {len(v1)}")

    b2 = vals[bvals][1]
    v2 = db.values(dti, bvecs, filter={bvals:b2, 'SliceLocation': z0})
    print(f"Number of b-vectors per slice for the second b-value: {len(v2)}")

    print(f"Expected number of images: {len(vals['SliceLocation']) * (len(v1) + len(v2))}")
    print(f"Number of images: {len(vals['InstanceNumber'])}")


def ibeat_dti():
    ibeat = r'C:\Users\md1spsx\Downloads\iBEAt'

    dti = db.series(ibeat)[0]
    manuf = db.unique('Manufacturer', dti)[0]

    print(f'Manufacturer: {manuf}')

    bvals = (0x19, 0x100c)
    bvecs = (0x19, 0x100e)

    dims = ['SliceLocation', bvals, bvecs, 'InstanceNumber']
    vals = db.unique(dims, dti)
    
    print(f"Number of slices: {len(vals['SliceLocation'])}")
    print(f"Number of unique b-values: {len(vals[bvals])}")
    print(f"Number of unique b-vectors: {len(vals[bvecs])}")
    print(f"Number of images: {len(vals['InstanceNumber'])}")

    z0 = vals['SliceLocation'][0]
    b1 = vals[bvals][0]
    v1 = db.values(dti, bvecs, filter={bvals:b1, 'SliceLocation': z0})
    print(f"Number of b-vectors per slice for the first b-value: {len(v1)}")

    b2 = vals[bvals][1]
    v2 = db.values(dti, bvecs, filter={bvals:b2, 'SliceLocation': z0})
    print(f"Number of b-vectors per slice for the second b-value: {len(v2)}")

    print(f"Expected number of images: {len(vals['SliceLocation']) * (len(v1) + len(v2))}")
    print(f"Number of images: {len(vals['InstanceNumber'])}")

if __name__ == '__main__':

    # ibeat_dti()
    dfoot_dti()
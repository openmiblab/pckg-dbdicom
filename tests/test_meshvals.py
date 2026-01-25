# pytest test_meshvals.py

from itertools import product
import pytest
import numpy as np

from dbdicom.utils.arrays import meshvals



def test_meshvals():

    z = [0  ,  0 ,  0 ,  0 ,  1 ,  1 ,  1 ,  1 ]
    p = ['A', 'A', 'B', 'B', 'A', 'A', 'B', 'B']
    q = ['X', 'Y', 'X', 'Y', 'X', 'Y', 'X', 'Y']
    coords, inds = meshvals([z,p,q])
    assert coords[0].shape == (2,2,2)
    assert np.array_equal(coords[2][0,0,:], ['X', 'Y'])

    z = [0,1,2,0,1,2]
    p = ['A','A','A','B','B','B']
    coords, inds = meshvals([z,p])
    assert coords[0].shape == (3,2)
    assert np.array_equal(coords[1][0,:], ['A','B'])
    assert np.array_equal(coords[0][0,:], [0,0])
    assert np.array_equal(coords[0][:,0], [0,1,2])
    assert np.array_equal(coords[1][:,0], ['A','A','A'])
    assert np.array_equal(inds, [0,3,1,4,2,5])

    z = [0,1,2,0,1,2]
    p = [['A','B'],['A','B'],['A','B'],['C','D'],['C','D'],['C','D']]
    coords, inds = meshvals([z,p])
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
        meshvals([z,p])
    except:
        assert True
    else:
        assert False

    z = [0,1,2,0,1,3]
    p = ['A','A','A','B','B','B']
    try:
        meshvals([z,p])
    except:
        assert True
    else:
        assert False

    z = [0,1,2,0,1,2]
    p = ['A','A','A','B','B','A']
    try:
        meshvals([z,p])
    except:
        assert True
    else:
        assert False

    z = [0,1]
    p = ['A','B']
    coords, inds = meshvals([z,p])
    assert coords[0].shape == (2,1)
    assert np.array_equal(coords[0], [[0],[1]])
    assert np.array_equal(coords[1], [['A'],['B']])

    z = [1,0]
    p = ['A','B']
    coords, inds = meshvals([z,p])
    assert coords[0].shape == (2,1)
    assert np.array_equal(coords[0], [[0],[1]])
    assert np.array_equal(coords[1], [['B'],['A']])

    z = [0,1]
    p = [['X', 'Y'], ['Z', 'W']]
    coords, inds = meshvals([z,p])
    assert coords[0].shape == (2,1)
    assert np.array_equal(coords[0], [[0],[1]])
    assert np.array_equal(coords[0][0,0], 0)
    assert np.array_equal(coords[0][1,0], 1)
    assert np.array_equal(coords[1][0,0], ['X', 'Y'])
    assert np.array_equal(coords[1][1,0], ['Z', 'W'])

    z = [0,1]
    p = ['A','B']
    q = ['X', 'Y']
    coords, inds = meshvals([z,p,q])
    assert coords[0].shape == (2,1,1)
    assert coords[2][0,0,0] == 'X'
    assert coords[2][1,0,0] == 'Y'

    z = [0,1,2,0,1,2]
    p = ['A','B','C','D','E','F']
    coords, inds = meshvals([z,p])
    assert coords[0].shape == (3,2)
    assert np.array_equal(coords[1][0,:], ['A','D'])
    assert np.array_equal(coords[1][1,:], ['B','E'])


test_meshvals()

print('Passed!')
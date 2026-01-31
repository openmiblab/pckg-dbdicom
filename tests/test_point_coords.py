# pytest test_point_coords.py
import numpy as np
import pytest
from dbdicom.utils.points import point_coords

def test_point_coords():

    # Case 1: 2x2x2 Cube
    z = [0  ,  0 ,  0 ,  0 ,  1 ,  1 ,  1 ,  1 ]
    p = ['A', 'A', 'B', 'B', 'A', 'A', 'B', 'B']
    q = ['X', 'Y', 'X', 'Y', 'X', 'Y', 'X', 'Y']
    coords, inds = point_coords(list(zip(z,p,q)))
    assert coords[0].shape == (2,2,2)
    assert np.array_equal(coords[2][0,0,:], ['X', 'Y'])

    # Case 2: 3x2 Grid (Standard)
    z = [0,1,2,0,1,2]
    p = ['A','A','A','B','B','B']
    coords, inds = point_coords(list(zip(z,p)))
    assert coords[0].shape == (3,2)
    assert np.array_equal(coords[1][0,:], ['A','B'])
    assert np.array_equal(coords[0][0,:], [0,0])
    assert np.array_equal(coords[0][:,0], [0,1,2])
    assert np.array_equal(coords[1][:,0], ['A','A','A'])
    assert np.array_equal(inds, [0,3,1,4,2,5])

    # Case 3: List entries in coordinates
    z = [0,1,2,0,1,2]
    p = [['A','B'],['A','B'],['A','B'],['C','D'],['C','D'],['C','D']]
    coords, inds = point_coords(list(zip(z,p)))
    assert coords[0].shape == (3,2)
    assert np.array_equal(coords[1][0,0], ['A','B'])
    assert np.array_equal(coords[1][0,1], ['C','D'])
    assert np.array_equal(coords[0][0,:], [0,0])
    assert np.array_equal(coords[0][:,0], [0,1,2])
    # Note: Accessing object arrays might require slightly different indexing depending on numpy version, 
    # but the logic remains the same.
    assert np.array_equal(coords[1][0,0], ['A','B']) 
    assert np.array_equal(coords[1][1,0], ['A','B'])
    assert np.array_equal(coords[1][2,0], ['A','B'])
    assert np.array_equal(inds, [0,3,1,4,2,5])

    # Case 4: 2x1 Grid
    z = [0,1]
    p = ['A','B']
    coords, inds = point_coords(list(zip(z,p)))
    assert coords[0].shape == (2,1)
    assert np.array_equal(coords[0], [[0],[1]])
    assert np.array_equal(coords[1], [['A'],['B']])

    # Case 5: 2x1 Grid (Reversed input order)
    z = [1,0]
    p = ['A','B']
    coords, inds = point_coords(list(zip(z,p)))
    assert coords[0].shape == (2,1)
    assert np.array_equal(coords[0], [[0],[1]])
    assert np.array_equal(coords[1], [['B'],['A']])

    # Case 6: 2x1 Grid with List Objects
    z = [0,1]
    p = [['X', 'Y'], ['Z', 'W']]
    coords, inds = point_coords(list(zip(z,p)))
    assert coords[0].shape == (2,1)
    assert np.array_equal(coords[0], [[0],[1]])
    assert np.array_equal(coords[0][0,0], 0)
    assert np.array_equal(coords[0][1,0], 1)
    assert np.array_equal(coords[1][0,0], ['X', 'Y'])
    assert np.array_equal(coords[1][1,0], ['Z', 'W'])

    # Case 7: 2x1x1 Grid with Tuples
    z = [0,1]
    p = ['A','B']
    q = [(0, 1), (2, 3)]
    coords, inds = point_coords(list(zip(z,p,q)))
    assert coords[0].shape == (2,1,1)
    assert coords[2][0,0,0] == (0, 1)
    assert coords[2][1,0,0] == (2, 3)

    # Case 8: 3x2 Grid (Mixed chars)
    z = [0,1,2,0,1,2]
    p = ['A','B','C','D','E','F']
    coords, inds = point_coords(list(zip(z,p)))
    assert coords[0].shape == (3,2)
    assert np.array_equal(coords[1][0,:], ['A','D'])
    assert np.array_equal(coords[1][1,:], ['B','E'])

    # Case 9: 2x1x3 Grid
    z = [0,0,0,1,1,1]
    p = ['B', 'B', 'B', 'C','C', 'C']
    q = [(0, 1), (2, 3), (3, 4), (0, 1), (2, 3), (3, 4)]
    coords, inds = point_coords(list(zip(z,p,q)))
    assert coords[0].shape == (2,1,3)
    assert coords[2][0,0,2] == (3, 4)
    assert coords[1][0,0,2] == 'B'
    assert coords[0][0,0,2] == 0

    # Case 10: 2x1x3 Grid (Shuffled q input)
    z = [0,0,0,1,1,1]
    p = ['B', 'B', 'B', 'C','C', 'C']
    q = [(3, 4), (2, 3), (0, 1), (0, 1), (2, 3), (3, 4)]
    coords, inds = point_coords(list(zip(z,p,q)))
    assert coords[0].shape == (2,1,3)
    assert coords[2][0,0,2] == (3, 4)
    assert coords[1][0,0,2] == 'B'
    assert coords[0][0,0,2] == 0

    # Case 11: Fail - Inconsistent Grid (Dimensions don't match)
    # p 'A' maps to z 0,1,2 (3 items)
    # p 'B' maps to z 0,1 (2 items) -> Inconsistent splits
    z = [0,0,0,1,1,1]
    p = ['A', 'B', 'B', 'A','B', 'C']
    q = [(0, 1), (2, 3), (3, 4), (0, 1), (2, 3), (3, 4)]
    
    # We must zip cautiously here because the previous logic relied on lists of specific values.
    # The structure here is actually:
    # (0,A), (0,B), (0,B)... this creates duplicates and irregularity.
    with pytest.raises(ValueError):
         point_coords(list(zip(z,p,q)))

    # Case 12: Fail - Inconsistent Rectangularity
    # Block A has 3 items (z=0,1,2)
    # Block B has 2 items (z=0,1) -> z=3 was the sort key break in previous, here we simulate structure break
    z = [0,1,2,0,1]
    p = ['A','A','A','B','B']
    with pytest.raises(ValueError):
        point_coords(list(zip(z,p)))

    # Case 13: Fail - Another Structural Inconsistency
    # Block 'A' has z=[0,1,2] (3 items)
    # Block 'B' has z=[0,1,2,0] -> Technically 4 items, or duplicate (0,B)
    # If we duplicate (0,B), it changes the count and breaks grid shape inference.
    z = [0,1,2,0,1,2]
    p = ['A','A','A','B','B','A'] 
    # This creates (2,A) twice in the input if we sort? No:
    # (0,A), (1,A), (2,A), (0,B), (1,B), (2,A)
    # Sorted: (0,A), (1,A), (2,A), (2,A), (0,B), (1,B)
    # Block A has 4 points. Block B has 2 points.
    # Block A splits into z=0,1,2 (3 unique). Block B splits into z=0,1 (2 unique).
    # 3 != 2. Fails.
    with pytest.raises(ValueError):
        point_coords(list(zip(z,p)))

if __name__ == "__main__":
    test_point_coords()
    print('Passed!')
import numpy as np

from typing import List, Tuple

def to_array_list(values):
    arrays = []
    for v in values:
        if np.isscalar(v[0]):
            v_arr = np.array(v)
        else:
            v_arr = np.empty(len(v), dtype=object)
            v_arr[:] = v
        arrays.append(v_arr)
    return arrays


def meshvals(values) -> Tuple[List[np.ndarray], np.ndarray]:
    """
    Lexicographically sort flattened N coordinate arrays and reshape back to inferred grid shape,
    preserving original type of each input array.
    
    Parameters
    ----------
    *arrays : array-like
        Flattened coordinate arrays of the same length. Can be numbers, strings, or list objects.
    
    Returns
    -------
    sorted_arrays : list[np.ndarray]
        Coordinate arrays reshaped to inferred N-D grid shape, dtype/type preserved.
    indices : np.ndarray
        Permutation indices applied to the flattened arrays.
    shape : tuple[int, ...]
        Inferred grid shape (number of unique values per axis).
    """
    # Ensure the list elements are arrays
    arrays = to_array_list(values)

    # Remember original type/dtype for each array
    orig_types = [a.dtype if isinstance(a[0], np.ndarray) else type(a[0]) for a in arrays]

    # Convert non arrays to object arrays
    arrs = []
    for a in arrays:
        arrs_a = np.empty(len(a), dtype=object)
        arrs_a[:] = a
        arrs.append(arrs_a)
    
    # Stack arrays as columns (M x N)
    coords = np.stack(arrs, axis=1)
    
    # Lexicographic sort using structured array
    indices = np.lexsort(coords.T[::-1])
    sorted_coords = coords[indices]

    # Check that all coordinates are unique
    points = [tuple(col) for col in sorted_coords]
    if not all_elements_unique(points):
        raise ValueError(
            f"Improper coordinates. Coordinate values are not unique."
        )
    
    # Infer shape from unique values per axis
    shape = tuple(len(np.unique(sorted_coords[:, i])) for i in range(sorted_coords.shape[1]))
    
    # Check perfect grid
    if np.prod(shape) != sorted_coords.shape[0]:
        raise ValueError(
            f"Coordinates do not form a perfect Cartesian grid: inferred shape {shape} "
            f"does not match number of points {sorted_coords.shape[0]}"
        )
    
    # Split back into individual arrays and cast to original type
    sorted_arrays = []
    for i, orig_type in enumerate(orig_types):
        arr = sorted_coords[:, i]
        arr = arr.astype(orig_type).reshape(shape)
        sorted_arrays.append(arr)
    
    return sorted_arrays, indices


def all_elements_unique(items):
    """
    The most general uniqueness check, but also the slowest (O(n^2)).
    
    It works for ANY type that supports equality checking (==), including
    lists, dicts, and custom objects, without requiring them to be hashable.
    """
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            if items[i] == items[j]:
                return False
    return True



# def NEWmeshvals(coords):
#     stack_coords = [np.array(c, dtype=object) for c in coords]
#     stack_coords = np.stack(stack_coords)
#     mesh_coords, sorted_indices = _meshvals(stack_coords)
#     mesh_coords = [mesh_coords[d,...] for d in range(mesh_coords.shape[0])]
#     return mesh_coords, sorted_indices


# def _meshvals(coords):
#     # Input array shape: (d, f) with d = nr of dims and f = nr of frames
#     # Output array shape: (d, f1,..., fd)
#     if coords.size == 0:
#         return np.array([])
#     # Sort by column
#     sorted_indices = np.lexsort(coords[::-1])
#     sorted_array = coords[:, sorted_indices]
#     # Find shape
#     shape = _mesh_shape(sorted_array)  
#     # Reshape
#     mesh_array = sorted_array.reshape(shape)
#     return mesh_array, sorted_indices


# def _mesh_shape(sorted_array):
    
#     nd = np.unique(sorted_array[0,:]).size
#     shape = (sorted_array.shape[0], nd)

#     for dim in range(1,shape[0]):
#         shape_dim = (shape[0], np.prod(shape[1:]), -1)
#         sorted_array = sorted_array.reshape(shape_dim)
#         nd = [np.unique(sorted_array[dim,d,:]).size for d in range(shape_dim[1])]
#         shape = shape + (max(nd),)

#     if np.prod(shape) != sorted_array.size:
#         raise ValueError(
#             'Improper dimensions for the series. This usually means '
#             'that there are multiple images at the same location, \n or that '
#             'there are no images at one or more locations. \n\n'
#             'Make sure to specify proper dimensions when reading a pixel array or volume. \n'
#             'If the default dimensions of pixel_array (InstanceNumber) generate this error, '
#             'the DICOM data may be corrupted.'
#         ) 
    
#     return shape
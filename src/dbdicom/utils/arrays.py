from typing import List, Tuple
import numpy as np
import itertools
from functools import cmp_to_key




def _unique_preserve_order(values):
    out = []
    for v in values:
        if not any(_equal(v, u) for u in out):
            out.append(v)
    return out


def _compare_points(p1, p2):
    """Lexicographic comparison for mixed-type points."""
    for a, b in zip(p1, p2):
        if isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
            for x, y in zip(a, b):
                if x < y:
                    return -1
                elif x > y:
                    return 1
            if len(a) < len(b):
                return -1
            elif len(a) > len(b):
                return 1
        elif isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
            res = _compare_points(a, b)
            if res != 0:
                return res
        else:
            if a < b:
                return -1
            elif a > b:
                return 1
    return 0


def _has_duplicate_points(points):
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            if all(a == b for a, b in zip(points[i], points[j])):
                return i, j
    return None, None



def _equal(a, b):
    """Compare any type of object recursively: scalars, lists, tuples, arrays."""
    
    # Both are numpy arrays
    if isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
        return np.array_equal(a, b)
    
    # Both are lists or tuples
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        if len(a) != len(b):
            return False
        return all(_equal(x, y) for x, y in zip(a, b))
    
    # Fallback to normal equality
    return a == b







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



def duplicates(values):
    """
    List indices of duplicate points.

    Works for ANY type of coordinate values (scalars, lists, tuples, arrays).
    
    Args:
        values (list of lists): coordinate axes of equal length

    Returns:
        list: indices of points that are duplicates
    """
    points = list(zip(*values))
    duplicates_idx = []

    for i, p in enumerate(points):
        # Compare with all previous points
        if any(_equal(p, points[j]) for j in range(i)):
            duplicates_idx.append(i)

    return duplicates_idx







# Original version
def _meshvals(values) -> Tuple[List[np.ndarray], np.ndarray]:
    """
    Lexicographically sort flattened N coordinate arrays and reshape back to inferred grid shape,
    preserving original type of each input array.
    
    Parameters
    ----------
    *arrays : array-like
        Flattened coordinate arrays of the same length. All elements of a coordinate array have 
        the same type but different coordinate arrays can have differen types. All python data types 
        can be coordinats, including lists or tuples.
    
    Returns
    -------
    sorted_arrays : list[np.ndarray]
        N Coordinate arrays reshaped to inferred N-D grid shape, dtype/type preserved.
    indices : np.ndarray
        Permutation indices applied to the flattened arrays.
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
    if not _all_elements_unique(points):
        raise ValueError(
            f"Improper coordinates. Coordinate values are not unique."
        )
    
    # Infer shape from unique values per axis
    shape = tuple(len(np.unique(sorted_coords[:, i])) for i in range(sorted_coords.shape[1]))
    
    # Check perfect grid
    if np.prod(shape) != sorted_coords.shape[0]:
        raise ValueError(
            f"Coordinates do not form a rectangular grid: inferred shape {shape} "
            f"does not match number of points {sorted_coords.shape[0]}"
        )
    
    # Split back into individual arrays and cast to original type
    sorted_arrays = []
    for i, orig_type in enumerate(orig_types):
        arr = sorted_coords[:, i]
        arr = arr.astype(orig_type).reshape(shape)
        sorted_arrays.append(arr)
    
    return sorted_arrays, indices



def _all_elements_unique(items):
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




import numpy as np
import itertools
from typing import Tuple, List, Any

def meshvals(values: List[Any]) -> Tuple[List[np.ndarray], np.ndarray]:
    """
    Lexicographically sort flattened N coordinate arrays and reshape back to 
    inferred grid shape.
    
    Handles coupled coordinates (non-Cartesian) by inferring shape hierarchically.
    Example: z=[0,1], p=['A','B'] -> shape (2, 1)
    """
    n_points = len(values[0])
    arrays = []
    orig_dtypes = []
    
    # 1. Normalize inputs to object arrays safely
    for v in values:
        if len(v) != n_points:
            raise ValueError("All input arrays must have the same length.")
            
        if isinstance(v, np.ndarray):
            orig_dtypes.append(v.dtype)
            arr_obj = v.astype(object)
        else:
            orig_dtypes.append(type(v[0])) 
            arr_obj = np.empty(n_points, dtype=object)
            arr_obj[:] = list(v)
            
        arrays.append(arr_obj)

    coords = np.stack(arrays, axis=1)

    # 2. Lexicographic sort (Cols: 0 -> N)
    # We transpose and reverse so lexsort uses col 0 as the primary key
    try:
        indices = np.lexsort(coords.T[::-1])
    except TypeError as e:
        raise TypeError("Coordinate elements must be comparable") from e
        
    sorted_coords = coords[indices]

    # 3. Hierarchical Shape Inference
    # We trace down the columns. At each level, we check how the current
    # column subdivides the blocks defined by previous columns.
    
    # Start with one block covering the whole array [0, N]
    block_boundaries = [0, n_points] 
    shape = []
    
    for i in range(sorted_coords.shape[1]):
        col = sorted_coords[:, i]
        
        new_boundaries = [0]
        splits_per_block = []
        
        # Iterate over the blocks defined by the PREVIOUS dimension
        for j in range(len(block_boundaries) - 1):
            start, end = block_boundaries[j], block_boundaries[j+1]
            segment = col[start:end]
            
            # Group by value in this segment (effectively finding unique values)
            # groupby is safe for unhashable types
            segment_splits = 0
            seg_idx = start
            for _, group in itertools.groupby(segment):
                count = sum(1 for _ in group)
                seg_idx += count
                new_boundaries.append(seg_idx)
                segment_splits += 1
                
            splits_per_block.append(segment_splits)
            
        # 4. Check Rectangularity
        # For a valid grid, every block from dim (i-1) must split into 
        # the SAME number of sub-blocks in dim (i).
        if len(set(splits_per_block)) != 1:
            raise ValueError(
                f"Coordinates do not form a rectangular grid. "
                f"Dimension {i} has inconsistent sizes across the grid."
            )
            
        dim_size = splits_per_block[0]
        shape.append(dim_size)
        
        # Update boundaries for the next dimension
        block_boundaries = new_boundaries

    shape = tuple(shape)
    
    # Final sanity check: Product of shape must equal number of points
    if np.prod(shape) != n_points:
         raise ValueError(f"Inferred shape {shape} does not match N points {n_points}")

    # 5. Reshape and Restore
    sorted_arrays = []
    for i, dtype_or_type in enumerate(orig_dtypes):
        arr = sorted_coords[:, i]
        reshaped_arr = arr.reshape(shape)
        
        if isinstance(dtype_or_type, np.dtype):
            reshaped_arr = reshaped_arr.astype(dtype_or_type)
        
        sorted_arrays.append(reshaped_arr)

    return sorted_arrays, indices



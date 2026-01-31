from typing import List, Tuple, Any, Union
import numpy as np
import itertools



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


def to_array(v) -> np.ndarray:
    if np.isscalar(v[0]):
        v_arr = np.array(v)
    else:
        v_arr = np.empty(len(v), dtype=object)
        v_arr[:] = v
    return v_arr




def duplicate_points(points):
    """
    List indices of duplicate points by sorting first.
    
    Handles generic types (numpy arrays, lists, tuples) by converting
    them to a sortable 'key' format internally.
    """
    
    # Helper: Recursively convert arrays/lists to hashable/sortable tuples.
    def make_sortable(obj):
        if isinstance(obj, np.ndarray):
            # .tolist() converts numpy hierarchy to pure python lists/scalars
            return make_sortable(obj.tolist())
        if isinstance(obj, (list, tuple)):
            return tuple(make_sortable(x) for x in obj)
        return obj

    # 1. Build a list of (sort_key, original_index)
    #    We pre-compute keys to avoid expensive recursion during the sort itself.
    indexed_points = []
    for i, p in enumerate(points):
        indexed_points.append((make_sortable(p), i))

    # 2. Sort. Python's sort is stable, so original order is preserved for ties.
    #    We sort primarily by the key (item[0]).
    indexed_points.sort(key=lambda item: item[0])

    duplicates_idx = []

    # 3. Linear scan to find duplicates
    for i in range(1, len(indexed_points)):
        curr_key, curr_idx = indexed_points[i]
        prev_key, _ = indexed_points[i-1]

        # Since we sorted by the key, identical items are now adjacent
        if curr_key == prev_key:
            duplicates_idx.append(curr_idx)

    return sorted(duplicates_idx)




def point_coords(points: Union[List[Any], np.ndarray]) -> Tuple[List[np.ndarray], np.ndarray]:
    """
    Takes a list of N points (each of dimension D) and reconstructs the 
    underlying orthogonal grid structure.
    
    Args:
        points: List or array of shape (N, D). 
                Example: [(x1, y1), (x2, y2), ...]

    Returns:
        sorted_arrays: List of D arrays (one per dimension), reshaped to the grid.
        indices: The indices that sort the input points into the grid order.
    """
    if len(points) == 0:
        return [], np.array([])

    # 0. Transpose Points to Columns (Axes)
    # We convert (N points x D dims) -> (D dims x N points) to handle types per axis
    if isinstance(points, np.ndarray):
        n_points, n_dim = points.shape
        # Split numpy array into list of columns to preserve distinct dtypes per column
        values = [points[:, i] for i in range(n_dim)]
    else:
        n_points = len(points)
        # zip(*points) effectively transposes a list of tuples/lists
        values = list(zip(*points))

    # --- Below is the original logic (adapted to use 'values') ---

    arrays = []
    orig_dtypes = []
    
    # 1. Normalize inputs to object arrays safely
    # This ensures we can sort mixed types (e.g. Ints and Strings) together
    for v in values:
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
    # Transpose and reverse so lexsort uses col 0 as the primary key
    try:
        indices = np.lexsort(coords.T[::-1])
    except TypeError as e:
        raise TypeError("Coordinate elements must be comparable") from e
        
    sorted_coords = coords[indices]

    # 3. Hierarchical Shape Inference
    block_boundaries = [0, n_points] 
    shape = []
    
    for i in range(sorted_coords.shape[1]):
        col = sorted_coords[:, i]
        new_boundaries = [0]
        splits_per_block = []
        
        # Iterate over blocks defined by the PREVIOUS dimension
        for j in range(len(block_boundaries) - 1):
            start, end = block_boundaries[j], block_boundaries[j+1]
            segment = col[start:end]
            
            # Count unique values in this segment to find sub-blocks
            segment_splits = 0
            seg_idx = start
            for _, group in itertools.groupby(segment):
                count = sum(1 for _ in group)
                seg_idx += count
                new_boundaries.append(seg_idx)
                segment_splits += 1
                
            splits_per_block.append(segment_splits)
            
        # 4. Check Rectangularity
        if len(set(splits_per_block)) != 1:
            raise ValueError(
                f"Coordinates do not form a rectangular grid. "
                f"Dimension {i} has inconsistent sizes across the grid."
            )
            
        dim_size = splits_per_block[0]
        shape.append(dim_size)
        block_boundaries = new_boundaries

    shape = tuple(shape)
    
    if np.prod(shape) != n_points:
         raise ValueError(f"Inferred shape {shape} does not match N points {n_points}")

    # 5. Reshape and Restore Types
    sorted_arrays = []
    for i, dtype_or_type in enumerate(orig_dtypes):
        arr = sorted_coords[:, i]
        reshaped_arr = arr.reshape(shape)
        
        if isinstance(dtype_or_type, np.dtype):
            reshaped_arr = reshaped_arr.astype(dtype_or_type)
        
        sorted_arrays.append(reshaped_arr)

    return sorted_arrays, indices




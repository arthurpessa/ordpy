import numpy as np
from typing import List

def lbb(series, b:int, B:float, size:int) -> List[List[float]]:
    """
    Local block bootstrap (lbb) method \\ [#martins2023]_ .
    It is guaranteed to estimate the mean well if the series is
    locally stationary.

    Parameters
    ----------
        series (np.ndarray): original series.
        b (int): block size. Must be between 0 and `len(series)`.
        B (float): size of the neighborhood. Must be between 0 and 1.
        size (int): how many bootstrap replicas you want.

    Raises:
        ValueError: whenever b or B are out of bounds.

    Returns:
        list[list[float]]: representing an array of bootstrap samples from the original.
    """
    n = len(series)
    if not (0 <= b <= n):
        raise ValueError("'b' argument not valid.")
    if not (0 <= B <= 1):
        raise ValueError("'B' argument not valid.")
    
    nB = int(n * B)
    q = int(n / b) - 1
    k_range = np.arange(-nB, nB + 1)

    result_boot = []
    for _ in range(size):
        new_sample = []
        k = np.random.choice(k_range, size=q+1, replace=True, p=np.repeat(1 / (2 * nB + 1), 2 * nB + 1))
        for i in range(q+1):
            for j in range(1, b+1):
                k_i = k[i]
                neighbor_index = j + i*b + k_i
                if not (1 <= neighbor_index <= n):
                    neighbor_index = j + i*b - k_i
                new_sample.append(series[neighbor_index - 1])
        result_boot.append(new_sample)
    return result_boot
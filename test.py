import numpy as np

def load_dem(path, width, height, endian="<"):
    dtype = np.dtype(endian + "i2")  # 16-bit signed
    n = width * height
    data = np.fromfile(path, dtype=dtype, count=n)
    if data.size != n:
        raise ValueError(f"Expected {n} samples, got {data.size}")
    return data.reshape(height, width)

# Try loading your 5x5 DEM
arr = load_dem("W100N40.bin", 5, 5, "<")

print("DEM shape:", arr.shape)
print("DEM values with coordinates:")
for y in range(arr.shape[0]):
    for x in range(arr.shape[1]):
        print(f"(row={y}, col={x}) -> {arr[y,x]}")

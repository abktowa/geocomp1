import numpy as np
import sys

def load_dem(path, width, height, endian="<"):
    dtype = np.dtype(endian + "i2")  # 16-bit signed
    n = width * height
    data = np.fromfile(path, dtype=dtype, count=n)
    if data.size != n:
        raise ValueError(f"Expected {n} samples, got {data.size}. Check width/height.")
    return data.reshape(height, width)

def main():
    if len(sys.argv) < 4:
        print("Usage: python check_dem.py <binfile> <width> <height> [endian]")
        return
    
    path   = sys.argv[1]
    width  = int(sys.argv[2])
    height = int(sys.argv[3])
    endian = sys.argv[4] if len(sys.argv) > 4 else "<"

    arr = load_dem(path, width, height, endian)
    print(f"DEM shape: {arr.shape} (rows, cols)")

    # print first few values (row 0)
    print("First row values:", arr[0, :10])

    # print a sample in the middle
    mid_y, mid_x = arr.shape[0] // 2, arr.shape[1] // 2
    print(f"Middle value at (row={mid_y}, col={mid_x}):", arr[mid_y, mid_x])

if __name__ == "__main__":
    main()

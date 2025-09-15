import argparse
import numpy as np
import sys

#input = dem grid of 6000 x 4800
#output = list of 100 peaks (prominence, peak_x, peak_y, peak_height, col_x, col_y, col_height)
global mountain_ranges
prominences = []


# Simplifies 2D coordinates from DEM grid to 1D index
def build_peaks(array, W, H):
    global mountain_ranges
    for i in range(W):
        for j in range(H):
            try:
                mountain_ranges.make_set((W * j) + i, Peak(i,j, array[i][j]))
            except:
                print("Failure at " + str(i) + "," + str(j))


# Array to check surrounding cells
neighbors = [(-1,-1), (-1,0), (-1,1), (0,-1), (0,1), (1,-1), (1,0), (1,1)]

# Container for coordinates and elevation of each peak
class Peak:
    def __init__(self, x, y, height):
        self.x = x
        self.y = y
        self.height = height

# Union-Find Data Structure
class UniFi:
    # Sets up DSU, each cell is its own island
    def __init__(self, n):
        self.parent = list(range(n)) # Each cell is own parent
        self.rank = [0] * n # Keep trees balanced when merging, attaches shorter tree to taller one
        self.peak = [None] * n # Stores highest point for each mountain

    # Finds mountain each cell belongs to
    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])  # path compression
        return self.parent[x]
    
    # Creates a set for a new cell
    def make_set(self, x, peak):
        self.parent[x] = x
        self.rank[x] = 0
        self.peak[x] = peak

    # Combines two cells/islands when they "touch", use rank to decide which one is superior
    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return rx

        # Union by rank
        if self.rank[rx] < self.rank[ry]:
            rx, ry = ry, rx
        elif self.rank[rx] == self.rank[ry]:
            self.rank[rx] += 1

        self.parent[ry] = rx

        # Update peak for the merged set
        peak_rx = self.peak[rx]
        peak_ry = self.peak[ry]
        if peak_ry and (not peak_rx or peak_ry.height > peak_rx.height):
            self.peak[rx] = peak_ry

        return rx
"""
# Algorithm
# sweeping line from high to low 
for (height, x, y) in cells:

    idx = simplify(x, y,6000)   # unique index for DSU
    active[y][x] = True

    # Create a new set for this cell
    DSU.make_set(idx)
    DSU.set_peak(idx, peak_x=x, peak_y=y, peak_height=height)

    for each (nx, ny) in 8_neighbors(x, y):
        if inside_grid(nx, ny) and active[ny][nx]:
            n_idx = simplify(nx, ny, 6000)
            root1 = DSU.find(idx)
            root2 = DSU.find(n_idx)

            if root1 != root2:
                # Each component knows its highest peak
                peak1 = DSU.peak[root1]
                peak2 = DSU.peak[root2]

                if peak1.height > peak2.height:
                    higher, lower = root1, root2
                else:
                    higher, lower = root2, root1

                # Col height is current sweep elevation
                col_height = height

                # Prominence of the lower peak
                prominence = DSU.peak[lower].height - col_height

                record = (prominence,
                          DSU.peak[lower].x,
                          DSU.peak[lower].y,
                          DSU.peak[lower].height,
                          x, y, col_height)

                results.append(record)

                DSU.union(higher, lower)

# parse tallest peak
# After all unions, the highest peak never got merged.
# Its prominence is simply its height, col info is N/A.
global_peak = DSU.peak[ DSU.find(any_active_idx) ]
results.append( (global_peak.height,
                 global_peak.x,
                 global_peak.y,
                 global_peak.height,
                 None, None, None) )

# sort and then output
sort results by prominence DESC
return top 100
"""

def load_dem(path, width, height, endian="<"):
    """
    Load raw 16-bit signed integers into a (H,W) array.
    endian: '<' little-endian, '>' big-endian
    """
    dtype = np.dtype(endian + "i2")  # 16-bit signed
    n = width * height
    data = np.fromfile(path, dtype=dtype, count=n)
    if data.size != n:
        raise ValueError(f"Expected {n} samples, got {data.size}. Check width/height.")
    return data.reshape(height, width)

def local_maxima_8(arr):
    """
    Simple 8-neighborhood peak detector.
    Returns a boolean mask where True = strict local maximum.
    (Equal neighbors do *not* count as higher.)
    """
    H, W = arr.shape
    # pad with -inf so edges can still be compared
    neg_inf = np.iinfo(arr.dtype).min
    P = np.pad(arr, 1, mode="constant", constant_values=neg_inf)

    c  = P[1:H+1, 1:W+1]
    n8 = [
        P[0:H,   0:W  ], P[0:H,   1:W+1], P[0:H,   2:W+2],
        P[1:H+1, 0:W  ],                   P[1:H+1, 2:W+2],
        P[2:H+2, 0:W  ], P[2:H+2, 1:W+1], P[2:H+2, 2:W+2],
    ]
    higher_neighbor = np.zeros_like(c, dtype=bool)
    for nb in n8:
        higher_neighbor |= (nb > c)
    peaks = ~higher_neighbor
    return peaks


def main():
    global xs, ys, sorted_vals, mountain_ranges
    ap = argparse.ArgumentParser(description="Read and inspect a raw 16-bit DEM tile.")
    ap.add_argument("binfile", type=str)
    ap.add_argument("--width",  type=int, required=True)
    ap.add_argument("--height", type=int, required=True)
    ap.add_argument("--endian", choices=["<", ">"], default="<",
                    help="Byte order: '<' little-endian, '>' big-endian (default: <)")
    ap.add_argument("--preview", type=str, default="preview.png",
                    help="Path to save a quick preview image.")
    ap.add_argument("--list-peaks", action="store_true",
                    help="Detect local maxima and print top 10 by elevation.")
    args = ap.parse_args()

    arr = load_dem(args.binfile, args.width, args.height, endian=args.endian)
    mountain_ranges = UniFi(args.width * args.height)
    build_peaks(arr, args.width, args.height)

    if args.list_peaks:
        peaks_mask = local_maxima_8(arr)
        ys, xs = np.where(peaks_mask)
        # Take top 10 peaks by elevation
        if ys.size > 0:
            vals = arr[ys, xs]
            sorted_vals = np.argsort(vals)[0:][::-1]
            print("\nTop 1000 local peaks (y, x, elevation):")
            for i in sorted_vals:
                y, x, h = int(ys[i]), int(xs[i]), int(vals[i])
                print(f"  ({y}, {x}) -> {h}")
        else:
            print("No local maxima found (unexpected).")

if __name__ == "__main__":
    main()


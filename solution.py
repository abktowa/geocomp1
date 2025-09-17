"""
########################################
Program Name: Mountain Prominence Calculator
Class: CSC615 Computational Geometry
Program Authors: Aban Khan and Matthew Glennon
Program Description: This program is a prominence calculator for geographical data, being able to convert bin file coordinates and heights to 
a cell grid, and use a sweep-line algorithm to first sort the peaks by height, then slowly activate them from top to bottom.  As new points are revealed,
the program checks surrounding cells to look for other mountain ranges, and union them together.  Upon hitting two mountain ranges at once in this way, we 
calculate prominence for the mountain range with the smaller peak height, then union them together with the taller peak as their overall peak.  Finally,
the program displays the top 100 prominences
########################################
"""
import argparse
import numpy as np
import sys

global mountain_ranges, active_cells
prominences = []
neighbors = [(-1,-1), (-1,0), (-1,1), (0,-1), (0,1), (1,-1), (1,0), (1,1)]

# Simplifies 2D coordinates from DEM grid to 1D index
def build_peaks(array, W, H):
    global mountain_ranges
    for y in range(H):
        for x in range(W):
            try:
                idc = simplify_to_1d(x, y, W)
                h   = array[y][x]
                mountain_ranges.make_set(idc, Peak(x, y, W, h))
            except:
                print("Failure at " + str(x) + "," + str(y))


#Method which can take a 2D coordinate, and transfer it to a 1D index
def simplify_to_1d(x, y, w):
    return ( y * w ) + x


#Checks to make sure if a coordinate is within a given grid or not
def in_bounds(x, y, W, H):
    return 0 <= x < W and 0 <= y < H


#Method that checks surrounding neighbors to a given cell, in this case a newly activated cell within our sweep line algorithm, to see if there are other active mountain ranges in the surrounding 8 neighbors
def check_neighbors(start_cell, W, H, current_height):
    global mountain_ranges, active_cells

    x, y = start_cell.x, start_cell.y

    # Gather roots of all active neighbor components
    neighbor_roots = []
    for dx, dy in neighbors:
        nx, ny = x + dx, y + dy
        if in_bounds(nx, ny, W, H) and active_cells[ny][nx]:   # row-major
            n_idc = simplify_to_1d(nx, ny, W)
            root = mountain_ranges.find(n_idc)
            neighbor_roots.append(root)

    # Deduplicate roots while preserving order
    if neighbor_roots:
        seen = set()
        uniq = []
        for r in neighbor_roots:
            if r not in seen:
                seen.add(r)
                uniq.append(r)
        neighbor_roots = uniq

    # No active neighbors → nothing to merge
    if not neighbor_roots:
        return

    # Exactly one neighbor component → just attach this cell to it (no col event)
    if len(neighbor_roots) == 1:
        mountain_ranges.union(start_cell.set_index, neighbor_roots[0], None)
        return

    # ≥2 distinct neighbor components → this cell is a col
    col_peak = Peak(x, y, W, current_height)  # col location & elevation

    # First, merge all neighbor components together *at this col*
    merged_root = neighbor_roots[0]
    for other in neighbor_roots[1:]:
        if mountain_ranges.find(merged_root) != mountain_ranges.find(other):
            merged_root = mountain_ranges.union(merged_root, other, col_peak)

    # Finally, attach the col cell's own component into the merged mountain
    mountain_ranges.union(merged_root, start_cell.set_index, None)


# Container for coordinates and elevation of each peak
class Peak:
    def __init__(self, x, y, width, height):
        self.x = x
        self.y = y
        self.set_index = simplify_to_1d(x, y, width)
        self.height = height

# Union-Find Data Structure
class UniFi:
    # Sets up DSU, each cell is its own island/parent
    def __init__(self, n):
        self.size = n
        self.parent = list(range(n))
        self.rank = [0] * n
        self.peak = [None] * n
        self.records = []

    # Finds mountain each cell belongs to
    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]
    
    # Creates a set for a new cell
    def make_set(self, x, peak):
        self.parent[x] = x
        self.rank[x] = 0
        self.peak[x] = peak

    # Combines two cells/islands when they "touch", use rank to decide which one is superior
    def union(self, x, y, col):
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return rx

        # Read peaks BEFORE changing parents/ranks
        peak_rx = self.peak[rx]
        peak_ry = self.peak[ry]

        # Decide higher vs lower by peak height
        if peak_rx.height >= peak_ry.height:
            higher_root, lower_root = rx, ry
            higher_peak, lower_peak = peak_rx, peak_ry
        else:
            higher_root, lower_root = ry, rx
            higher_peak, lower_peak = peak_ry, peak_rx

        # If this merge happens at a col, lock the lower peak's prominence
        if col is not None:
            prom = lower_peak.height - col.height
            self.records.append(
                (prom,
                lower_peak.x, lower_peak.y, lower_peak.height,
                col.x,        col.y,        col.height)
            )

        # Union by rank
        rx, ry = higher_root, lower_root
        if self.rank[rx] < self.rank[ry]:
            rx, ry = ry, rx
        elif self.rank[rx] == self.rank[ry]:
            self.rank[rx] += 1
        self.parent[ry] = rx

        # Winner peak is the higher one we computed
        self.peak[rx] = higher_peak

        return rx


# Testing on differing data sets required using endian modifications, in-case of variation of data.
def load_dem(path, width, height, endian="<"):
    dtype = np.dtype(endian + "i2")
    n = width * height
    data = np.fromfile(path, dtype=dtype, count=n)
    if data.size != n:
        raise ValueError(f"Expected {n} samples, got {data.size}. Check width/height.")
    return data.reshape(height, width)


# Detects neighborhood peaks. True = local maximum.
def local_maxima_neighbor(arr):
    H, W = arr.shape
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

#Main method
def main():
    global xs, ys, sorted_vals, mountain_ranges, active_cells
    active_cells = []
    
    #First collect args variables
    ap = argparse.ArgumentParser(description="Read and inspect DEM tile.")
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

    #Load in our data, set it up as a UniFi variable called mountain_ranges
    arr = load_dem(args.binfile, args.width, args.height, endian=args.endian)
    mountain_ranges = UniFi(args.width * args.height)
    print("Setting up peaks")
    build_peaks(arr, args.width, args.height)

    H, W = args.height, args.width

    #Set up our active cells, to be used when checking neighbors
    cells = []
    for y in range(H):
        for x in range(W):
            h   = arr[y][x]
            idc = simplify_to_1d(x, y, W)
            cells.append((h, x, y, idc))

    # Sorting tallest to shortest (one at a time)
    cells.sort(key=lambda t: t[0], reverse=True)

    active_cells = [[False] * W for _ in range(H)]

    # (MINE)SWEEP(ER), ONE CELL AT A TIME.  This joke was brought by Aban, thank you Aban.
    i = 0
    while i < len(cells):
        h = cells[i][0]
        j = i
        while j < len(cells) and cells[j][0] == h:
            j += 1

        # "Look" at all cells at this height
        for _, x, y, idc in cells[i:j]:
            active_cells[y][x] = True

        # Union for neighbors of equal-height (we merge this into one component)
        for _, x, y, idc in cells[i:j]:
            for dx, dy in neighbors:
                nx, ny = x + dx, y + dy
                if in_bounds(nx, ny, W, H) and active_cells[ny][nx]:
                    if arr[ny][nx] == h:  # only same-height neighbors
                        n_idc = simplify_to_1d(nx, ny, W)
                        mountain_ranges.union(idc, n_idc, None)

        # Handling for merges across different components
        for _, x, y, idc in cells[i:j]:
            start_cell = mountain_ranges.peak[idc]
            check_neighbors(start_cell, W, H, h)

        i = j

    # Add one record for the global tallest peak
    h_max, x_max, y_max, idc_max = cells[0]
    root_max = mountain_ranges.find(idc_max)
    gpeak = mountain_ranges.peak[root_max]

    mountain_ranges.records.append(
        (gpeak.height, gpeak.x, gpeak.y, gpeak.height, None, None, None)
    )

    nonzero_records = [rec for rec in mountain_ranges.records if rec[0] > 0]
    # Sort and print top 100
    top = sorted(nonzero_records, key=lambda r: r[0], reverse=True)[:100]

    print("Peaks by prominence:")
    print("  prom    row    col   elev   crow   ccol  celev")
    print("--------------------------------------------------")
    for prom, px, py, pz, cx, cy, cz in top:
        row, col = py, px            # row = y, col = x
        if cx is None:
            # Tallest peak: no col
            print(f"{prom:6d} {row:6d} {col:6d} {pz:6d} {'NA':>6} {'NA':>6} {'NA':>6}")
        else:
            # Note: crow = cy (row), ccol = cx (col)
            print(f"{prom:6d} {row:6d} {col:6d} {pz:6d} {cy:6d} {cx:6d} {cz:6d}")


if __name__ == "__main__":
    main()


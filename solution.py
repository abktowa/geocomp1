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

def simplify_to_1d(x, y, w):
    return ( y * w ) + x

def in_bounds(x, y, W, H):
    return 0 <= x < W and 0 <= y < H

# def check_neighbors(start_cell, grid_width, grid_height):
#     global mountain_ranges, active_cells
#     # Array to check surrounding cells
#     hits = []
#     for coordinates in neighbors:
#         if not ((start_cell.x == 0 and neighbors[coordinates][0] == -1) or (start_cell.x == (grid_width - 1) and neighbors[coordinates][0] == 1) or (start_cell.y == 0 and neighbors[coordinates][1] == -1) or (start_cell.y == (grid_height - 1) and neighbors[coordinates][1] == 1)):
#             if(active_cells[start_cell.x + neighbors[coordinates][0]][start_cell.y + neighbors[coordinates][1]] == True):
#                 hits.append(simplify_to_1d(start_cell.x + neighbors[coordinates][0], start_cell.y + neighbors[coordinates][1], grid_width))

#     if len(hits) > 0:
#         #Only one mountain of either equal or greater height detected. Combine those mountain ranges in UniFi
#         if len(hits) == 1:
#             mountain_ranges.union(start_cell.set_index, hits[0], None)
#         #If more than one mountain fulfills that role, that means this cell is a col, for at least 2 separate mountain ranges.  First, we unionize each surrounding mountain, then add in the col last.
#         if len(hits) >= 2:
#             new_mountain = None
#             new_mountain = mountain_ranges.union(hits[0],hits[1], start_cell)
#             for i in range(len(hits) - 1):
#                 if new_mountain == None:
#                     if mountain_ranges.find(hits[i]) != mountain_ranges.find(hits[i+1]):
#                         new_mountain = mountain_ranges.union(hits[i],hits[i+1], start_cell)
#                 else:
#                     if mountain_ranges.find(new_mountain) != mountain_ranges.find(hits[i+1]):
#                         new_mountain = mountain_ranges.union(new_mountain,hits[i+1], start_cell)
#             mountain_ranges.union(new_mountain, start_cell, None)

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


def main():
    global xs, ys, sorted_vals, mountain_ranges, active_cells
    active_cells = []
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

    arr = load_dem(args.binfile, args.width, args.height, endian=args.endian)
    mountain_ranges = UniFi(args.width * args.height)
    print("Setting up peaks")
    build_peaks(arr, args.width, args.height)

    #Make cell block to track if a given cell has been activated or not, which will procedurally turned on as height decreases
    # for i in range(args.width):
        # if i % 100 == 0:
            # print("Now writing at x " + str(i))
        # active_cells.append([])
        # for j in range(args.height):
            # active_cells[i].append(False)

    #Slowly 'lower the water' to reveal
    # for i in range(6500):
        # if i % 100 == 0:
            # print("Now checking at depth " + str(6500 - i))
        # height_req = 6500 - i
        # for j in range(mountain_ranges.size):
            # if mountain_ranges.peak[j].height >= height_req:
                # active_cells[mountain_ranges.peak[j].x][mountain_ranges.peak[j].y] = True
                # check_neighbors(mountain_ranges.peak[j], args.width, args.height)
        
    H, W = args.height, args.width

    cells = []
    for y in range(H):
        for x in range(W):
            h   = arr[y][x]
            idc = simplify_to_1d(x, y, W)
            cells.append((h, x, y, idc))

    # Sorting tallest to shortest (one at a time)
    cells.sort(key=lambda t: t[0], reverse=True)

    active_cells = [[False] * W for _ in range(H)]

    # (MINE)SWEEP(ER), ONE CELL AT A TIME.
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

    # Sort and print top 100
    top = sorted(mountain_ranges.records, key=lambda r: r[0], reverse=True)[:100]

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


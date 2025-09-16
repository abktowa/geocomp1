import argparse
import numpy as np
import sys

#input = dem grid of 6000 x 4800
#output = list of 100 peaks (prominence, peak_x, peak_y, peak_height, col_x, col_y, col_height)
global mountain_ranges, active_cells
prominences = []


# Simplifies 2D coordinates from DEM grid to 1D index
def build_peaks(array, W, H):
    global mountain_ranges
    for i in range(W):
        for j in range(H):
            try:
                mountain_ranges.make_set((W * j) + i, Peak(i,j, W, array[i][j]))
            except:
                print("Failure at " + str(i) + "," + str(j))


def convert_to_1d(x, y, w):
    return (y*w) + x

def check_neighbors(start_cell, grid_width, grid_height):
    global mountain_ranges, active_cells
    # Array to check surrounding cells
    neighbors = [(-1,-1), (-1,0), (-1,1), (0,-1), (0,1), (1,-1), (1,0), (1,1)]
    hits = []
    for coordinates in neighbors:
        if not ((start_cell.x == 0 and neighbors[coordinates][0] == -1) or (start_cell.x == (grid_width - 1) and neighbors[coordinates][0] == 1) or (start_cell.y == 0 and neighbors[coordinates][1] == -1) or (start_cell.y == (grid_height - 1) and neighbors[coordinates][1] == 1)):
            if(active_cells[start_cell.x + neighbors[coordinates][0]][start_cell.y + neighbors[coordinates][1]] == True):
                hits.append(convert_to_1d(start_cell.x + neighbors[coordinates][0], start_cell.y + neighbors[coordinates][1], grid_width))

    if len(hits) > 0:
        #Only one mountain of either equal or greater height detected. Combine those mountain ranges in UniFi
        if len(hits) == 1:
            mountain_ranges.union(start_cell.set_index, hits[0], None)
        #If more than one mountain fulfills that role, that means this cell is a col, for at least 2 separate mountain ranges.  First, we unionize each surrounding mountain, then add in the col last.
        if len(hits) >= 2:
            new_mountain = None
            new_mountain = mountain_ranges.union(hits[0],hits[1], start_cell)
            for i in range(len(hits) - 1):
                if new_mountain == None:
                    if mountain_ranges.find(hits[i]) != mountain_ranges.find(hits[i+1]):
                        new_mountain = mountain_ranges.union(hits[i],hits[i+1], start_cell)
                else:
                    if mountain_ranges.find(new_mountain) != mountain_ranges.find(hits[i+1]):
                        new_mountain = mountain_ranges.union(new_mountain,hits[i+1], start_cell)
            mountain_ranges.union(new_mountain, start_cell, None)


# Container for coordinates and elevation of each peak
class Peak:
    def __init__(self, x, y, width, height):
        self.x = x
        self.y = y
        self.set_index = convert_to_1d(x, y, width)
        self.height = height



# Union-Find Data Structure
class UniFi:
    # Sets up DSU, each cell is its own island
    def __init__(self, n):
        self.size = n
        self.parent = list(range(n)) # Each cell is own parent
        self.rank = [0] * n # Keep trees balanced when merging, attaches shorter tree to taller one
        self.peak = [Peak(-1,-1,-1,-1)] * n # Stores highest point for each mountain
        self.prominences = {}

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
    def union(self, x, y, col):
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
            if peak_rx and col:
                self.prominences[(peak_rx.x, peak_rx.y)] = peak_rx.height - col.height

            self.peak[rx] = peak_ry

        elif peak_ry and col: 
            self.prominences[(peak_ry.x, peak_ry.y)] = peak_ry.height - col.height

        return rx


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
    return data.reshape(width, height)

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
    global xs, ys, sorted_vals, mountain_ranges, active_cells
    active_cells = []
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
    print("Setting up peaks")
    build_peaks(arr, args.width, args.height)

    #Make cell block to track if a given cell has been activated or not, which will procedurally turned on as height decreases
    for i in range(args.width):
        if i % 100 == 0:
            print("Now writing at x " + str(i))
        active_cells.append([])
        for j in range(args.height):
            active_cells[i].append(False)

    #Slowly 'lower the water' to reveal
    for i in range(6500):
        if i % 100 == 0:
            print("Now checking at depth " + str(6500 - i))
        height_req = 6500 - i
        for j in range(mountain_ranges.size):
            if mountain_ranges.peak[j].height >= height_req:
                active_cells[mountain_ranges.peak[j].x][mountain_ranges.peak[j].y] = True
                check_neighbors(mountain_ranges.peak[j], args.width, args.height)
        
    print("Now checking from sea level")
    for j in range(mountain_ranges.size):
        if (mountain_ranges.peak[mountain_ranges.find(j)].x, mountain_ranges.peak[mountain_ranges.find(j)].y) not in mountain_ranges.prominences and mountain_ranges.peak[mountain_ranges.find(j)].height > 0:
            mountain_ranges.prominences[(mountain_ranges.peak[mountain_ranges.find(j)].x, mountain_ranges.peak[mountain_ranges.find(j)].y)] = mountain_ranges.peak[mountain_ranges.find(j)].height
    

    sorted_prominences = sorted(mountain_ranges.prominences.items(), key=lambda x: x[1], reverse=True)


    for key, value in sorted_prominences[:100]:
        print(f"{key}: {value}")

if __name__ == "__main__":
    main()


#input = dem grid of 6000 x 4800
#output = list of 100 peaks (prominence, peak_x, peak_y, peak_height, col_x, col_y, col_height)

# Simplifies 2D coordinates from DEM grid to 1D index
def simplify(x, y, W):
    return y * W + x

# Array to check surrounding cells
neighbors = [(-1,-1), (-1,0), (-1,1), (0,-1), (0,1), (1,-1), (1,0), (1,1)]

# Container for coordinates and elevation of each peak
class Peak:
    def __init__(self, x, y, height):
        self.x = x
        self.y = y
        self.height + height

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

# Algorithm
# sweeping line from high to low 
for (height, x, y) in cells:

    idx = flatten(x, y)   # unique index for DSU
    active[y][x] = True

    # Create a new set for this cell
    DSU.make_set(idx)
    DSU.set_peak(idx, peak_x=x, peak_y=y, peak_height=height)

    for each (nx, ny) in 8_neighbors(x, y):
        if inside_grid(nx, ny) and active[ny][nx]:
            n_idx = flatten(nx, ny)
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

#input = dem grid of 6000 x 4800
#output = list of 100 peaks (prominence, peak_x, peak_y, peak_height, col_x, col_y, col_height)

# set up
cells = []
for y in 0..H-1:
    for x in 0..W-1:
        cells.append( (grid[y][x], x, y) )

sort cells by height DESCENDING

DSU = new UnionFind(H * W)
active[H][W] = all False
results = []

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
                 NA, NA, NA) )

# sort and then outpu
sort results by prominence DESC
return top 100

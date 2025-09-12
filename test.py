#!/usr/bin/env python3
import argparse
import numpy as np
import sys
from pathlib import Path

global xs, ys, sorted_vals

mountain_ranges = []
prominences = []

def getProm(peak, col_height):  #peak represents the second highest peak in a NEW mountain range, made by the combination of two smaller mountain ranges.
    global xs, ys, sorted_vals
    peak


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

def print_quick_checks(arr, n=12):
    flat = arr.ravel()
    print(f"Array shape: {arr.shape}, dtype: {arr.dtype}")
    print(f"Min/Max: {int(arr.min())} / {int(arr.max())}")
    print(f"First {n} values:", flat[:n].tolist())

def save_preview(arr, out_png, step=10):
    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        print("matplotlib not available; skipping preview.", file=sys.stderr)
        return
    # Light downsample for speed
    thumb = arr[::step, ::step]
    plt.figure()
    plt.imshow(thumb)  # default colormap (no explicit colors)
    plt.title(f"Preview (every {step}th pixel)")
    plt.colorbar(label="elevation")
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    print("Saved preview to:", out_png)

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
    global xs, ys, sorted_vals
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
    print_quick_checks(arr)

    # sanity check for the assignment's early numbers:
    print("First 4 (should start 691, 692, 699, 696,... if endian is correct):",
          arr.ravel()[:4].tolist())

    save_preview(arr, args.preview)

    if args.list_peaks:
        peaks_mask = local_maxima_8(arr)
        ys, xs = np.where(peaks_mask)
        # Take top 10 peaks by elevation
        if ys.size > 0:
            vals = arr[ys, xs]
            sorted_vals = np.argsort(vals)[0:][::-1]
            print("\nTop 1000 local peaks (y, x, elevation):")
            for i in sorted_vals[0:1000]:
                y, x, h = int(ys[i]), int(xs[i]), int(vals[i])
                print(f"  ({y}, {x}) -> {h}")
        else:
            print("No local maxima found (unexpected).")

if __name__ == "__main__":
    main()

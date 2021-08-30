#!/usr/bin/env python

import yt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str, help="Path to input file to slice.")
parser.add_argument("-f", "--field", type=str, default="u", help="Name of field to plot. Default is u.")
parser.add_argument("-ax", "--axis", type=str, default="z", help="Axis to slice across (x, y, or z). Default is z.")
parser.add_argument("-min", "--field_min", type=float, default=1.0, help="Minimum value of field to plot.")
parser.add_argument("-max", "--field_max", type=float, default=-1.0, help="Maximum value of field to plot.")

args = parser.parse_args()

if __name__ == "__main__":
    ds = yt.load(args.infile)
    
    field_to_plot = args.field
    for cat, name in ds.field_list:
        if name == args.field:
            field_to_plot = (cat, name)

    slice = yt.SlicePlot(ds, args.axis, field_to_plot)

    slice.set_log(field_to_plot, False)
    if args.field_min < args.field_max:
        slice.set_zlim(args.field, args.field_min, args.field_max)

    slice.save("{}.png".format(args.infile))

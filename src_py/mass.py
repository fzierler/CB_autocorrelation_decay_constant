#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
import numpy as np

from bootstrap import get_rng, sample_bootstrap_1d

def get_args():
    parser = ArgumentParser(
        description="Compute the mass and matrix element from correlators in an HDF5 file"
    )
    parser.add_argument("h5file", help="The file to read")
    parser.add_argument(
        "--ensemble_name",
        default=None,
        help="Name of the ensemble to analyse.",
    )
    parser.add_argument(
        "--plateau_start",
        type=int,
        default=None,
        help="Time slice at which plateau starts",
    )
    parser.add_argument(
        "--plateau_end", type=int, default=None, help="Time slice at which plateau ends"
    )
    parser.add_argument(
        "--output_file_mean",
        type=FileType("w"),
        default="-",
        help="Where to output the mean and uncertainty of meson mass. (Defaults to stdout.)",
    )
    parser.add_argument(
        "--channel",
        choices=["ps", "v", "t", "av", "at", "s"],
        default=None,
        help="Measuring channel",
    )
    return parser.parse_args()


def get_correlator_samples(
    ensemble,
    measurement,
):
    C = ensemble[measurement][:, :]

    return sample_bootstrap_1d(C.T, get_rng(ensemble.name))


def get_channel_tags(ch):
    return {
        "ps": ["g5"],
        "v": ["g1", "g2", "g3"],
        "t": ["g0g1", "g0g2", "g0g3"],
        "av": ["g5g1", "g5g2", "g5g3"],
        "at": ["g0g5g1", "g0g5g2", "g0g5g3"],
        "s": ["id"],
    }.get(ch, ch)


def fold_correlators(C):
    return (C + np.roll(np.flip(C, axis=1), 1, axis=1)) / 2


def fold_correlators_cross(C):
    C_fold = (C - np.roll(np.flip(C, axis=1), 1, axis=1)) / 2

    C_fold[:, 0] = C[:, 0]

    return C_fold

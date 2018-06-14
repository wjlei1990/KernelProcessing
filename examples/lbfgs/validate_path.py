from __future__ import print_function, division
import os
import sys


def load_txt(fn):
    with open(fn) as fh:
        return [line.rstrip() for line in fh]


def main(path_file):
    print("Check path_file:", path_file)
    content = load_txt(path_file)

    niter = int(content[0])
    if len(content) != (2 * niter + 2):
        raise ValueError("niter and number of files are not consistent")

    idx = 0
    for i in range(0, niter):
        idx += 1
        grad_file = content[idx]
        if not os.path.exists(grad_file):
            raise ValueError("Missing grad_file: ", grad_file)

        idx += 1
        dkernel_file = content[idx]
        if not os.path.exists(dkernel_file):
            raise ValueError("Missing dkernel_file: ", dkernel_file)

    idx += 1
    grad_file = content[idx]
    if not os.path.exists(grad_file):
        raise ValueError("Missing current iter grad_file: ", grad_file)

    print("Pass Check...")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise ValueError("Usage: python validate_path.py path_file")
    main(sys.argv[1])

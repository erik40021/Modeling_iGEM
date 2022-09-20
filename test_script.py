#!/usr/bin/env python

import sys


def run_script():
    pass


def main():
    args = sys.argv

    print('Number of arguments:', len(args), 'arguments.')
    print ('Argument List:', str(args))

    num_required_args = 3

    if len(args) < num_required_args + 1:
        print('[Error]: number of arguments required: ', num_required_args + 1)
        return

    run_script()    # the actual functionality



if __name__ == "__main__":
    main()
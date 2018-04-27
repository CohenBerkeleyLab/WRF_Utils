from __future__ import print_function

import argparse
import datetime as dt
import os
import sys


def shell_error(msg, exit_code=1):
    print(msg, file=sys.stderr)
    exit(exit_code)


def parse_wrf_date(datestr):
    return dt.datetime.strptime(datestr, '%Y-%m-%d_%H:%M:%S')


def parse_timestep(stepstr):
    vals = [int(x) for x in stepstr.split(':')]
    while len(vals) < 0:
        vals.append(0)

    return dt.timedelta(hours=vals[0], minutes=vals[1], seconds=vals[2])


def iter_date(start, end, step):
    curr_date = start
    while curr_date <= end:
        yield curr_date
        curr_date += step


def parse_args():
    parser = argparse.ArgumentParser(description='List missing WRF output files',
                                     epilog='This program will exit with code 0 if all files were found and code 1 if '
                                            'any files were missing. Other errors will be indicated by higher exit codes.')
    parser.add_argument('start_date', type=parse_wrf_date, help='First expected date in yyyy-mm-dd_HH:MM:SS format')
    parser.add_argument('end_date', type=parse_wrf_date, help='Last expected date in yyyy-mm-dd_HH:MM:SS format')
    parser.add_argument('-t', '--time-step', type=parse_timestep, default=dt.timedelta(hours=1),
                        help='Increment between files in HH, HH:MM, or HH:MM:SS format. Default is %(default)s')
    parser.add_argument('-d', '--dir', default='.', help='Directory to search. Default is the current directory.')
    parser.add_argument('-n', '--domain', default=1, type=int, help='Domain number to check')
    parser.add_argument('-f', '--filename', default='wrfout_d{domain:02}_%Y-%m-%d_%H:%M:%S',
                        help='The format of the WRF output filenames. For each file, it is first formatted with the '
                             '.format(domain=d) call (to insert the domain number) then passed through '
                             'datetime.strftime to insert the date. Default is %(default)s.')

    args = parser.parse_args()

    if not os.path.isdir(args.dir):
        shell_error('Argument for "--dir" is not a valid directory', exit_code=2)

    return args


def main():
    args = parse_args()
    files_missing = False
    for date in iter_date(args.start_date, args.end_date, args.time_step):
        wrf_fname = os.path.join(args.dir, date.strftime(args.filename.format(domain=args.domain)))
        if not os.path.isfile(wrf_fname):
            files_missing = True
            print(wrf_fname, 'MISSING')

    if files_missing:
        exit(1)
    else:
        exit(0)


if __name__ == '__main__':
    main()

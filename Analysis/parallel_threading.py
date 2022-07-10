#!/usr/bin/env python
from __future__ import print_function

import sys
import shlex
import subprocess
import multiprocessing as mp
from multiprocessing.pool import ThreadPool
from tqdm import tqdm
# import numpy as np
# import time

def get_list_of_commands_from_file(file_path):
    result = list()
    with open(file_path) as file:
        for line in file:
            if line.empty(): continue
            result.append(line.strip('\n'))
    return result

def call_process(command, pbar):
    p = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    pbar.update(1)
    return (out, err)

def call_function(func, args_tuple, pbar):
    result = func(**args_tuple)
    pbar.update(1)
    return result

def run_with_pool(commands=None, function=None, sets_of_args=None, callback=None, max_workers=None): # commands: list(); sets_of_args: list of dicts (keys = args of function)
    if max_workers == None:
        max_workers = min(mp.cpu_count()-1, 12)
    print('Used number of threads:', max_workers)
    pool = ThreadPool(max_workers)
    pbar = None # tqdm progress bar
    results = list()
    if commands != None:
        print('Jobs to be spawned:', len(commands))
        pbar = tqdm(total=len(commands), desc='Jobs done', dynamic_ncols=True, leave=False)
        for cmd in commands:
            results.append(pool.apply_async(func=call_process, args=(cmd, pbar,)))
    elif function != None and sets_of_args != None:
        print('Jobs to be spawned:', len(sets_of_args))
        pbar = tqdm(total=len(sets_of_args), desc='Jobs done', dynamic_ncols=True, leave=False)
        for args_tuple in sets_of_args:
            results.append(pool.apply_async(func=call_function, args=(function, args_tuple, pbar,)))
    else:
        sys.exit('run_with_pool(): You need to specify "commands" or "function"+"sets_of_args" arguments.')
    pool.close()
    pool.join()
    pbar.close()
    # for result in results:
    #     out, err = result.get()
    #     print('out: {} err: {}'.format(out, err))
    print('Finished.')
    return results

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage:', 'python parallel_threading.py <file-with-one-command-per-line.txt>')
        sys.exit(0)
    commands = get_list_of_commands_from_file(sys.argv[1])
    run_with_pool(commands=commands)

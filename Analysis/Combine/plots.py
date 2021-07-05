#!/usr/bin/env python
import sys
import shlex
import subprocess
import multiprocessing as mp
from multiprocessing.pool import ThreadPool
from tqdm import tqdm
import numpy as np
import time
from tqdm import tqdm

def get_list_of_commands_from_file(file_path):
    result = list()
    with open(file_path) as file:
        for line in file:
            result.append(line.strip('\n'))
    return result

def call_process(command, pbar):
    p = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    pbar.update(1)
    return (out, err)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage:', 'python plots.py <file-with-plotting-commands.txt>'
        sys.exit(0)
    max_workers = min(mp.cpu_count()-1, 12)
    print 'Used number of threads:', max_workers
    pool = ThreadPool(max_workers)
    commands = get_list_of_commands_from_file(sys.argv[1])
    print 'Jobs to be spawned:', len(commands)
    results = list()
    pbar = tqdm(total=len(commands), desc='Jobs done', dynamic_ncols=True, leave=False)
    for cmd in commands:
        results.append(pool.apply_async(call_process, (cmd, pbar,)))
    pool.close()
    pool.join()
    pbar.close()
    # for result in results:
    #     out, err = result.get()
    #     print 'out: {} err: {}'.format(out, err)
    print 'Finished.'
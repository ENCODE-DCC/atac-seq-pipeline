from multiprocessing import Pool
import subprocess
import random
import time
import shlex
import os, sys
from encode_common import *
from encode_common_genomic import *

# def func(i1,i2):
#     # print('* started target {} {}'.format(i1,i2))
#     # time.sleep(5)
#     # print('* done target {} {}'.format(i1,i2))
#     # return 'target {} {}'.format(i1,i2)
#     cmd = 'echo -e "{} {}"; sleep 5'.format(i1,i2)
#     subprocess.Popen(shlex.split(cmd),stdout=PIPE,stderr=PIPE)

# pool = Pool(2) # maximum two processes at time.
# ret1 = pool.apply_async(func,(1,2))
# ret2 = pool.apply_async(func,(3,4))
# ret3 = pool.apply_async(func,(5,6))
# ret4 = pool.apply_async(func,(7,8))
# print("hello")
# pool.close()
# pool.join()
# print("hello2")

# def test_log():
#     log.setLevel('DEBUG')
#     log.warning('WARN')
#     log.info('INFO')
#     log.debug('X')

# def test1():
#     log.setLevel('DEBUG')
    
#     pool = Pool(2)
#     # cmd1 = 'zcat -f /srv/scratch/shared/wotan/leepc12/run/ENCODE3_NEW2/data/ENCSR176BYZ/in_progress/reads/fastq/rep1/pair1/ENCFF071AXB.fastq.gz | gzip -nc > test1.fastq.gz'
#     # cmd2 = 'zcat -f /srv/scratch/shared/wotan/leepc12/run/ENCODE3_NEW2/data/ENCSR176BYZ/in_progress/reads/fastq/rep1/pair2/ENCFF628CBN.fastq.gz | gzip -nc > test2.fastq.gz'

#     cmd1 = 'echo 2; sleep 1; echo 2; sleep 1; echo 2; sleep 1; echo 2; sleep 1; echo 2; sleep 1; echo 2; sleep 1;'
#     cmd2 = 'echo 2; sleep 1; echo 2; sleep 1; echo 2; sleep 1; echo 2; sleep 1; echo 2; sleep 1; echo 2; sleep 1;'

#     ret_val1 = pool.apply_async(run_shell_cmd,(cmd1,))
#     ret_val2 = pool.apply_async(run_shell_cmd,(cmd2,))
#     print(ret_val1.get())
#     print(ret_val2.get())

#     pool.close()
#     pool.join()

# test1()

#def test3():
#    path = 'out_rep2_macs2/ENCFF463QCX.trim.merged.nodup.tn5.narrowPeak.gz'
#    print(get_ext(path))
#    print(strip_ext(path))
#    print(strip_ext(path,'narrowP'))
#    print(strip_ext(path,'narrowPeak'))

#test3()

def test4():
    log.setLevel('DEBUG')
    log.info("INFO DA")
    log.error("ERROR DA")
    run_shell_cmd('ls -l')
    run_shell_cmd('asdfasdf')
    
    
    
    

test4()

# run_shell_cmd("echo 2; echo 2; echo 3; sleep 10")

# subprocess.Popen()
#     non blocking
#         .communicate()

# multiprocessig.Pool()
#     non blocking
#         .get()

# python function
#     single-threaded
#         detect_adapter

# shell command
#     single-threaded
#         idr
#         macs2
#     multi-threaded    
#         bowtie2
#         samtools sort
#         spp(xcor)

# temp.reserve_threads(nth=1)

# TaskManager
#     dict[task_name,Task]
# Task
#     subprocess.Popen
#     multiprocessing.Process

# TaskOutput
#     how to wait
#         subprocess.Popen.communicate()


# trimmed_fastq = 


# EncodeFile

# class EncodeMultiThreadedTask:
#     def __init__(self):
#         pass
#     def 


# samtools_flagstat(

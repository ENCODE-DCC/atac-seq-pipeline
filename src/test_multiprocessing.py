from multiprocessing import Pool
import random
import time


def func(i1,i2):
    print('* started target {} {}'.format(i1,i2))
    time.sleep(5)
    print('* done target {} {}'.format(i1,i2))
    return 'target {} {}'.format(i1,i2)

pool = Pool(2) # maximum two processes at time.
ret1 = pool.apply_async(func,(1,2))
ret2 = pool.apply_async(func,(3,4))
ret3 = pool.apply_async(func,(5,6))
ret4 = pool.apply_async(func,(7,8))
print("hello")
pool.close()
pool.join()
print("hello2")
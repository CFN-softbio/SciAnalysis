#!/usr/bin/python3

#import sys
#code_PATH='./'
#code_PATH in sys.path or sys.path.append(code_PATH)
from CustomS3 import Queue_measure as queue
q = queue()

def loop(wait_time=300, max_iterations=None):
    import time
    
    iterations = 0
    while max_iterations is None or iterations<max_iterations:
        q.msg('Looping update iteration {}'.format(iterations), 2, 0)
        iterations += 1
        
        q.get_status_files()
        
        time.sleep(wait_time)


if __name__=='__main__':
    q.get_status_files()
    #loop()
    

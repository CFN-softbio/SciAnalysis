#!/usr/bin/python3
from pathlib import Path

#import sys
#code_PATH='../../'
#code_PATH in sys.path or sys.path.append(code_PATH)
from CustomS3 import Queue_decision as queue
q = queue()

def send(decision=True, measure=False, analyze=True):

    if decision:
        p = Path('../gpcam/')
        files = sorted( p.glob('Data_*.npy'), key=lambda x: x.stat().st_mtime )
        q.send = 'decision'
        q.publish_status_file(files[-1])
        q.publish_status_file(p/'decisionS3-received.npy')
        q.publish_status_file(p/'decisionS3-sent.npy')

    if measure:
        p = Path('../')
        q.send = 'measure'
        q.publish_status_file(p/'measureS3-received.npy')
        q.publish_status_file(p/'measureS3-sent.npy')

    if analyze:
        p = Path('../saxs/analysis/')
        q.send = 'analyze'
        q.publish_status_file(p/'analyzeS3-received.npy')
        q.publish_status_file(p/'analyzeS3-sent.npy')



def loop(wait_time=300, max_iterations=None):
    import time
    
    iterations = 0
    while max_iterations is None or iterations<max_iterations:
        q.msg('Looping update iteration {}'.format(iterations), 2, 0)
        iterations += 1

        send()

        time.sleep(wait_time)


if __name__=='__main__':
    send()
    #loop()


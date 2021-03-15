import os
import time
import queue
import threading
import subprocess


def funcThread(work_queue, tid=0, returnCodeQueue=None):
    '''Runs python function jobs received from work_queue'''
    # TODO - switch to runinng the functions using the multiprocessing library
    #  to avoid GIL global thread interpreter lock contention
    print('start thread')
    try:
        while job := work_queue.get(block=False):
            args = job.get('args', ())
            kwargs = job.get('kwargs', {})
            func = job.get('func', None)
            if func is None:
                raise ValueError('Function not defined')
            print(f'Thread {tid}: func: {func.__name__}, args: {args}, kwargs: {kwargs}')
            func(*args, **kwargs)
    except queue.Empty:
        pass
    return


def jobThread(work_queue, tid=0, returnCodeQueue=None, outputQueue=None):
    '''Runs command line jobs received from the work_queue'''
    print('start thread')
    captureOutput = False
    if outputQueue is not None:
        captureOutput = True
    try:
        while cmd := work_queue.get(block=False):
            useShell = False
            if type(cmd) is str:
                useShell = True
            print(f'Thread {tid}: cmd: {cmd}')
            ret = subprocess.run(cmd, shell=useShell, capture_output=captureOutput)
            if returnCodeQueue is not None:
                returnCodeQueue.put(ret.returncode)
            if outputQueue is not None:
                outputQueue.put(ret.stdout.decode('utf-8'))
    except queue.Empty:
        pass
    return
    

def runParallelJobs(cmdlist, concurrency=2, isPyFunction=False):
    '''
    Start up concurrency num threads all sharing a work queue.
    Each thread dequeues an item from the list and runs it in a subprocess
    Ehen the queue is empty the threads exit
    Join all the threads here and exit
    '''
    if type(cmdlist[0]) is dict and isPyFunction is False:
        # cmdlist will be a set of dictionaries when using python functions
        # cmdlist will be a set of lists when using bash commands
        raise ValueError("cmdlist contains dictionaries when isPyFunction is False")

    # put the commands in a queue
    work_queue = queue.Queue()
    for cmd in cmdlist:
        work_queue.put(cmd)

    target = None
    checkReturnCode = False
    if isPyFunction:
        target = funcThread
    else:
        target = jobThread
        checkReturnCode = True

    returnCodeQueue = None
    if checkReturnCode is True:
        returnCodeQueue = queue.Queue()

    threads = []
    for idx in range(concurrency):
        jthread = threading.Thread(name=f'thread_{idx}',
                                   target=target,
                                   args=(work_queue,),
                                   kwargs={'tid': idx,
                                           'returnCodeQueue': returnCodeQueue,})
        jthread.setDaemon(True)
        jthread.start()
        threads.append(jthread)
    
    for jthread in threads:
        jthread.join()
    
    if checkReturnCode is True:
        for i in range(len(cmdlist)):
            code = returnCodeQueue.get(block=False)
            if code != 0:
                raise RuntimeError(f"Exit code {code} for cmd: {cmdlist[i]}")

    print("Complete")
    return
     

def test_func(a, b=0):
    print(f'a {a}, b {b}')
    time.sleep(b)


if __name__=="__main__":
    # test the parallel jobs
    cmdlist = [['sleep', '2'], 'sleep 2', ['sleep', '2'], 'sleep 2', ['sleep', '2'], 'sleep 2']
    startTime = time.time()
    runParallelJobs(cmdlist, concurrency=2)
    print('parallel jobs time: {}s'.format(time.time() - startTime))

    # test parallel funcs
    cmdlist = [{'func': test_func, 'args': ('hello',), 'kwargs': {'b': 2}},
               {'func': test_func, 'args': ('world',), 'kwargs': {'b': 2}},
               {'func': test_func, 'args': ('test',), 'kwargs': {'b': 2}}]
    # cmdlist = [{'func': test_func, 'args': ('hello',), 'kwargs': {'b': 2}}]

    startTime = time.time()
    runParallelJobs(cmdlist, concurrency=len(cmdlist), isPyFunction=True)
    print('parallel funcs time: {}s'.format(time.time() - startTime))

    print("Test Complete")

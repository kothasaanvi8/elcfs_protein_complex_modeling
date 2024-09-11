#Generates notifications (e.g., iPhone, etc.) based on Pushover, refer to support.pushover.net/
import requests, time
from multiprocessing import Process, Event


class statusCheck(object):
    #example call: cmdReport = alertMe.statusCheck(pushoverAPI, pushoverKey_user)
    #must set up an account and generate an application token, api, and use account's associated user key
    def __init__(self, api, key,
                 url='https://api.pushover.net/1/messages.json'):

        self.api = api
        self.key = key
        self.url = url
        self.waitTime = None
        self.func = None
        self.funcInputs = None

    #sends notification
    def release(self, title, msg):
        cPigeon = \
            {'token': self.api, 'user': self.key,
             'title': title, 'message': msg}
        requests.post(self.url, data=cPigeon)

    #generates notification upon completion of job
    def finishPush(self):
        #example call: cmdReport.finishPush()
        title = 'executionSuccess'
        msg = 'done'
        self.release(title, msg)

    #generates notification that command is still running periodically
    # as defined by wait time until target function completes
    def runningPush(self, event):
        title = 'statusUpdate'
        msg = 'command running...'

        exeTime = time.time() + self.waitTime
        while True:
            if time.time() > exeTime:
                self.release(title, msg)
                exeTime+=self.waitTime
            if event.is_set():
                break

    #wrapper function for target function
    def envelopFcn(self, event):
        self.func(self.funcInputs)
        event.set()

    #uses multiprocessing to run notification and target function concurrently
    def progress(self, waitTime, func, funcInputs):
        #example call: cmdReport.progress(10, fcn, [arg1, arg2, ..., argN])
        self.waitTime = waitTime
        self.func = func
        self.funcInputs = funcInputs
        event = Event()

        p1 = Process(target=self.runningPush, args=(event,))
        p2 = Process(target=self.envelopFcn, args=(event,))
        p1.start()
        p2.start()
        p1.join()
        p2.join()

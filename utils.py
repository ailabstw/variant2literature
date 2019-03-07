import threading
import functools


class FuncThread(threading.Thread):
    def __init__(self, target, args, kwargs):
        threading.Thread.__init__(self)
        self.target = target
        self.args = args
        self.kwargs = kwargs
        self.result = None

    def run(self):
        self.result = self.target(*self.args, **self.kwargs)


def timeout(time):
    def _timeout(fun):
        @functools.wraps(fun)
        def wrapper(*args, **kwargs):
            t = FuncThread(target=fun, args=args, kwargs=kwargs)
            t.start()
            t.join(time)
            if t.isAlive():
                raise TimeoutError
            return t.result
        return wrapper
    return _timeout

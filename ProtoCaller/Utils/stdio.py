import os as _os
import sys as _sys
import warnings as _warnings

class stdout_stderr_cls:
    def __init__(self, stdout=_os.devnull, stderr=_os.devnull):
        _sys.stdout.flush()
        _sys.stderr.flush()
        self.stdout_prev, self.stderr_prev = _os.dup(1), _os.dup(2)
        self.stdout, self.stderr = stdout, stderr

    def __enter__(self):
        self.out = open(self.stdout, "w")
        self.err = open(self.stderr, "w")
        _os.dup2(self.out.fileno(), 1)
        _os.dup2(self.err.fileno(), 2)

    def __exit__(self, type, value, traceback):
        self.out.close()
        self.err.close()
        _sys.stdout.flush()
        _sys.stderr.flush()
        _os.dup2(self.stdout_prev, 1)
        _os.dup2(self.stderr_prev, 2)

def stdout_stderr(stdout=_os.devnull, stderr=_os.devnull):
    """A decorator which redirects the stdout and stderr to a custom stream. Default is suppression of all output."""
    def decorator(f):
        def f_mod(*args, **kwargs):
            _sys.stdout.flush()
            _sys.stderr.flush()
            stdout_prev, stderr_prev = _os.dup(1), _os.dup(2)
            with open(stdout, "w") as out, open(stderr, "w") as err:
                _os.dup2(out.fileno(), 1)
                _os.dup2(err.fileno(), 2)
                try:
                    result = f(*args, **kwargs)
                finally:
                    _sys.stdout.flush()
                    _sys.stderr.flush()
                    _os.dup2(stdout_prev, 1)
                    _os.dup2(stderr_prev, 2)
            return result
        return f_mod
    return decorator

def warnings_as_errors(f):
    """A decorator which casts all warnings as errors."""
    def f_mod(*args, **kwargs):
        with _warnings.catch_warnings():
            _warnings.filterwarnings("error")
            result = f(*args, **kwargs)
        return result
    return f_mod

def ignore_warnings(f):
    """A decorator which suppresses all warnings."""
    def f_mod(*args, **kwargs):
        with _warnings.catch_warnings():
            _warnings.filterwarnings("ignore")
            result = f(*args, **kwargs)
        return result
    return f_mod

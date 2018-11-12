import os as _os
import sys as _sys
import warnings as _warnings

def stdout_stderr(stdout=_os.devnull, stderr=_os.devnull):
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
    def f_mod(*args, **kwargs):
        _warnings.filterwarnings("error")
        try:
            result = f(*args, **kwargs)
        finally:
            _warnings.filterwarnings("always")
        return result
    return f_mod

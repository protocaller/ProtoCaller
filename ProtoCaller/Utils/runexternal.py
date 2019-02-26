import re as _re
import subprocess as _subprocess
import sys as _sys


def runExternal(*commands, procname=None, output_filebase=None):
    """
    A wrapper for subprocess.check_call.

    Parameters
    ----------
    commands
        Positional arguments of type str. These are the shell commands to be executed.
    procname : str or None, optional
        The name of the process. None is the first word of the first executed command.
    output_filebase : str or None, optional
        The file base of the output logs. None is equal to the procname.
    """
    try:
        if(procname in ["", None]): procname = _re.search(r"([\w]*)\s", commands[0]).group(1)
    except:
        if(procname in ["", None]): procname = "UNDEFINED PROCESS"
    if(output_filebase is None): output_filebase = procname

    stdout = open("%s.out" % output_filebase, "a")
    stderr = open("%s.err" % output_filebase, "a")
    try:
        print("Running %s... " % procname, end="")
        _sys.stdout.flush()
        _subprocess.check_call(commands, shell=True, stdout=stdout, stderr=stderr)
        print("OK")
    except _subprocess.CalledProcessError:
        print("FAILED")
        raise OSError("Error while calling command '%s'" % procname)
    stdout.close()
    stderr.close()

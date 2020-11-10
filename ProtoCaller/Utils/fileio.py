import atexit as _atexit
import os as _os
import shutil as _shutil


class Dir:
    """
    A class that handles directory changes in Python in a more intuitive way.

    Parameters
    ----------
    dirname : str
        Initialises dirname.
    copydirname : str or None
        Initialises copydirname.
    overwrite : bool
        Initialises overwrite.
    temp : bool
        Initialises temp.
    purge_immediately : bool
        Initialise purge_immediately.


    Attributes
    ----------
    workdirname : str
        The absolute path to the parent directory of Dir.
    dirname : str
        The directory name.
    copydirname : str or None
        If not None, copies the directory in copydirname as a basis for Dir.
    path : str
        Full absolute path of the directory.
    overwrite : bool
        Whether to overwrite any directories with the same name as Dir.
    temp : bool
        If True, the directory is deleted upon __exit__.
    purge_immediately : bool
        If True, the directory is only deleted upon program exit. Only valid if temp is True.
    initialdirnames : [str]
        Keeps track of different places from which __enter__ has been called in order to avoid mistakes.
    """
    def __init__(self, dirname, copydirname=None, overwrite=False, temp=False, purge_immediately=True):
        self.path = _os.path.normpath(dirname)
        if not _os.path.isabs(self.path):
            self.path = "%s/%s" % (_os.getcwd(), self.path)
        self.workdirname = _os.path.dirname(self.path)
        self.dirname = _os.path.basename(self.path)

        self.copydirname = _os.path.normpath(copydirname)
        if self.copydirname and not _os.path.isabs(self.copydirname):
            self.copydirname = "%s/%s" % (_os.getcwd(), self.copydirname)

        self.overwrite = overwrite
        self.temp = temp
        self.purge_immediately = purge_immediately
        self.initialdirnames = []

    def __enter__(self):
        self.initialdirnames += [_os.getcwd()]
        if not _os.path.exists(self.workdirname):
            _os.makedirs(self.workdirname)
        _os.chdir(self.workdirname)
        if self.overwrite and _os.path.exists(self.dirname):
            _shutil.rmtree("%s/%s" % (self.workdirname, self.dirname))
        if self.overwrite and self.copydirname and _os.path.exists(self.copydirname):
            _shutil.copytree("%s/%s" % (self.workdirname, self.copydirname), "%s/%s" % (self.workdirname, self.dirname))
        elif not _os.path.exists(self.dirname):
            _os.makedirs(self.dirname)
        _os.chdir(self.dirname)
        self.path = "%s/%s" % (self.workdirname, self.dirname)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        _os.chdir(self.workdirname)
        # removes temporary directory at the end of execution
        if self.temp:
            def delete():
                try:
                    _shutil.rmtree(self.path)
                except:
                    pass

            if self.purge_immediately:
                delete()
            else:
                _atexit.register(delete)
        _os.chdir(self.initialdirnames.pop())


def checkFileExists(files):
    """A simple wrapper around os.path.exists which throws an error if False."""
    if isinstance(files, str):
        single_file = True
        files = [files]
    else:
        single_file = False

    if all(_os.path.exists(x) for x in files):
        files = [_os.path.abspath(x) for x in files]
    else:
        raise ValueError("File(s) {} do(es) not exist. Please provide valid filename(s).".format(files))

    return files[0] if single_file else files

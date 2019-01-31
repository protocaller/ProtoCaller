import atexit as _atexit
import os as _os
import shutil as _shutil


class Dir:
    def __init__(self, dirname, copydirname=None, overwrite=False, temp=False, purge_immediately=True):
        if _os.path.isabs(dirname):
            self.workdirname = _os.path.dirname(dirname)
        else:
            self.workdirname = _os.getcwd()
        self.dirname = _os.path.basename(dirname)
        self.copydirname = copydirname
        self.overwrite = overwrite
        self.temp = temp
        self.purge_immediately = purge_immediately

    def __enter__(self):
        self.initialdirname = _os.getcwd()
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
        _os.chdir(self.initialdirname)


def checkFileExists(file):
    if not _os.path.exists(file):
        raise ValueError("File %s does not exist. Please provide a valid filename." % file)
    return _os.path.abspath(file)

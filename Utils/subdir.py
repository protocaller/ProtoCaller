#makes the generation of subdirectories a bit easier, more readable and less prone to errors
import os as _os
import shutil as _shutil

class Subdir():
    def __init__(self, dirname, copydirname=None, overwrite=False, temp=False):
        self.dirname = dirname
        self.copydirname = copydirname
        self.overwrite = overwrite
        self.temp = temp
    def __enter__(self):
        self.workdirname = _os.getcwd()
        if(self.overwrite and _os.path.exists(self.dirname)):
            _shutil.rmtree("%s/%s" % (self.workdirname, self.dirname))
        if(self.overwrite and self.copydirname and _os.path.exists(self.copydirname)):
            _shutil.copytree("%s/%s" % (self.workdirname, self.copydirname), "%s/%s" % (self.workdirname, self.dirname))
        elif(not _os.path.exists(self.dirname)):
            _os.makedirs(self.dirname)
        _os.chdir(self.dirname)
        self.path = "%s/%s" % (self.workdirname, self.dirname)
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        _os.chdir(self.workdirname)
        if(self.temp):
            _shutil.rmtree(self.path)
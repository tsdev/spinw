__author__ = "github.com/wardsimon"
__version__ = "0.1.0"

import os
import json
import libpymcr

_VERSION_DIR = os.path.join(os.path.dirname(__file__), 'ctfs')
_VERSIONS = []
for file in os.scandir(_VERSION_DIR):
    if file.is_file() and file.name.endswith('.ctf'):
        _VERSIONS.append({'file':    os.path.join(_VERSION_DIR, file.name),
                          'version': 'R' + file.name.split('.')[0].split('SpinW_')[1]
                          })


class Matlab(libpymcr.Matlab):
    def __init__(self, *args, **kwargs):
        initialized = False
        for version in _VERSIONS:
            try:
                print(f"Trying MATLAB version: {version['version']} ({version['file']}))")
                super().__init__(version['file'], *args, **kwargs)
                initialized = True
            except RuntimeError:
                continue
        if not initialized:
            raise RuntimeError(
                f"No MATLAB version found. Please use: [{', '.join([version['version'] for version in _VERSIONS])}]\n "
                f"If installed, please specify the root directory of the MATLAB installation.")

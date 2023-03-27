__author__ = "github.com/wardsimon"
__version__ = "0.1.0"

import os
import json
import libpymcr

_VERSION_INFO = os.path.join(os.path.dirname(__file__), 'ctfs', 'versions')
_VERSIONS = None
with open(_VERSION_INFO, 'r') as f:
    _VERSIONS = json.load(f)


class Matlab(libpymcr.Matlab):
    def __init__(self, *args, **kwargs):
        for version in _VERSIONS:
            initialized = False
            try:
                super().__init__(os.path.join(os.path.dirname(__file__), 'ctfs', version['file']), *args, **kwargs)
                initialized = True
            except RuntimeError:
                continue
        if not initialized:
            raise RuntimeError(f"No MATLAB version found. Please use: [{','.join([version['version'] for version in _VERSIONS])}]\n If installed, please specify the root directory of the MATLAB installation.")

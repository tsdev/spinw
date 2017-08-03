from transplant import Matlab
from transplant import MatlabStruct
import numpy as np

m = Matlab('/Users/sandortoth/spinw_git/dev/standalone/MacOS/pyspinw.app/pyspinw.sh')

tri = m.sw_model('triAF',1.)

spec = tri.powspec(np.linspace(0.,4.,51),'Evect',np.linspace(0.,6.,501))

spec = MatlabStruct(spec)
spec['param'] = MatlabStruct(spec['param'])

m.figure()
m.sw_plotspec(spec)
m.waitforgui()
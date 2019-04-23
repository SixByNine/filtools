from __future__ import absolute_import


from .sigproc import FilterbankIO
from ._data_scaling import *

import filtools.inject_pulsar

#FilterbankIO.__module__ = "filtools"
FluxScale.__module__ = "filtools"
SimpleIntegerise.__module__ = "filtools"

__all__ = ['FilterbankIO','FluxScale','SimpleIntegerise']

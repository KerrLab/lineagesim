# -*- coding: utf-8 -*-

VERSION = (0, 2, 2)
__version__ = ".".join(map(str, VERSION[0:3])) + "".join(VERSION[3:])
__license__ = "BSD"

from lineagesim.population import create_population, get_lineage_abundances, mutate, prune_frequency, reproduce
import lineagesim.output

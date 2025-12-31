Install environment with python 3.9 and Ocean libraries
See: https://support.dwavesys.com/hc/en-us/articles/360003718553-Install-the-Ocean-SDK-Locally

Arrange for the following: 
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import dimod
import csv
import os
from dwave.system import EmbeddingComposite
from dwave.system.samplers import DWaveSampler
from dwave.system import LeapHybridSampler
from dwave.cloud.exceptions import *
import time
import sys

#!/usr/bin/env python
import os
import sys

os.system( 'g++ example_variacional.cpp rayleighritz/linearrayleighritz/linearrayleighritz.cpp rayleighritz/csrayleighritz/csrayleighritz.cpp rayleighritz/csrayleighritz/bspline.cpp')
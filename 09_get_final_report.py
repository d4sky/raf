import sys, os

from _utils import *

#got_pMatches.Align()
    
pMatches = Matches.from_json("matches.json")
#print(pMatches.As_aligned())

pMatches.Populate_beg()
pMatches.Write('aligned', "matches.txt")


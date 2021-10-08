#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
# Copyright or © or Copr. Anthony Mathelier and Alessandra Carbone (02-05-2010) 
# anthony.mathelier@gmail.com, alessandra.carbone@lip6.fr
# 
# This software is a computer program whose purpose is to provide a platform to
# check several questions around miRNAs and pre-miRNAs prediction and play
# between sensitivity and specificity by using parameters variability.
# 
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
# 
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author, the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
# 
# In this respect, the user's attention is drawn to the risks associated
# with loading, using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that  also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
# 
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
################################################################################


###############################################################################
# script permettant de filter les précurseurs chevauchants des gènes définis
# dans le fichier GenBank passé en paramètre
###############################################################################

import sys, string, re
from Bio import GenBank

def readGbk(cur):
  featList = []
  for feat in cur.features:
    if (feat.type == "CDS" or feat.type == "exon" or feat.type == "ncRNA" or
        feat.type == "rRNA" or feat.type == "scRNA" or feat.type == "snRNA" or
        feat.type == "snoRNA" or feat.type == "tRNA" or feat.type == "21U-RNA"):
      if feat.type == "ncRNA":
        classe = feat.qualifiers['ncRNA_class'][0]
        if (classe <> "rRNA" and classe <> "scRNA" and classe <> "snRNA" and
            classe <> "snoRNA" and classe <> "tRNA" and classe <> "21U-RNA"):
          continue
      if feat.sub_features:
        for subfeat in feat.sub_features:
          featList.append((subfeat.location.nofuzzy_start,
            subfeat.location.nofuzzy_end))
      else:
        featList.append((feat.location.nofuzzy_start, feat.location.nofuzzy_end))
  return featList

def find(begin, end, feat):
  for fstart, fend in feat:
    if fstart > end:
      return 0
    if ((begin <= fend and begin > fstart) or
        (end <= fend and end > fstart) or
        (begin <= fstart and end >= fend)):
      return 1
  return 0

################################################################################
# Lecture des précurseurs et test pour savoir si il chevauche un gène ou non
################################################################################
def filter(input, feat):
  query = ""
  prec = ""
  seq = ""
  struct = ""
  for line in input:
    grp = re.match('>.*_(\d+):(\d+)_(\d+)_(\d+)$', line)
    if grp:
      description = line
      begin = eval(grp.group(1)) - eval(grp.group(3))
      end = eval(grp.group(1)) + eval(grp.group(2)) + eval(grp.group(4)) - 1
    elif re.search('[AGCTUacgtu]+', line):
      seq = line
      if not find(begin, end, feat):
        print "%s"%(description),
        print seq,
  input.close()

################################################################################
#           MAIN
################################################################################
if __name__ == "__main__": 
  if (len (sys.argv) <> 3):
    sys.exit ("\nUsage : " + sys.argv[0]+ " <predictions> <GenBank>\n")
  
  input = open(sys.argv[1], 'r')
  gbk = open(sys.argv[2], 'r')
  
  feature_parser = GenBank.FeatureParser()
  
  gb_iterator = GenBank.Iterator(gbk, feature_parser)
  gb_cur = gb_iterator.next()

  feat = readGbk(gb_cur)
  feat.sort()
  filter(input, feat)
  
  gbk.close()

#!/usr/bin/python
#-*- coding:utf-8 *-*

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



################################################################################
# script permettant d'éliminer les doublons
################################################################################

import re
from string import *
import sys, getopt

# Calcul du mfei étant données une séquence et sa mfe associée
def mfei (seq, mfe):
  length = len (seq)
  gc = seq.lower().count ('g')
  gc += seq.lower().count ('c')
  gc_content = float(gc) * 100 / float(length)
  if gc_content == 0.0:
    return float (-sys.maxint - 1)
  amfe = eval (mfe) * 100.0 / float(length)
  return amfe / gc_content

# classe représentant les pre-miRNAs
class PreMirna:
  def __init__ (self, name="", begin=0, end=0, chr="", seq="", struct="",
      mfe="", mirbegin=0, mirend=0, strand=0):
    self.name = name
    self.begin = begin
    self.end = end
    self.chr = chr
    self.seq = seq.upper()
    self.struct = struct
    self.mfe = mfe
    self.mfei = mfei (seq, mfe)
    self.mirbegin = mirbegin
    self.mirend = mirend
    self.strand = strand

  def strictly_included (self, preMi):
    return ((self.begin >= preMi.begin and self.begin <= preMi.end and
      self.end >= preMi.begin and self.end <= preMi.end) or
      (preMi.begin >= self.begin and preMi.begin <= self.end and
        preMi.end >= self.begin and preMi.end <= self.end))

  def overlap (self, preMi):
    return ((self.strictly_included (preMi)) or 
        (self.begin >= preMi.begin and self.begin <= preMi.end) or 
        (self.end >= preMi.begin and self.end <= preMi.end))

  def overlapMir(self, other):
    return ((self.mirbegin >= other.mirbegin and self.mirbegin <= other.mirend)
        or (self.mirend >= other.mirbegin and self.mirend <= other.mirend) or
        (self.mirbegin <= other.mirbegin and self.mirend >= other.mirend))

  def __str__ (self):
    return self.name+"\n"+self.seq+"\n"#+self.struct+" "+self.mfe+"\n"

  def miRNAfromPre(self):
    if self.strand == 1:
      beginFrom = self.mirbegin - self.begin + 1
      endFrom = self.mirend - self.begin + 1
    elif self.strand == -1:
      beginFrom = self.end - self.mirend + 1
      endFrom = self.end - self.mirbegin + 1
    else:
      sys.exit('pre-miRNA:\n%s\nHAS NO STRAND\n'%(self))
    return beginFrom, endFrom

  def miRNAmatching(self):
    beginFrom, endFrom = self.miRNAfromPre()
    miRNAstruct = self.struct[beginFrom-1:endFrom]
    nbPoint = float(miRNAstruct.count('.'))
    return (len(miRNAstruct) - nbPoint) / len(miRNAstruct)

  def structMatches(self):
    matchesOpen = {}
    matchesClose = {}
    index = 0
    stack = []
    for char in self.struct:
      index += 1
      if char == '(':
        stack.append(index)
      elif char == ')':
        i = stack.pop()
        matchesOpen[i] = index
        matchesClose[index] = i
      elif char == '.':
        continue
      else:
        sys.exit("\nError parsing structure:\n%s\nnear character %c"%(self.struct, char))
    if stack:
      sys.exit("\nError parsing structure:\n%s\nStack not empty\n"%(struct))
    return matchesOpen, matchesClose

  def miRNAstarMatch(self):
    matchesOpen, matchesClose = self.structMatches()
    beginFrom, endFrom = self.miRNAfromPre()
    miRNAstruct = self.struct[beginFrom-1:endFrom]
    if re.search('\(', miRNAstruct) and re.search('\)', miRNAstruct):
      return 1.
    deb = 0
    while miRNAstruct[deb] == '.':
      deb += 1
    deb += beginFrom
    fin = len(miRNAstruct) - 1
    while miRNAstruct[fin] == '.':
      fin -= 1
    fin += beginFrom
    if re.search('\(', miRNAstruct):
      deb = matchesOpen[deb]
      fin = matchesOpen[fin]
      return (abs(deb - fin) + 1.) / len(miRNAstruct)
    elif re.search('\)', miRNAstruct):
      deb = matchesClose[deb]
      fin = matchesClose[fin]
      return (abs(fin - deb) + 1.) / len(miRNAstruct)
    else:
      sys.exit('\nCannot find miRNA* in preMirna %s\n'%(self))

# Fonction pour la comparaison des mfei de deux pre-miRNAs
def cmp_mfei (p1, p2):
  if p1.mfei > p2.mfei:
    return 1
  elif p1.mfei == p2.mfei:
    if p1.miRNAlength() < p2.miRNAlength():
      return 1
    elif p1.miRNAlength() > p2.miRNAlength():
      return -1
    else:
      return 0
  else:
    return -1

def cmp_mfei_match(p1, p2):
#  print "\ncmp_mfei_match:\n%s\n%s\n"%(p1,p2)
  if p1.mfei > p2.mfei:
    return 1
  elif p1.mfei == p2.mfei:
    if p1.miRNAmatching() < p2.miRNAmatching():
      return 1
    elif p1.miRNAmatching() == p2.miRNAmatching():
      if abs(1-p1.miRNAstarMatch()) > abs(1-p2.miRNAstarMatch()):
        return 1
      elif abs(1-p1.miRNAstarMatch()) == abs(1-p2.miRNAstarMatch()):
        if p1.mirend - p1.mirbegin < p2.mirend - p2.mirbegin:
          return 1
        elif p1.mirend - p1.mirbegin > p2.mirend - p2.mirbegin:
          return -1
        else:
          if p1.mirbegin > p2.mirbegin:
            return 1
          elif p1.mirbegin < p2.mirbegin:
            return -1
          else:
            return 0
      else:
        return -1
    else:
      return -1
  else:
    return -1

# Fonction auxiliaire à remove_overlaps pour les miRNA avec effet de bord
def auxMirna (overlaps, non_overlap):
  if len (overlaps) == 0:
    return

  overlaps.sort(cmp=cmp_mfei_match)

  preMi_in = [overlaps[0]]
  for preMi_overlaps in overlaps[1:]:
    is_overlapping = 0
    for preMin in preMi_in:
      if preMin.overlapMir(preMi_overlaps):
          is_overlapping = 1

################################################################################
# Modifications pour conserver les miRNAs chevauchants qui ont les mêmes
# caractéristiques mais dont les pre-miRNAs sont différents avec le même MFEI
################################################################################
          if ( not cmp_mfei_match(preMin, preMi_overlaps) 
              and preMin.begin <> preMi_overlaps.begin
              and preMin.end <> preMi_overlaps.end):
            sys.stderr.write("\nequality between:\n%s\nand\n%s\n"%(preMin,
              preMi_overlaps))

          break
    if not is_overlapping:
      preMi_in.append (preMi_overlaps)
  non_overlap.extend (preMi_in)

# Retire les pre-miRNAs chevauchants par minimisation du MFEI
# La liste contient les pre-miRNAs par ordre croissant de leurs positions
def remove_overlaps (list):
  non_overlap = []
  overlaps = []
  for preMi_list in list:
    bool = 1
    for preMi_overlaps in overlaps:
      if preMi_list.overlap (preMi_overlaps):
        overlaps.append (preMi_list)
        bool = 0
        break
    if bool:
      auxMirna(overlaps, non_overlap)
      overlaps = [preMi_list]
  auxMirna(overlaps, non_overlap)
  return non_overlap

################################################################################
#                               MAIN
################################################################################
if __name__ == "__main__":
  usage = "\n%s -f <input file> -o <output file>\n"%(sys.argv[0])
  try:
    opts, args = getopt.getopt(sys.argv[1:], "f:o:")
  except getopt.GetoptError, err:
    sys.stdout.write(str(err))
    sys.exit(usage)
  if len(opts) <> 2:
    sys.exit(usage)
  for o, a in opts:
    if o == "-f":
      input = a
    elif o == "-o":
      output = a
  
  infile = open (input, "r")
  outfile = open (output, "w")
  
  dic = {}
  line = infile.readline()
  while (line <> ''):
    l = re.search (">.*_(.*)_(\d+):(\d+)_(\d+)_(\d+)", line)
    if (l <> None):
      name = rstrip (line)
      strand = 0
      begin = eval(l.group(2)) - eval(l.group(4))
      end = eval(l.group(2)) + eval(l.group(3)) + eval(l.group(5)) - 1
      strand = 1
      mirbegin = eval(l.group(2))
      mirend = eval(l.group(2)) + eval(l.group(3)) - 1
      chr = l.group(1)
      seq = rstrip (infile.readline())
      struct = rstrip (infile.readline())
      partition = re.search ("(.*) (\(.*\))", struct)
      struct = partition.group (1)
      mfe = partition.group (2)
      premirna = PreMirna (name, begin, end, chr, seq, struct, mfe, mirbegin,
          mirend, strand)
      if dic.has_key(premirna.chr):
        dic[premirna.chr].append (premirna)
      else:
        dic[premirna.chr] = [premirna]
    else:
      sys.exit ("Error in reading the file, not in right format\nError with line: %s" %(line))
    line = infile.readline()
  infile.close ()

  for geno in dic.keys ():
    #On trie la liste en fonction des positions des pre-miRNAs par ordre
    #croissant
    dic[geno].sort (cmp=lambda p1,p2: p1.begin - p2.begin)
    dic[geno] = remove_overlaps (dic[geno])

  for geno in dic.keys ():
    dic[geno].sort (cmp=lambda p1,p2: p1.begin - p2.begin)
    for premirna in dic[geno]:
      outfile.write (premirna.__str__())
  outfile.close ()

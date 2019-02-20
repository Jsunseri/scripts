#!/usr/bin/env python
import argparse,os
import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
import pymol
pymol.finish_launching()

parser = argparse.ArgumentParser(description="Extract pocket residues from a"
        " protein structure using a reference molecule to identify the binding"
        " site")
parser.add_argument("-r", "--receptor", help="Protein structure to extract"
        " residues from", required=True)
parser.add_argument("-l", "--ligand", help="Ligand used to identify the"
        " pocket", required=True)
parser.add_argument("-d", "--distance", default=8, help="Max distance from"
        " ligand to be considered part of the pocket")
parser.add_argument("-o", "--outname", default='', help="Prefix for output"
        " filename, by default it is [receptorname]_pocket")
args = parser.parse_args()

assert os.path.isfile(args.receptor), "File %s not found" %args.receptor
assert os.path.isfile(args.ligand), "File %s not found" %args.ligand
assert args.distance > 0, "Distance has to be positive"

pymol.cmd.load(args.receptor, "rec")
pymol.cmd.load(args.ligand, "lig")
pymol.cmd.select("sele", "byres rec within %s of lig" %args.distance)

if not args.outname:
    args.outname = os.path.splitext(args.receptor)[0] + "_pocket"

pymol.cmd.save(args.outname + ".pdb", "sele")

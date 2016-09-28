#!/usr/bin/python

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import is_aa
from Bio import pairwise2
import clustering
import os.path
import parsepdb
import numpy as np
import subprocess as sp
import argparse

def readinSeqs(parser, fname_array, seq_array, name_array):
    '''
    Read in aa string using PDBParser
    '''
    for fname in fname_array:
        id = os.path.basename(fname).split('.')[0]
        structure = parser.get_structure(id,fname)
        seq = clustering.getResidueString(structure)
        seq_array.append(seq)
    return

def calcSim(test_sequences, test_names, exclusion_sequences, exclusion_names,
        method, alignment):
    '''
    Returns nonsymmetric n X m matrix of dissimilarities for proteins in the
    proposed test set to those in the exclusion set. Choose between global and
    local alignment methods and whether to output the final alignment.
    '''
    n = len(test_sequences)
    m = len(exclusion_sequences)
    distanceMatrix = np.zeros((n,m))
    for r in range(n):
    	for c in range(m):
            if method == "global":
                if alignment:
    	            result = pairwise2.align.globalxx(test_sequences[r],
                            exclusion_sequences[c], one_alignment_only=True)
                    score = result.score
                    assert len(seqA) == len(seqB), "sequences not equal"
                    with open(test_names[r] + '_to_' + exclusion_names[c] + '_global.txt') as f:
                        for i in range(len(seqA)):
                            f.write('%s\t%s\n' (result.seqA[i], result.seqB[i]))
                else:
    	            score = pairwise2.align.globalxx(test_sequences[r],
                            exclusion_sequences[c], score_only=True)
            else:
                if alignment:
    	            result = pairwise2.align.localxx(test_sequences[r],
                            exclusion_sequences[c], one_alignment_only=True)
                    score = result.score
                    assert len(seqA) == len(seqB), "sequences not equal"
                    with open(test_names[r] + '_to_' + exclusion_names[c] + '_local.txt') as f:
                        for i in range(len(seqA)):
                            f.write('%s\t%s\n' (result.seqA[i], result.seqB[i]))
                else:
    	            score = pairwise2.align.localxx(test_sequences[r],
                            exclusion_sequences[c], score_only=True)
    	    length = max(len(test_sequences[r]), len(exclusion_sequences[c]))
    	    distance = (length-score)/length
    	    distanceMatrix[n,m] = distance
    return distanceMatrix

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generates a set of test\
    compounds by removing any whose similarity to an exclusion list exceeds a\
    chosen threshold')
    parser.add_argument('-t', '--test', nargs='+', required=True,
            help='Proposed test set which will be screened. Given as a list of\
            file paths')
    parser.add_argument('-e', '--exclude', nargs='+', required=True,
            help='Exclusion set; for supervised learning this is at least the\
            training set')
    parser.add_argument('--threshold', type=float, default=.2, help='Dissimilarity\
            threshold for exclusion')
    parser.add_argument('--no-sequence', dest='sequence', action='store_false',
            help='Do not compute sequence similarity')
    parser.add_argument('--no-structure', dest='structure',
            action='store_false', help='Do not compute structural similarity')
    parser.add_argument('--no-alignment', dest='alignment',
            action='store_false', help='Do not store alignments')
    parser.set_defaults(sequence=True, structure=True, alignment=True)
    args = parser.parse_args()

    if args.sequence:
        p = PDBParser(PERMISSIVE=1, QUIET=1)
        test_sequences = []
        test_names = []
        exclusion_sequences = []
        exclusion_names = []

        # Get strings for proposed test set and exclusion set
        readinSeqs(p, args.test, test_sequences, test_names)
        readinSeqs(p, args.exclude, exclusion_sequences, exclusion_names)

        # Generate global sequence similarity matrix
        globalDist = calcSim(test_sequences, exclusion_sequences, test_names, 
                exclusion_names, "global", args.alignment)
    
        # Generate local sequence similarity matrix
        localDist = calcSim(test_sequences, exclusion_sequences, test_names, 
                exclusion_names, "local", args.alignment)

    if args.structure:
        # Prep for ProBiS - if ligands, do structural alignment to binding
        # site, inferred from ligand position. Otherwise do a full
        # superposition. 
        for struct in args.test:
            mol_data,outfile = parsepdb.fix_unlabeled_chains(parsepdb.parsepdb(struct))
            rec_chains,lig_chains = parsepdb.get_chains(mol_data)
            resnums,resnames,natoms = parsepdb.get_nonprot(mol_data)
            possible_ligand_indices = [i for i,val in enumerate(natoms) if val > 1]
            possible_ligands = zip([resnames[i] for i in
                possible_ligand_indices], [resnums[i] for i in
                possible_ligand_indices])
            
            if outfile:
                struct = os.path.dirname(struct) + '/' + outfile

            if not possible_ligands:
                # No ligands, so we can't identify the binding site. 
                print('Doing full structural superposition. Provide complex if\
                        you want to base alignment on binding site.\n')
                for exc in args.exclude:
                    exc_data = parsepdb.parsepdb(exc)
                    exc_chains,_ = parsepdb.get_chains(exc_data)
                    cmd = ['probis','-compare','-super','-f1',struct,'-c1', 
                            ''.join(rec_chains),'-f2',exc,'-c2',''.join(exc_chains),
                            '-out',os.path.dirname(struct)]
                    p = sp.Popen(cmd, stdout = sp.PIPE, stderr = sp.PIPE)
                    out,err = p.communicate()
            else:
                # If there's only one plausible ligand after excluding ions, go
                # with that. Otherwise, prompt.
                if len(possible_ligands) > 1:
                    while lig is None:
                        print "Multiple ligands found; please choose one"
                        print ' '.join([':'.join([str(i[0]),i[1]]) for i in
                            possible_ligands])
                        choice = int(raw_input())
                        if choice > len(possible_ligands):
                            print "Invalid selection"
                            lig = None
                        lig = possible_ligands[choice-1]
                else:
                    lig = possible_ligands[0]
                for exc in args.exclude:
                    exc_data = parsepdb.parsepdb(exc)
                    exc_chains,_ = parsepdb.get_chains(exc_data)
                    cmd = ['probis','-compare','-super','-dist','3.0','-f1',struct,'-c1',
                            ''.join(rec_chains),'-bsite1',
                            '.'.join([lig[0],lig[1],lig_chains[lig[1]]]),
                            '-f2',exc,'-c2',''.join(exc_chains),'-o',
                            os.path.dirname(struct)]
                    p = sp.Popen(cmd, stdout = sp.PIPE, stderr = sp.PIPE)
                    out,err = p.communicate()


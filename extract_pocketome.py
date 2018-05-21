#!/usr/bin/python
import prody
import os, tarfile
from xml.dom import minidom
import argparse

def parse_pocketome(xml):
    """
    Takes XML from Pocketome with information about related proteins, returns
    dictionary with those proteins' accession IDs, main chains according to
    Pocketome, binding site residues according to Pocketome, and ligands.
    """
    struct_list = xml.getElementsByTagName("object")
    #TODO: ehhhh? need to decide whether the Pocketome 'site' or 'pocket' is
    #what we want, and whether we actually need the specific contact residues
    #for each protein
    binding_site =\
    xml.getElementsByTagName('residueNumbersPocket')[0].firstChild.nodeValue.split(',')
    binding_site = 'resnum ' + ' '.join(binding_site)
    output = {}
    output['pdb'] = []
    output['lig'] = []
    output['chain'] = []
    output['site'] = binding_site

    for i,entry in enumerate(struct_list):
        # we don't mess around with covalently bound ligands
        if entry.getElementsByTagName('covlig'): continue
        ligs = entry.getElementsByTagName('ligand')
        # might want to change this
        if not ligs: continue
        id = entry.attributes["id"].value
        output['pdb'].append(id[1:-1])
        output['chain'].append(str(id[-1].upper()))
        ligs = ligs[0].firstChild.nodeValue.split(',')
        output['lig'].append(ligs)
    return output

def protein_align(pdb, ref, **kwargs):
    """
    Align proteins with prody
    """
    mobile = prody.parsePDB(pdb)
    static = prody.parsePDB(ref)
    if kwargs['pdb_chain']:
        mobile = mobile.select('chain %s ' %
        (kwargs['pdb_chain']))
    if kwargs['ref_chain']:
        static = static.select('chain %s ' % (kwargs['ref_chain']))
    try:
        rmatch, tmatch, seqid, overlap = prody.matchChains(static, mobile,
                pwalign=True, overlap=70)[0]
        moved, transformation = prody.superpose(tmatch, rmatch)
        return mobile
    except:
        print "Failed to match %s, continuing...\n" %pdb
        return prody.parsePDB(pdb)

def write_output(structure, name, outdir, chain, ligs):
    outname = outdir + '/' + name
    prody.writePDB(outname + '_rec.pdb', structure.select('protein and chain \
%s or ion and chain %s' %(chain, chain)))
    for lig in ligs:
        if structure.select('resname %s and chain %s' %(lig.upper(), chain)):
            prody.writePDB('%s_%s_lig.pdb' %(outname, lig),
            structure.select('resname %s and chain %s' %(lig.upper(), chain)))
    return

def extract_output(xml, reference, no_align, refchain):
    outdir = xml.getElementsByTagName("entry")[0].attributes["name"].value
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    #TODO: sanity check input?
    output = parse_pocketome(xml)
    if not reference:
        reference = output['pdb'][0]
        refchain = output['chain'][0]
    for i,pdb in enumerate(output['pdb']):
        if not no_align:
            structure = protein_align(pdb, reference,
            pdb_chain=output['chain'][i],
            ref_chain=refchain)
        else:
            structure = prody.parsePDB(pdb)
        write_output(structure, pdb, outdir, output['chain'][i],
                output['lig'][i])
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generates set of related \
    PDBs for cross-docking from Pocketome datafile. Optionally aligns them, \
    either to the first usable PDB specified by Pocketome or to a user-\
    provided reference.")

    parser.add_argument('-i', '--input', required=True, nargs='+', help="Input \
    data from Pocketome. This can be one or more XML files, or it can be a tar \
    archive containing those files.")

    parser.add_argument('-na', '--no_align', required=False, default=False,
            action='store_true',
    help="Indicate if you don't want alignment to be performed on the \
    resulting structures.")

    parser.add_argument('-r', '--reference', required=False, help="Optionally \
    specify a reference structure for alignment. If you don't specify anything \
    and alignment is on, the first structure returned from Pocketome is used.")

    parser.add_argument('-rc', '--refchain', required=False, help="Optionally \
    specify a reference chain to be used for alignment with the reference \
    structure.")

    args = parser.parse_args()

    for file in args.input:
        input_ext = os.path.splitext(file)[-1]
        if input_ext == '.tbz':
            dataset = tarfile.open(file)
            for datum in dataset:
                xml = minidom.parse(dataset.extractfile(datum))
                extract_output(xml, args.reference, args.no_align,
                        args.refchain)
        else:
            assert input_ext == '.xml', 'Only XML files or a tar archive of XML \
files are acceptable inputs. Got %s.\n' %(input_ext)
            xml = minidom.parse(file)
            extract_output(xml, args.reference, args.no_align, args.refchain)

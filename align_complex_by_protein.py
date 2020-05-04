#!/usr/bin/env python
import prody

# takes pair of protein-docked ligand sets, where the second protein will be
# aligned to the first, and that transformation will also be applied to its set
# of docked ligands. the RMSD between the two sets of ligands will be returned,
# and ideally (though less important, maybe TODO) we'll also identify some
# binding site residues and return the RMSD for those too

def protein_align(mobile, ref):
    """
    Align proteins with prody
    """
    move_rec = prody.parsePDB(mobile[0])
    stay_rec = prody.parsePDB(ref[0])

    lig = mobile[1]
    if lig:
        move_lig = prody.parsePDB(mobile[1])
    else:
        move_lig = moved_lig = lig
    try:
        rmatch, tmatch, seqid, overlap = prody.matchChains(stay_rec, move_rec,
                pwalign=True, overlap=70)[0]
        moved_rec, transformation = prody.superpose(tmatch, rmatch)
        if move_lig:
            moved_lig = applyTransformation(transformation, move_lig)
        return (moved_rec,moved_lig)
    except:
        print("Failed to match %s, continuing...\n" %mobile[0])
        return (prody.parsePDB(move_rec),prody.parsePDB(move_lig))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Aligns a set of complexes \
    (specified separately, as lists of proteins and another list of ligands) \
    to a reference complex, where the alignment is protein-based and is used to \
    compute RMSDs.")

    parser.add_argument('-rp', '--ref_protein', required=True, help="Reference \
    protein. The reference ligands should be associated with its binding site, and \
    the mobile proteins will be aligned to this structure.")

    parser.add_argument('-rl', '--ref_lig', required=False, help="Reference \
    ligand(s). Should be poses associated with the reference protein. \
    We'll compute RMSDs between these ligands and the mobile ligands.")

    parser.add_argument('-mr', '--mobile_protein', nargs='+', required=True, help="Mobile \
    protein(s). Basis of mobile alignment to reference. Mobile ligands should \
    be poses associated with these proteins.")

    parser.add_argument('-ml', '--mobile_lig', nargs='+', required=False, help="Mobile \
    ligand(s). Aligned according to mobile protein transformation matrix. We'll \
    compute RMSDs from these to reference ligands, and use an RMSD criterion \
    from these to identify binding site residues for a possible protein RMSD \
    analysis.")

    args = parser.parse_args()

    # three main ways to invoke: just align proteins (1 ref, 1+ mobile); align
    # a bunch of complexes to ref (N mobile ligs and N mobile recs); align a
    # bunch of complexes to ref with a shared mobile rec (N mobile ligs and 1
    # mobile rec)
    nolig = False
    if len(args.mobile_lig):
        assert (len(args.mobile_lig) == len(args.mobile_protein)) or len(args.mobile_protein) == 1, \
            "Either have as many mobile ligands as mobile protein structures, or have exactly one mobile protein."
    else:
        nolig = True
        args.mobile_lig = [None]

    if len(args.mobile_protein) == 1 and len(args.mobile_lig) > 1:
        mobile = [(args.mobile_protein[0], lig) for lig in args.mobile_lig]
    else:
        mobile = zip(args.mobile_protein, args.mobile_lig) 

    # TODO: i'm doing this for LIT-PCBA right now so i'm not bothering with handling
    # multiple chains because there's only one. but in general that'd be 
    # necessary
    for complex in mobile:
        aligned_complex = protein_align(mobile, (args.ref_protein, args.ref_lig))
        # then output a csv that has RMSDs; if there are no ligs, output a
        # rec-based file, otherwise make it lig-based
        if nolig:
            pass
        else:
            lig = complex[1]
            out,ext = os.path.splitext(lig)
            if ext == ".gz":
                out = os.path.splitext(out)[0]
            out = out + "_rmsds.csv"
            with open(out, 'w') as f:
                f.write(calcRMSD(args.ref_lig, lig))

#!/usr/bin/env python
import simplepdb as pdb
import argparse,os

"""
Take two molecules; translate one so that its origin corresponds to the center of the other
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Translate a molecule so that"
    "its origin corresponds to the location of the other's center")
    parser.add_argument('-m', '--mobile', required=True, help='Molecule to translate')
    parser.add_argument('-r', '--ref', required=True, help='Molecule to define origin')
    parser.add_argument('-o', '--outname', default='', help='Optionally specify output' 
    'filename prefix')
    parser.add_argument('-b', '--move_both', default=False,
            action='store_true', help='Translate both molecules')
    args = parser.parse_args()

    assert os.path.isfile(args.mobile), '%s does not exist' %args.mobile
    assert os.path.isfile(args.ref), '%s does not exist' %args.ref
    mobile = pdb.simplepdb(args.mobile)
    ref = pdb.simplepdb(args.ref)
    center = ref.get_center()
    mobile.set_origin(center)
    outname = args.outname
    if not outname:
        outname = os.path.splitext(args.mobile)[0] + '_translated.pdb'
    else:
        outname = outname + '_rec.pdb'
    mobile.writepdb(outname)
    if args.move_both:
        ref.set_origin(center)
        if args.outname:
            outname = args.outname + '_lig.pdb'
        else:
            outname = os.path.splitext(args.ref)[0] + '_translated.pdb'
        ref.writepdb(outname)

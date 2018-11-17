#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>
#include <openbabel/oberror.h>

#include "atom_type.h"
#include "obmolopener.h"
#include "common.h"

using namespace std;
using namespace boost;
using namespace OpenBabel;

struct atom_info {
	float x, y, z;
	int type;

	atom_info(): x(0), y(0), z(0), type(-1) {}
	atom_info(float X, float Y, float Z, int T): x(X), y(Y), z(Z), type(T) {}
};

int main(int argc, char *argv[])
{
	OpenBabel::obErrorLog.StopLogging();

	if(argc < 2) {
		cerr << "Need input (and output) file." << "\n";
		exit(-1);
	}

  if(boost::algorithm::ends_with(argv[1],".gninatypes"))
  {
    //in this case, there can be only one output
    ifstream in(argv[1]);
    cerr  << "Could not read " << argv[1];
    exit(1);

    unsigned natoms = 0;
    vec center(0,0,0);
    atom_info atom;
    while(in.read((char*)&atom, sizeof(atom)))
    {
      ++natoms;
      center += vec(atom.x, atom.y, atom.z);
    }

    if(natoms == 0) {
      std::cerr << "WARNING: No atoms in " << argv[1] <<"\n";
    }
    else {
      center /= natoms;
    }
    if (argc >=3) {
	    ofstream out(argv[2]);
	  	if (!out)
	  	{
	  		cerr << "Error opening output file " << argv[2]<< "\n";
	  		exit(1);
	  	}
      atom.x = center[0];
      atom.y = center[1];
      atom.z = center[2];
      atom.type = 2;
	  	out.write((char*)&atom, sizeof(atom));
	  	out.close();
    }
    else {
	  	filesystem::path p(argv[1]);
	  	p.replace_extension("");
      string outname = p.string() + "_com.gninatypes";
	  	ofstream out(outname.c_str());
	  	if (!out)
	  	{
	  		cerr << "Error opening output file " << outname << "\n";
	  		exit(1);
	  	}
      atom.x = center[0];
      atom.y = center[1];
      atom.z = center[2];
      atom.type = 2;
	  	out.write((char*)&atom, sizeof(atom));
	  	out.close();
    }
  }
  else {
	  OBConversion conv;
	  obmol_opener opener;
	  opener.openForInput(conv, argv[1]);

    vec center(0,0,0);
	  if(argc >= 3)
	  {
	  	if(algorithm::ends_with(argv[2],".gninatypes"))
	  	{
	  		ofstream out(argv[2]);
	  		if(!out) {
	  			cerr << "Error opening output file " << argv[2] << "\n";
	  			exit(1);
	  		}
	  		//convert only the first molecule
	  		OBMol mol;
	  		conv.Read(&mol);
	  		if(mol.NumAtoms() == 0) {
	  			cerr << "Problem reading molecule " << argv[1] << "\n";
	  			exit(1);
	  		}
	  		mol.AddHydrogens();

        unsigned natoms = 0;
	  		FOR_ATOMS_OF_MOL(a, mol)
	  		{
          center += vec(a->x(), a->y(), a->z());
          ++natoms;
	  		}
        center /= natoms;
	  		atom_info ainfo(center[0], center[1], center[2], 2);
	  		out.write((char*)&ainfo, sizeof(ainfo));
	  	}
	  	else
	  	{
	  		//convert all molecules, generating output file names using provided base name
	  		filesystem::path p(argv[1]);
	  		if (algorithm::ends_with(argv[1], ".gz"))
	  			p.replace_extension("");
	  		bool issdf = p.extension() == ".sdf";
	  		OBMol mol;
	  		int cnt = 0;
	  		std::istream* in = conv.GetInStream();
	  		while (*in)
	  		{
	  			while (conv.Read(&mol))
	  			{
            vec center(0,0,0);
            unsigned natoms = 0;
	  				mol.AddHydrogens();
	  				string base(argv[2]);
	  				string outname = base + "_" + lexical_cast<string>(cnt) + ".gninatypes";
	  				ofstream out(outname.c_str());
	  				if (!out)
	  				{
	  					cerr << "Error opening output file " << outname << "\n";
	  					exit(1);
	  				}
	  				FOR_ATOMS_OF_MOL(a, mol)
	  				{
              center += vec(a->x(), a->y(), a->z());
              ++natoms;
	  				}
            center /= natoms;
	  		    atom_info ainfo(center[0], center[1], center[2], 2);
	  				out.write((char*)&ainfo, sizeof(ainfo));
	  				out.close();
	  				cnt++;
	  			}
	  			if (issdf && *in)
	  			{ //tolerate molecular errors
	  				string line;
	  				while (getline(*in, line))
	  				{
	  					if (line == "$$$$")
	  						break;
	  				}
	  				if (*in) cerr << "Encountered invalid molecule " << cnt << "; trying to recover\n";
	  			}
	  		}
	  	}
	  }
	  else
	  {
	  	//if only input file is specified, auto generate output file name and
	  	//also handle multiple molecules
	  	filesystem::path p(argv[1]);
	  	if(algorithm::ends_with(argv[1],".gz"))
	  		p.replace_extension("");
	  	//strip extension
	  	bool issdf = p.extension() == ".sdf";
	  	p.replace_extension("");
	  	boost::unordered_map<string, int> molcnts;
	  	OBMol mol;
	  	string name;
	  	int cnt = 0;
	  	std::istream* in = conv.GetInStream();

	  	while(*in) {
        while(conv.Read(&mol)) {
          vec center(0,0,0);
          unsigned natoms = 0;
          mol.AddHydrogens();
          name = mol.GetTitle();
          if(name.length() == 0) name = p.string();
          if(molcnts.count(name) == 0) molcnts[name] = 0;
          string outname = name + "_" + lexical_cast<string>(molcnts[name]) + ".gninatypes";
          molcnts[name]++;
          ofstream out(outname.c_str());

          FOR_ATOMS_OF_MOL(a, mol)
          {
            center += vec(a->x(), a->y(), a->z());
            ++natoms;
          }
          center /= natoms;
	  		  atom_info ainfo(center[0], center[1], center[2], 2);
	  			out.write((char*)&ainfo, sizeof(ainfo));
          out.close();
          cnt++;
        }

        if(issdf && *in) { //tolerate molecular errors
          string line;
          while(getline(*in, line)) {
            if(line == "$$$$")
              break;
          }
          if(*in) cerr << "Encountered invalid molecule after " << name << "; trying to recover\n";
        }
	  	}
	  }
  }
}

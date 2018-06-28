#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <map>
#include <cstdlib>
#include <stdio.h>
#include <cstring>
#include <cctype>
#include <iomanip>
#include <algorithm>
#include "binafile.h"
#include "chck_jbnd.h"
#include "leng_stat.h"
#include "formfactor.h"
#include "pcf_calc.h"
#include "Sofk.h"
#include "ConfigFile.hpp"
#define MXSZ 500

using namespace std;

bool str_to_bool(string str){
	std::transform(str.begin(), str.end(),str.begin(), ::tolower);
	istringstream is(str);
	bool b;
	is >> std::boolalpha >> b;
	return b;
}

int main()
{


  string line;
  char base_[MXSZ];     //  INPUT: DATA DIRECTORY
  char cord_[MXSZ];     //  INPUT: SITE COORDINATES
  char fram_[MXSZ];     //  INPUT: SIMULATION TIME STEPS
  bool __binary__=false;

	ConfigFile cfg("config.cfg");

  // Check for model parameters
  int Nmono = cfg.getValueOfKey<int>("Nmono",50); // Number of monomers per chain. Default 50
	cout << "Number of monomers per chain: " << Nmono << endl;
	int poly = cfg.getValueOfKey<int>("nPoly",300); // Number of chains. Default 300
	cout << "Number of polymer chains: " << poly << endl;
	double boxs = cfg.getValueOfKey<double>("LBox",40.0); // Simulation box size. Default 40
	cout << "Simulations box length: " << boxs << endl;
	int atmtypes = cfg.getValueOfKey<double>("ATOM_TYPES",2); // Number of atom types. Default 2
	cout << "Number of atom types: " << atmtypes << endl;
  vector<int> charges;
	int tmp_charge;
	for(int i=0; i<atmtypes; i++)
	{
		tmp_charge = cfg.getValueOfKey<double>("CHARGE",0,i); // Charges on each atom type. Default 0
		charges.push_back(tmp_charge);
	}
	for(int i=0; i<atmtypes; i++)
	{
			cout << "Atom type " << i+1 << " Charge: " << charges[i] << endl;
	}
  // Get coordinate file and frame file
  string stval = cfg.getValueOfKey<string>("FILE_PATH",""); // path to coordinates. Default ""
	cout << stval << endl;
	strcpy(base_, stval.c_str());

	stval = cfg.getValueOfKey<string>("COORD_FILE","coords.dat"); // Name of coordinate file. Default "coords.dat"
	cout << "Reading coordinates in " << base_ << stval << endl;
	strcpy(cord_, stval.c_str());

	stval = cfg.getValueOfKey<string>("FRAME_FILE","frames.dat"); // Name of frame file. Default "frames.dat"
	cout << "Reading frames in " << base_ << stval << endl;
	strcpy(fram_, stval.c_str());

  // Reading what to calculate

  bool __fixbonds__=false;
	bool exists = cfg.keyExists("FIX_BONDS"); // Correct for pbc. Default is false
	if(exists)
	{
		string stval = cfg.getValueOfKey<string>("FIX_BONDS");
		__fixbonds__ = str_to_bool(stval);
		std::cout << "Coordinates will be corrected for pbc" << endl;
	}

	bool __computelength__=false;
	exists = cfg.keyExists("LENGTH_STATISTICS"); // Compute polymer length statistics. Default is false
	if(exists)
	{
		string stval = cfg.getValueOfKey<string>("LENGTH_STATISTICS");
		__computelength__ = str_to_bool(stval);
		std::cout << "Polymer length statistics will be computed" << endl;
	}

	bool __computewmmk__=false;
	exists = cfg.keyExists("FORMFACTOR"); // Compute monomer form factor. Default is false
	if(exists)
	{
		string stval = cfg.getValueOfKey<string>("FORMFACTOR");
		__computewmmk__ = str_to_bool(stval);
		std::cout << "Intramolecular form factor will be computed" << endl;
    __binary__=true;
	}

	bool __computewcck__=false;
	exists = cfg.keyExists("CHARGE_FORMFACTOR"); // Compute charge form factor. Default is false
	if(exists)
	{
		string stval = cfg.getValueOfKey<string>("CHARGE_FORMFACTOR");
		__computewcck__ = str_to_bool(stval);
		std::cout << "Charge form factor will be computed" << endl;
    __binary__=true;
	}

	bool __computehmmr__=false;
	exists = cfg.keyExists("PAIR_CORRELATION_FUNC"); // Compute pair correlation function. Default is false
	if(exists)
	{
		string stval = cfg.getValueOfKey<string>("PAIR_CORRELATION_FUNC");
	  __computehmmr__ = str_to_bool(stval);
		std::cout << "Pair correlation function will be computed " << endl;
    __binary__=true;
	}

	bool __computehcmr__=false;
	exists = cfg.keyExists("COM_PAIR_CORRELATION_FUNC"); // Compute center of mass pair correlations. Default is false
	if(exists)
	{
		string stval = cfg.getValueOfKey<string>("COM_PAIR_CORRELATION_FUNC");
		__computehcmr__ = str_to_bool(stval);
		std::cout << "center of mass pair correlation function will be computed" << endl;
    __binary__=true;
	}

	bool __computeSmmk__=false;
	exists = cfg.keyExists("STRUCTURE_FACTOR"); // Compute monomer structure factor. Default is false
	if(exists)
	{
		string stval = cfg.getValueOfKey<string>("STRUCTURE_FACTOR");
		__computeSmmk__ = str_to_bool(stval);
		std::cout << "Structure factor will be computed" << endl;
    __binary__=true;
	}

	bool __computeScck__=false;
	exists = cfg.keyExists("CHARGE_STRUCTURE_FACTOR"); // Compute charge structure factor. Default is false
	if(exists)
	{
		string stval = cfg.getValueOfKey<string>("CHARGE_STRUCTURE_FACTOR");
		__computeScck__ = str_to_bool(stval);
		std::cout << "Charge structure factor will be computed" << endl;
    __binary__=true;
	}
  // Done parsing inputs

  if(__fixbonds__)
  {
    cout << "fixing for pbc" << endl;
    checkbond bond;
    bond.fixbrokenbonds(Nmono,poly,boxs,base_,cord_,fram_);
    stval="jbnds.dat";
    strcpy(cord_, stval.c_str());
    cout << "New coordinates will be read from " << base_ << stval << endl;
  }
  if(__computelength__)
  {
    cout << "computing polymer length statistics" << endl;
    lengthstats length;
    length.computelengthstats(Nmono,poly,base_,cord_,fram_);
  }
  if(__binary__)
  {
    cout << "generating binary coordinates" << endl;
    binafile bina;
    bina.convertbinary(Nmono,poly,base_,cord_,fram_);
  }
  if(__computewmmk__)
  {
    cout << "computing monomer form factor" << endl;
    if(!__binary__)
    {
    cerr << "form factor must use binary coordinates" << endl;
    exit(1);
    }
    FORMFACTOR formfact;
    formfact.calculate_wmmk(Nmono,poly,base_,fram_);
  }
  if(__computewcck__)
  {
    cout << "computing charge-charge form factor" << endl;
    if(!__binary__)
    {
    cerr << "form factor must use binary coordinates" << endl;
    exit(1);
    }
    FORMFACTOR formfactcharge;
    formfactcharge.calculate_wmmk_charge(Nmono,poly,charges,base_,fram_);
  }
  if(__computehmmr__)
  {
    cout << "computing monomer pair correlation function" << endl;
    if(!__binary__)
    {
    cerr << "pair correlation must use binary coordinates" << endl;
    exit(1);
    }
    PAIRCORRELATION pcfhmmr;
    pcfhmmr.calculate_hmmr(Nmono,poly,boxs,base_,fram_);
  }
  if(__computehcmr__)
  {
    cout << "computing center of mass pair correlation function" << endl;
    if(!__binary__)
    {
    cerr << "pair correlation must use binary coordinates" << endl;
    exit(1);
    }
    PAIRCORRELATION pcfhcomr;
    pcfhcomr.calculate_hcomr(Nmono,poly,boxs,base_,fram_);
  }

  if(__computeSmmk__)
  {
    cout << "computing monomer total structure factor" << endl;
    if(!__binary__)
    {
    cerr << "structure factor must use binary coordinates" << endl;
    exit(1);
    }
    STRUCTUREFACTOR smonok;
    smonok.calculate_Smmk(Nmono,poly,boxs,base_,fram_);
  }
  if(__computeScck__)
  {
    cout << "computing charge charge structure factors" << endl;
    if(!__binary__)
    {
    cerr << "structure factor must use binary coordinates" << endl;
    exit(1);
    }
    STRUCTUREFACTOR sccok;
    sccok.calculate_Scck(Nmono,poly,boxs,charges,base_,fram_);
  }
  return 0;
}

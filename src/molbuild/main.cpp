/**main.cpp
 * Author: Nicholas Browning
 * **/


#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/detail/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <openbabel/atom.h>
#include <openbabel/base.h>
#include <openbabel/bond.h>
#include <openbabel/builder.h>
#include <openbabel/generic.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/residue.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "MoleculeBuilder.h"

using namespace boost;
using namespace boost::filesystem;
using namespace std;

string fragmentLibraryPath;
string inputFilePath;
string outFilePath;
string outFileExtension;

std::map<std::string, std::string> fragmentPathMap;
std::map<std::string, std::string> fragmentExtensionMap;

OBConversion conv;
OBMol molecule;
MoleculeBuilder build;

vector<int> psiphiResidueIndex;
vector<double> psiAngle;
vector<double> phiAngle;

vector<string> actions;
vector<int> molatom1Index;
vector<int> molatom2Index;
vector<int> fragatom1Index;
vector<int> fragatom2Index;
vector<bool> referenceBondTorsionFromFragment;
vector<bool> updateResidueInfo;
vector<int> bondO;
vector<string> fragmentSubstitution;
string loadMolPath;

vector<OBMol *> fragments;

MoleculeBuilder mb;

string molPath;
string molType;

string fragmentPath;
string fragmentType;

string finalPath;
string finalType;

int molatom1index = -1;
int molatom2index = -1;

int fragatom1index = -1;
int fragatom2index = -1;

void performSubstitutions(string fraglibpath, string xyzpath,
		string outputpath) {
//	cout << "Start" << endl;
	OBConversion conv;
	conv.SetInAndOutFormats("xyz", "xyz");

	string xyzs =
			"/Users/roethlisbergergroup/EPFL/GA-QMML/22k/modoptclasses/structs";
	string fraglib =
			"/Users/roethlisbergergroup/EPFL/GA-QMML/22k/modoptclasses/fraglib";

	for (boost::filesystem::recursive_directory_iterator end, dir(xyzs);
			dir != end; ++dir) {

		vector < OBMol > frags;
		for (boost::filesystem::recursive_directory_iterator end, dir(fraglib);
				dir != end; ++dir) {
			OBMol frag;

			path currentPath = dir->path();
			path localFilenamePath = currentPath.filename();
			string localFileame = localFilenamePath.string();

			vector < string > strs;
			boost::split(strs, localFileame, boost::is_any_of("."));
			string filename = strs[0];
			string ext = strs[1];

			if (!boost::iequals(ext, "xyz")) {
				continue;
			}
			cout << "Reading Frag: " << currentPath.string() << endl;
			if (!conv.ReadFile(&frag, currentPath.string())) {
				cout << "ERROR READING FRAG" << endl;
			}

			frags.push_back(frag);

		}

		OBMol mol;

		path currentPath = dir->path();
		path localFilenamePath = currentPath.filename();
		string localFileame = localFilenamePath.string();

		vector < string > strs;
		boost::split(strs, localFileame, boost::is_any_of("."));
		string filename = strs[0];
		string ext = strs[1];

		if (boost::iequals(ext, "dat")) {
			continue;
		}

		conv.ReadFile(&mol, currentPath.string());

		cout << "Performing Substitution on: " << filename << endl;

		vector<OBBond*>::iterator j;
		OBBond * bond;

		vector<OBBond *> availableConnectionPoints;
		for (bond = mol.BeginBond(j); bond; bond = mol.NextBond(j)) {
			OBAtom * beginAtom = bond->GetBeginAtom();
			OBAtom * endAtom = bond->GetEndAtom();
			if (beginAtom->GetType()[0] == 'H'
					|| endAtom->GetType()[0] == 'H') {
				//cout << beginAtom->GetType() << " " << endAtom->GetType() << endl;
				availableConnectionPoints.push_back(bond);
			}
		}

		string outputFile = xyzs + "/" + filename + "_connections.dat";

		cout << outputFile << endl;
		ofstream s(outputFile);

		for (int i = 0; i < availableConnectionPoints.size(); i++) {

			OBBond * bond = availableConnectionPoints[i];

			OBAtom * beginAtom = bond->GetBeginAtom();
			OBAtom * endAtom = bond->GetEndAtom();
			cout << "Substituting mol atoms: " << beginAtom->GetType() << " "
					<< endAtom->GetType() << endl;

			int molatom1 = -1;
			int molatom2 = -1;
			if (beginAtom->GetType()[0] == 'H') {
				s << endAtom->GetIdx() << " " << beginAtom->GetIdx() << "\n";
				molatom1 = endAtom->GetIdx();
				molatom2 = beginAtom->GetIdx();
			} else if (endAtom->GetType()[0] == 'H') {
				s << beginAtom->GetIdx() << " " << endAtom->GetIdx() << "\n";
				molatom1 = beginAtom->GetIdx();
				molatom2 = endAtom->GetIdx();
			}

			int arr[3][2] = { { 1, 2 }, { 2, 1 }, { 1, 2 } };
			string fragnames[3] = { "OH", "NH2", "CH3" };

			for (int j = 0; j < frags.size(); j++) {

				cout << "--Frag: " << fragnames[j] << " "
						<< std::to_string(molatom1) << "-"
						<< std::to_string(molatom2) << endl;

				OBMol frag = frags[j];
				OBMol copyMol = mol;

				OBAtom* molatom1P = copyMol.GetAtom(molatom1);
				OBPairInteger *connectionPointM1 = new OBPairInteger;
				connectionPointM1->SetAttribute("molatom1");
				connectionPointM1->SetValue(0);
				molatom1P->SetData(connectionPointM1);

				OBAtom* molatom2P = copyMol.GetAtom(molatom2);
				OBPairInteger *connectionPointM2 = new OBPairInteger;
				connectionPointM2->SetAttribute("molatom2");
				connectionPointM2->SetValue(0);
				molatom2P->SetData(connectionPointM2);

				copyMol.BeginModify();

				OBAtom* fragatom1P = frag.GetAtom(arr[j][0]);
				OBPairInteger *connectionPointF1 = new OBPairInteger;
				connectionPointF1->SetAttribute("fragatom1");
				connectionPointF1->SetValue(0);
				fragatom1P->SetData(connectionPointF1);

				OBAtom* fragatom2P = frag.GetAtom(arr[j][1]);
				OBPairInteger *connectionPointF2 = new OBPairInteger;
				connectionPointF2->SetAttribute("fragatom2");
				connectionPointF2->SetValue(0);
				fragatom2P->SetData(connectionPointF2);

				build.addFragmentToMol(copyMol, frag);

				int ma1index = build.getAtomIndexFromData(copyMol, "molatom1",
						0);
				int ma2index = build.getAtomIndexFromData(copyMol, "molatom2",
						0);
				int fa1index = build.getAtomIndexFromData(copyMol, "fragatom1",
						0);
				int fa2index = build.getAtomIndexFromData(copyMol, "fragatom2",
						0);

				build.connect(copyMol, frag, ma1index, ma2index, fa1index,
						fa2index, 0, 1, false);

				copyMol.EndModify();

				string finalName = filename + "_" + std::to_string(molatom1)
						+ "-" + std::to_string(molatom2) + "_" + fragnames[j];

				if (!conv.WriteFile(&copyMol,
						outputpath + finalName + ".xyz")) {
					cout << "strange error" << endl;
				} else {
					cout << "next loop" << endl;
				}
			}
		}

		s.close();
	}
	cout << "DONE" << endl;
}

void writeHydrogenContainingBonds(string XYZDIR, string OUTDIR) {
//	cout << "Start" << endl;
	OBConversion conv;
	conv.SetInAndOutFormats("xyz", "xyz");

	for (boost::filesystem::recursive_directory_iterator end, dir(XYZDIR);
			dir != end; ++dir) {

		OBMol mol;

		path currentPath = dir->path();
		path localFilenamePath = currentPath.filename();
		string localFileame = localFilenamePath.string();

		vector < string > strs;
		boost::split(strs, localFileame, boost::is_any_of("."));
		string filename = strs[0];
		string ext = strs[1];

		if (boost::iequals(ext, "xyz")) {

			conv.ReadFile(&mol, currentPath.string());

			vector<OBBond*>::iterator j;
			OBBond * bond;

			vector<OBBond *> availableConnectionPoints;
			for (bond = mol.BeginBond(j); bond; bond = mol.NextBond(j)) {
				OBAtom * beginAtom = bond->GetBeginAtom();
				OBAtom * endAtom = bond->GetEndAtom();
				if (beginAtom->GetType()[0] == 'H'
						|| endAtom->GetType()[0] == 'H') {
					//cout << beginAtom->GetType() << " " << endAtom->GetType() << endl;
					availableConnectionPoints.push_back(bond);
				}
			}

			string outputFile = OUTDIR + "/" + filename + "_c.dat";

			cout << outputFile << endl;

			ofstream s(outputFile);

			for (int i = 0; i < availableConnectionPoints.size(); i++) {

				OBBond * bond = availableConnectionPoints[i];

				OBAtom * beginAtom = bond->GetBeginAtom();
				OBAtom * endAtom = bond->GetEndAtom();

				int molatom1 = -1;
				int molatom2 = -1;
				if (beginAtom->GetType()[0] == 'H') {
					s << endAtom->GetIdx() << " " << beginAtom->GetIdx()
							<< "\n";
					molatom1 = endAtom->GetIdx();
					molatom2 = beginAtom->GetIdx();
				} else if (endAtom->GetType()[0] == 'H') {
					s << beginAtom->GetIdx() << " " << endAtom->GetIdx()
							<< "\n";
					molatom1 = beginAtom->GetIdx();
					molatom2 = endAtom->GetIdx();
				}
			}

			s.close();
		}
	}
	cout << "DONE" << endl;
}

void testMolAddition() {
	OBConversion conv;
	OBMol mol;
	OBBuilder build;

	conv.SetInAndOutFormats("smi", "mol2");
	conv.ReadString(&mol, "c1cc2ccccc2cc1");

	mol.AddHydrogens(false, false, 7.0);

	build.Build(mol);

	if (!conv.WriteFile(&mol, "napth.mol2")) {
		cout << "Error when writing finalPath" << endl;
		exit(0);
	}

	OBMol frag;
	conv.ReadString(&frag, "C");
	frag.AddHydrogens(false, false, 7.0);
	build.Build(frag);

	OBAtom * molAtomConnect1 = mol.GetAtom(1);
	OBAtom * molAtomConnect2 = mol.GetAtom(11);
	molAtomConnect1->SetId(999);
	molAtomConnect2->SetId(1000);

	OBAtom * fragAtomConnect1 = frag.GetAtom(1);
	OBAtom * fragAtomConnect2 = frag.GetAtom(2);

//fragAtomConnect1->SetId(1001);
//	fragAtomConnect2->SetId(1002);

	molAtomConnect1->GetResidue()->SetName("MOL");
	fragAtomConnect2->GetResidue()->SetName("FRAG");

	if (!conv.WriteFile(&frag, "CH4.mol2")) {
		cout << "Error when writing finalPath" << endl;
		exit(0);
	}

	OBMol newMol(mol);

	mb.addFragmentToMol(newMol, frag);

	for (unsigned int i = 1; i < mol.NumAtoms(); i++) {
		OBAtom * cAtom = newMol.GetAtom(i);
		//cout << cAtom->GetIdx() << ", " << cAtom->GetCIdx() << ", "
		//		<< cAtom->GetType() << " " << cAtom->GetId() << endl;
	}

	if (!conv.WriteFile(&newMol, "napth_CH4.mol2")) {
		cout << "Error when writing finalPath" << endl;
		exit(0);
	}
}

void testPsiPhiTorsion() {
	OBConversion conv;
	OBMol mol;
	OBBuilder build;

	conv.SetInAndOutFormats("smi", "mol2");
	conv.ReadString(&mol, "OC([C@@H](N)CCSC)=O");

	mol.AddHydrogens(false, false, 7.0);

	build.Build(mol);

	if (!conv.WriteFile(&mol, "M.mol2")) {
		cout << "Error when writing finalPath" << endl;
		exit(0);
	}

	mb.setPhiPsi(mol, 0, -175, 180);

	if (!conv.WriteFile(&mol, "M_rot.mol2")) {
		cout << "Error when writing finalPath" << endl;
		exit(0);
	}
}

int main3(int argc, char * argv[]) {
//testPsiPhiTorsion();
	//testConnect();

	/*OBConversion conv;
	 OBBuilder build;
	 conv.SetInAndOutFormats("smi", "mol2");

	 OBMol mol;
	 conv.ReadString(&mol, "O=C(O)[C@@H](N)CCCNC(N)=N");

	 build.Build(mol);

	 conv.WriteFile(&mol, "arginine.mol2"); */
}

void loadFragmentLibrary() {
	for (boost::filesystem::recursive_directory_iterator end, dir(
			fragmentLibraryPath); dir != end; ++dir) {

		path currentPath = dir->path();

		path fileNamePath = currentPath.filename();

		string pathString = currentPath.string();
		string fileName = fileNamePath.string();

		if (currentPath.has_extension()) {

			string extension = currentPath.extension().string();
			boost::replace_first(extension, ".", "");

			if (conv.FindFormat(extension)) {

				fragmentPathMap.insert(
						std::pair<string, string>(fileName, pathString));
				fragmentExtensionMap.insert(
						std::pair<string, string>(fileName, extension));
			}
		}
	}

	cout << "Fragment library size: " << fragmentPathMap.size() << endl;

}

bool isSubstring(string str, string sub) {
	return (str.find(sub) != std::string::npos);
}

bool readInputFile(string inputFilePath) {

	std::ifstream infile(inputFilePath.c_str());
	std::string line;

	vector < string > fileData;

	while (std::getline(infile, line)) {

		if (isSubstring(line, "#") || isSubstring(line, "//")) {
			continue;
		}

		if (line.find_first_not_of(' ') != std::string::npos) {
			fileData.push_back(line);
			cout << line << endl;
		}
	}

	for (unsigned int i = 0; i < fileData.size(); i++) {
		string line = fileData[i];
		boost::algorithm::trim(line);

		std::vector < std::string > results;

		if (isSubstring(line, "%SET_PSIPHI")) {

			for (unsigned int j = i + 1; j < fileData.size(); j++) {
				string newLine = fileData[j];

				boost::split(results, newLine, boost::is_any_of(" "));

				if (isSubstring(newLine, "%")) {
					break;
				}

				psiphiResidueIndex.push_back(stoi(results[0]));
				psiAngle.push_back(stof(results[1]));
				phiAngle.push_back(stof(results[2]));
			}

		} else if (isSubstring(line, "%ACTION")) {

			for (unsigned int j = i + 1; j < fileData.size(); j++) {

				string newLine = fileData[j];

				boost::split(results, newLine, boost::is_any_of(" "));

				if (isSubstring(newLine, "%")) {
					break;
				}

				string action = results[0];

				if (action == "SUBST") {
					actions.push_back(action);
					molatom1Index.push_back(stoi(results[1]));
					molatom2Index.push_back(stoi(results[2]));
					fragatom1Index.push_back(stoi(results[3]));
					fragatom2Index.push_back(stoi(results[4]));
					fragmentSubstitution.push_back(results[5]);
					bondO.push_back(stoi(results[6]));
					referenceBondTorsionFromFragment.push_back(
							stoi(results[7]));
					updateResidueInfo.push_back(stoi(results[8]));
				} else {
					cout << "action not recognised  \"" << results[0]
							<< "\", exiting program." << endl;
					return false;
				}

			}
		} else if (isSubstring(line, "%LOAD_MOL")) {

			if (isSubstring(fileData[i + 1], "%") || i + 1 >= fileData.size()) {
				break;
			}

			loadMolPath = fileData[i + 1];
		} else if (isSubstring(line, "%FRAGMENT_LIBRARY")) {

			if (isSubstring(fileData[i + 1], "%") || i + 1 >= fileData.size()) {
				break;
			}

			fragmentLibraryPath = fileData[i + 1];
		} else if (isSubstring(line, "%OUT_FILE")) {

			if (isSubstring(fileData[i + 1], "%") || i + 1 >= fileData.size()) {
				break;
			}

			outFilePath = fileData[i + 1];

			std::vector < std::string > results;
			boost::split(results, outFilePath, boost::is_any_of("."));
			outFileExtension = results[1];
		}
	}
	infile.close();
	return true;
}

void doActions() {

}

void cleanup() {
	for (unsigned int i = 0; i < fragments.size(); i++) {
		delete &fragments[i];
	}

}

int normalstuff(int argc, char * argv[]) {

	inputFilePath = string(argv[1]);

	readInputFile(inputFilePath);

	if (!loadMolPath.empty()) {

		string extension = loadMolPath.substr(
				loadMolPath.find_last_of(".") + 1);

		if (!conv.SetInFormat(extension.c_str())) {
			cout << "OBConversion failed to set extension \"" << extension
					<< "\"" << endl;
			return -1;
		}

		if (!conv.ReadFile(&molecule, loadMolPath)) {
			cout << "OBConversion failed to read in molecule from path \""
					<< loadMolPath << "\"" << endl;
			return -1;
		}
	}

	cout << "Num Residues: " << molecule.NumResidues() << endl;
	build.fixProlines(molecule);
	for (int a = 0; a < molecule.NumResidues(); a++) {
		OBResidue * res = molecule.GetResidue(a);

		if (res == NULL) {
			cout << "res idx: " << a << " is null." << endl;
		} else {
			/*
			 cout << a << ": " << res->GetName() << ", " << res->GetNumAtoms()
			 << endl;

			 vector<OBBond*> resbonds = res->GetBonds(false);

			 cout << "---BONDS---" << endl;

			 for (int b = 0; b < resbonds.size(); b++) {
			 OBAtom* beginAtom = resbonds[b]->GetBeginAtom();
			 OBAtom* endAtom = resbonds[b]->GetEndAtom();

			 cout << res->GetAtomID(beginAtom) << "---"
			 << res->GetAtomID(endAtom) << endl;
			 }

			 OBAtom * CA = build.getAlphaCarbon(res);
			 OBAtom * N = build.getBackboneNitrogen(res);
			 OBAtom * C = build.getBackboneCarboxylCarbon(res);
			 OBAtom * N_C = build.getBackboneNitrogen_CarboxylCarbon(res);
			 OBAtom * CB = build.getBetaAtom(res);
			 OBAtom * C_N = build.getCarboxylCarbon_Nitrogen(res);

			 vector<OBAtom*> atoms;

			 molecule.FindChildren(atoms, C, CA);

			 cout << "psi update size: " << atoms.size() << endl;

			 for (int b = 0; b < atoms.size(); b++) {
			 OBResidue *atomRes = atoms[b]->GetResidue();
			 cout << atomRes->GetAtomID(atoms[b]) << ", "
			 << atoms[b]->GetIdx() << ", " << atomRes->GetName()
			 << endl;
			 }
			 atoms.clear();
			 molecule.FindChildren(atoms, N, CA);

			 cout << "phi update size: " << atoms.size() << endl;

			 for (int b = 0; b < atoms.size(); b++) {
			 OBResidue *atomRes = atoms[b]->GetResidue();
			 cout << atomRes->GetAtomID(atoms[b]) << ", "
			 << atoms[b]->GetIdx() << ", " << atomRes->GetName()
			 << endl;
			 }

			 build.printPhiPsi(molecule, a);
			 } */

		}

		if (!fragmentLibraryPath.empty()) {
			loadFragmentLibrary();

			if (fragmentPathMap.size() == 0) {
				cout << "Fragment library failed to load, exiting program"
						<< endl;
				return -1;
			}
		}

		for (unsigned int i = 0; i < fragmentSubstitution.size(); i++) {

			cout << fragmentSubstitution[i] << ", "
					<< fragmentPathMap[fragmentSubstitution[i]] << ","
					<< fragmentExtensionMap[fragmentSubstitution[i]] << endl;
			OBMol * frag = new OBMol;

			if (!conv.SetInFormat(
					fragmentExtensionMap[fragmentSubstitution[i]].c_str())) {
				cout << "OBConversion failed to set extension \""
						<< fragmentExtensionMap[fragmentSubstitution[i]]
						<< "\" for fragment \"" << fragmentSubstitution[i]
						<< "\"" << endl;
				return -1;
			}

			if (!conv.ReadFile(frag,
					fragmentPathMap[fragmentSubstitution[i]])) {
				cout << "OBConversion failed to read in fragment from path \""
						<< fragmentPathMap[fragmentSubstitution[i]] << "\""
						<< endl;
				return -1;
			}

			cout << "fragment loaded: " << frag->GetTitle() << endl;

			fragments.push_back(frag);
		}

		doActions();

		if (!conv.SetOutFormat(outFileExtension.c_str())) {
			cout << "OBConversion failed to set out extension \""
					<< outFileExtension << "\" for molecule" << endl;
			return -1;
		}

		if (!conv.WriteFile(&molecule, outFilePath)) {
			cout << "OBConversion failed to write to outfile \"" << outFilePath
					<< endl;
			return -1;
		}

		int molStartNumAtoms = molecule.NumAtoms();

//cout << "ACTION SIZE:" << actions.size() << endl;

		for (unsigned int a = 0; a < actions.size(); a++) {
			//cout << "IDXS: " << molatom1Index[a] << " " << molatom2Index[a] << " " << fragatom1Index[a] << " " << fragatom2Index[a] << endl;

			OBMol * frag = fragments[a];

			if (molatom1Index[a] != -1) {
				OBAtom* molatom1 = molecule.GetAtom(molatom1Index[a]);

				OBPairInteger *connectionPoint = new OBPairInteger;
				connectionPoint->SetAttribute("molatom1");
				connectionPoint->SetValue(a);
				molatom1->SetData(connectionPoint);

			}

			if (molatom2Index[a] != -1) {

				OBAtom* molatom2 = molecule.GetAtom(molatom2Index[a]);

				OBPairInteger *connectionPoint = new OBPairInteger;
				connectionPoint->SetAttribute("molatom2");
				connectionPoint->SetValue(a);
				molatom2->SetData(connectionPoint);
			}

			if (fragatom1Index[a] != -1) {

				OBAtom* fragatom1 = frag->GetAtom(fragatom1Index[a]);
				OBPairInteger *connectionPoint = new OBPairInteger;
				connectionPoint->SetAttribute("fragatom1");
				connectionPoint->SetValue(a);
				fragatom1->SetData(connectionPoint);
			}

			if (fragatom2Index[a] != -1) {

				OBAtom* fragatom2 = frag->GetAtom(fragatom2Index[a]);

				OBPairInteger *connectionPoint = new OBPairInteger;
				connectionPoint->SetAttribute("fragatom2");
				connectionPoint->SetValue(a);
				fragatom2->SetData(connectionPoint);
				//cout << fragatom2->GetData().size() << endl;
			}
		}

//cout << fragments.size() << " " << molatom1Index.size() << " " << molatom2Index.size() << " " << fragatom1Index.size() << " "
//		<< fragatom2Index.size() << " " << bondO.size() << " " << referenceBondTorsionFromFragment.size() << " "
//		<< updateResidueInfo.size() << endl;

//cout << "Starting building." << endl;
		for (unsigned int a = 0; a < fragments.size(); a++) {
			OBMol * frag = fragments[a];
			molecule.BeginModify();
			build.addFragmentToMol(molecule, *frag);

			int ma1index = build.getAtomIndexFromData(molecule, "molatom1", a);
			int ma2index = build.getAtomIndexFromData(molecule, "molatom2", a);
			int fa1index = build.getAtomIndexFromData(molecule, "fragatom1", a);
			int fa2index = build.getAtomIndexFromData(molecule, "fragatom2", a);
			int bo = bondO[a];
			bool refTorsion = referenceBondTorsionFromFragment[a];
			bool updateRes = updateResidueInfo[a];

			//cout << " " << ma1index << " " << ma2index << " " << fa1index << " " << fa2index << " " << bo << " " << refTorsion << " "
			//		<< updateRes << endl;

			build.connect(molecule, *frag, ma1index, ma2index, fa1index,
					fa2index, a, bo, updateRes);
			//cout << "Finished connection" << endl;
			molecule.EndModify();

			if (refTorsion) {
				ma1index = build.getAtomIndexFromData(molecule, "molatom1", a);
				fa2index = build.getAtomIndexFromData(molecule, "fragatom2", a);
				//cout << "Attempting to reference torsion by AtomID profile. " << ma1index << " " << fa2index << endl;
				build.referenceTorsionByAtomIDProfile(molecule, *frag, ma1index,
						fragatom1Index[a], fragatom2Index[a], fa2index);
			}

		}

		for (int a = 0; a < psiAngle.size(); a++) {
			build.setPhiPsi(molecule, psiphiResidueIndex[a], phiAngle[a],
					psiAngle[a]);
		}
//cout << "Finished building." << endl;
		if (!conv.SetOutFormat(outFileExtension.c_str())) {
			cout << "Could not set output file format" << endl;
			return -1;
		}

		if (!conv.WriteFile(&molecule, outFilePath)) {
			cout << "Could not write out file" << endl;
			return -1;
		}

		cout << "FINISHED" << endl;
		return 0;
	}

int main(int argc, char * argv[]) {
	normalstuff(argc, argv);
//testConnect();
//writeHydrogenContainingBonds("/Users/roethlisbergergroup/EPFL/GA-QMML/overfitData/direct/OptClasses/DM/10",
//			"/Users/roethlisbergergroup/EPFL/GA-QMML/overfitData/direct/OptClasses/DM/10/closestgdb9");
}


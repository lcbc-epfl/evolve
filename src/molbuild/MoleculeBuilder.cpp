/**MoleculeBuilder.cpp
 * Author: Nicholas Browning
 * **/

#include "MoleculeBuilder.h"

#include <openbabel/builder.h>
#include <openbabel/oberror.h>
#include <openbabel/mol.h>
#include <string>
#include <vector>
#include <cmath>
#include <openbabel/math/vector3.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>
#include <boost/algorithm/string/trim.hpp>

using namespace OpenBabel;

bool find_str(string parent, string tofind) {
	std::size_t found = parent.find(tofind);
	return found != std::string::npos;
}

OBAtom * MoleculeBuilder::getAlphaCarbon(OBResidue * res) {
	vector<OBAtom *>::iterator i;
	OBAtom * atom;
	for (atom = res->BeginAtom(i); atom; atom = res->NextAtom(i)) {
		string atomID = res->GetAtomID(atom);
		//cout << "AtomID:" << atomID << "|" << endl;
		boost::trim(atomID);
		if (atomID == "CA") {
			return atom;
		}
	}
	cout << "FAIL (getAlphaCarbon) -> " << "residue " << res->GetIdx() << " does not seem to have a C_alpha" << endl;
	return NULL;
}

/** returns true if atom is nitrogen and bonded to alpha carbon**/
bool MoleculeBuilder::isBackboneNitrogen(OBResidue * res, OBAtom * atom) {

	OBAtom *nbratom;
	OBBond *bond;

	OBBondIterator i;

	if (!atom->IsNitrogen()) {
		return false;
	}
	for (bond = atom->BeginBond(i); bond; bond = atom->NextBond(i)) {

		nbratom = bond->GetNbrAtom(atom);
		string atomID = res->GetAtomID(nbratom);
		boost::trim(atomID);
		if (atomID == "CA") {
			return true;
		}
	}
	return false;
}

OBAtom * MoleculeBuilder::getBackboneNitrogen(OBResidue * res) {

	OBAtom *nbratom, *atom;
	OBBond *bond;

	OBBondIterator i;

	atom = getAlphaCarbon(res);

	if (atom == NULL) {
		return NULL;
	}

	for (bond = atom->BeginBond(i); bond; bond = atom->NextBond(i)) {

		nbratom = bond->GetNbrAtom(atom);

		if (!nbratom->IsNitrogen()) {
			continue;
		}

		return nbratom;
	}
	cout << "FAIL (getAmideNitrogen) -> " << "residue " << res->GetIdx() << " does not seem to have a backbone nitrogen" << endl;
	return NULL;
}

OBAtom * MoleculeBuilder::getBackboneNitrogen_CarboxylCarbon(OBResidue * res) {

	OBAtom *nbratom, *natom, *catom;
	OBBond *bond, *abbond;

	OBBondIterator i, j;

	natom = getBackboneNitrogen(res);
	catom = getAlphaCarbon(res);

	if (natom == NULL) {
		return NULL;
	}

	for (bond = natom->BeginBond(i); bond; bond = natom->NextBond(i)) {

		nbratom = bond->GetNbrAtom(natom);

		if (nbratom == catom) {
			continue;
		}

		if (nbratom->IsHydrogen()) {
			continue;
		}

		for (abbond = nbratom->BeginBond(j); abbond; abbond = nbratom->NextBond(j))
			if (abbond->GetBO() == 2 && (abbond->GetNbrAtom(nbratom))->GetAtomicNum() == 8) {
				return nbratom;
			} else {
				continue;
			}
	}
	cout << "FAIL (getBackboneNitrogen_CarboxylCarbon) -> " << "residue " << res->GetIdx() << " does not seem to have a backbone (N)-C=(O)"
			<< endl;
	return NULL;
}

bool MoleculeBuilder::isBackboneCarboxylCarbon(OBAtom * atom) {

	OBAtom *nbratom;
	OBBond *abbond, *bond;

	OBBondIterator i, j;

	for (bond = atom->BeginBond(i); bond; bond = atom->NextBond(i)) {

		nbratom = bond->GetNbrAtom(atom);

		if (!nbratom->IsCarbon()) {
			continue;
		}

		for (abbond = nbratom->BeginBond(j); abbond; abbond = nbratom->NextBond(j))
			if (abbond->GetBO() == 2 && (abbond->GetNbrAtom(nbratom))->GetAtomicNum() == 8) {
				return true;
			}

	}
	return false;
}

OBAtom * MoleculeBuilder::getBackboneCarboxylCarbon(OBResidue * res) {

	OBAtom * nbratom, *atom;
	OBBond *abbond, *bond;

	OBBondIterator i, j;

	atom = getAlphaCarbon(res);

	if (atom == NULL) {
		return NULL;
	}

	for (bond = atom->BeginBond(i); bond; bond = atom->NextBond(i)) {

		nbratom = bond->GetNbrAtom(atom);

		if (!nbratom->IsCarbon()) {
			continue;
		}

		for (abbond = nbratom->BeginBond(j); abbond; abbond = nbratom->NextBond(j)) {

			if (abbond->GetBO() == 2 && (abbond->GetNbrAtom(nbratom))->GetAtomicNum() == 8) {
				return nbratom;
			}
		}

	}
	cout << "FAIL (getBackboneCarboxylCarbon) -> " << "residue " << res->GetIdx() << " does not seem to have a backbone Carbon" << endl;
	return NULL;
}

OBAtom * MoleculeBuilder::getCarboxylCarbon_Nitrogen(OBResidue * res) {

	OBAtom * nbratom, *atom;
	OBBond * bond;

	OBBondIterator i, j;

	atom = getBackboneCarboxylCarbon(res);

	if (atom == NULL) {
		return NULL;
	}

	for (bond = atom->BeginBond(i); bond; bond = atom->NextBond(i)) {

		nbratom = bond->GetNbrAtom(atom);

		if (!nbratom->IsNitrogen()) {
			continue;
		}

		return nbratom;
	}

	cout << "FAIL (getCarboxylOxygen) -> " << "residue " << res->GetIdx() << " does not seem to have a backbone (O=C)-N-(C_a)" << endl;
	return NULL;
}

OBAtom * MoleculeBuilder::getBetaAtom(OBResidue * res) {

	OBAtom * nbratom, *atom;
	OBBond * bond;

	OBBondIterator i, j;

	atom = getAlphaCarbon(res);

	OBAtom * Natom = getBackboneNitrogen(res);
	OBAtom * CAtom = getBackboneCarboxylCarbon(res);

	if (atom == NULL) {
		return NULL;
	}

	OBAtom * betaAtom = NULL;
	for (bond = atom->BeginBond(i); bond; bond = atom->NextBond(i)) {

		nbratom = bond->GetNbrAtom(atom);

		if (nbratom == Natom || nbratom == CAtom) {
			continue;
		}

		if (nbratom->IsHydrogen()) {
			betaAtom = nbratom;
		} else {
			return nbratom;
		}
	}
	if (betaAtom != NULL) {
		//residue is glycine, 2H @ alphaC
		return betaAtom;
	}

	cout << "FAIL (getBetaAtom) -> " << "residue " << res->GetIdx() << " does not seem to have an atom in beta position" << endl;
	return NULL;

}

void MoleculeBuilder::printPhiPsi(OBMol &mol, int resIndex) {

	OBResidue* res = mol.GetResidue(resIndex);

	OBAtom * CA = getAlphaCarbon(res);

	if (CA == NULL) {
		return;
	}

	OBAtom * N = getBackboneNitrogen(res);
	OBAtom * C = getBackboneCarboxylCarbon(res);
	OBAtom * N_C = getBackboneNitrogen_CarboxylCarbon(res);
	OBAtom * CB = getBetaAtom(res);
	OBAtom * C_N = getCarboxylCarbon_Nitrogen(res);
	double phi = INFINITY, psi = INFINITY;
	if (CA == NULL || CB == NULL) {
		cout << "Residue " << resIndex << ": Cannot read psi or phi torsion: C_alpha or (C|H)_beta were not found" << endl;
		return;
	}

	if (N_C != NULL && N != NULL) {
		//mol.SetTorsion(NH, N, CA, CB, psi_deg * (M_PI / 180));
		//want side chain to remain fixed, so do it the other way
		phi = mol.GetTorsion(N_C, N, CA, C);
	} else {
		cout << "Residue " << resIndex << ": Cannot read phi torsion: N-C=(O) was not completely found" << endl;
	}

	if (C_N != NULL && C != NULL) {
		//mol.SetTorsion(CO, C, CA, CB, phi_deg * (M_PI / 180));
		psi = mol.GetTorsion(N, CA, C, C_N);
	} else {
		cout << "Residue " << resIndex << ": Cannot read psi torsion: C-N-(Ca) was not completely found" << endl;
	}

	cout << "Residue: " << resIndex << "(" << mol.GetResidue(resIndex)->GetName() << "): " << "PSI: " << psi << " PHI: " << phi << endl;
}

void MoleculeBuilder::setPhiPsi(OBMol &mol, int resIndex, double phi_deg, double psi_deg) {

//TODO need to add proline support
//TODO dihedrals should be set ALONG main chain:
//phi: C'-N-----C_a-C'
//psi: N-C_a-----C'-N
//http://images.slideplayer.com/25/7968678/slides/slide_39.jpg
	OBResidue* res = mol.GetResidue(resIndex);

	OBAtom * CA = getAlphaCarbon(res);

	if (CA == NULL) {
		return;
	}

	OBAtom * N = getBackboneNitrogen(res);
	OBAtom * C = getBackboneCarboxylCarbon(res);
	OBAtom * N_C = getBackboneNitrogen_CarboxylCarbon(res);
	OBAtom * CB = getBetaAtom(res);
	OBAtom * C_N = getCarboxylCarbon_Nitrogen(res);

	if (CA == NULL || CB == NULL) {
		cout << "Residue " << resIndex << ": Cannot set psi or phi torsion: C_alpha or (C|H)_beta were not found" << endl;
		return;
	}

	if (phi_deg != INFINITY) {
		if (N_C != NULL && N != NULL) {
			//mol.SetTorsion(NH, N, CA, CB, psi_deg * (M_PI / 180));
			//want side chain to remain fixed, so do it the other way
			mol.SetTorsion(N_C, N, CA, C, phi_deg * (M_PI / 180));
		} else {
			cout << "Residue " << resIndex << ": Cannot set phi torsion: N-C=(O) was not completely found" << endl;
		}

	}
	if (psi_deg != INFINITY) {
		if (C_N != NULL && C != NULL) {
			//mol.SetTorsion(CO, C, CA, CB, phi_deg * (M_PI / 180));
			mol.SetTorsion(N, CA, C, C_N, psi_deg * (M_PI / 180));
		} else {
			cout << "Residue " << resIndex << ": Cannot set psi torsion: C-N-(Ca) was not completely found" << endl;
		}
	}
}

/**Updated version of overloaded += operator in OBMol class**/
void MoleculeBuilder::addFragmentToMol(OBMol &mol, OBMol &fragment) {

	if (mol.GetDimension() != fragment.GetDimension()) {
		cout << "FAIL (addFragmentToMol) -> Cannot add two mols with different dimensions " << mol.GetDimension() + "D + "
				<< fragment.GetDimension() << "D" << endl;
		return;
	}

	vector<OBAtom*>::iterator i;
	vector<OBBond*>::iterator j;
	vector<OBResidue*>::iterator k;
	OBAtom *atom;
	OBBond *bond;
	OBResidue *residue;
	mol.BeginModify();

	int prevatms = mol.NumAtoms();

	string fragTitle(fragment.GetTitle());
	string srcTitle(mol.GetTitle());
	string newTitle = srcTitle + "_" + fragTitle;

	if (!fragTitle.empty())
		mol.SetTitle(newTitle);
	map<unsigned long int, unsigned long int> correspondingId;
// First, handle atoms and bonds
	int id = mol.NumAtoms();
	for (atom = fragment.BeginAtom(i); atom; atom = fragment.NextAtom(i)) {
//	atom->SetId(NoId); //Need to remove ID which relates to source mol rather than this mol
		atom->SetId(id);
		mol.AddAtom(*atom);

		OBAtom *addedAtom = mol.GetAtom(mol.NumAtoms());
		for (OBGenericData* data : atom->GetData()) {
			addedAtom->SetData(data);

		}

		correspondingId[atom->GetId()] = addedAtom->GetId();
		id++;
	}

	for (bond = fragment.BeginBond(j); bond; bond = fragment.NextBond(j)) {
		bond->SetId(NoId); //Need to remove ID which relates to source mol rather than this mol
		mol.AddBond(bond->GetBeginAtomIdx() + prevatms, bond->GetEndAtomIdx() + prevatms, bond->GetBO(), bond->GetFlags());
	}

// Now update all copied residues too
	for (residue = fragment.BeginResidue(k); residue; residue = fragment.NextResidue(k)) {

		mol.AddResidue(*residue);
		OBResidue * molRes = mol.GetResidue(mol.NumResidues() - 1);

		FOR_ATOMS_OF_RESIDUE(resAtom, residue)
		{
			// This is the equivalent atom in our combined molecule
			atom = mol.GetAtom(resAtom->GetIdx() + prevatms);
			// So we add this to the last-added residue
			// (i.e., what we just copied)
			molRes->AddAtom(atom);
		}
	}

	std::vector<OBGenericData*> vdata = fragment.GetAllData(OBGenericDataType::StereoData);
	for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data) {
		OBStereo::Type datatype = ((OBStereoBase*) *data)->GetType();
		if (datatype == OBStereo::CisTrans) {
			OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
			OBCisTransStereo *nct = new OBCisTransStereo(&mol);
			OBCisTransStereo::Config ct_cfg = ct->GetConfig();
			ct_cfg.begin = correspondingId[ct_cfg.begin];
			ct_cfg.end = correspondingId[ct_cfg.end];
			for (OBStereo::RefIter ri = ct_cfg.refs.begin(); ri != ct_cfg.refs.end(); ++ri)
				*ri = correspondingId[*ri];
			nct->SetConfig(ct_cfg);
			mol.SetData(nct);
		} else if (datatype == OBStereo::Tetrahedral) {
			OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
			OBTetrahedralStereo *nts = new OBTetrahedralStereo(&mol);
			OBTetrahedralStereo::Config ts_cfg = ts->GetConfig();
			ts_cfg.center = correspondingId[ts_cfg.center];
			ts_cfg.from = correspondingId[ts_cfg.from];
			for (OBStereo::RefIter ri = ts_cfg.refs.begin(); ri != ts_cfg.refs.end(); ++ri)
				*ri = correspondingId[*ri];
			nts->SetConfig(ts_cfg);
			mol.SetData(nts);
		}
	}

	mol.EndModify();
}

OBAtom * MoleculeBuilder::getAtomByAtomID(OBResidue *res, string type) {
	for (OBAtom * atom : res->GetAtoms()) {
		string atomID = res->GetAtomID(atom);
		boost::trim(atomID);
		if (atomID == type) {
			return atom;
		}
	}
	return NULL;
}

void MoleculeBuilder::fixProlines(OBMol & mol) {
	for (int a = 0; a < mol.NumResidues(); a++) {
		OBResidue* res = mol.GetResidue(a);

		if (res->GetName() == "PRO") {
			std::vector<OBBond *> bonds = res->GetBonds(false);

			OBAtom * CA = getAtomByAtomID(res, "CA");
			OBAtom * CB = getAtomByAtomID(res, "CB");
			OBAtom * CY = getAtomByAtomID(res, "CG");
			OBAtom * CD = getAtomByAtomID(res, "CD");

			if (CA == NULL || CB == NULL || CY == NULL || CD == NULL) {
				cout << "Cannot find proline ring carbons! Torsion updates will not work." << endl;
				cout << CA << " " << CB << " " << CY << " " << CD << endl;
				return;
			}

			if (!CA->IsConnected(CB)) {
				cout << "CA and CB are not connected in OBMol, connecting." << endl;
				OBBond *bond = mol.NewBond();
				bond->SetBegin(CA);
				bond->SetEnd(CB);

				bond->SetBondOrder(1);

				CA->AddBond(bond);
				CB->AddBond(bond);
				bond->SetLength(CA->GetDistance(CB));
				cout << "New bond length: " << CA->GetDistance(CB) << endl;
			}

			if (!CB->IsConnected(CY)) {
				cout << "CB and CY are not connected in OBMol, connecting." << endl;
				OBBond *bond = mol.NewBond();
				bond->SetBegin(CB);
				bond->SetEnd(CY);

				bond->SetBondOrder(1);

				CB->AddBond(bond);
				CY->AddBond(bond);
				bond->SetLength(CB->GetDistance(CY));
				cout << "New bond length: " << CB->GetDistance(CY) << endl;
			}
			if (!CY->IsConnected(CD)) {
				cout << "CY and CD are not connected in OBMol, connecting." << endl;
				OBBond *bond = mol.NewBond();
				bond->SetBegin(CY);
				bond->SetEnd(CD);

				bond->SetBondOrder(1);

				CY->AddBond(bond);
				CD->AddBond(bond);
				bond->SetLength(CY->GetDistance(CD));
				cout << "New bond length: " << CY->GetDistance(CD) << endl;
			}

		}
	}
}

int MoleculeBuilder::getAtomIdFromData(OBMol &mol, string dataVal, int value) {
	OBAtom * atom;
	vector<OBAtom*>::iterator i;
	for (atom = mol.BeginAtom(i); atom; atom = mol.NextAtom(i)) {
		if (atom->HasData(dataVal)) {
			OBPairInteger *data = dynamic_cast<OBPairInteger*>(atom->GetData(dataVal.c_str()));
			if (data->GetGenericValue() == value) {
				return atom->GetId();
			}
		}
	}
	return -1;
}

int MoleculeBuilder::getAtomIndexFromData(OBMol &mol, string dataVal, int value) {
	OBAtom * atom;
	vector<OBAtom*>::iterator i;
	for (atom = mol.BeginAtom(i); atom; atom = mol.NextAtom(i)) {
		if (atom->HasData(dataVal)) {
			OBPairInteger *data = dynamic_cast<OBPairInteger*>(atom->GetData(dataVal.c_str()));
			if (data->GetGenericValue() == value) {
				return atom->GetIdx();
			}
		}
	}
	return -1;
}
/** Substitutes group or atom at matom2 (matom1 - matom2) with atom or group of fatom2.
 * matom2 and anything attached to it (except matom1) is deleted from OBMol mol. Anything attached to fatom1 (except fatom2) is also deleted from OBMol fragment.
 *
 *useFragmentTorsion updates new matom1 - fatom2 bond torsion to correspond to that of fragment. R1 - matom1 - matom2 - R2 molecular signature must therefore match
 * R1 - fatom1 - fatom2 -R2
 * **/

void MoleculeBuilder::connect(OBMol &mol, OBMol &fragment, int ma1index = -1, int ma2index = -1, int fa1index = -1, int fa2index = -1,
		int index = -1, int bondOrder = 1, bool updateResidueInfo = false) {

	cout << ma1index << " " << ma2index << " " << fa1index << " " << fa2index << endl;
	OBBuilder build;

	OBAtom * molatom1 = mol.GetAtom(ma1index);

	if (molatom1 == NULL) {
		cout << "mol atom 1 should not be null." << endl;
		return;
	}

	OBAtom * molatom2 = NULL;

	if (ma2index != -1) {
		molatom2 = mol.GetAtom(ma2index);
	}

	OBAtom * fragatom1 = NULL;
	OBAtom * fragatom2 = mol.GetAtom(fa2index);

	if (fa1index != -1) {
		fragatom1 = mol.GetAtom(fa1index);
	}

	vector<OBAtom*> deleteAtoms;

	vector<OBAtom*> molDelete;
	vector<OBAtom*> fragDelete;

	cout << fragatom1->GetType() << " " << fragatom2->GetType() << endl;

	if (fragatom1 != NULL) {
		mol.FindChildren(fragDelete, fragatom2, fragatom1);
		fragDelete.push_back(fragatom1);
	}

	if (molatom2 != NULL) {
		mol.FindChildren(molDelete, molatom1, molatom2);
		molDelete.push_back(molatom2);
	}

	deleteAtoms.insert(deleteAtoms.end(), fragDelete.begin(), fragDelete.end());
	deleteAtoms.insert(deleteAtoms.end(), molDelete.begin(), molDelete.end());

	OBBitVec fragBitVec = build.GetFragment(fragatom2);

	vector<OBAtom *> fragAtoms;

	int bit;

	vector3 molatom1vec = molatom1->GetVector();

	vector3 fragatom1vec = fragatom1->GetVector();

	bit = fragBitVec.NextBit(0);

	while (bit != fragBitVec.EndBit()) {

		OBAtom * fatom = mol.GetAtom(bit);

//cout << fatom->GetIdx() << " " << fatom->GetType() << endl;

		fragAtoms.push_back(fatom);

		vector3 current = fatom->GetVector();

		vector3 diff = current - fragatom1vec;

		fatom->SetVector(molatom1vec + diff);

		bit = fragBitVec.NextBit(bit);
	}

	vector3 molatom2vec = molatom2->GetVector();
	vector3 molbondvec = molatom2vec - molatom1vec;

	fragatom1vec = fragatom1->GetVector();
	vector3 fragatom2vec = fragatom2->GetVector();
	vector3 fragbondvec = fragatom2vec - fragatom1vec;

//cout << molbondvec << " " << fragbondvec << endl;

// molbondvec and fragbondvec could be equal - no rotation needed
	if (molbondvec != fragbondvec) {
		vector3 rotationAxisU = OBAPI::cross(molbondvec, fragbondvec);

		vector3 rotationAxisN = rotationAxisU.normalize();

//cout << rotationAxisN << endl;

		double dp;
		dp = OBAPI::dot(molbondvec, fragbondvec) / (molbondvec.length() * fragbondvec.length());
		double angle = acos(dp);
		double diffAng = -angle;
//cout << "rotating" << endl;
		rotateAtoms(fragAtoms, fragatom1vec, rotationAxisN, diffAng);

	}

	if (deleteAtoms.size() > 0) {
		for (unsigned int aa = 0; aa < deleteAtoms.size(); aa++) {
			OBResidue* res = deleteAtoms[aa]->GetResidue();
			if (res != NULL) {
				res->RemoveAtom(deleteAtoms[aa]);
			}
			mol.DeleteAtom(deleteAtoms[aa], true);
		}
	}

	fragAtoms.clear();

	fragBitVec = build.GetFragment(fragatom2);

	bit = fragBitVec.NextBit(0);

	while (bit != fragBitVec.EndBit()) {

		OBAtom * fatom = mol.GetAtom(bit);

		fragAtoms.push_back(fatom);

		bit = fragBitVec.NextBit(bit);
	}

	OBBitVec molBitVec = build.GetFragment(molatom1);
	vector<OBAtom*> molAtoms;
	int molbit = molBitVec.NextBit(0);

	while (molbit != molBitVec.EndBit()) {

		OBAtom * fatom = mol.GetAtom(molbit);

		molAtoms.push_back(fatom);

		molbit = molBitVec.NextBit(molbit);
	}

	if (ma2index != -1) {
//insert at C_b position
//cout << molAtoms.size() << " " << fragAtoms.size() << " " << ma2index << endl;
		molAtoms.insert(molAtoms.begin() + ma2index - 1, fragAtoms.begin(), fragAtoms.end());
	} else {
//insert after  C_a position
		molAtoms.insert(molAtoms.begin() + molatom1->GetIdx(), fragAtoms.begin(), fragAtoms.end());
	}



	mol.RenumberAtoms(molAtoms);

	OBBond *bond = mol.NewBond();
	bond->SetBegin(molatom1);
	bond->SetEnd(fragatom2);

	bond->SetBondOrder(bondOrder);

	fragatom2->AddBond(bond);
	molatom1->AddBond(bond);
	bond->SetLength(fragbondvec.length());

	if (updateResidueInfo) {

		cout << "Updating res info" << endl;
//TODO -- DeleteResidue will delete atom names, or can keep them if we delete after doing switch, but then somewhere a perception routine changes
//molRes resName from e.g K24 to LYS

		OBResidue* molRes = molatom1->GetResidue();
		OBResidue* fragRes = fragatom2->GetResidue();

//cout << molRes->GetName() << endl;
//cout << fragRes->GetName() << endl;
		string fragName = fragRes->GetName();

		vector<OBAtom*> fragAtoms = fragRes->GetAtoms();

		vector < string > atomIDs;
		for (unsigned int a = 0; a < fragAtoms.size(); a++) {
			//cout << fragAtoms[a]->GetType() << endl;
			string atomID = fragRes->GetAtomID(fragAtoms[a]);
			atomIDs.push_back(atomID);

			molRes->AddAtom(fragAtoms[a]);
			molRes->SetAtomID(fragAtoms[a], atomIDs[a]);
			fragAtoms[a]->SetResidue(molRes);

		}
		mol.DeleteResidue(fragRes, false);

		molRes->SetName(fragName);
	}
}

/**Redefinition of OBBond::SetLength to remove annoying cerr output: "cerr << "v3: " << v3 << " v4: " << v4 << endl;" **/
void OBBond::SetLength(OBAtom *fixed, double length) {
	unsigned int i;
	OBMol *mol = (OBMol*) fixed->GetParent();
	vector3 v1, v2, v3, v4, v5;
	vector<int> children;

	obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::SetBondLength", obAuditMsg);

	int a = fixed->GetIdx();
	int b = GetNbrAtom(fixed)->GetIdx();

	if (a == b)
		return; // this would be a problem...

	mol->FindChildren(children, a, b);
	children.push_back(b);

	v1 = GetNbrAtom(fixed)->GetVector();
	v2 = fixed->GetVector();
	v3 = v1 - v2;

	if (IsNearZero(v3.length_2())) { // too small to normalize, move the atoms apart
		obErrorLog.ThrowError(__FUNCTION__, "Atoms are both at the same location, moving out of the way.", obWarning);
		v3.randomUnitVector();
	} else {
		v3.normalize();
	}

	v3 *= length;
	v3 += v2;
	v4 = v3 - v1;

	for (i = 0; i < children.size(); i++) {
		v1 = mol->GetAtom(children[i])->GetVector();
		v1 += v4;
		mol->GetAtom(children[i])->SetVector(v1);
	}
}

void MoleculeBuilder::referenceTorsionByAtomIDProfile(OBMol &mol, OBMol &frag, int molatom1index, int fragatom1index, int fragatom2index,
		int molfragatom2index) {

	OBAtom * fragatom1 = frag.GetAtom(fragatom1index);
	OBAtom * fragatom2 = frag.GetAtom(fragatom2index);

	OBAtom * molatom1 = mol.GetAtom(molatom1index);
	OBAtom * molfragatom2 = mol.GetAtom(molfragatom2index);

//cout << molatom1->GetId() << " " << molfragatom2->GetId() << endl;

	OBAtom*ma, *mb, *mc, *md;
	double moltor;

//cout << "sear" << endl;

	outer: FOR_TORSIONS_OF_MOL(t, mol)
	{
		ma = mol.GetAtom((*t)[0] + 1);
		mb = mol.GetAtom((*t)[1] + 1);
		mc = mol.GetAtom((*t)[2] + 1);
		md = mol.GetAtom((*t)[3] + 1);

		moltor = mol.GetTorsion(ma->GetIdx(), mb->GetIdx(), mc->GetIdx(), md->GetIdx());
//cout << mb->GetId() << " " <<mc->GetId() << endl;
		if (mb->GetId() == molatom1->GetId() && mc->GetId() == molfragatom2->GetId()) {

			bool foundMatch = false;

			OBResidue * molresa = ma->GetResidue();
			OBResidue * molresb = mb->GetResidue();
			OBResidue * molresc = mc->GetResidue();
			OBResidue * molresd = md->GetResidue();

			string molatom1ID = molresa->GetAtomID(ma);
			string molatom2ID = molresb->GetAtomID(mb);
			string molatom3ID = molresc->GetAtomID(mc);
			string molatom4ID = molresd->GetAtomID(md);

			boost::trim(molatom1ID);
			boost::trim(molatom2ID);
			boost::trim(molatom3ID);
			boost::trim(molatom4ID);

			cout << molatom1ID << " (" << molresa->GetName() << ") " << molatom2ID << " (" << molresb->GetName() << ") " << molatom3ID
					<< " (" << molresc->GetName() << ") " << molatom4ID << " (" << molresd->GetName() << ") " << moltor << endl;

			OBAtom*fa, *fb, *fc, *fd;

			double fragtor;

			FOR_TORSIONS_OF_MOL(t, frag)
			{
				fa = frag.GetAtom((*t)[0] + 1);
				fb = frag.GetAtom((*t)[1] + 1);
				fc = frag.GetAtom((*t)[2] + 1);
				fd = frag.GetAtom((*t)[3] + 1);

				fragtor = frag.GetTorsion(fa->GetIdx(), fb->GetIdx(), fc->GetIdx(), fd->GetIdx());

				OBResidue * fragresa = fa->GetResidue();
				OBResidue * fragresb = fb->GetResidue();
				OBResidue * fragresc = fc->GetResidue();
				OBResidue * fragresd = fd->GetResidue();

				string fragatom1ID = fragresa->GetAtomID(fa);
				string fragatom2ID = fragresa->GetAtomID(fb);
				string fragatom3ID = fragresa->GetAtomID(fc);
				string fragatom4ID = fragresa->GetAtomID(fd);

				boost::trim(fragatom1ID);
				boost::trim(fragatom2ID);
				boost::trim(fragatom3ID);
				boost::trim(fragatom4ID);

				if (fb->GetId() == fragatom1->GetId() && fc->GetId() == fragatom2->GetId()) {
					cout << "potential match: " << fragatom1ID << " " << fragatom2ID << " " << fragatom3ID << " " << fragatom4ID << " "
							<< fragtor << endl;

					if (molatom1ID.compare(fragatom1ID) == 0 && molatom2ID.compare(fragatom2ID) == 0 && molatom3ID.compare(fragatom3ID) == 0
							&& molatom4ID.compare(fragatom4ID) == 0) {
						cout << "found match:" << fragatom1ID << " (" << fragresa->GetName() << ") " << fragatom2ID << " ("
								<< fragresb->GetName() << ") " << fragatom3ID << " (" << fragresc->GetName() << ") " << fragatom4ID << " ("
								<< fragresd->GetName() << ") " << fragtor << endl;
						mol.SetTorsion(ma, mb, mc, md, fragtor * DEG_TO_RAD);
						foundMatch = true;
						break;
					}

				}
			}

			if (foundMatch) {
				break;
			}

		}

	}
}

void referenceTorsionByElementProfile() {

}
/**finds the torsion in mol which should be updated with the one in frag - used after connect() (and ... addFragmentToMol)**/
void MoleculeBuilder::referenceTorsion(OBMol &mol, OBMol &frag, int molatom1id, int molatom2id, int fragatom1id, int fragatom2id) {

	OBAtom*ma, *mb, *mc, *md;
	double moltor;

	outer: FOR_TORSIONS_OF_MOL(t, mol)
	{

		ma = mol.GetAtom((*t)[0] + 1);
		mb = mol.GetAtom((*t)[1] + 1);
		mc = mol.GetAtom((*t)[2] + 1);
		md = mol.GetAtom((*t)[3] + 1);
		moltor = mol.GetTorsion(ma->GetIdx(), mb->GetIdx(), mc->GetIdx(), md->GetIdx());

		if (mb->GetId() == molatom1id && mc->GetId() == fragatom2id) {

			bool foundMatch = false;

			OBResidue * molresa = ma->GetResidue();
			OBResidue * molresb = mb->GetResidue();
			OBResidue * molresc = mc->GetResidue();
			OBResidue * molresd = md->GetResidue();

			string molatom1ID = molresa->GetAtomID(ma);
			string molatom2ID = molresb->GetAtomID(mb);
			string molatom3ID = molresc->GetAtomID(mc);
			string molatom4ID = molresd->GetAtomID(md);

			//cout << molresa->GetAtomID(ma) << " (" << molresa->GetName() << ") " << molresb->GetAtomID(mb) << " (" << molresb->GetName() << ") " << molresc->GetAtomID(mc) << " (" << molresc->GetName() << ") " << molresd->GetAtomID(md) << " (" << molresd->GetName() << ") "<< moltor<< endl;

			OBAtom*fa, *fb, *fc, *fd;

			double fragtor;

			FOR_TORSIONS_OF_MOL(t, frag)
			{
				fa = frag.GetAtom((*t)[0] + 1);
				fb = frag.GetAtom((*t)[1] + 1);
				fc = frag.GetAtom((*t)[2] + 1);
				fd = frag.GetAtom((*t)[3] + 1);

				fragtor = frag.GetTorsion(fa->GetIdx(), fb->GetIdx(), fc->GetIdx(), fd->GetIdx());

				OBResidue * fragresa = fa->GetResidue();
				OBResidue * fragresb = fb->GetResidue();
				OBResidue * fragresc = fc->GetResidue();
				OBResidue * fragresd = fd->GetResidue();

				string fragatom1ID = fragresa->GetAtomID(fa);
				string fragatom2ID = fragresa->GetAtomID(fb);
				string fragatom3ID = fragresa->GetAtomID(fc);
				string fragatom4ID = fragresa->GetAtomID(fd);

				if (fb->GetId() == fragatom1id && fc->GetId() == fragatom2id) {
					if (molatom1ID == fragatom1ID && molatom2ID == fragatom2ID && molatom3ID == fragatom3ID && molatom4ID == fragatom4ID) {
						cout << "found match:" << fragatom1ID << " (" << fragresa->GetName() << ") " << fragatom2ID << " ("
								<< fragresb->GetName() << ") " << fragatom3ID << " (" << fragresc->GetName() << ") " << fragatom4ID << " ("
								<< fragresd->GetName() << ") " << fragtor << endl;
						mol.SetTorsion(ma, mb, mc, md, fragtor * DEG_TO_RAD);
						foundMatch = true;
						break;
					}
				}
			}

			if (foundMatch) {
				break;
			}

		}

	}
}

void MoleculeBuilder::rotateAtoms(vector<OBAtom *> atoms, vector3 &origin, vector3 &rotationAxis, double angle) {
	double rotMatrix[9];

	double sn = sin(angle);
	double cs = cos(angle);
	double t = 1 - cs;

	double x = rotationAxis.x();
	double y = rotationAxis.y();
	double z = rotationAxis.z();

	rotMatrix[0] = t * x * x + cs;
	rotMatrix[1] = t * x * y + sn * z;
	rotMatrix[2] = t * x * z - sn * y;
	rotMatrix[3] = t * x * y - sn * z;
	rotMatrix[4] = t * y * y + cs;
	rotMatrix[5] = t * y * z + sn * x;
	rotMatrix[6] = t * x * z + sn * y;
	rotMatrix[7] = t * y * z - sn * x;
	rotMatrix[8] = t * z * z + cs;

	for (unsigned int i = 0; i < atoms.size(); i++) {

		OBAtom * atom = atoms.at(i);

		vector3 tempvec = atom->GetVector();

		double xtemp = tempvec.x();
		double ytemp = tempvec.y();
		double ztemp = tempvec.z();

		xtemp -= origin.x();
		ytemp -= origin.y();
		ztemp -= origin.z();

		x = xtemp * rotMatrix[0] + ytemp * rotMatrix[3] + ztemp * rotMatrix[6];
		y = xtemp * rotMatrix[1] + ytemp * rotMatrix[4] + ztemp * rotMatrix[7];
		z = xtemp * rotMatrix[2] + ytemp * rotMatrix[5] + ztemp * rotMatrix[8];

		xtemp = x;
		ytemp = y;
		ztemp = z;

		xtemp += origin.x();
		ytemp += origin.y();
		ztemp += origin.z();

		vector3 newvec(xtemp, ytemp, ztemp);

		atom->SetVector(newvec);
	}
}

void OBMol::FindChildren(vector<OBAtom*> &children, OBAtom *bgn, OBAtom *end) {

	cout << "HAHAHA" << endl;
	OBBitVec visited, curr, next;

	OBResidue *beginRes = bgn->GetResidue();
	OBResidue *endRes = end->GetResidue();

	visited |= bgn->GetIdx();
	visited |= end->GetIdx();
	curr |= end->GetIdx();
	children.clear();

	int i;
	OBAtom *atom;
	vector<OBBond*>::iterator j;

	OBBond * bond;

	OBResidue *res;

	for (;;) {
		next.Clear();
		for (i = curr.NextBit(-1); i != curr.EndBit(); i = curr.NextBit(i)) {

			atom = GetAtom(i);
			//IsConnected
			for (bond = atom->BeginBond(j); bond; bond = atom->NextBond(j)) {
				OBAtom* bondBegin = bond->GetBeginAtom();
				OBAtom* bondEnd = bond->GetEndAtom();

				if (!visited[bondBegin->GetIdx()]) {
					children.push_back(bondBegin);
					next |= bondBegin->GetIdx();
					visited |= bondBegin->GetIdx();
				}

				if (!visited[bondEnd->GetIdx()]) {
					children.push_back(bondEnd);
					next |= bondEnd->GetIdx();
					visited |= bondEnd->GetIdx();
				}

			}
		}
		if (next.Empty())
			break;
		curr = next;
	}
}

//Reimplementation of find children - causes problems for Proline/cycles
void FindChildrenOld(vector<OBAtom*> &children, OBAtom *bgn, OBAtom *end) {
	OBBitVec used, curr, next;

	used |= bgn->GetIdx();
	used |= end->GetIdx();
	curr |= end->GetIdx();
	children.clear();

	int i;
	OBAtom *atom, *nbr;
	vector<OBBond*>::iterator j;

	OBMol * mol = bgn->GetParent();

	for (;;) {
		next.Clear();
		for (i = curr.NextBit(-1); i != curr.EndBit(); i = curr.NextBit(i)) {

			atom = mol->GetAtom(i);

			OBResidue* res = atom->GetResidue();
			string atomID = res->GetAtomID(atom);
			boost::trim(atomID);
			if (atomID == "CD" && res->GetName() == "PRO") {
				cout << "WA!" << endl;
			}

			for (nbr = atom->BeginNbrAtom(j); nbr; nbr = atom->NextNbrAtom(j))
				if (!used[nbr->GetIdx()]) {
					children.push_back(nbr);
					next |= nbr->GetIdx();
					used |= nbr->GetIdx();
				}
		}
		if (next.Empty())
			break;
		curr = next;
	}
}

//Custom Implementation of SetTorsion to account for Proline - FindChildren does not return all the atoms it should
void OBMol::SetTorsion(OBAtom *a, OBAtom *b, OBAtom *c, OBAtom *d, double ang) {

	vector<int> tor;
	vector<int> atoms;

	obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::SetTorsion", obAuditMsg);

	tor.push_back(a->GetCIdx());
	tor.push_back(b->GetCIdx());
	tor.push_back(c->GetCIdx());
	tor.push_back(d->GetCIdx());

	OBResidue* res_a = a->GetResidue();
	OBResidue* res_b = b->GetResidue();
	OBResidue* res_c = c->GetResidue();
	OBResidue* res_d = d->GetResidue();

	if (res_b->GetName() == "PRO" && res_c->GetName() == "PRO") {

	}

	FindChildren(atoms, b->GetIdx(), c->GetIdx());
	int j;
	for (j = 0; (unsigned) j < atoms.size(); j++)
		atoms[j] = (atoms[j] - 1) * 3;

	double v2x, v2y, v2z;
	double radang, m[9];
	double x, y, z, mag, rotang, sn, cs, t, tx, ty, tz;

//calculate the torsion angle
	radang = CalcTorsionAngle(a->GetVector(), b->GetVector(), c->GetVector(), d->GetVector()) / RAD_TO_DEG;
//
// now we have the torsion angle (radang) - set up the rot matrix
//

//find the difference between current and requested
	rotang = ang - radang;

	sn = sin(rotang);
	cs = cos(rotang);
	t = 1 - cs;

	v2x = _c[tor[1]] - _c[tor[2]];
	v2y = _c[tor[1] + 1] - _c[tor[2] + 1];
	v2z = _c[tor[1] + 2] - _c[tor[2] + 2];

//normalize the rotation vector
	mag = sqrt(SQUARE(v2x) + SQUARE(v2y) + SQUARE(v2z));
	x = v2x / mag;
	y = v2y / mag;
	z = v2z / mag;

//set up the rotation matrix
	m[0] = t * x * x + cs;
	m[1] = t * x * y + sn * z;
	m[2] = t * x * z - sn * y;
	m[3] = t * x * y - sn * z;
	m[4] = t * y * y + cs;
	m[5] = t * y * z + sn * x;
	m[6] = t * x * z + sn * y;
	m[7] = t * y * z - sn * x;
	m[8] = t * z * z + cs;

//
//now the matrix is set - time to rotate the atoms
//
	tx = _c[tor[1]];
	ty = _c[tor[1] + 1];
	tz = _c[tor[1] + 2];
	vector<int>::iterator i;
	for (i = atoms.begin(); i != atoms.end(); ++i) {
		j = *i;

		_c[j] -= tx;
		_c[j + 1] -= ty;
		_c[j + 2] -= tz;
		x = _c[j] * m[0] + _c[j + 1] * m[1] + _c[j + 2] * m[2];
		y = _c[j] * m[3] + _c[j + 1] * m[4] + _c[j + 2] * m[5];
		z = _c[j] * m[6] + _c[j + 1] * m[7] + _c[j + 2] * m[8];
		_c[j] = x;
		_c[j + 1] = y;
		_c[j + 2] = z;
		_c[j] += tx;
		_c[j + 1] += ty;
		_c[j + 2] += tz;
	}
}

void MoleculeBuilder::rotateMolecule(OBMol &mol, vector3 &origin, vector3 &rotationAxis, double angle) {
	double rotMatrix[9];

	double sn = sin(angle);
	double cs = cos(angle);
	double t = 1 - cs;

	double x = rotationAxis.x();
	double y = rotationAxis.y();
	double z = rotationAxis.z();

	rotMatrix[0] = t * x * x + cs;
	rotMatrix[1] = t * x * y + sn * z;
	rotMatrix[2] = t * x * z - sn * y;
	rotMatrix[3] = t * x * y - sn * z;
	rotMatrix[4] = t * y * y + cs;
	rotMatrix[5] = t * y * z + sn * x;
	rotMatrix[6] = t * x * z + sn * y;
	rotMatrix[7] = t * y * z - sn * x;
	rotMatrix[8] = t * z * z + cs;

	for (unsigned int i = 0; i < mol.NumAtoms(); i++) {

		OBAtom * atom = mol.GetAtom(i);

		vector3 tempvec = atom->GetVector();

		double xtemp = tempvec.x();
		double ytemp = tempvec.y();
		double ztemp = tempvec.z();

		xtemp -= origin.x();
		ytemp -= origin.y();
		ztemp -= origin.z();

		x = xtemp * rotMatrix[0] + ytemp * rotMatrix[3] + ztemp * rotMatrix[6];
		y = xtemp * rotMatrix[1] + ytemp * rotMatrix[4] + ztemp * rotMatrix[7];
		z = xtemp * rotMatrix[2] + ytemp * rotMatrix[5] + ztemp * rotMatrix[8];

		xtemp = x;
		ytemp = y;
		ztemp = z;

		xtemp += origin.x();
		ytemp += origin.y();
		ztemp += origin.z();

		vector3 newvec(xtemp, ytemp, ztemp);

		atom->SetVector(newvec);
	}
}


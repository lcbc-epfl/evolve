/**MoleculeBuilder.h
 * Author: Nicholas Browning
 * **/

#ifndef MOLECULEBUILDER_H_
#define MOLECULEBUILDER_H_

#include <openbabel/mol.h>
#include <vector>

using namespace OpenBabel;
using namespace std;

class MoleculeBuilder {

public:
	OBAtom * getAlphaCarbon(OBResidue * res);
	OBAtom * getBackboneNitrogen(OBResidue * res);
	OBAtom * getBackboneNitrogen_CarboxylCarbon(OBResidue * res);
	OBAtom * getBackboneCarboxylCarbon(OBResidue * res);
	OBAtom * getCarboxylCarbon_Nitrogen(OBResidue * res);
	OBAtom * getBetaAtom(OBResidue * res);

	OBAtom * getAtomByAtomID(OBResidue * res, string type);

	bool isBackboneNitrogen(OBResidue * res, OBAtom * atom);
	bool isBackboneCarboxylCarbon(OBAtom * atom);

	void setPhiPsi(OBMol &mol, int resIndex, double phi_deg, double psi_deg);
	void printPhiPsi(OBMol &mol, int resIndex);
	void addFragmentToMol(OBMol &mol, OBMol &fragment);
	void connect(OBMol &mol, OBMol &fragment, int matom1Idx, int matom2Idx,
			int fatom1Idx, int fatom2Idx, int a, int bondOrder,
			bool updateResidueInfo);
	void referenceTorsion(OBMol &mol, OBMol &frag, int molatom1id,
			int molatom2id, int fragatom1id, int fragatom2id);
	void referenceTorsionByAtomIDProfile(OBMol &mol, OBMol &frag,
			int molatom1id, int fragatom1id, int fragatom2id,
			int molfragatom2id);
	void rotateMolecule(OBMol &mol, vector3 &origin, vector3 &rotationAxis,
			double angle);
	void rotateAtoms(vector<OBAtom *> atoms, vector3 &origin,
			vector3 &rotationAxis, double angle);

	void fixProlines(OBMol &mol);

	int getAtomIdFromData(OBMol &mol, string dataVal, int value);
	int getAtomIndexFromData(OBMol &mol, string dataVal, int value);

}
;

#endif /* MOLECULEBUILDER_H_ */

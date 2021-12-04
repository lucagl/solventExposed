# solventExposed
A script based on NanoShaper protein surface software returning the list of atoms and residues solvent exposed (or according to a user defined probe radius)

## Requirements
- Install NanoShaper https://gitlab.iit.it/SDecherchi/nanoshaper
- or directly use executable provided in the "NS" folder:
  - To correctly link the libraries install *patchelf* (sudo apt get install patchelf or https://gist.github.com/ruario/80fefd174b3395d34c14)
  then run the "install_script"

## Usage
**Important**:
1. in the working folder the provided *temp* folder must be present with a working NanoShaper executable inside (run the install_script before, see above)
2. The input file must be in PQR format


python3 expRes.py <options> <structure_name>
For help type python3 expRes.py -h


**Default parameters**: Probe_radius=1.4 (water), Thread_number=1 (for NanoShaper triangulation)

**Note**: The probe radius can be changed.

### Example:
Hen Egg-White Lisozyme 1hew.pqr:

python3 expRes.py --savePQR --saveRES 1hew

--> 1hew_exposed.pqr file containing the subset of exposed atoms of the original PQR

--> 1hew_resList.txt file containing the list of exposed residues

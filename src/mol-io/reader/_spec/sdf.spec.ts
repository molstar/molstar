
import { parseSdf } from '../sdf/parser';

const SdfString = `
  Mrv1718007121815122D          

  5  4  0  0  0  0            999 V2000
    0.0000    0.8250    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -0.8250    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.8250    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.8250    0.0000    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
  4  1  1  0  0  0  0
  4  2  2  0  0  0  0
  4  3  1  0  0  0  0
  4  5  1  0  0  0  0
M  CHG  3   1  -1   3  -1   5  -1
M  END
> <DATABASE_ID>
DB14523

> <DATABASE_NAME>
drugbank

> <SMILES>
[O-]P([O-])([O-])=O

> <INCHI_IDENTIFIER>
InChI=1S/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)/p-3

> <INCHI_KEY>
NBIIXXVUZAFLBC-UHFFFAOYSA-K

> <FORMULA>
O4P

> <MOLECULAR_WEIGHT>
94.9714

> <EXACT_MASS>
94.95342

> <JCHEM_ACCEPTOR_COUNT>
4

> <JCHEM_ATOM_COUNT>
5

> <JCHEM_AVERAGE_POLARIZABILITY>
4.932162910070488

> <JCHEM_BIOAVAILABILITY>
1

> <JCHEM_DONOR_COUNT>
0

> <JCHEM_FORMAL_CHARGE>
-3

> <JCHEM_GHOSE_FILTER>
0

> <JCHEM_IUPAC>
phosphate

> <JCHEM_LOGP>
-1.0201038226666665

> <JCHEM_MDDR_LIKE_RULE>
0

> <JCHEM_NUMBER_OF_RINGS>
0

> <JCHEM_PHYSIOLOGICAL_CHARGE>
-2

> <JCHEM_PKA>
6.951626889535468

> <JCHEM_PKA_STRONGEST_ACIDIC>
1.7961261340181292

> <JCHEM_POLAR_SURFACE_AREA>
86.25

> <JCHEM_REFRACTIVITY>
11.2868

> <JCHEM_ROTATABLE_BOND_COUNT>
0

> <JCHEM_RULE_OF_FIVE>
1

> <JCHEM_TRADITIONAL_IUPAC>
phosphate

> <JCHEM_VEBER_RULE>
0

> <DRUGBANK_ID>
DB14523

> <DRUG_GROUPS>
experimental

> <GENERIC_NAME>
Phosphate ion

> <SYNONYMS>
Orthophosphate; Phosphate

$$$$

Comp 2

5  4  0  0  0  0            999 V2000
  0.0000    0.8250    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
 -0.8250    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  0.0000   -0.8250    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
  0.0000    0.0000    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
  0.8250    0.0000    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
4  1  1  0  0  0  0
4  2  2  0  0  0  0
4  3  1  0  0  0  0
4  5  1  0  0  0  0
M  CHG  3   1  -1   3  -1   5  -1
M  END`;

describe('sdf reader', () => {
    it('basic', async () => {
        const parsed =  await parseSdf(SdfString).run();
        if (parsed.isError) {
            throw new Error(parsed.message);
        }
        const compound1 = parsed.result.compounds[0];
        const compound2 = parsed.result.compounds[1];
        const { molFile, dataItems } = compound1;
        const { atoms, bonds } = molFile;

        expect(parsed.result.compounds.length).toBe(2);

        // number of structures
        expect(atoms.count).toBe(5);
        expect(bonds.count).toBe(4);

        expect(compound2.molFile.atoms.count).toBe(5);
        expect(compound2.molFile.bonds.count).toBe(4);

        expect(atoms.x.value(0)).toBeCloseTo(0, 0.001);
        expect(atoms.y.value(0)).toBeCloseTo(0.8250, 0.0001);
        expect(atoms.z.value(0)).toBeCloseTo(0, 0.0001);
        expect(atoms.type_symbol.value(0)).toBe('O');

        expect(bonds.atomIdxA.value(3)).toBe(4);
        expect(bonds.atomIdxB.value(3)).toBe(5);
        expect(bonds.order.value(3)).toBe(1);

        expect(dataItems.dataHeader.value(0)).toBe('DATABASE_ID');
        expect(dataItems.data.value(0)).toBe('DB14523');

        expect(dataItems.dataHeader.value(1)).toBe('DATABASE_NAME');
        expect(dataItems.data.value(1)).toBe('drugbank');

        expect(dataItems.dataHeader.value(31)).toBe('SYNONYMS');
        expect(dataItems.data.value(31)).toBe('Orthophosphate; Phosphate');
    });
});

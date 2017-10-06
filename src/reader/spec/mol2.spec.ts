

import Mol2 from '../mol2/parser'

const Mol2String = `@<TRIPOS>MOLECULE
5816
 26 26 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 O           1.7394   -2.1169   -1.0894 O.3     1  LIG1       -0.3859
      2 O          -2.2941    1.0781   -1.7979 O.3     1  LIG1       -0.5033
@<TRIPOS>BOND
     1     1     5    1
     2     1    21    1`


////////// nothing works until add async and await and promise to parser.mol2 file.
describe('mol2 reader', () => {
    it('basic', async () => {
        const parsed = await Mol2(Mol2String)();

        if (parsed.isError) {
            console.log(parsed)
            return;
        }

        const mol2File = parsed.result;
        const data = mol2File.structures[0];

        const { molecule, atoms, bonds } = data;

        expect(molecule.mol_name).toBe(5816)
        expect(molecule.num_atoms).toBe(26)
        expect(molecule.num_bonds).toBe(26);
        expect(molecule.num_subst).toBe(0);
        expect(molecule.num_feat).toBe(0);
        expect(molecule.num_sets).toBe(0);
        expect(molecule.mol_type).toBe("")
        expect(molecule.charge_type).toBe("");
        expect(molecule.status_bits).toBe("");
        expect(molecule.mol_comment).toBe("");

        expect(atoms.count).toBe(2);
        expect(atoms.atom_id.value(0)).toBe(1);
        expect(atoms.atom_name.value(0)).toBe('o');
        expect(atoms.x.value(0)).toBeCloseTo(1.7394, 0.001);
        expect(atoms.y.value(0)).toBeCloseTo(-2.1169, 0.0001);
        expect(atoms.z.value(0)).toBeCloseTo(-1.0893, 0.0001);
        expect(atoms.atom_type.value(0)).toBe('');
        ///// optionals
        expect(atoms.subst_id.value(0)).toBe(0);
        expect(atoms.subst_name.value(0)).toBe('');
        expect(atoms.charge.value(0)).toBeCloseTo(0.000);
        expect(atoms.status_bits.value(0)).toBe('');

        expect(bonds.count).toBe(2);
        expect(bonds.bond_id.value(0)).toBe(1);
        expect(bonds.origin_atom_id.value(0)).toBe(1);
        expect(bonds.target_atom_id.value(0)).toBe(5);
        expect(bonds.bond_type.value(0)).toBe('1');
        /////// optional
        expect(bonds.status_bits.value(0)).toBe('');

    });
});

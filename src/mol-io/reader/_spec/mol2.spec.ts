
import { parseMol2 } from '../mol2/parser';

const Mol2String = `@<TRIPOS>MOLECULE
5816
 26 26 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 O           1.7394   -2.1169   -1.0894 O.3     1  LIG1       -0.3859
      2 O          -2.2941    1.0781   -1.7979 O.3     1  LIG1       -0.5033
      3 O          -3.6584    0.5842    0.5722 O.3     1  LIG1       -0.5033
      4 N           2.6359    1.0243    0.7030 N.3     1  LIG1       -0.3162
      5 C           1.6787   -1.1447   -0.0373 C.3     1  LIG1        0.0927
      6 C           0.2684   -0.6866    0.1208 C.ar    1  LIG1       -0.0143
      7 C           2.6376    0.0193   -0.3576 C.3     1  LIG1        0.0258
      8 C          -0.3658   -0.0099   -0.9212 C.ar    1  LIG1       -0.0109
      9 C          -0.4164   -0.9343    1.3105 C.ar    1  LIG1       -0.0524
     10 C          -1.6849    0.4191   -0.7732 C.ar    1  LIG1        0.1586
     11 C          -1.7353   -0.5053    1.4585 C.ar    1  LIG1       -0.0162
     12 C          -2.3696    0.1713    0.4166 C.ar    1  LIG1        0.1582
     13 C           3.5645    2.1013    0.3950 C.3     1  LIG1       -0.0157
     14 H           2.0210   -1.6511    0.8741 H       1  LIG1        0.0656
     15 H           2.3808    0.4742   -1.3225 H       1  LIG1        0.0453
     16 H           3.6478   -0.3931   -0.4831 H       1  LIG1        0.0453
     17 H           0.1501    0.1801   -1.8589 H       1  LIG1        0.0659
     18 H           0.0640   -1.4598    2.1315 H       1  LIG1        0.0622
     19 H           2.9013    0.5888    1.5858 H       1  LIG1        0.1217
     20 H          -2.2571   -0.7050    2.3907 H       1  LIG1        0.0655
     21 H           2.6646   -2.4067   -1.1652 H       1  LIG1        0.2103
     22 H           3.2862    2.6124   -0.5325 H       1  LIG1        0.0388
     23 H           4.5925    1.7346    0.3078 H       1  LIG1        0.0388
     24 H           3.5401    2.8441    1.1985 H       1  LIG1        0.0388
     25 H          -3.2008    1.2997   -1.5231 H       1  LIG1        0.2923
     26 H          -3.9690    0.3259    1.4570 H       1  LIG1        0.2923
@<TRIPOS>BOND
     1     1     5    1
     2     1    21    1
     3     2    10    1
     4     2    25    1
     5     3    12    1
     6     3    26    1
     7     4     7    1
     8     4    13    1
     9     4    19    1
    10     5     6    1
    11     5     7    1
    12     5    14    1
    13     6     8   ar
    14     6     9   ar
    15     7    15    1
    16     7    16    1
    17     8    10   ar
    18     8    17    1
    19     9    11   ar
    20     9    18    1
    21    10    12   ar
    22    11    12   ar
    23    11    20    1
    24    13    22    1
    25    13    23    1
    26    13    24    1`;

const Mol2StringMultiBlocks = `@<TRIPOS>MOLECULE
5816
 26 26 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 O           1.7394   -2.1169   -1.0894 O.3     1  LIG1       -0.3859
      2 O          -2.2941    1.0781   -1.7979 O.3     1  LIG1       -0.5033
      3 O          -3.6584    0.5842    0.5722 O.3     1  LIG1       -0.5033
      4 N           2.6359    1.0243    0.7030 N.3     1  LIG1       -0.3162
      5 C           1.6787   -1.1447   -0.0373 C.3     1  LIG1        0.0927
      6 C           0.2684   -0.6866    0.1208 C.ar    1  LIG1       -0.0143
      7 C           2.6376    0.0193   -0.3576 C.3     1  LIG1        0.0258
      8 C          -0.3658   -0.0099   -0.9212 C.ar    1  LIG1       -0.0109
      9 C          -0.4164   -0.9343    1.3105 C.ar    1  LIG1       -0.0524
     10 C          -1.6849    0.4191   -0.7732 C.ar    1  LIG1        0.1586
     11 C          -1.7353   -0.5053    1.4585 C.ar    1  LIG1       -0.0162
     12 C          -2.3696    0.1713    0.4166 C.ar    1  LIG1        0.1582
     13 C           3.5645    2.1013    0.3950 C.3     1  LIG1       -0.0157
     14 H           2.0210   -1.6511    0.8741 H       1  LIG1        0.0656
     15 H           2.3808    0.4742   -1.3225 H       1  LIG1        0.0453
     16 H           3.6478   -0.3931   -0.4831 H       1  LIG1        0.0453
     17 H           0.1501    0.1801   -1.8589 H       1  LIG1        0.0659
     18 H           0.0640   -1.4598    2.1315 H       1  LIG1        0.0622
     19 H           2.9013    0.5888    1.5858 H       1  LIG1        0.1217
     20 H          -2.2571   -0.7050    2.3907 H       1  LIG1        0.0655
     21 H           2.6646   -2.4067   -1.1652 H       1  LIG1        0.2103
     22 H           3.2862    2.6124   -0.5325 H       1  LIG1        0.0388
     23 H           4.5925    1.7346    0.3078 H       1  LIG1        0.0388
     24 H           3.5401    2.8441    1.1985 H       1  LIG1        0.0388
     25 H          -3.2008    1.2997   -1.5231 H       1  LIG1        0.2923
     26 H          -3.9690    0.3259    1.4570 H       1  LIG1        0.2923
@<TRIPOS>BOND
     1     1     5    1
     2     1    21    1
     3     2    10    1
     4     2    25    1
     5     3    12    1
     6     3    26    1
     7     4     7    1
     8     4    13    1
     9     4    19    1
    10     5     6    1
    11     5     7    1
    12     5    14    1
    13     6     8   ar
    14     6     9   ar
    15     7    15    1
    16     7    16    1
    17     8    10   ar
    18     8    17    1
    19     9    11   ar
    20     9    18    1
    21    10    12   ar
    22    11    12   ar
    23    11    20    1
    24    13    22    1
    25    13    23    1
    26    13    24    1
@<TRIPOS>MOLECULE
5816
 26 26 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 O           1.7394   -2.1169   -1.0894 O.3     1  LIG1       -0.3859
      2 O          -2.2941    1.0781   -1.7979 O.3     1  LIG1       -0.5033
      3 O          -3.6584    0.5842    0.5722 O.3     1  LIG1       -0.5033
      4 N           2.6359    1.0243    0.7030 N.3     1  LIG1       -0.3162
      5 C           1.6787   -1.1447   -0.0373 C.3     1  LIG1        0.0927
      6 C           0.2684   -0.6866    0.1208 C.ar    1  LIG1       -0.0143
      7 C           2.6376    0.0193   -0.3576 C.3     1  LIG1        0.0258
      8 C          -0.3658   -0.0099   -0.9212 C.ar    1  LIG1       -0.0109
      9 C          -0.4164   -0.9343    1.3105 C.ar    1  LIG1       -0.0524
     10 C          -1.6849    0.4191   -0.7732 C.ar    1  LIG1        0.1586
     11 C          -1.7353   -0.5053    1.4585 C.ar    1  LIG1       -0.0162
     12 C          -2.3696    0.1713    0.4166 C.ar    1  LIG1        0.1582
     13 C           3.5645    2.1013    0.3950 C.3     1  LIG1       -0.0157
     14 H           2.0210   -1.6511    0.8741 H       1  LIG1        0.0656
     15 H           2.3808    0.4742   -1.3225 H       1  LIG1        0.0453
     16 H           3.6478   -0.3931   -0.4831 H       1  LIG1        0.0453
     17 H           0.1501    0.1801   -1.8589 H       1  LIG1        0.0659
     18 H           0.0640   -1.4598    2.1315 H       1  LIG1        0.0622
     19 H           2.9013    0.5888    1.5858 H       1  LIG1        0.1217
     20 H          -2.2571   -0.7050    2.3907 H       1  LIG1        0.0655
     21 H           2.6646   -2.4067   -1.1652 H       1  LIG1        0.2103
     22 H           3.2862    2.6124   -0.5325 H       1  LIG1        0.0388
     23 H           4.5925    1.7346    0.3078 H       1  LIG1        0.0388
     24 H           3.5401    2.8441    1.1985 H       1  LIG1        0.0388
     25 H          -3.2008    1.2997   -1.5231 H       1  LIG1        0.2923
     26 H          -3.9690    0.3259    1.4570 H       1  LIG1        0.2923
@<TRIPOS>BOND
     1     1     5    1
     2     1    21    1
     3     2    10    1
     4     2    25    1
     5     3    12    1
     6     3    26    1
     7     4     7    1
     8     4    13    1
     9     4    19    1
    10     5     6    1
    11     5     7    1
    12     5    14    1
    13     6     8   ar
    14     6     9   ar
    15     7    15    1
    16     7    16    1
    17     8    10   ar
    18     8    17    1
    19     9    11   ar
    20     9    18    1
    21    10    12   ar
    22    11    12   ar
    23    11    20    1
    24    13    22    1
    25    13    23    1
    26    13    24    1`;

const Mol2StringMinimal = `@<TRIPOS>MOLECULE
5816
 26 26 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 O           1.7394   -2.1169   -1.0894 O.3
      2 O          -2.2941    1.0781   -1.7979 O.3
      3 O          -3.6584    0.5842    0.5722 O.3
      4 N           2.6359    1.0243    0.7030 N.3
      5 C           1.6787   -1.1447   -0.0373 C.3
      6 C           0.2684   -0.6866    0.1208 C.ar
      7 C           2.6376    0.0193   -0.3576 C.3
      8 C          -0.3658   -0.0099   -0.9212 C.ar
      9 C          -0.4164   -0.9343    1.3105 C.ar
     10 C          -1.6849    0.4191   -0.7732 C.ar
     11 C          -1.7353   -0.5053    1.4585 C.ar
     12 C          -2.3696    0.1713    0.4166 C.ar
     13 C           3.5645    2.1013    0.3950 C.3
     14 H           2.0210   -1.6511    0.8741 H
     15 H           2.3808    0.4742   -1.3225 H
     16 H           3.6478   -0.3931   -0.4831 H
     17 H           0.1501    0.1801   -1.8589 H
     18 H           0.0640   -1.4598    2.1315 H
     19 H           2.9013    0.5888    1.5858 H
     20 H          -2.2571   -0.7050    2.3907 H
     21 H           2.6646   -2.4067   -1.1652 H
     22 H           3.2862    2.6124   -0.5325 H
     23 H           4.5925    1.7346    0.3078 H
     24 H           3.5401    2.8441    1.1985 H
     25 H          -3.2008    1.2997   -1.5231 H
     26 H          -3.9690    0.3259    1.4570 H
@<TRIPOS>BOND
     1     1     5    1
     2     1    21    1
     3     2    10    1
     4     2    25    1
     5     3    12    1
     6     3    26    1
     7     4     7    1
     8     4    13    1
     9     4    19    1
    10     5     6    1
    11     5     7    1
    12     5    14    1
    13     6     8   ar
    14     6     9   ar
    15     7    15    1
    16     7    16    1
    17     8    10   ar
    18     8    17    1
    19     9    11   ar
    20     9    18    1
    21    10    12   ar
    22    11    12   ar
    23    11    20    1
    24    13    22    1
    25    13    23    1
    26    13    24    1`;

const Mol2StringCrysin = `@<TRIPOS>MOLECULE
1144204
    12    11     2     0     0
SMALL
USER_CHARGES
****
Generated from the CSD

@<TRIPOS>ATOM
     1 Cl1      0.0925   3.6184   1.9845   Cl        1 RES1  -1.0000
     2 C1      -4.7391   0.3350   0.4215   C.ar      2 RES2   0.0000
     3 C2      -3.4121   0.2604   0.9351   C.ar      2 RES2   0.0000
     4 C3      -2.9169   1.2555   1.7726   C.ar      2 RES2   0.0000
     5 C4      -3.7118   2.3440   2.1099   C.ar      2 RES2   0.0000
     6 C5      -5.0314   2.4052   1.6209   C.ar      2 RES2   0.0000
     7 C6      -5.5372   1.4057   0.7962   C.ar      2 RES2   0.0000
     8 C7      -6.9925   1.4547   0.3334   C.3       2 RES2   0.0000
     9 C8      -7.8537   0.5554   1.1859   C.3       2 RES2   0.0000
    10 N1      -9.3089   0.7134   0.8192   N.3       2 RES2   1.0000
    11 O1      -2.6613  -0.8147   0.5707   O.3       2 RES2   0.0000
    12 O2      -1.6204   1.0919   2.2584   O.3       2 RES2   0.0000
@<TRIPOS>BOND
     1     2     3   ar
     2     3     4   ar
     3     4     5   ar
     4     5     6   ar
     5     6     7   ar
     6     7     2   ar
     7     8     7    1
     8     9     8    1
     9    10     9    1
    10    11     3    1
    11    12     4    1
@<TRIPOS>SUBSTRUCTURE
     1 RES1        1 GROUP             0 ****  ****    0
     2 RES2        2 GROUP             0 ****  ****    0
@<TRIPOS>CRYSIN
   10.5150   11.1300    7.9380   90.0000   90.0000   90.0000    29     5
@<TRIPOS>MOLECULE
   1144204
    12    11     2     0     0
SMALL
USER_CHARGES
****
Generated from the CSD

@<TRIPOS>ATOM
     1 Cl1      0.0925   3.6184   1.9845   Cl        1 RES1  -1.0000
     2 C1      -4.7391   0.3350   0.4215   C.ar      2 RES2   0.0000
     3 C2      -3.4121   0.2604   0.9351   C.ar      2 RES2   0.0000
     4 C3      -2.9169   1.2555   1.7726   C.ar      2 RES2   0.0000
     5 C4      -3.7118   2.3440   2.1099   C.ar      2 RES2   0.0000
     6 C5      -5.0314   2.4052   1.6209   C.ar      2 RES2   0.0000
     7 C6      -5.5372   1.4057   0.7962   C.ar      2 RES2   0.0000
     8 C7      -6.9925   1.4547   0.3334   C.3       2 RES2   0.0000
     9 C8      -7.8537   0.5554   1.1859   C.3       2 RES2   0.0000
    10 N1      -9.3089   0.7134   0.8192   N.3       2 RES2   1.0000
    11 O1      -2.6613  -0.8147   0.5707   O.3       2 RES2   0.0000
    12 O2      -1.6204   1.0919   2.2584   O.3       2 RES2   0.0000
@<TRIPOS>BOND
     1     2     3   ar
     2     3     4   ar
     3     4     5   ar
     4     5     6   ar
     5     6     7   ar
     6     7     2   ar
     7     8     7    1
     8     9     8    1
     9    10     9    1
    10    11     3    1
    11    12     4    1
@<TRIPOS>SUBSTRUCTURE
     1 RES1        1 GROUP             0 ****  ****    0
     2 RES2        2 GROUP             0 ****  ****    0
@<TRIPOS>CRYSIN
   10.5150   11.1300    7.9380   90.0000   90.0000   90.0000    29     5
`;

describe('mol2 reader', () => {
    it('basic', async () => {
        const parsed = await parseMol2(Mol2String, '').run();
        if (parsed.isError) {
            throw new Error(parsed.message);
        }
        const mol2File = parsed.result;

        // number of structures
        expect(mol2File.structures.length).toBe(1);

        const data = mol2File.structures[0];
        const { molecule, atoms, bonds } = data;

        // molecule fields
        expect(molecule.mol_name).toBe('5816');
        expect(molecule.num_atoms).toBe(26);
        expect(molecule.num_bonds).toBe(26);
        expect(molecule.num_subst).toBe(0);
        expect(molecule.num_feat).toBe(0);
        expect(molecule.num_sets).toBe(0);
        expect(molecule.mol_type).toBe('SMALL');
        expect(molecule.charge_type).toBe('GASTEIGER');
        expect(molecule.status_bits).toBe('');
        expect(molecule.mol_comment).toBe('');

        // required atom fields
        expect(atoms.count).toBe(26);
        expect(atoms.atom_id.value(0)).toBe(1);
        expect(atoms.atom_name.value(0)).toBe('O');
        expect(atoms.x.value(0)).toBeCloseTo(1.7394, 0.001);
        expect(atoms.y.value(0)).toBeCloseTo(-2.1169, 0.0001);
        expect(atoms.z.value(0)).toBeCloseTo(-1.0893, 0.0001);
        expect(atoms.atom_type.value(0)).toBe('O.3');

        // optional atom fields
        expect(atoms.subst_id.value(0)).toBe(1);
        expect(atoms.subst_name.value(0)).toBe('LIG1');
        expect(atoms.charge.value(0)).toBeCloseTo(-0.3859);
        expect(atoms.status_bit.value(0)).toBe('');

        // required bond fields
        expect(bonds.count).toBe(26);
        expect(bonds.bond_id.value(0)).toBe(1);
        expect(bonds.origin_atom_id.value(0)).toBe(1);
        expect(bonds.target_atom_id.value(0)).toBe(5);
        expect(bonds.bond_type.value(0)).toBe('1');

        // optional bond fields
        expect(bonds.status_bits.value(0)).toBe('');
    });

    it('multiblocks', async () => {
        const parsed = await parseMol2(Mol2StringMultiBlocks, '').run();
        if (parsed.isError) {
            throw new Error(parsed.message);
        }
        const mol2File = parsed.result;

        // number of structures
        expect(mol2File.structures.length).toBe(2);

        const data = mol2File.structures[1];
        const { molecule, atoms, bonds } = data;

        // molecule fields
        expect(molecule.mol_name).toBe('5816');
        expect(molecule.num_atoms).toBe(26);
        expect(molecule.num_bonds).toBe(26);
        expect(molecule.num_subst).toBe(0);
        expect(molecule.num_feat).toBe(0);
        expect(molecule.num_sets).toBe(0);
        expect(molecule.mol_type).toBe('SMALL');
        expect(molecule.charge_type).toBe('GASTEIGER');
        expect(molecule.status_bits).toBe('');
        expect(molecule.mol_comment).toBe('');

        // required atom fields
        expect(atoms.count).toBe(26);
        expect(atoms.atom_id.value(0)).toBe(1);
        expect(atoms.atom_name.value(0)).toBe('O');
        expect(atoms.x.value(0)).toBeCloseTo(1.7394, 0.001);
        expect(atoms.y.value(0)).toBeCloseTo(-2.1169, 0.0001);
        expect(atoms.z.value(0)).toBeCloseTo(-1.0893, 0.0001);
        expect(atoms.atom_type.value(0)).toBe('O.3');

        // optional atom fields
        expect(atoms.subst_id.value(0)).toBe(1);
        expect(atoms.subst_name.value(0)).toBe('LIG1');
        expect(atoms.charge.value(0)).toBeCloseTo(-0.3859);
        expect(atoms.status_bit.value(0)).toBe('');

        // required bond fields
        expect(bonds.count).toBe(26);
        expect(bonds.bond_id.value(0)).toBe(1);
        expect(bonds.origin_atom_id.value(0)).toBe(1);
        expect(bonds.target_atom_id.value(0)).toBe(5);
        expect(bonds.bond_type.value(0)).toBe('1');

        // optional bond fields
        expect(bonds.status_bits.value(0)).toBe('');
    });

    it('minimal', async () => {
        const parsed = await parseMol2(Mol2StringMinimal, '').run();
        if (parsed.isError) {
            throw new Error(parsed.message);
        }
        const mol2File = parsed.result;

        // number of structures
        expect(mol2File.structures.length).toBe(1);

        const data = mol2File.structures[0];
        const { molecule, atoms, bonds } = data;

        // molecule fields
        expect(molecule.mol_name).toBe('5816');
        expect(molecule.num_atoms).toBe(26);
        expect(molecule.num_bonds).toBe(26);
        expect(molecule.num_subst).toBe(0);
        expect(molecule.num_feat).toBe(0);
        expect(molecule.num_sets).toBe(0);
        expect(molecule.mol_type).toBe('SMALL');
        expect(molecule.charge_type).toBe('GASTEIGER');
        expect(molecule.status_bits).toBe('');
        expect(molecule.mol_comment).toBe('');

        // required atom fields
        expect(atoms.count).toBe(26);
        expect(atoms.atom_id.value(0)).toBe(1);
        expect(atoms.atom_name.value(0)).toBe('O');
        expect(atoms.x.value(0)).toBeCloseTo(1.7394, 0.001);
        expect(atoms.y.value(0)).toBeCloseTo(-2.1169, 0.0001);
        expect(atoms.z.value(0)).toBeCloseTo(-1.0893, 0.0001);
        expect(atoms.atom_type.value(0)).toBe('O.3');

        // optional atom fields
        expect(atoms.subst_id.value(0)).toBe(0);
        expect(atoms.subst_name.value(0)).toBe('');
        expect(atoms.charge.value(0)).toBeCloseTo(0);
        expect(atoms.status_bit.value(0)).toBe('');

        // required bond fields
        expect(bonds.count).toBe(26);
        expect(bonds.bond_id.value(0)).toBe(1);
        expect(bonds.origin_atom_id.value(0)).toBe(1);
        expect(bonds.target_atom_id.value(0)).toBe(5);
        expect(bonds.bond_type.value(0)).toBe('1');

        // optional bond fields
        expect(bonds.status_bits.value(0)).toBe('');
    });

    it('crysin', async () => {
        const parsed = await parseMol2(Mol2StringCrysin, '').run();
        if (parsed.isError) {
            throw new Error(parsed.message);
        }
        const mol2File = parsed.result;

        // number of structures
        expect(mol2File.structures.length).toBe(2);

        // crysin fields
        for (const data of mol2File.structures) {
            expect(data.crysin).toEqual({
                a: 10.5150,
                b: 11.1300,
                c: 7.9380,
                alpha: 90.0,
                beta: 90.0,
                gamma: 90.0,
                spaceGroup: 29,
                setting: 5
            });
        }
    });
});

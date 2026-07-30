/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Himanshu Raj <himanshuraj6771@gmail.com>
 */

import { parseLammpsData } from '../lammps/data/parser';
import { trajectoryFromLammpsData } from '../../../mol-model-formats/structure/lammps-data';

async function parse(data: string) {
    const result = await parseLammpsData(data).run();
    if (result.isError) throw new Error(result.message);
    return result.result;
}

describe('lammps data parser', () => {

    describe('atom_style variants', () => {
        it('parses atomic style (atomId, atomType, x, y, z)', async () => {
            const data = `LAMMPS data file

2 atoms
1 atom types

Atoms # atomic

1 1 0.0 0.0 0.0
2 1 1.5 0.0 0.0
`;
            const file = await parse(data);
            expect(file.atoms?.count).toBe(2);
            expect(file.atoms?.atomId.value(0)).toBe(1);
            expect(file.atoms?.atomType.value(1)).toBe(1);
        });

        it('parses full style (atomId, moleculeId, atomType, charge, x, y, z)', async () => {
            const data = `LAMMPS data file

1 atoms
1 atom types

Atoms # full

1 1 1 -0.5 0.0 0.0 0.0
`;
            const file = await parse(data);
            expect(file.atoms?.count).toBe(1);
            expect(file.atoms?.moleculeId?.value(0)).toBe(1);
            expect(file.atoms?.charge?.value(0)).toBe(-0.5);
        });

        it('parses bond style (atomId, moleculeId, atomType, x, y, z)', async () => {
            const data = `LAMMPS data file

1 atoms
1 atom types

Atoms # bond

1 1 1 0.0 0.0 0.0
`;
            const file = await parse(data);
            expect(file.atoms?.count).toBe(1);
            expect(file.atoms?.moleculeId?.value(0)).toBe(1);
        });

        it('parses molecular style (atomId, moleculeId, atomType, x, y, z)', async () => {
            const data = `LAMMPS data file

1 atoms
1 atom types

Atoms # molecular

1 2 1 0.0 0.0 0.0
`;
            const file = await parse(data);
            expect(file.atoms?.count).toBe(1);
            expect(file.atoms?.moleculeId?.value(0)).toBe(2);
        });

        it('defaults to full style when no style comment is present', async () => {
            const data = `LAMMPS data file

1 atoms
1 atom types

Atoms

1 1 1 0.0 0.0 0.0 0.0
`;
            const file = await parse(data);
            expect(file.atoms?.count).toBe(1);
        });
    });

    describe('Masses section', () => {
        it('parses element symbol from a trailing "# <symbol>" comment', async () => {
            const data = `LAMMPS data file

2 atoms
2 atom types

Masses

1 12.011  # C
2 1.008   # H

Atoms # atomic

1 1 0.0 0.0 0.0
2 2 1.0 0.0 0.0
`;
            const file = await parse(data);
            expect(file.masses?.count).toBe(2);
            expect(file.masses?.symbol.value(0)).toBe('C');
            expect(file.masses?.symbol.value(1)).toBe('H');
        });

        it('skips a whole-line "#" comment before the first data row', async () => {
            const data = `LAMMPS data file

3 atoms
3 atom types

Masses

# masses for each atom type
1 12.011 # C
2 1.008  # H
3 15.999 # O

Atoms # atomic

1 1 0.0 0.0 0.0
2 2 1.0 0.0 0.0
3 3 2.0 0.0 0.0
`;
            const file = await parse(data);
            expect(file.masses?.count).toBe(3);
            expect(file.masses?.atomType.value(0)).toBe(1);
            expect(file.masses?.mass.value(0)).toBe(12.011);
            expect(file.masses?.symbol.value(2)).toBe('O');
        });

        it('skips multiple consecutive comment/blank lines before data rows', async () => {
            const data = `LAMMPS data file

1 atoms
1 atom types

Masses

# first comment line
# second comment line

1 12.011 # C

Atoms # atomic

1 1 0.0 0.0 0.0
`;
            const file = await parse(data);
            expect(file.masses?.count).toBe(1);
            expect(file.masses?.atomType.value(0)).toBe(1);
            expect(file.masses?.symbol.value(0)).toBe('C');
        });

        it('handles a mass line with no trailing comment (empty symbol token, not a crash)', async () => {
            const data = `LAMMPS data file

1 atoms
1 atom types

Masses

1 12.011

Atoms # atomic

1 1 0.0 0.0 0.0
`;
            const file = await parse(data);
            expect(file.masses?.count).toBe(1);
            expect(file.masses?.mass.value(0)).toBe(12.011);
            expect(file.masses?.symbol.value(0)).toBe('');
        });

        it('leaves masses undefined when the section is absent', async () => {
            const data = `LAMMPS data file

1 atoms
1 atom types

Atoms # atomic

1 1 0.0 0.0 0.0
`;
            const file = await parse(data);
            expect(file.masses).toBeUndefined();
        });

        it('infers element symbol from mass when comment is absent', async () => {
            const data = `LAMMPS data file

1 atoms
1 atom types

Masses

1 12.011

Atoms # atomic

1 1 0 0 0
`;

            const file = await parse(data);

            expect(file.masses?.symbol.value(0)).toBe('');
        });

        it('returns empty symbol when mass cannot be matched', async () => {
            const data = `LAMMPS data file

1 atoms
1 atom types

Masses

1 999.999

Atoms # atomic

1 1 0 0 0
`;

            const file = await parse(data);

            expect(file.masses?.symbol.value(0)).toBe('');
        });

        it('infers element symbol from mass', async () => {
            const data = `LAMMPS data file

1 atoms
1 atom types

Masses

1 12 # Hello Mol* Community !!!

Atoms # atomic

1 1 0 0 0
`;
            const file = await parse(data);

            const traj = await trajectoryFromLammpsData(file).run();
            const model = traj.representative;

            expect(model.atomicHierarchy.atoms.type_symbol.value(0)).toBe('C');
        });
    });

    describe('box and header counts', () => {
        it('parses atom/atom-type/bond counts and the box from the header', async () => {
            const data = `LAMMPS data file

4 atoms
2 atom types
1 bonds

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 10.0 zlo zhi

Atoms # atomic

1 1 0.0 0.0 0.0
2 1 1.0 0.0 0.0
3 2 0.0 1.0 0.0
4 2 1.0 1.0 0.0

Bonds

1 1 1 2
`;
            const file = await parse(data);
            expect(file.atoms?.count).toBe(4);
            expect(file.bonds?.count).toBe(1);
            expect(file.box).toBeDefined();
            expect(file.box?.lower).toEqual([0.0, 0.0, 0.0]);
            expect(file.box?.length).toEqual([10.0, 10.0, 10.0]);
        });
    });
});
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Aniruddha Adak <aniruddhaadak80@gmail.com>
 */

import { Column, Table } from '../../../../mol-data/db';
import { RuntimeContext } from '../../../../mol-task';
import { createModels } from '../../../../mol-model-formats/structure/basic/parser';
import { BasicSchema, createBasic } from '../../../../mol-model-formats/structure/basic/schema';
import { EntityBuilder } from '../../../../mol-model-formats/structure/common/entity';
import { MoleculeType } from '../../../../mol-model/structure/model/types';
import { Structure } from '../../../../mol-model/structure/structure';
import { Unit } from '../../../../mol-model/structure/structure/unit';
import {
    createNucleicIndices, getNucleotideBaseType,
    setPurinIndices, hasPurinIndices, setPyrimidineIndices, hasPyrimidineIndices,
    setSugarIndices, hasSugarIndices
} from '../nucleotide';

type Atom = { name: string, element: string, xyz: [number, number, number] };

function offset(atoms: Atom[], dx: number) {
    return atoms.map(a => ({ ...a, xyz: [a.xyz[0] + dx, a.xyz[1], a.xyz[2]] as [number, number, number] }));
}

/** phosphate and sugar atoms, identical for all of the nucleotides below */
const SugarAtoms: Atom[] = [
    { name: 'P', element: 'P', xyz: [11.617, 3.631, -0.459] },
    { name: 'OP1', element: 'O', xyz: [11.814, 5.331, -0.636] },
    { name: 'OP2', element: 'O', xyz: [11.779, 3.143, 1.183] },
    { name: 'OP3', element: 'O', xyz: [13.097, 2.756, -0.413] },
    { name: "O5'", element: 'O', xyz: [10.121, 3.216, -1.198] },
    { name: "C5'", element: 'C', xyz: [9.058, 2.676, -0.408] },
    { name: "C4'", element: 'C', xyz: [7.924, 2.396, -1.411] },
    { name: "O4'", element: 'O', xyz: [6.695, 2.912, -0.970] },
    { name: "C3'", element: 'C', xyz: [7.754, 0.916, -1.767] },
    { name: "O3'", element: 'O', xyz: [8.076, 1.303, -3.105] },
    { name: "C2'", element: 'C', xyz: [6.253, 0.801, -1.985] },
    { name: "C1'", element: 'C', xyz: [5.689, 1.952, -1.155] },
];

/** purine ring of 0DA, atom names and model coordinates from the wwPDB CCD entry for 0DA */
const ZeroDARingAtoms: Atom[] = [
    { name: 'C8A', element: 'C', xyz: [5.903, 1.463, 1.322] },
    { name: 'N9A', element: 'N', xyz: [5.175, 1.605, 0.177] },
    { name: 'C4A', element: 'C', xyz: [3.859, 1.627, 0.567] },
    { name: 'C5A', element: 'C', xyz: [3.875, 1.565, 1.933] },
    { name: 'N7A', element: 'N', xyz: [5.187, 1.468, 2.416] },
    { name: 'N3A', element: 'N', xyz: [2.766, 1.690, -0.240] },
    { name: 'C2A', element: 'C', xyz: [1.634, 1.668, 0.464] },
    { name: 'N1A', element: 'N', xyz: [1.483, 1.626, 1.788] },
    { name: 'C6A', element: 'C', xyz: [2.598, 1.586, 2.535] },
    { name: 'N6A', element: 'N', xyz: [2.427, 1.540, 3.844] },
];

/** the same purine ring with the conventional, unsuffixed names */
const PurinePlainAtoms: Atom[] = ZeroDARingAtoms.map(a => ({ ...a, name: a.name.slice(0, -1) }));

/** pyrimidine ring, laid out as a regular hexagon */
const PyrimidineRingAtoms: Atom[] = [
    { name: 'N1', element: 'N', xyz: [31.40, 0.00, 0.00] },
    { name: 'C2', element: 'C', xyz: [30.70, 1.21, 0.00] },
    { name: 'N3', element: 'N', xyz: [29.30, 1.21, 0.00] },
    { name: 'C4', element: 'C', xyz: [28.60, 0.00, 0.00] },
    { name: 'C5', element: 'C', xyz: [29.30, -1.21, 0.00] },
    { name: 'C6', element: 'C', xyz: [30.70, -1.21, 0.00] },
];

/**
 * Build a single-residue nucleotide from explicit chem_comp type, so that the residue
 * molecule type is derived the same way as for a deposited structure. 0DA/0DC/0DG/0DT
 * are typed `l-dna linking` in the wwPDB CCD.
 */
async function buildNucleotideUnit(compId: string, compType: string, atoms: Atom[]) {
    const count = atoms.length;
    const atomNames = atoms.map(a => a.name);

    const auth_asym_id = Column.ofConst('A', count, Column.Schema.str);
    const auth_seq_id = Column.ofConst(1, count, Column.Schema.int);
    const auth_comp_id = Column.ofConst(compId, count, Column.Schema.str);

    const atom_site = Table.ofPartialColumns(BasicSchema.atom_site, {
        auth_asym_id,
        auth_atom_id: Column.ofStringArray(atomNames),
        auth_comp_id,
        auth_seq_id,
        Cartn_x: Column.ofFloatArray(Float32Array.from(atoms.map(a => a.xyz[0]))),
        Cartn_y: Column.ofFloatArray(Float32Array.from(atoms.map(a => a.xyz[1]))),
        Cartn_z: Column.ofFloatArray(Float32Array.from(atoms.map(a => a.xyz[2]))),
        id: Column.ofIntArray(atoms.map((_, i) => i + 1)),
        label_asym_id: auth_asym_id,
        label_atom_id: Column.ofStringArray(atomNames),
        label_comp_id: auth_comp_id,
        label_seq_id: auth_seq_id,
        label_entity_id: Column.ofConst('1', count, Column.Schema.str),
        occupancy: Column.ofConst(1, count, Column.Schema.float),
        type_symbol: Column.ofStringArray(atoms.map(a => a.element)),
        pdbx_PDB_model_num: Column.ofConst(1, count, Column.Schema.int),
    }, count);

    const entityBuilder = new EntityBuilder();
    entityBuilder.getEntityId(compId, MoleculeType.DNA, 'A');

    const basic = createBasic({
        entity: entityBuilder.getEntityTable(),
        chem_comp: Table.ofPartialColumns(BasicSchema.chem_comp, {
            id: Column.ofStringArray([compId]),
            name: Column.ofStringArray([compId]),
            type: Column.ofStringAliasArray([compType]),
            mon_nstd_flag: Column.ofStringAliasArray(['y']),
        }, 1),
        atom_site
    });

    const trajectory = await createModels(basic, { kind: 'test', name: 'synthetic-nucleotide', data: undefined }, RuntimeContext.Synchronous);
    const structure = Structure.ofModel(trajectory.representative);
    const atomicUnits = structure.units.filter(Unit.isAtomic);
    if (atomicUnits.length === 0) throw new Error('expected at least one atomic unit');
    return atomicUnits[0];
}

function atomName(unit: Unit.Atomic, elementIndex: number) {
    return unit.model.atomicHierarchy.atoms.label_atom_id.value(elementIndex);
}

describe('nucleotide ring atom lookup', () => {
    it('resolves the A-suffixed purine ring atoms of 0DA', async () => {
        const unit = await buildNucleotideUnit('0DA', 'l-dna linking', [...SugarAtoms, ...ZeroDARingAtoms]);
        const residueIndex = unit.residueIndex[0];

        const { isPurine, isPyrimidine } = getNucleotideBaseType(unit, residueIndex);
        expect(isPurine).toBe(true);
        expect(isPyrimidine).toBe(false);

        const idx = setPurinIndices(createNucleicIndices(), unit, residueIndex);
        expect(hasPurinIndices(idx)).toBe(true);
        expect(atomName(unit, idx.N1)).toBe('N1A');
        expect(atomName(unit, idx.C4)).toBe('C4A');
        expect(atomName(unit, idx.C5)).toBe('C5A');
        expect(atomName(unit, idx.N7)).toBe('N7A');
        expect(atomName(unit, idx.C8)).toBe('C8A');
        expect(atomName(unit, idx.N9)).toBe('N9A');

        expect(hasSugarIndices(setSugarIndices(createNucleicIndices(), unit, residueIndex))).toBe(true);
    });

    it('still resolves the plain purine ring atoms of DA', async () => {
        const unit = await buildNucleotideUnit('DA', 'dna linking', [...offset(SugarAtoms, 20), ...offset(PurinePlainAtoms, 20)]);
        const residueIndex = unit.residueIndex[0];

        const { isPurine } = getNucleotideBaseType(unit, residueIndex);
        expect(isPurine).toBe(true);

        const idx = setPurinIndices(createNucleicIndices(), unit, residueIndex);
        expect(hasPurinIndices(idx)).toBe(true);
        expect(atomName(unit, idx.N1)).toBe('N1');
        expect(atomName(unit, idx.C4)).toBe('C4');
        expect(atomName(unit, idx.N9)).toBe('N9');
    });

    it('still resolves the plain pyrimidine ring atoms of 0DC', async () => {
        const unit = await buildNucleotideUnit('0DC', 'l-dna linking', [...offset(SugarAtoms, 25), ...PyrimidineRingAtoms]);
        const residueIndex = unit.residueIndex[0];

        const { isPurine, isPyrimidine } = getNucleotideBaseType(unit, residueIndex);
        expect(isPurine).toBe(false);
        expect(isPyrimidine).toBe(true);

        const idx = setPyrimidineIndices(createNucleicIndices(), unit, residueIndex);
        expect(hasPyrimidineIndices(idx)).toBe(true);
        expect(atomName(unit, idx.N1)).toBe('N1');
        expect(atomName(unit, idx.C6)).toBe('C6');
    });
});

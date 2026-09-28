/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { Column } from '../../../mol-data/db';
import { SortedArray } from '../../../mol-data/int';
import { CIF } from '../../../mol-io/reader/cif';
import { parseMol } from '../../../mol-io/reader/mol/parser';
import { parsePDB } from '../../../mol-io/reader/pdb/parser';
import { IntAdjacencyGraph } from '../../../mol-math/graph';
import { trajectoryFromMmCIF } from '../../../mol-model-formats/structure/mmcif';
import { trajectoryFromMol } from '../../../mol-model-formats/structure/mol';
import { IndexPairBonds } from '../../../mol-model-formats/structure/property/bonds/index-pair';
import { trajectoryFromPDB } from '../../../mol-model-formats/structure/pdb';
import { Structure, Unit } from '../../../mol-model/structure';
import { ElementIndex } from '../../../mol-model/structure/model';
import { BondType } from '../../../mol-model/structure/model/types';
import { computeIntraUnitBonds } from '../../../mol-model/structure/structure/unit/bonds/intra-compute';
import { IntraUnitBonds } from '../../../mol-model/structure/structure/unit/bonds/data';
import { ElementSetIntraBondCache } from '../../../mol-model/structure/structure/unit/bonds/element-set-intra-bond-cache';
import { BondOrderProviderName, registerBondOrderProviders, unregisterBondOrderProviders } from '../provider';
import { BondProvider, BondProviderRegistry } from '../../../mol-model/structure/structure/unit/bonds/bond-provider';

const PdbWithConectLigand = [
    'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  ',
    'ATOM      2  CA  ALA A   1       1.460   0.000   0.000  1.00  0.00           C  ',
    'ATOM      3  C   ALA A   1       1.960   1.420   0.000  1.00  0.00           C  ',
    'ATOM      4  O   ALA A   1       1.220   2.350   0.000  1.00  0.00           O  ',
    'HETATM    5  C1  LIG A   2       5.000   0.000   0.000  1.00  0.00           C  ',
    'HETATM    6  C2  LIG A   2       6.400   0.000   0.000  1.00  0.00           C  ',
    'HETATM    7  C3  LIG A   2       7.100   1.210   0.000  1.00  0.00           C  ',
    'CONECT    5    6                                                                 ',
    'CONECT    6    5    7                                                            ',
    'CONECT    7    6                                                                 ',
    'END                                                                             ',
].join('\n');

const PdbWithKnownHis = [
    'ATOM      1  CG  HIS A   1       0.000   0.000   0.000  1.00  0.00           C  ',
    'ATOM      2  CD2 HIS A   1       1.400   0.000   0.000  1.00  0.00           C  ',
    'END                                                                             ',
].join('\n');

const MolWithSingleCC = `test
  mol

  2  1  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
`;

const MmcifWithChemCompBond = `data_TST
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_atom_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
HETATM 1 C C1 . LIG A 1 . ? 0.000 0.000 0.000 1.00 0.00 1 LIG C1 A 1
HETATM 2 C C2 . LIG A 1 . ? 1.400 0.000 0.000 1.00 0.00 1 LIG C2 A 1
#
loop_
_chem_comp.id
_chem_comp.type
LIG non-polymer
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_ordinal
LIG C1 C2 sing 1
#
`;

async function structureFromPdb(pdbText: string): Promise<Structure> {
    const parsed = await parsePDB(pdbText, 'TST').run();
    if (parsed.isError) throw new Error(parsed.message);
    const trajectory = await trajectoryFromPDB(parsed.result).run();
    return Structure.ofModel(trajectory.representative);
}

function atomicUnits(structure: Structure): Unit.Atomic[] {
    return structure.units.filter(Unit.isAtomic);
}

function unitWithComp(structure: Structure, compId: string): Unit.Atomic {
    const { label_comp_id } = structure.model.atomicHierarchy.atoms;
    for (const unit of atomicUnits(structure)) {
        for (let i = 0, il = unit.elements.length; i < il; i++) {
            if (label_comp_id.value(unit.elements[i]) === compId) return unit;
        }
    }
    throw new Error(`no unit with ${compId}`);
}

function countCompCCOrder(unit: Unit.Atomic, compId: string, wantOrder: number) {
    const { offset, b, edgeProps: { order } } = unit.bonds;
    const { elements } = unit;
    const { type_symbol, label_comp_id } = unit.model.atomicHierarchy.atoms;
    let n = 0;
    for (let u = 0, ul = elements.length; u < ul; u++) {
        if (label_comp_id.value(elements[u]) !== compId) continue;
        for (let i = offset[u], il = offset[u + 1]; i < il; i++) {
            if (u >= b[i]) continue;
            if (order[i] === wantOrder && type_symbol.value(elements[u]) === 'C' && type_symbol.value(elements[b[i]]) === 'C') n++;
        }
    }
    return n;
}

function collectCompElements(unit: Unit.Atomic, compId: string): ElementIndex[] {
    const { label_comp_id } = unit.model.atomicHierarchy.atoms;
    const els: ElementIndex[] = [];
    for (let i = 0, il = unit.elements.length; i < il; i++) {
        if (label_comp_id.value(unit.elements[i]) === compId) els.push(unit.elements[i]);
    }
    return els;
}

describe('bond provider registry', () => {
    it('perceives CONECT ligand C-C lazily via unit.bonds and getChild subsets', async () => {
        const structure = await structureFromPdb(PdbWithConectLigand);
        for (const unit of atomicUnits(structure)) {
            expect(unit.props.bonds).toBeUndefined();
        }

        const registered = registerBondOrderProviders(structure, 'model');
        expect(registered.length).toBe(1);
        for (const unit of atomicUnits(structure)) {
            expect(unit.props.bonds).toBeUndefined();
        }
        expect(IndexPairBonds.Provider.get(structure.model)).toBeUndefined();

        const ligUnit = unitWithComp(structure, 'LIG');
        expect(countCompCCOrder(ligUnit, 'LIG', 2)).toBeGreaterThan(0);
        const perceived = ligUnit.bonds;
        expect(ligUnit.bonds).toBe(perceived);
        expect(ligUnit.props.bonds).toBeUndefined();
        expect(ElementSetIntraBondCache.get(ligUnit.model).get(ligUnit.elements)).toBeUndefined();
        expect(ligUnit.props.rings).toBeDefined();
        expect((ligUnit.props.rings as any)._aromaticRings).toBeUndefined();
        expect((ligUnit.props.rings as any)._index).toBeUndefined();

        const child = ligUnit.getChild(SortedArray.ofUnsortedArray(collectCompElements(ligUnit, 'LIG').slice(0, 2))) as Unit.Atomic;
        expect(countCompCCOrder(child, 'LIG', 2)).toBeGreaterThan(0);

        const alaUnit = unitWithComp(structure, 'ALA');
        expect(countCompCCOrder(alaUnit, 'ALA', 2)).toBeGreaterThan(0);
    });

    it('does not customize models that already have IndexPairBonds', async () => {
        const structure = await structureFromPdb(PdbWithConectLigand);
        const model = structure.model;
        const existing = IndexPairBonds.fromData({
            pairs: { indexA: Column.ofIntArray([]), indexB: Column.ofIntArray([]) },
            count: model.atomicHierarchy.atoms._rowCount,
        });
        IndexPairBonds.Provider.set(model, existing);

        const registered = registerBondOrderProviders(structure, 'model');
        expect(registered.length).toBe(1);
        expect(IndexPairBonds.Provider.get(model)).toBe(existing);
        void (structure.units[0] as Unit.Atomic).bonds;
    });

    it('does not perceive components covered by IntraBondOrderTable', async () => {
        const structure = await structureFromPdb(PdbWithKnownHis);
        const registered = registerBondOrderProviders(structure, 'model');
        const unit = unitWithComp(structure, 'HIS');

        expect(countCompCCOrder(unit, 'HIS', 2)).toBe(1);
        expect(registered.length).toBe(1);
    });

    it('skips mol/SDF graphs and leaves file orders unchanged', async () => {
        const parsed = await parseMol(MolWithSingleCC).run();
        if (parsed.isError) throw new Error(parsed.message);
        const trajectory = await trajectoryFromMol(parsed.result).run();
        const structure = Structure.ofModel(trajectory.representative);
        const existing = IndexPairBonds.Provider.get(structure.model);
        expect(existing).toBeDefined();

        const unit = structure.units[0] as Unit.Atomic;
        const before = Array.from(unit.bonds.edgeProps.order);

        const registered = registerBondOrderProviders(structure, 'force');
        expect(registered.length).toBe(1);
        expect(IndexPairBonds.Provider.get(structure.model)).toBe(existing);
        expect(Array.from(unit.bonds.edgeProps.order)).toEqual(before);
        expect(countCompCCOrder(unit, 'MOL', 1)).toBeGreaterThan(0);
        expect(countCompCCOrder(unit, 'MOL', 2)).toBe(0);
    });

    it('leaves chem_comp_bond dictionary orders alone in model mode', async () => {
        const parsed = await CIF.parseText(MmcifWithChemCompBond).run();
        if (parsed.isError) throw new Error(parsed.message);
        const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
        const structure = Structure.ofModel(trajectory.representative);

        const registered = registerBondOrderProviders(structure, 'model');
        expect(registered.length).toBe(1);

        const ligUnit = unitWithComp(structure, 'LIG');
        expect(countCompCCOrder(ligUnit, 'LIG', 1)).toBeGreaterThan(0);
        expect(countCompCCOrder(ligUnit, 'LIG', 2)).toBe(0);

        const { offset, b, edgeProps: { flags } } = ligUnit.bonds;
        let sawComputed = false;
        for (let u = 0; u < ligUnit.elements.length; u++) {
            for (let i = offset[u]; i < offset[u + 1]; i++) {
                if (u >= b[i]) continue;
                if (BondType.is(flags[i], BondType.Flag.Computed)) sawComputed = true;
            }
        }
        expect(sawComputed).toBe(false);
    });

    it('looks up providers by name and explicitly selects the active provider', async () => {
        const structure = await structureFromPdb(PdbWithConectLigand);
        const unit = unitWithComp(structure, 'LIG');
        const calls: string[] = [];

        const double: BondProvider = {
            name: 'double',
            getBonds: providerUnit => {
                calls.push('double');
                return withOrder(computeIntraUnitBonds(providerUnit), 2);
            }
        };
        const triple: BondProvider = {
            name: 'triple',
            getBonds: providerUnit => {
                calls.push('triple');
                return withOrder(computeIntraUnitBonds(providerUnit), 3);
            }
        };

        const defaultBonds = unit.bonds;
        expect(defaultBonds.edgeProps.order[0]).toBe(1);
        const registry = BondProviderRegistry.get(structure.model);
        registry.add(double);
        registry.add({ ...triple, name: double.name });
        registry.add(triple);
        expect(registry.get('double')).toBe(double);
        expect(registry.get('triple')).toBe(triple);
        expect(unit.bonds.edgeProps.order[0]).toBe(1);
        expect(calls).toEqual([]);

        registry.select('double');
        expect(unit.bonds.edgeProps.order[0]).toBe(2);
        expect(calls).toEqual(['double']);

        registry.select('triple');
        expect(unit.bonds.edgeProps.order[0]).toBe(3);
        expect(calls).toEqual(['double', 'triple']);

        registry.remove(triple);
        expect(unit.bonds).toBe(defaultBonds);
        expect(computeIntraUnitBonds(unit).edgeProps.order[0]).toBe(1);
        expect(calls).toEqual(['double', 'triple']);
    });

    it('returns to default computation when the active provider detaches', async () => {
        const structure = await structureFromPdb(PdbWithConectLigand);
        const unit = unitWithComp(structure, 'LIG');
        const previous: BondProvider = {
            name: 'previous',
            getBonds: providerUnit => withOrder(computeIntraUnitBonds(providerUnit), 3)
        };

        const registry = BondProviderRegistry.get(structure.model);
        registry.add(previous);
        registry.select(previous.name);

        const registered = registerBondOrderProviders(structure, 'model');
        expect(registry.active).toBe(registered[0].provider);

        unregisterBondOrderProviders(registered);
        expect(registry.active).toBeUndefined();
        expect(registry.get(previous.name)).toBe(previous);
        expect(unit.bonds.edgeProps.order[0]).toBe(1);
    });

    it('keeps duplicate provider registration idempotent', async () => {
        const structure = await structureFromPdb(PdbWithConectLigand);
        const registry = BondProviderRegistry.get(structure.model);
        const first = registerBondOrderProviders(structure, 'model');
        const duplicate = registerBondOrderProviders(structure, 'force');

        expect(registry.active).toBe(first[0].provider);
        expect(registry.get(BondOrderProviderName)).toBe(first[0].provider);

        unregisterBondOrderProviders(duplicate);
        expect(registry.active).toBe(first[0].provider);

        unregisterBondOrderProviders(first);
        expect(registry.active).toBeUndefined();
    });

    it('falls back without registering on multi-model structures', async () => {
        const structureA = await structureFromPdb(PdbWithConectLigand);
        const structureB = await structureFromPdb(PdbWithKnownHis);
        const source = structureB.units[0];
        const unitB = source.getCopy(
            structureA.units.length,
            source.invariantId + structureA.units.length,
            source.chainGroupId + structureA.units.length
        );
        const structure = Structure.create([...structureA.units, unitB]);

        expect(structure.models.length).toBe(2);
        expect(registerBondOrderProviders(structure, 'model')).toEqual([]);

        const unit = unitWithComp(structureA, 'LIG');
        expect(countCompCCOrder(unit, 'LIG', 1)).toBeGreaterThan(0);
        expect(countCompCCOrder(unit, 'LIG', 2)).toBe(0);
    });
});

function withOrder(bonds: IntraUnitBonds, value: number): IntraUnitBonds {
    const order = new Int8Array(bonds.edgeProps.order.length);
    order.fill(value);
    return IntAdjacencyGraph.create(bonds.offset, bonds.a, bonds.b, bonds.edgeCount, {
        ...bonds.edgeProps,
        order,
    }, bonds.props) as IntraUnitBonds;
}

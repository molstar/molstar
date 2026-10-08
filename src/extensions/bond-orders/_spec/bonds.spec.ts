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
import { trajectoryFromMmCIF } from '../../../mol-model-formats/structure/mmcif';
import { trajectoryFromMol } from '../../../mol-model-formats/structure/mol';
import { IndexPairBonds } from '../../../mol-model-formats/structure/property/bonds/index-pair';
import { trajectoryFromPDB } from '../../../mol-model-formats/structure/pdb';
import { Structure, Unit } from '../../../mol-model/structure';
import { ElementIndex, Model } from '../../../mol-model/structure/model';
import { BondType } from '../../../mol-model/structure/model/types';
import { ElementSetIntraBondCache } from '../../../mol-model/structure/structure/unit/bonds/element-set-intra-bond-cache';
import { BondProvider, BondProviderRegistry, ModelBondProvider } from '../../../mol-model/structure/structure/unit/bonds/bond-provider';
import { createModelBondProviderProperty } from '../../../mol-model-props/common/model-bond-provider';
import { AssetManager } from '../../../mol-util/assets';
import { SyncRuntimeContext } from '../../../mol-task/execution/synchronous';
import { BondOrderProvider, BondOrderProviderName } from '../provider';
import { BondOrdersMode } from '../perceiver';

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

async function modelFromPdb(pdbText: string): Promise<Model> {
    const parsed = await parsePDB(pdbText, 'TST').run();
    if (parsed.isError) throw new Error(parsed.message);
    const trajectory = await trajectoryFromPDB(parsed.result).run();
    return trajectory.representative;
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

const bondProviderRegistry = new BondProviderRegistry();
bondProviderRegistry.add(BondOrderProvider);
const modelBondProviderProperty = createModelBondProviderProperty(bondProviderRegistry);
const propertyContext = {
    runtime: SyncRuntimeContext,
    assetManager: new AssetManager(),
};

async function selectBondOrderProvider(model: Model, mode: BondOrdersMode) {
    const props = {
        provider: {
            name: BondOrderProviderName,
            params: { mode },
        },
    };
    const isAttached = model.customProperties.hasReference(ModelBondProvider.Descriptor);
    await modelBondProviderProperty.attach(propertyContext, model, props, !isAttached);
}

function structureWithSelectedBondProvider(model: Model): Structure {
    const bondProviderProps = ModelBondProvider.get(model);
    const context: BondProvider.Context = {};
    const provider = bondProviderProps ? bondProviderRegistry.get(bondProviderProps.name) : undefined;
    const bondProvider = provider?.isApplicable(model)
        ? provider.factory(model, bondProviderProps!.params, context)
        : undefined;
    const structure = Structure.ofModel(model, { bondProvider });
    if (bondProvider) context.structure = structure;
    return structure;
}

async function structureFromPdb(pdbText: string, mode: BondOrdersMode): Promise<Structure> {
    const model = await modelFromPdb(pdbText);
    await selectBondOrderProvider(model, mode);
    return structureWithSelectedBondProvider(model);
}

describe('bond provider custom model property', () => {
    it('perceives CONECT ligand C-C lazily via unit.bonds and getChild subsets', async () => {
        const structure = await structureFromPdb(PdbWithConectLigand, 'auto');
        for (const unit of atomicUnits(structure)) {
            expect(unit.props.bonds).toBeUndefined();
        }

        for (const unit of atomicUnits(structure)) {
            expect(unit.props.bonds).toBeUndefined();
        }
        expect(IndexPairBonds.Provider.get(structure.model)).toBeUndefined();

        const ligUnit = unitWithComp(structure, 'LIG');
        expect(countCompCCOrder(ligUnit, 'LIG', 2)).toBeGreaterThan(0);
        const perceived = ligUnit.bonds;
        expect(ligUnit.bonds).toBe(perceived);
        expect(ligUnit.props.bonds).toBe(perceived);
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
        const model = await modelFromPdb(PdbWithConectLigand);
        const existing = IndexPairBonds.fromData({
            pairs: { indexA: Column.ofIntArray([]), indexB: Column.ofIntArray([]) },
            count: model.atomicHierarchy.atoms._rowCount,
        });
        IndexPairBonds.Provider.set(model, existing);

        await selectBondOrderProvider(model, 'auto');
        const structure = structureWithSelectedBondProvider(model);
        expect(IndexPairBonds.Provider.get(model)).toBe(existing);
        void (structure.units[0] as Unit.Atomic).bonds;
    });

    it('does not perceive components covered by IntraBondOrderTable', async () => {
        const structure = await structureFromPdb(PdbWithKnownHis, 'auto');
        const unit = unitWithComp(structure, 'HIS');

        expect(countCompCCOrder(unit, 'HIS', 2)).toBe(1);
    });

    it('skips mol/SDF graphs and leaves file orders unchanged', async () => {
        const parsed = await parseMol(MolWithSingleCC).run();
        if (parsed.isError) throw new Error(parsed.message);
        const trajectory = await trajectoryFromMol(parsed.result).run();
        const model = trajectory.representative;
        const existing = IndexPairBonds.Provider.get(model);
        expect(existing).toBeDefined();

        const before = Array.from((Structure.ofModel(model).units[0] as Unit.Atomic).bonds.edgeProps.order);
        await selectBondOrderProvider(model, 'forceCompute');
        const structure = structureWithSelectedBondProvider(model);
        const unit = structure.units[0] as Unit.Atomic;
        expect(IndexPairBonds.Provider.get(model)).toBe(existing);
        expect(Array.from(unit.bonds.edgeProps.order)).toEqual(before);
        expect(countCompCCOrder(unit, 'MOL', 1)).toBeGreaterThan(0);
        expect(countCompCCOrder(unit, 'MOL', 2)).toBe(0);
    });

    it('leaves chem_comp_bond dictionary orders alone in auto mode', async () => {
        const parsed = await CIF.parseText(MmcifWithChemCompBond).run();
        if (parsed.isError) throw new Error(parsed.message);
        const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
        const model = trajectory.representative;
        await selectBondOrderProvider(model, 'auto');
        const structure = structureWithSelectedBondProvider(model);

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

    it('supports none, auto, and forceCompute modes through the model property', async () => {
        const parsed = await CIF.parseText(MmcifWithChemCompBond).run();
        if (parsed.isError) throw new Error(parsed.message);
        const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
        const model = trajectory.representative;

        await selectBondOrderProvider(model, 'auto');
        const autoStructure = structureWithSelectedBondProvider(model);
        const autoUnit = unitWithComp(autoStructure, 'LIG');
        const autoProvider = autoUnit.props.bondProvider;
        expect(autoProvider).toBeDefined();
        expect(countCompCCOrder(autoUnit, 'LIG', 1)).toBeGreaterThan(0);
        expect(countCompCCOrder(autoUnit, 'LIG', 2)).toBe(0);

        await selectBondOrderProvider(model, 'forceCompute');
        expect(countCompCCOrder(autoUnit, 'LIG', 2)).toBe(0);
        const forceStructure = structureWithSelectedBondProvider(model);
        const forceUnit = unitWithComp(forceStructure, 'LIG');
        expect(forceUnit.props.bondProvider).not.toBe(autoProvider);
        expect(countCompCCOrder(forceUnit, 'LIG', 2)).toBeGreaterThan(0);

        await selectBondOrderProvider(model, 'none');
        const noneStructure = structureWithSelectedBondProvider(model);
        const noneUnit = unitWithComp(noneStructure, 'LIG');
        expect(countCompCCOrder(noneUnit, 'LIG', 1)).toBeGreaterThan(0);
        expect(countCompCCOrder(noneUnit, 'LIG', 2)).toBe(0);
    });
});

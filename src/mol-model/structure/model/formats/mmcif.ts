/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, Table } from 'mol-data/db';
import { Interval, Segmentation } from 'mol-data/int';
import { Spacegroup, SpacegroupCell } from 'mol-math/geometry';
import { Vec3 } from 'mol-math/linear-algebra';
import UUID from 'mol-util/uuid';
import Format from '../format';
import Model from '../model';
import { AtomicConformation, AtomicData, AtomicSegments, AtomsSchema, ChainsSchema, ResiduesSchema } from '../properties/atomic';
import { Entities } from '../properties/common';
import { ModelSymmetry } from '../properties/symmetry';
import { getAtomicKeys } from '../properties/utils/atomic-keys';
import { ElementSymbol } from '../types';
import { createAssemblies } from './mmcif/assembly';
import { getIHMCoarse } from './mmcif/ihm';
import { getSequence } from './mmcif/sequence';

import mmCIF_Format = Format.mmCIF
import { Task } from 'mol-task';
import { getSecondaryStructureMmCif } from './mmcif/secondary-structure';

function findModelBounds({ data }: mmCIF_Format, startIndex: number) {
    const num = data.atom_site.pdbx_PDB_model_num;
    const atomCount = num.rowCount;
    if (!num.isDefined) return Interval.ofBounds(startIndex, atomCount);
    let endIndex = startIndex + 1;
    while (endIndex < atomCount && num.areValuesEqual(startIndex, endIndex)) endIndex++;
    return Interval.ofBounds(startIndex, endIndex);
}

function findHierarchyOffsets({ data }: mmCIF_Format, bounds: Interval) {
    if (Interval.size(bounds) === 0) return { residues: [], chains: [] };

    const start = Interval.start(bounds), end = Interval.end(bounds);
    const residues = [start], chains = [start];

    const { label_entity_id, auth_asym_id, auth_seq_id, pdbx_PDB_ins_code, label_comp_id } = data.atom_site;

    for (let i = start + 1; i < end; i++) {
        const newChain = !label_entity_id.areValuesEqual(i - 1, i) || !auth_asym_id.areValuesEqual(i - 1, i);
        const newResidue = newChain
            || !auth_seq_id.areValuesEqual(i - 1, i)
            || !pdbx_PDB_ins_code.areValuesEqual(i - 1, i)
            || !label_comp_id.areValuesEqual(i - 1, i);

        if (newResidue) residues[residues.length] = i;
        if (newChain) chains[chains.length] = i;
    }
    return { residues, chains };
}

function createHierarchyData({ data }: mmCIF_Format, bounds: Interval, offsets: { residues: ArrayLike<number>, chains: ArrayLike<number> }): AtomicData {
    const { atom_site } = data;
    const start = Interval.start(bounds), end = Interval.end(bounds);
    const atoms = Table.ofColumns(AtomsSchema, {
        type_symbol: Column.ofArray({ array: Column.mapToArray(Column.window(atom_site.type_symbol, start, end), ElementSymbol), schema: Column.Schema.Aliased<ElementSymbol>(Column.Schema.str) }),
        label_atom_id: Column.window(atom_site.label_atom_id, start, end),
        auth_atom_id: Column.window(atom_site.auth_atom_id, start, end),
        label_alt_id: Column.window(atom_site.label_alt_id, start, end),
        pdbx_formal_charge: Column.window(atom_site.pdbx_formal_charge, start, end)
    });
    const residues = Table.view(atom_site, ResiduesSchema, offsets.residues);
    // Optimize the numeric columns
    Table.columnToArray(residues, 'label_seq_id', Int32Array);
    Table.columnToArray(residues, 'auth_seq_id', Int32Array);
    const chains = Table.view(atom_site, ChainsSchema, offsets.chains);
    return { atoms, residues, chains };
}

function getConformation({ data }: mmCIF_Format, bounds: Interval): AtomicConformation {
    const start = Interval.start(bounds), end = Interval.end(bounds);
    const { atom_site } = data;
    return {
        id: UUID.create(),
        atomId: Column.window(atom_site.id, start, end),
        occupancy: Column.window(atom_site.occupancy, start, end),
        B_iso_or_equiv: Column.window(atom_site.B_iso_or_equiv, start, end),
        x: atom_site.Cartn_x.toArray({ array: Float32Array, start, end }),
        y: atom_site.Cartn_y.toArray({ array: Float32Array, start, end }),
        z: atom_site.Cartn_z.toArray({ array: Float32Array, start, end }),
    }
}

function getSymmetry(format: mmCIF_Format): ModelSymmetry {
    const assemblies = createAssemblies(format);
    const spacegroup = getSpacegroup(format);
    const isNonStandardCrytalFrame = checkNonStandardCrystalFrame(format, spacegroup);
    return { assemblies, spacegroup, isNonStandardCrytalFrame };
}

function checkNonStandardCrystalFrame(format: mmCIF_Format, spacegroup: Spacegroup) {
    const { atom_sites } = format.data;
    if (atom_sites._rowCount === 0) return false;
    // TODO: parse atom_sites transform and check if it corresponds to the toFractional matrix
    return false;
}

function getSpacegroup(format: mmCIF_Format): Spacegroup {
    const { symmetry, cell } = format.data;
    if (symmetry._rowCount === 0 || cell._rowCount === 0) return Spacegroup.ZeroP1;
    const groupName = symmetry['space_group_name_H-M'].value(0);
    const spaceCell = SpacegroupCell.create(groupName,
        Vec3.create(cell.length_a.value(0), cell.length_b.value(0), cell.length_c.value(0)),
        Vec3.scale(Vec3.zero(), Vec3.create(cell.angle_alpha.value(0), cell.angle_beta.value(0), cell.angle_gamma.value(0)), Math.PI / 180));

    return Spacegroup.create(spaceCell);
}

function isHierarchyDataEqual(a: AtomicData, b: AtomicData) {
    // need to cast because of how TS handles type resolution for interfaces https://github.com/Microsoft/TypeScript/issues/15300
    return Table.areEqual(a.chains as Table<ChainsSchema>, b.chains as Table<ChainsSchema>)
        && Table.areEqual(a.residues as Table<ResiduesSchema>, b.residues as Table<ResiduesSchema>)
        && Table.areEqual(a.atoms as Table<AtomsSchema>, b.atoms as Table<AtomsSchema>)
}

function createModel(format: mmCIF_Format, bounds: Interval, previous?: Model): Model {
    const hierarchyOffsets = findHierarchyOffsets(format, bounds);
    const hierarchyData = createHierarchyData(format, bounds, hierarchyOffsets);

    if (previous && isHierarchyDataEqual(previous.atomicHierarchy, hierarchyData)) {
        return {
            ...previous,
            atomicConformation: getConformation(format, bounds)
        };
    }

    const hierarchySegments: AtomicSegments = {
        residueSegments: Segmentation.ofOffsets(hierarchyOffsets.residues, bounds),
        chainSegments: Segmentation.ofOffsets(hierarchyOffsets.chains, bounds),
    }

    const entities: Entities = { data: format.data.entity, getEntityIndex: Column.createIndexer(format.data.entity.id) };

    const hierarchyKeys = getAtomicKeys(hierarchyData, entities, hierarchySegments);

    const atomicHierarchy = { ...hierarchyData, ...hierarchyKeys, ...hierarchySegments };

    const coarse = getIHMCoarse(format.data, entities);

    return {
        id: UUID.create(),
        label: format.data.entry.id.value(0),
        sourceData: format,
        modelNum: format.data.atom_site.pdbx_PDB_model_num.value(Interval.start(bounds)),
        entities,
        atomicHierarchy,
        sequence: getSequence(format.data, entities, atomicHierarchy),
        atomicConformation: getConformation(format, bounds),
        coarseHierarchy: coarse.hierarchy,
        coarseConformation: coarse.conformation,
        properties: {
            secondaryStructure: getSecondaryStructureMmCif(format.data, atomicHierarchy)
        },
        symmetry: getSymmetry(format)
    };
}

function buildModels(format: mmCIF_Format): Task<ReadonlyArray<Model>> {
    return Task.create('Create mmCIF Model', async ctx => {
        const atomCount = format.data.atom_site._rowCount;
        const isIHM = format.data.ihm_model_list._rowCount > 0;

        if (atomCount === 0) {
            return isIHM
                ? [createModel(format, Interval.Empty, void 0)]
                : [];
        }

        const models: Model[] = [];
        let modelStart = 0;
        while (modelStart < atomCount) {
            const bounds = findModelBounds(format, modelStart);
            const model = createModel(format, bounds, models.length > 0 ? models[models.length - 1] : void 0);
            models.push(model);
            modelStart = Interval.end(bounds);
        }
        return models;
    });
}

export default buildModels;
/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import UUID from 'mol-util/uuid'
import { Column, Table } from 'mol-data/db'
import { Interval, Segmentation } from 'mol-data/int'
import Format from '../format'
import Model from '../model'
import * as Hierarchy from '../properties/hierarchy'
import AtomSiteConformation from '../properties/atom-site-conformation'
import CoarseGrained from '../properties/coarse-grained'
import Symmetry from '../properties/symmetry'
import findHierarchyKeys from '../utils/hierarchy-keys'
import { ElementSymbol} from '../types'
import createAssemblies from './mmcif/assembly'

import mmCIF_Format = Format.mmCIF
import { getSequence } from './mmcif/sequence';
import { Entities } from '../properties/common';

function findModelBounds({ data }: mmCIF_Format, startIndex: number) {
    const num = data.atom_site.pdbx_PDB_model_num;
    const atomCount = num.rowCount;
    if (!num.isDefined) return Interval.ofBounds(startIndex, atomCount);
    let endIndex = startIndex + 1;
    while (endIndex < atomCount && num.areValuesEqual(startIndex, endIndex)) endIndex++;
    return Interval.ofBounds(startIndex, endIndex);
}

function findHierarchyOffsets({ data }: mmCIF_Format, bounds: Interval) {
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

function createHierarchyData({ data }: mmCIF_Format, bounds: Interval, offsets: { residues: ArrayLike<number>, chains: ArrayLike<number> }): Hierarchy.Data {
    const { atom_site } = data;
    const start = Interval.start(bounds), end = Interval.end(bounds);
    const atoms = Table.ofColumns(Hierarchy.AtomsSchema, {
        type_symbol: Column.ofArray({ array: Column.mapToArray(Column.window(atom_site.type_symbol, start, end), ElementSymbol), schema: Column.Schema.Aliased<ElementSymbol>(Column.Schema.str) }),
        label_atom_id: Column.window(atom_site.label_atom_id, start, end),
        auth_atom_id: Column.window(atom_site.auth_atom_id, start, end),
        label_alt_id: Column.window(atom_site.label_alt_id, start, end),
        pdbx_formal_charge: Column.window(atom_site.pdbx_formal_charge, start, end)
    });
    const residues = Table.view(atom_site, Hierarchy.ResiduesSchema, offsets.residues);
    // Optimize the numeric columns
    Table.columnToArray(residues, 'label_seq_id', Int32Array);
    Table.columnToArray(residues, 'auth_seq_id', Int32Array);
    const chains = Table.view(atom_site, Hierarchy.ChainsSchema, offsets.chains);
    return { atoms, residues, chains };
}

function getConformation({ data }: mmCIF_Format, bounds: Interval): AtomSiteConformation {
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

function getSymmetry(format: mmCIF_Format): Symmetry {
    return { assemblies: createAssemblies(format) };
}

function isHierarchyDataEqual(a: Hierarchy.Hierarchy, b: Hierarchy.Data) {
    // need to cast because of how TS handles type resolution for interfaces https://github.com/Microsoft/TypeScript/issues/15300
    return Table.areEqual(a.chains as Table<Hierarchy.ChainsSchema>, b.chains as Table<Hierarchy.ChainsSchema>)
        && Table.areEqual(a.residues as Table<Hierarchy.ResiduesSchema>, b.residues as Table<Hierarchy.ResiduesSchema>)
        && Table.areEqual(a.atoms as Table<Hierarchy.AtomsSchema>, b.atoms as Table<Hierarchy.AtomsSchema>)
}

function createModel(format: mmCIF_Format, bounds: Interval, previous?: Model): Model {
    const hierarchyOffsets = findHierarchyOffsets(format, bounds);
    const hierarchyData = createHierarchyData(format, bounds, hierarchyOffsets);

    if (previous && isHierarchyDataEqual(previous.hierarchy, hierarchyData)) {
        return {
            ...previous,
            atomSiteConformation: getConformation(format, bounds)
        };
    }

    const hierarchySegments: Hierarchy.Segments = {
        residueSegments: Segmentation.ofOffsets(hierarchyOffsets.residues, bounds),
        chainSegments: Segmentation.ofOffsets(hierarchyOffsets.chains, bounds),
    }

    const entities: Entities = { data: format.data.entity, getEntityIndex: Column.createIndexer(format.data.entity.id) };

    const hierarchyKeys = findHierarchyKeys(hierarchyData, entities, hierarchySegments);

    const hierarchy = { ...hierarchyData, ...hierarchyKeys, ...hierarchySegments };

    return {
        id: UUID.create(),
        sourceData: format,
        modelNum: format.data.atom_site.pdbx_PDB_model_num.value(Interval.start(bounds)),
        entities,
        hierarchy,
        sequence: getSequence(format.data, entities, hierarchy),
        atomSiteConformation: getConformation(format, bounds),
        coarseGrained: CoarseGrained.Empty,
        symmetry: getSymmetry(format),
        atomCount: Interval.size(bounds)
    };
}

function buildModels(format: mmCIF_Format): ReadonlyArray<Model> {
    const models: Model[] = [];
    const atomCount = format.data.atom_site._rowCount;

    if (atomCount === 0) return models;

    let modelStart = 0;
    while (modelStart < atomCount) {
        const bounds = findModelBounds(format, modelStart);
        const model = createModel(format, bounds, models.length > 0 ? models[models.length - 1] : void 0);
        models.push(model);
        modelStart = Interval.end(bounds);
    }
    return models;
}

export default buildModels;
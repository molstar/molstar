/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import UUID from 'mol-util/uuid'
import { Column, Table } from 'mol-data/db'
import { Interval, Segmentation } from 'mol-data/int'
import { Atoms } from 'mol-io/reader/gro/schema'
import Format from '../format'
import Model from '../model'
import * as Hierarchy from '../properties/hierarchy'
import Conformation from '../properties/conformation'
import findHierarchyKeys from '../utils/hierarchy-keys'
import { guessElement } from '../utils/guess-element'
import { ElementSymbol} from '../types'

import gro_Format = Format.gro

type HierarchyOffsets = { residues: ArrayLike<number>, chains: ArrayLike<number> }

function findHierarchyOffsets(atomsData: Atoms, bounds: Interval) {
    const start = Interval.start(bounds), end = Interval.end(bounds);
    const residues = [start], chains = [start];

    const { residueName, residueNumber } = atomsData;

    for (let i = start + 1; i < end; i++) {
        const newResidue = !residueNumber.areValuesEqual(i - 1, i)
            || !residueName.areValuesEqual(i - 1, i);
        console.log(residueName.value(i - 1), residueName.value(i), residueNumber.value(i - 1), residueNumber.value(i), newResidue)
        if (newResidue) residues[residues.length] = i;
    }
    console.log(residues, residues.length)
    return { residues, chains };
}

function guessElementSymbol (value: string) {
    return ElementSymbol(guessElement(value));
}

function createHierarchyData(atomsData: Atoms, offsets: HierarchyOffsets): Hierarchy.Data {
    console.log(atomsData.atomName)
    const atoms = Table.ofColumns(Hierarchy.AtomsSchema, {
        type_symbol: Column.ofArray({ array: Column.mapToArray(atomsData.atomName, guessElementSymbol), schema: Column.Schema.Aliased<ElementSymbol>(Column.Schema.str) }),
        label_atom_id: atomsData.atomName,
        auth_atom_id: atomsData.atomName,
        label_alt_id: Column.Undefined(atomsData.count, Column.Schema.str),
        pdbx_formal_charge: Column.Undefined(atomsData.count, Column.Schema.int)
    });

    const residues = Table.view(Table.ofColumns(Hierarchy.ResiduesSchema, {
        group_PDB: Column.Undefined(atomsData.count, Column.Schema.str),
        label_comp_id: atomsData.residueName,
        auth_comp_id: atomsData.residueName,
        label_seq_id: atomsData.residueNumber,
        auth_seq_id: atomsData.residueNumber,
        pdbx_PDB_ins_code: Column.Undefined(atomsData.count, Column.Schema.str)
    }), Hierarchy.ResiduesSchema, offsets.residues);
    // Optimize the numeric columns
    Table.columnToArray(residues, 'label_seq_id', Int32Array);
    Table.columnToArray(residues, 'auth_seq_id', Int32Array);

    // const chains = Table.ofColumns(Hierarchy.ChainsSchema, {
    //     label_asym_id: Column.ofConst('A', atomsData.count, Column.Schema.str),
    //     auth_asym_id: Column.ofConst('A', atomsData.count, Column.Schema.str),
    //     label_entity_id: Column.Undefined(atomsData.count, Column.Schema.str)
    // });

    const chains = Table.ofUndefinedColumns(Hierarchy.ChainsSchema, 0);
    const entities = Table.ofUndefinedColumns(Hierarchy.EntitySchema, 0);

    return { atoms, residues, chains, entities };
}

function getConformation(atoms: Atoms): Conformation {
    return {
        id: UUID.create(),
        atomId: atoms.atomNumber,
        occupancy: Column.Undefined(atoms.count, Column.Schema.int),
        B_iso_or_equiv: Column.Undefined(atoms.count, Column.Schema.float),
        x: Column.mapToArray(atoms.x, x => x * 10, Float32Array),
        y: Column.mapToArray(atoms.y, y => y * 10, Float32Array),
        z: Column.mapToArray(atoms.z, z => z * 10, Float32Array)
    }
}

function isHierarchyDataEqual(a: Hierarchy.Hierarchy, b: Hierarchy.Data) {
    // need to cast because of how TS handles type resolution for interfaces https://github.com/Microsoft/TypeScript/issues/15300
    return Table.areEqual(a.residues as Table<Hierarchy.ResiduesSchema>, b.residues as Table<Hierarchy.ResiduesSchema>)
        && Table.areEqual(a.atoms as Table<Hierarchy.AtomsSchema>, b.atoms as Table<Hierarchy.AtomsSchema>)
}

function createModel(format: gro_Format, modelNum: number, previous?: Model): Model {
    const structure = format.data.structures[modelNum];
    const bounds = Interval.ofBounds(0, structure.atoms.count);

    const hierarchyOffsets = findHierarchyOffsets(structure.atoms, bounds);
    const hierarchyData = createHierarchyData(structure.atoms, hierarchyOffsets);

    if (previous && isHierarchyDataEqual(previous.hierarchy, hierarchyData)) {
        return {
            ...previous,
            conformation: getConformation(structure.atoms)
        };
    }

    const hierarchySegments: Hierarchy.Segments = {
        residueSegments: Segmentation.ofOffsets(hierarchyOffsets.residues, bounds),
        chainSegments: Segmentation.ofOffsets(hierarchyOffsets.chains, bounds),
    }
    const hierarchyKeys = findHierarchyKeys(hierarchyData, hierarchySegments);
    return {
        id: UUID.create(),
        sourceData: format,
        modelNum,
        hierarchy: { ...hierarchyData, ...hierarchyKeys, ...hierarchySegments },
        conformation: getConformation(structure.atoms),
        symmetry: { assemblies: [] },
        atomCount: structure.atoms.count
    };
}

function buildModels(format: gro_Format): ReadonlyArray<Model> {
    const models: Model[] = [];

    format.data.structures.forEach((_, i) => {
        const model = createModel(format, i, models.length > 0 ? models[models.length - 1] : void 0);
        models.push(model);
    });
    return models;
}

export default buildModels;
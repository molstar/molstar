/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, Table } from 'mol-data/db';
import { Interval, Segmentation } from 'mol-data/int';
import { Spacegroup, SpacegroupCell } from 'mol-math/geometry';
import { Vec3 } from 'mol-math/linear-algebra';
import { Task } from 'mol-task';
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
import { getSecondaryStructureMmCif } from './mmcif/secondary-structure';
import { getSequence } from './mmcif/sequence';
import { sortAtomSite } from './mmcif/sort';
import { mmCIF_Database } from 'mol-io/reader/cif/schema/mmcif';

import mmCIF_Format = Format.mmCIF
type AtomSite = mmCIF_Database['atom_site']

function findModelBounds({ data }: mmCIF_Format, startIndex: number) {
    const num = data.atom_site.pdbx_PDB_model_num;
    const atomCount = num.rowCount;
    if (!num.isDefined) return Interval.ofBounds(startIndex, atomCount);
    let endIndex = startIndex + 1;
    while (endIndex < atomCount && num.areValuesEqual(startIndex, endIndex)) endIndex++;
    return Interval.ofBounds(startIndex, endIndex);
}

function findHierarchyOffsets(atom_site: AtomSite) {
    if (atom_site._rowCount === 0) return { residues: [], chains: [] };

    const start = 0, end = atom_site._rowCount;
    const residues = [start], chains = [start];

    const { label_entity_id, label_asym_id, label_seq_id, auth_seq_id, pdbx_PDB_ins_code, label_comp_id } = atom_site;

    for (let i = start + 1; i < end; i++) {
        const newChain = !label_entity_id.areValuesEqual(i - 1, i) || !label_asym_id.areValuesEqual(i - 1, i);
        const newResidue = newChain
            || !label_seq_id.areValuesEqual(i - 1, i)
            || !auth_seq_id.areValuesEqual(i - 1, i)
            || !pdbx_PDB_ins_code.areValuesEqual(i - 1, i)
            || !label_comp_id.areValuesEqual(i - 1, i);

        if (newResidue) residues[residues.length] = i;
        if (newChain) chains[chains.length] = i;
    }
    return { residues, chains };
}

function createHierarchyData(atom_site: AtomSite, offsets: { residues: ArrayLike<number>, chains: ArrayLike<number> }): AtomicData {
    const atoms = Table.ofColumns(AtomsSchema, {
        type_symbol: Column.ofArray({ array: Column.mapToArray(atom_site.type_symbol, ElementSymbol), schema: Column.Schema.Aliased<ElementSymbol>(Column.Schema.str) }),
        label_atom_id: atom_site.label_atom_id,
        auth_atom_id: atom_site.auth_atom_id,
        label_alt_id: atom_site.label_alt_id,
        pdbx_formal_charge: atom_site.pdbx_formal_charge
    });
    const residues = Table.view(atom_site, ResiduesSchema, offsets.residues);
    // Optimize the numeric columns
    Table.columnToArray(residues, 'label_seq_id', Int32Array);
    Table.columnToArray(residues, 'auth_seq_id', Int32Array);
    const chains = Table.view(atom_site, ChainsSchema, offsets.chains);
    return { atoms, residues, chains };
}

function getConformation(atom_site: AtomSite): AtomicConformation {
    return {
        id: UUID.create(),
        atomId: atom_site.id,
        occupancy: atom_site.occupancy,
        B_iso_or_equiv: atom_site.B_iso_or_equiv,
        x: atom_site.Cartn_x.toArray({ array: Float32Array }),
        y: atom_site.Cartn_y.toArray({ array: Float32Array }),
        z: atom_site.Cartn_z.toArray({ array: Float32Array }),
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

function modResMap(format: mmCIF_Format) {
    const data = format.data.pdbx_struct_mod_residue;
    const map = new Map<string, string>();
    const comp_id = data.label_comp_id.isDefined ? data.label_comp_id : data.auth_comp_id;
    const parent_id = data.parent_comp_id;

    for (let i = 0; i < data._rowCount; i++) {
        map.set(comp_id.value(i), parent_id.value(i));
    }

    return map;
}

function createModel(format: mmCIF_Format, atom_site: AtomSite, previous?: Model): Model {
    const hierarchyOffsets = findHierarchyOffsets(atom_site);
    const hierarchyData = createHierarchyData(atom_site, hierarchyOffsets);

    if (previous && isHierarchyDataEqual(previous.atomicHierarchy, hierarchyData)) {
        return {
            ...previous,
            atomicConformation: getConformation(atom_site)
        };
    }

    const hierarchySegments: AtomicSegments = {
        residueSegments: Segmentation.ofOffsets(hierarchyOffsets.residues, Interval.ofBounds(0, atom_site._rowCount)),
        chainSegments: Segmentation.ofOffsets(hierarchyOffsets.chains, Interval.ofBounds(0, atom_site._rowCount)),
    }

    const entities: Entities = { data: format.data.entity, getEntityIndex: Column.createIndexer(format.data.entity.id) };

    const hierarchyKeys = getAtomicKeys(hierarchyData, entities, hierarchySegments);

    const atomicHierarchy = { ...hierarchyData, ...hierarchyKeys, ...hierarchySegments };

    const coarse = getIHMCoarse(format.data, entities);

    const label = format.data.entry.id.valueKind(0) === Column.ValueKind.Present
        ? format.data.entry.id.value(0)
        : format.data._name;

    const modifiedResidueNameMap = modResMap(format);

    return {
        id: UUID.create(),
        label,
        sourceData: format,
        modelNum: format.data.atom_site.pdbx_PDB_model_num.value(0),
        entities,
        atomicHierarchy,
        sequence: getSequence(format.data, entities, atomicHierarchy, modifiedResidueNameMap),
        atomicConformation: getConformation(atom_site),
        coarseHierarchy: coarse.hierarchy,
        coarseConformation: coarse.conformation,
        properties: {
            secondaryStructure: getSecondaryStructureMmCif(format.data, atomicHierarchy),
            modifiedResidueNameMap
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
                ? [createModel(format, format.data.atom_site, void 0)]
                : [];
        }

        const models: Model[] = [];
        let modelStart = 0;
        while (modelStart < atomCount) {
            const bounds = findModelBounds(format, modelStart);

            const atom_site = await sortAtomSite(ctx, format.data.atom_site, 0, Interval.end(bounds));
            const model = createModel(format, atom_site, models.length > 0 ? models[models.length - 1] : void 0);
            models.push(model);
            modelStart = Interval.end(bounds);
        }
        return models;
    });
}

export default buildModels;
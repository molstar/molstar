/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, Table } from 'mol-data/db';
import { mmCIF_Database, mmCIF_Schema } from 'mol-io/reader/cif/schema/mmcif';
import { Spacegroup, SpacegroupCell, SymmetryOperator } from 'mol-math/geometry';
import { Tensor, Vec3 } from 'mol-math/linear-algebra';
import { Task, RuntimeContext } from 'mol-task';
import UUID from 'mol-util/uuid';
import Format from '../format';
import Model from '../model';
import { Entities } from '../properties/common';
import { CustomProperties } from '../properties/custom';
import { ModelSymmetry } from '../properties/symmetry';
import { createAssemblies } from './mmcif/assembly';
import { getAtomicHierarchyAndConformation } from './mmcif/atomic';
import { ComponentBond } from './mmcif/bonds';
import { getIHMCoarse, EmptyIHMCoarse, IHMData } from './mmcif/ihm';
import { getSecondaryStructureMmCif } from './mmcif/secondary-structure';
import { getSequence } from './mmcif/sequence';
import { sortAtomSite } from './mmcif/sort';

import mmCIF_Format = Format.mmCIF

type AtomSite = mmCIF_Database['atom_site']

function getSymmetry(format: mmCIF_Format): ModelSymmetry {
    const assemblies = createAssemblies(format);
    const spacegroup = getSpacegroup(format);
    const isNonStandardCrytalFrame = checkNonStandardCrystalFrame(format, spacegroup);
    return { assemblies, spacegroup, isNonStandardCrytalFrame, ncsOperators: getNcsOperators(format) };
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

function getNcsOperators(format: mmCIF_Format) {
    const { struct_ncs_oper } = format.data;
    if (struct_ncs_oper._rowCount === 0) return void 0;
    const { id, matrix, vector } = struct_ncs_oper;

    const matrixSpace = mmCIF_Schema.struct_ncs_oper.matrix.space, vectorSpace = mmCIF_Schema.struct_ncs_oper.vector.space;

    const opers: SymmetryOperator[] = [];
    for (let i = 0; i < struct_ncs_oper._rowCount; i++) {
        const m = Tensor.toMat3(matrixSpace, matrix.value(i));
        const v = Tensor.toVec3(vectorSpace, vector.value(i));
        opers[i] = SymmetryOperator.ofRotationAndOffset(`ncs_${id.value(i)}`, m, v);
    }
    return opers;
}
function getModifiedResidueNameMap(format: mmCIF_Format) {
    const data = format.data.pdbx_struct_mod_residue;
    const map = new Map<string, string>();
    const comp_id = data.label_comp_id.isDefined ? data.label_comp_id : data.auth_comp_id;
    const parent_id = data.parent_comp_id;

    for (let i = 0; i < data._rowCount; i++) {
        map.set(comp_id.value(i), parent_id.value(i));
    }

    return map;
}

function getAsymIdSerialMap(format: mmCIF_Format) {
    const data = format.data.struct_asym;
    const map = new Map<string, number>();
    let serial = 0

    const id = data.id
    const count = data._rowCount
    for (let i = 0; i < count; ++i) {
        const _id = id.value(i)
        if (!map.has(_id)) {
            map.set(_id, serial)
            serial += 1
        }
    }

    return map;
}

function createStandardModel(format: mmCIF_Format, atom_site: AtomSite, entities: Entities, previous?: Model): Model {
    const atomic = getAtomicHierarchyAndConformation(format, atom_site, entities, previous);
    if (previous && atomic.sameAsPrevious) {
        return { ...previous, atomicConformation: atomic.conformation };
    }

    const coarse = EmptyIHMCoarse;
    const label = format.data.entry.id.valueKind(0) === Column.ValueKind.Present
        ? format.data.entry.id.value(0)
        : format.data._name;

    const modifiedResidueNameMap = getModifiedResidueNameMap(format);
    const asymIdSerialMap = getAsymIdSerialMap(format)

    return {
        id: UUID.create(),
        label,
        sourceData: format,
        modelNum: atom_site.pdbx_PDB_model_num.value(0),
        entities,
        symmetry: getSymmetry(format),
        sequence: getSequence(format.data, entities, atomic.hierarchy, modifiedResidueNameMap),
        atomicHierarchy: atomic.hierarchy,
        atomicConformation: atomic.conformation,
        coarseHierarchy: coarse.hierarchy,
        coarseConformation: coarse.conformation,
        properties: {
            secondaryStructure: getSecondaryStructureMmCif(format.data, atomic.hierarchy),
            modifiedResidueNameMap,
            asymIdSerialMap
        },
        customProperties: new CustomProperties(),
        _staticPropertyData: Object.create(null),
        _dynamicPropertyData: Object.create(null)
    };
}

function createModelIHM(format: mmCIF_Format, data: IHMData): Model {
    const atomic = getAtomicHierarchyAndConformation(format, data.atom_site, data.entities);
    const coarse = getIHMCoarse(data);
    const modifiedResidueNameMap = getModifiedResidueNameMap(format);
    const asymIdSerialMap = getAsymIdSerialMap(format)

    return {
        id: UUID.create(),
        label: data.model_name,
        sourceData: format,
        modelNum: data.model_id,
        entities: data.entities,
        symmetry: getSymmetry(format),
        sequence: getSequence(format.data, data.entities, atomic.hierarchy, modifiedResidueNameMap),
        atomicHierarchy: atomic.hierarchy,
        atomicConformation: atomic.conformation,
        coarseHierarchy: coarse.hierarchy,
        coarseConformation: coarse.conformation,
        properties: {
            secondaryStructure: getSecondaryStructureMmCif(format.data, atomic.hierarchy),
            modifiedResidueNameMap,
            asymIdSerialMap
        },
        customProperties: new CustomProperties(),
        _staticPropertyData: Object.create(null),
        _dynamicPropertyData: Object.create(null)
    };
}

function attachProps(model: Model) {
    ComponentBond.attachFromMmCif(model);
}

function findModelEnd(num: Column<number>, startIndex: number) {
    const rowCount = num.rowCount;
    if (!num.isDefined) return rowCount;
    let endIndex = startIndex + 1;
    while (endIndex < rowCount && num.areValuesEqual(startIndex, endIndex)) endIndex++;
    return endIndex;
}

async function readStandard(ctx: RuntimeContext, format: mmCIF_Format) {
    const atomCount = format.data.atom_site._rowCount;
    const entities: Entities = { data: format.data.entity, getEntityIndex: Column.createIndexer(format.data.entity.id) };

    const models: Model[] = [];
    let modelStart = 0;
    while (modelStart < atomCount) {
        const modelEnd = findModelEnd(format.data.atom_site.pdbx_PDB_model_num, modelStart);
        const atom_site = await sortAtomSite(ctx, format.data.atom_site, modelStart, modelEnd);
        const model = createStandardModel(format, atom_site, entities, models.length > 0 ? models[models.length - 1] : void 0);
        attachProps(model);
        models.push(model);
        modelStart = modelEnd;
    }
    return models;
}

function splitTable<T extends Table<any>>(table: T, col: Column<number>) {
    const ret = new Map<number, T>()
    const rowCount = table._rowCount;
    let modelStart = 0;
    while (modelStart < rowCount) {
        const modelEnd = findModelEnd(col, modelStart);
        const id = col.value(modelStart);
        const window = Table.window(table, table._schema, modelStart, modelEnd) as T;
        ret.set(id, window);
        modelStart = modelEnd;
    }
    return ret;
}

async function readIHM(ctx: RuntimeContext, format: mmCIF_Format) {
    const { ihm_model_list } = format.data;
    const entities: Entities = { data: format.data.entity, getEntityIndex: Column.createIndexer(format.data.entity.id) };

    // TODO: will IHM require sorting or will we trust it?
    const atom_sites = splitTable(format.data.atom_site, format.data.atom_site.ihm_model_id);
    const sphere_sites = splitTable(format.data.ihm_sphere_obj_site, format.data.ihm_sphere_obj_site.model_id);
    const gauss_sites = splitTable(format.data.ihm_gaussian_obj_site, format.data.ihm_gaussian_obj_site.model_id);

    const models: Model[] = [];

    const { model_id, model_name } = ihm_model_list;
    for (let i = 0; i < ihm_model_list._rowCount; i++) {
        const id = model_id.value(i);
        const data: IHMData = {
            model_id: id,
            model_name: model_name.value(i),
            entities: entities,
            atom_site: atom_sites.has(id) ? atom_sites.get(id)! : Table.window(format.data.atom_site, format.data.atom_site._schema, 0, 0),
            ihm_sphere_obj_site: sphere_sites.has(id) ? sphere_sites.get(id)! : Table.window(format.data.ihm_sphere_obj_site, format.data.ihm_sphere_obj_site._schema, 0, 0),
            ihm_gaussian_obj_site: gauss_sites.has(id) ? gauss_sites.get(id)! : Table.window(format.data.ihm_gaussian_obj_site, format.data.ihm_gaussian_obj_site._schema, 0, 0)
        };
        const model = createModelIHM(format, data);
        attachProps(model);
        models.push(createModelIHM(format, data));
    }

    return models;
}

function buildModels(format: mmCIF_Format): Task<ReadonlyArray<Model>> {
    return Task.create('Create mmCIF Model', async ctx => {
        const isIHM = format.data.ihm_model_list._rowCount > 0;
        return isIHM ? await readIHM(ctx, format) : await readStandard(ctx, format);
    });
}

export default buildModels;
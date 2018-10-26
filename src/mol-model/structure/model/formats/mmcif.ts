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
import { Model } from '../model';
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
import { StructConn } from './mmcif/bonds/struct_conn';
import { ChemicalComponent, ChemicalComponentMap } from '../properties/chemical-component';
import { ComponentType, getMoleculeType } from '../types';

import mmCIF_Format = Format.mmCIF
import { SaccharideComponentMap, SaccharideComponent, SaccharidesSnfgMap, SaccharideCompIdMap } from 'mol-model/structure/structure/carbohydrates/constants';

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
        if (!SymmetryOperator.checkIfRotationAndTranslation(m, v)) continue;
        opers[opers.length] = SymmetryOperator.ofRotationAndOffset(`ncs_${id.value(i)}`, m, v);
    }
    return opers;
}
function getModifiedResidueNameMap(format: mmCIF_Format): Model['properties']['modifiedResidues'] {
    const data = format.data.pdbx_struct_mod_residue;
    const parentId = new Map<string, string>();
    const details = new Map<string, string>();
    const comp_id = data.label_comp_id.isDefined ? data.label_comp_id : data.auth_comp_id;
    const parent_id = data.parent_comp_id, details_data = data.details;

    for (let i = 0; i < data._rowCount; i++) {
        const id = comp_id.value(i);
        parentId.set(id, parent_id.value(i));
        details.set(id, details_data.value(i));
    }

    return { parentId, details };
}

function getAsymIdSerialMap(format: mmCIF_Format): ReadonlyMap<string, number> {
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

function getChemicalComponentMap(format: mmCIF_Format): ChemicalComponentMap {
    const map = new Map<string, ChemicalComponent>();
    const { id, type, name, pdbx_synonyms, formula, formula_weight } = format.data.chem_comp
    for (let i = 0, il = id.rowCount; i < il; ++i) {
        const _id = id.value(i)
        const _type = type.value(i)
        const cc: ChemicalComponent = {
            id: _id,
            type: ComponentType[_type],
            moleculeType: getMoleculeType(_type, _id),
            name: name.value(i),
            synonyms: pdbx_synonyms.value(i),
            formula: formula.value(i),
            formulaWeight: formula_weight.value(i),
        }
        map.set(_id, cc)
    }
    return map
}

function getSaccharideComponentMap(format: mmCIF_Format): SaccharideComponentMap {
    const map = new Map<string, SaccharideComponent>();
    const { pdbx_chem_comp_identifier } = format.data
    if (pdbx_chem_comp_identifier._rowCount > 0) {
        const { type, comp_id, identifier } = pdbx_chem_comp_identifier
        for (let i = 0, il = pdbx_chem_comp_identifier._rowCount; i < il; ++i) {
            if (type.value(i) === 'SNFG CARB SYMBOL') {
                const snfgName = identifier.value(i)
                const saccharideComp = SaccharidesSnfgMap.get(snfgName)
                if (saccharideComp) {
                    map.set(comp_id.value(i), saccharideComp)
                } else {
                    console.warn(`Unknown SNFG name '${snfgName}'`)
                }
            }
        }
        return map
    } else {
        return SaccharideCompIdMap
    }
}

export interface FormatData {
    modifiedResidues: Model['properties']['modifiedResidues']
    asymIdSerialMap: Model['properties']['asymIdSerialMap']
    chemicalComponentMap: Model['properties']['chemicalComponentMap']
    saccharideComponentMap: Model['properties']['saccharideComponentMap']
}

function getFormatData(format: mmCIF_Format): FormatData {
    return {
        modifiedResidues: getModifiedResidueNameMap(format),
        asymIdSerialMap: getAsymIdSerialMap(format),
        chemicalComponentMap: getChemicalComponentMap(format),
        saccharideComponentMap: getSaccharideComponentMap(format)
    }
}

function createStandardModel(format: mmCIF_Format, atom_site: AtomSite, entities: Entities, formatData: FormatData, previous?: Model): Model {
    const atomic = getAtomicHierarchyAndConformation(format, atom_site, entities, formatData, previous);
    if (previous && atomic.sameAsPrevious) {
        return {
            ...previous,
            id: UUID.create(),
            modelNum: atom_site.pdbx_PDB_model_num.value(0),
            atomicConformation: atomic.conformation,
            _dynamicPropertyData: Object.create(null)
        };
    }

    const coarse = EmptyIHMCoarse;
    const label = format.data.entry.id.valueKind(0) === Column.ValueKind.Present
        ? format.data.entry.id.value(0)
        : format.data._name;

    return {
        id: UUID.create(),
        label,
        sourceData: format,
        modelNum: atom_site.pdbx_PDB_model_num.value(0),
        entities,
        symmetry: getSymmetry(format),
        sequence: getSequence(format.data, entities, atomic.hierarchy, formatData.modifiedResidues.parentId),
        atomicHierarchy: atomic.hierarchy,
        atomicConformation: atomic.conformation,
        coarseHierarchy: coarse.hierarchy,
        coarseConformation: coarse.conformation,
        properties: {
            secondaryStructure: getSecondaryStructureMmCif(format.data, atomic.hierarchy),
            ...formatData
        },
        customProperties: new CustomProperties(),
        _staticPropertyData: Object.create(null),
        _dynamicPropertyData: Object.create(null)
    };
}

function createModelIHM(format: mmCIF_Format, data: IHMData, formatData: FormatData): Model {
    const atomic = getAtomicHierarchyAndConformation(format, data.atom_site, data.entities, formatData);
    const coarse = getIHMCoarse(data, formatData);

    return {
        id: UUID.create(),
        label: data.model_name,
        sourceData: format,
        modelNum: data.model_id,
        entities: data.entities,
        symmetry: getSymmetry(format),
        sequence: getSequence(format.data, data.entities, atomic.hierarchy, formatData.modifiedResidues.parentId),
        atomicHierarchy: atomic.hierarchy,
        atomicConformation: atomic.conformation,
        coarseHierarchy: coarse.hierarchy,
        coarseConformation: coarse.conformation,
        properties: {
            secondaryStructure: getSecondaryStructureMmCif(format.data, atomic.hierarchy),
            ...formatData
        },
        customProperties: new CustomProperties(),
        _staticPropertyData: Object.create(null),
        _dynamicPropertyData: Object.create(null)
    };
}

function attachProps(model: Model) {
    ComponentBond.attachFromMmCif(model);
    StructConn.attachFromMmCif(model);
}

function findModelEnd(num: Column<number>, startIndex: number) {
    const rowCount = num.rowCount;
    if (!num.isDefined) return rowCount;
    let endIndex = startIndex + 1;
    while (endIndex < rowCount && num.areValuesEqual(startIndex, endIndex)) endIndex++;
    return endIndex;
}

async function readStandard(ctx: RuntimeContext, format: mmCIF_Format, formatData: FormatData) {
    const atomCount = format.data.atom_site._rowCount;
    const entities: Entities = { data: format.data.entity, getEntityIndex: Column.createIndexer(format.data.entity.id) };

    const models: Model[] = [];
    let modelStart = 0;
    while (modelStart < atomCount) {
        const modelEnd = findModelEnd(format.data.atom_site.pdbx_PDB_model_num, modelStart);
        const atom_site = await sortAtomSite(ctx, format.data.atom_site, modelStart, modelEnd);
        const model = createStandardModel(format, atom_site, entities, formatData, models.length > 0 ? models[models.length - 1] : void 0);
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

async function readIHM(ctx: RuntimeContext, format: mmCIF_Format, formatData: FormatData) {
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
        const model = createModelIHM(format, data, formatData);
        attachProps(model);
        models.push(model);
    }

    return models;
}

function buildModels(format: mmCIF_Format): Task<ReadonlyArray<Model>> {
    const formatData = getFormatData(format)
    return Task.create('Create mmCIF Model', async ctx => {
        const isIHM = format.data.ihm_model_list._rowCount > 0;
        return isIHM ? await readIHM(ctx, format, formatData) : await readStandard(ctx, format, formatData);
    });
}

export default buildModels;
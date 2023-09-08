/**
 * Copyright (c) 2017-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { Model } from '../../mol-model/structure/model/model';
import { RuntimeContext, Task } from '../../mol-task';
import { ModelFormat } from '../format';
import { CifFrame, CIF, CifFile } from '../../mol-io/reader/cif';
import { mmCIF_Database } from '../../mol-io/reader/cif/schema/mmcif';
import { createModels } from './basic/parser';
import { ModelSymmetry } from './property/symmetry';
import { ModelSecondaryStructure } from './property/secondary-structure';
import { Column, Table } from '../../mol-data/db';
import { AtomSiteAnisotrop } from './property/anisotropic';
import { ComponentBond } from './property/bonds/chem_comp';
import { StructConn } from './property/bonds/struct_conn';
import { ArrayTrajectory, Trajectory } from '../../mol-model/structure';
import { GlobalModelTransformInfo } from '../../mol-model/structure/model/properties/global-transform';
import { BasicSchema, createBasic } from './basic/schema';
import { CCD_Database } from '../../mol-io/reader/cif/schema/ccd';
import { EntityBuilder } from './common/entity';
import { MoleculeType } from '../../mol-model/structure/model/types';
import { ComponentBuilder } from './common/component';

function modelSymmetryFromMmcif(model: Model) {
    if (!MmcifFormat.is(model.sourceData)) return;
    return ModelSymmetry.fromData(model.sourceData.data.db);
}
ModelSymmetry.Provider.formatRegistry.add('mmCIF', modelSymmetryFromMmcif);

function secondaryStructureFromMmcif(model: Model) {
    if (!MmcifFormat.is(model.sourceData)) return;
    const { struct_conf, struct_sheet_range } = model.sourceData.data.db;
    return ModelSecondaryStructure.fromStruct(struct_conf, struct_sheet_range, model.atomicHierarchy);
}
ModelSecondaryStructure.Provider.formatRegistry.add('mmCIF', secondaryStructureFromMmcif);

function atomSiteAnisotropFromMmcif(model: Model) {
    if (!MmcifFormat.is(model.sourceData)) return;
    const { atom_site_anisotrop } = model.sourceData.data.db;
    const data = Table.ofColumns(AtomSiteAnisotrop.Schema, atom_site_anisotrop);
    const elementToAnsiotrop = AtomSiteAnisotrop.getElementToAnsiotrop(model.atomicConformation.atomId, atom_site_anisotrop.id);
    return { data, elementToAnsiotrop };
}
function atomSiteAnisotropApplicableMmcif(model: Model) {
    if (!MmcifFormat.is(model.sourceData)) return false;
    return model.sourceData.data.db.atom_site_anisotrop.U.isDefined;
}
AtomSiteAnisotrop.Provider.formatRegistry.add('mmCIF', atomSiteAnisotropFromMmcif, atomSiteAnisotropApplicableMmcif);

function componentBondFromMmcif(model: Model) {
    if (!MmcifFormat.is(model.sourceData)) return;
    const { chem_comp_bond } = model.sourceData.data.db;
    if (chem_comp_bond._rowCount === 0) return;
    return {
        data: chem_comp_bond,
        entries: ComponentBond.getEntriesFromChemCompBond(chem_comp_bond)
    };
}
ComponentBond.Provider.formatRegistry.add('mmCIF', componentBondFromMmcif);

function structConnFromMmcif(model: Model) {
    if (!MmcifFormat.is(model.sourceData)) return;
    const { struct_conn } = model.sourceData.data.db;
    if (struct_conn._rowCount === 0) return;
    const entries = StructConn.getEntriesFromStructConn(struct_conn, model);
    return {
        data: struct_conn,
        byAtomIndex: StructConn.getAtomIndexFromEntries(entries),
        entries,
    };
}
StructConn.Provider.formatRegistry.add('mmCIF', structConnFromMmcif);

GlobalModelTransformInfo.Provider.formatRegistry.add('mmCIF', GlobalModelTransformInfo.fromMmCif, GlobalModelTransformInfo.hasData);

//

export { MmcifFormat };

type MmcifFormat = ModelFormat<MmcifFormat.Data>

namespace MmcifFormat {
    export type Data = {
        db: mmCIF_Database,
        frame: CifFrame,
        file?: CifFile,
        /**
         * Original source format. Some formats, including PDB, are converted
         * to mmCIF before further processing.
         */
        source?: ModelFormat
    }
    export function is(x?: ModelFormat): x is MmcifFormat {
        return x?.kind === 'mmCIF';
    }

    export function fromFrame(frame: CifFrame, db?: mmCIF_Database, source?: ModelFormat, file?: CifFile): MmcifFormat {
        if (!db) db = CIF.schema.mmCIF(frame);
        return { kind: 'mmCIF', name: db._name, data: { db, file, frame, source } };
    }
}

export function trajectoryFromMmCIF(frame: CifFrame, file?: CifFile): Task<Trajectory> {
    const format = MmcifFormat.fromFrame(frame, undefined, undefined, file);
    const basic = createBasic(format.data.db, true);
    return Task.create('Create mmCIF Model', ctx => createModels(basic, format, ctx));
}

export { CCDFormat };

type CCDFormat = ModelFormat<CCDFormat.Data>

namespace CCDFormat {
    export type Data = {
        db: CCD_Database,
        frame: CifFrame
    }

    const CoordinateTypeProp = '__CcdCoordinateType__';
    export type CoordinateType = 'ideal' | 'model'
    export const CoordinateType = {
        get(model: Model): CoordinateType | undefined {
            return model._staticPropertyData[CoordinateTypeProp];
        },
        set(model: Model, type: CoordinateType) {
            return model._staticPropertyData[CoordinateTypeProp] = type;
        }
    };

    export function is(x?: ModelFormat): x is CCDFormat {
        return x?.kind === 'CCD';
    }

    export function fromFrame(frame: CifFrame, db?: CCD_Database): CCDFormat {
        if (!db) db = CIF.schema.CCD(frame);
        return { kind: 'CCD', name: db._name, data: { db, frame } };
    }
}

export function trajectoryFromCCD(frame: CifFrame): Task<Trajectory> {
    const format = CCDFormat.fromFrame(frame);
    return Task.create('Create CCD Models', ctx => createCcdModels(format.data.db, CCDFormat.fromFrame(frame), ctx));
}

async function createCcdModels(data: CCD_Database, format: CCDFormat, ctx: RuntimeContext) {
    const ideal = await createCcdModel(data, format, { coordinateType: 'ideal', cartn_x: 'pdbx_model_Cartn_x_ideal', cartn_y: 'pdbx_model_Cartn_y_ideal', cartn_z: 'pdbx_model_Cartn_z_ideal' }, ctx);
    const model = await createCcdModel(data, format, { coordinateType: 'model', cartn_x: 'model_Cartn_x', cartn_y: 'model_Cartn_y', cartn_z: 'model_Cartn_z' }, ctx);

    const models = [];
    if (ideal) models.push(ideal);
    if (model) models.push(model);
    for (let i = 0, il = models.length; i < il; ++i) {
        Model.TrajectoryInfo.set(models[i], { index: i, size: models.length });
    }

    return new ArrayTrajectory(models);
}

type CCDProps = { coordinateType: CCDFormat.CoordinateType, cartn_x: 'model_Cartn_x' | 'pdbx_model_Cartn_x_ideal', cartn_y: 'model_Cartn_y' | 'pdbx_model_Cartn_y_ideal', cartn_z: 'model_Cartn_z' | 'pdbx_model_Cartn_z_ideal' };
async function createCcdModel(data: CCD_Database, format: CCDFormat, props: CCDProps, ctx: RuntimeContext) {
    const { chem_comp, chem_comp_atom, chem_comp_bond } = data;
    const { coordinateType, cartn_x, cartn_y, cartn_z } = props;

    const name = chem_comp.name.value(0);
    const id = chem_comp.id.value(0);

    const { atom_id, charge, comp_id, pdbx_ordinal, type_symbol } = chem_comp_atom;
    const atomCount = chem_comp_atom._rowCount;

    const filteredRows: number[] = [];
    for (let i = 0; i < atomCount; i++) {
        if (chem_comp_atom[cartn_x].valueKind(i) > 0) continue;
        filteredRows[filteredRows.length] = i;
    }
    const filteredRowCount = filteredRows.length;

    const A = Column.ofConst('A', filteredRowCount, Column.Schema.str);
    const seq_id = Column.ofConst(1, filteredRowCount, Column.Schema.int);
    const entity_id = Column.ofConst('1', filteredRowCount, Column.Schema.str);
    const occupancy = Column.ofConst(1, filteredRowCount, Column.Schema.float);
    const model_num = Column.ofConst(1, filteredRowCount, Column.Schema.int);

    const filteredAtomId = Column.view(atom_id, filteredRows);
    const filteredCompId = Column.view(comp_id, filteredRows);
    const filteredX = Column.view(chem_comp_atom[cartn_x], filteredRows);
    const filteredY = Column.view(chem_comp_atom[cartn_y], filteredRows);
    const filteredZ = Column.view(chem_comp_atom[cartn_z], filteredRows);
    const filteredId = Column.view(pdbx_ordinal, filteredRows);
    const filteredTypeSymbol = Column.view(type_symbol, filteredRows);
    const filteredCharge = Column.view(charge, filteredRows);

    const model_atom_site = Table.ofPartialColumns(BasicSchema.atom_site, {
        auth_asym_id: A,
        auth_atom_id: filteredAtomId,
        auth_comp_id: filteredCompId,
        auth_seq_id: seq_id,
        Cartn_x: filteredX,
        Cartn_y: filteredY,
        Cartn_z: filteredZ,
        id: filteredId,

        label_asym_id: A,
        label_atom_id: filteredAtomId,
        label_comp_id: filteredCompId,
        label_seq_id: seq_id,
        label_entity_id: entity_id,

        occupancy,
        type_symbol: filteredTypeSymbol,

        pdbx_PDB_model_num: model_num,
        pdbx_formal_charge: filteredCharge
    }, filteredRowCount);

    const entityBuilder = new EntityBuilder();
    entityBuilder.setNames([[id, `${name} (${coordinateType})`]]);
    entityBuilder.getEntityId(id, MoleculeType.Unknown, 'A');

    const componentBuilder = new ComponentBuilder(seq_id, type_symbol);
    componentBuilder.setNames([[id, `${name} (${coordinateType})`]]);
    componentBuilder.add(id, 0);

    const basicModel = createBasic({
        entity: entityBuilder.getEntityTable(),
        chem_comp: componentBuilder.getChemCompTable(),
        atom_site: model_atom_site
    });
    const models = await createModels(basicModel, format, ctx);

    // all ideal or model coordinates might be absent
    if (!models.representative) return;

    const first = models.representative;
    const entries = ComponentBond.getEntriesFromChemCompBond(chem_comp_bond);
    ComponentBond.Provider.set(first, { data: chem_comp_bond, entries });
    CCDFormat.CoordinateType.set(first, coordinateType);

    return models.representative;
}
/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model } from '../../mol-model/structure/model/model';
import { RuntimeContext, Task } from '../../mol-task';
import { ModelFormat } from '../format';
import { CifFrame, CIF } from '../../mol-io/reader/cif';
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
import { CCD_Database, CCD_Schema } from '../../mol-io/reader/cif/schema/ccd';
import { EntityBuilder } from './common/entity';
import { BondType, MoleculeType } from '../../mol-model/structure/model/types';
import { ComponentBuilder } from './common/component';
import { cantorPairing } from '../../mol-data/util';
import { IndexPairBonds } from './property/bonds/index-pair';

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
        /**
         * Original source format. Some formats, including PDB, are converted
         * to mmCIF before further processing.
         */
        source?: ModelFormat
    }
    export function is(x?: ModelFormat): x is MmcifFormat {
        return x?.kind === 'mmCIF';
    }

    export function fromFrame(frame: CifFrame, db?: mmCIF_Database, source?: ModelFormat): MmcifFormat {
        if (!db) db = CIF.schema.mmCIF(frame);
        return { kind: 'mmCIF', name: db._name, data: { db, frame, source } };
    }
}

export function trajectoryFromMmCIF(frame: CifFrame): Task<Trajectory> {
    const format = MmcifFormat.fromFrame(frame);
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
    const model = await createCcdModel(data, format, { suffix: '(model)', cartn_x: 'model_Cartn_x', cartn_y: 'model_Cartn_y', cartn_z: 'model_Cartn_z' }, ctx);
    const ideal = await createCcdModel(data, format, { suffix: '(ideal)', cartn_x: 'pdbx_model_Cartn_x_ideal', cartn_y: 'pdbx_model_Cartn_y_ideal', cartn_z: 'pdbx_model_Cartn_z_ideal' }, ctx);

    const models = [model.representative, ideal.representative];
    Model.TrajectoryInfo.set(models[0], { index: 0, size: models.length });
    Model.TrajectoryInfo.set(models[1], { index: 1, size: models.length });

    return new ArrayTrajectory(models);
}

type x = keyof Pick<CCD_Schema['chem_comp_atom'], 'model_Cartn_x'> | keyof Pick<CCD_Schema['chem_comp_atom'], 'pdbx_model_Cartn_x_ideal'>;
type y = keyof Pick<CCD_Schema['chem_comp_atom'], 'model_Cartn_y'> | keyof Pick<CCD_Schema['chem_comp_atom'], 'pdbx_model_Cartn_y_ideal'>;
type z = keyof Pick<CCD_Schema['chem_comp_atom'], 'model_Cartn_z'> | keyof Pick<CCD_Schema['chem_comp_atom'], 'pdbx_model_Cartn_z_ideal'>;
type CCDProps = { suffix: string, cartn_x: x, cartn_y: y, cartn_z: z };
async function createCcdModel(data: CCD_Database, format: CCDFormat, props: CCDProps, ctx: RuntimeContext) {
    const { chem_comp, chem_comp_atom, chem_comp_bond } = data;
    const { suffix, cartn_x, cartn_y, cartn_z } = props;

    const name = chem_comp.name.value(0);

    const { atom_id, charge, comp_id, pdbx_ordinal, type_symbol } = chem_comp_atom;
    const atomCount = chem_comp_atom._rowCount;
    const A = Column.ofConst('A', atomCount, Column.Schema.str);
    const seq_id = Column.ofConst(1, atomCount, Column.Schema.int);
    const entity_id = Column.ofConst('1', atomCount, Column.Schema.str);
    const occupancy = Column.ofConst(1, atomCount, Column.Schema.float);
    const model_num = Column.ofConst(1, atomCount, Column.Schema.int);

    const model_atom_site = Table.ofPartialColumns(BasicSchema.atom_site, {
        auth_asym_id: A,
        auth_atom_id: atom_id,
        auth_comp_id: comp_id,
        auth_seq_id: seq_id,
        Cartn_x: chem_comp_atom[cartn_x],
        Cartn_y: chem_comp_atom[cartn_y],
        Cartn_z: chem_comp_atom[cartn_z],
        id: pdbx_ordinal,

        label_asym_id: A,
        label_atom_id: type_symbol,
        label_comp_id: comp_id,
        label_seq_id: seq_id,
        label_entity_id: entity_id,

        occupancy,
        type_symbol,

        pdbx_PDB_model_num: model_num,
        pdbx_formal_charge: charge
    }, atomCount);

    const entityBuilder = new EntityBuilder();
    entityBuilder.setNames([['MOL', `${(name || 'Unknown Entity')} ${suffix}`]]);
    entityBuilder.getEntityId('MOL', MoleculeType.Unknown, 'A');

    const componentBuilder = new ComponentBuilder(seq_id, type_symbol);
    componentBuilder.setNames([['MOL', `${(name || 'Unknown Molecule')} ${suffix}`]]);
    componentBuilder.add('MOL', 0);

    const basicModel = createBasic({
        entity: entityBuilder.getEntityTable(),
        atom_site: model_atom_site
    });
    const modelsModel = await createModels(basicModel, format, ctx);

    if (modelsModel.frameCount > 0) {
        const first = modelsModel.representative;

        const bondCount = chem_comp_bond._rowCount;
        if (bondCount > 0) {
            const labelIndexMap: { [label: string]: number } = {};
            const { atom_id } = chem_comp_atom;
            for (let i = 0, il = atom_id.rowCount; i < il; ++i) {
                labelIndexMap[atom_id.value(i)] = i;
            }

            const indexA: number[] = [];
            const indexB: number[] = [];
            const order: number[] = [];
            const flag: number[] = [];

            const included = new Set<number>();
            let j = 0;

            const { atom_id_1, atom_id_2, pdbx_aromatic_flag, value_order } = chem_comp_bond;
            for (let i = 0; i < bondCount; ++i) {
                const iA = labelIndexMap[atom_id_1.value(i)];
                const iB = labelIndexMap[atom_id_2.value(i)];
                const id = iA < iB ? cantorPairing(iA, iB) : cantorPairing(iB, iA);
                if (included.has(id)) continue;
                included.add(id);

                indexA[j] = iA;
                indexB[j] = iB;

                let flags: number = BondType.Flag.Covalent;
                let ord = 1;
                if (pdbx_aromatic_flag.value(i) === 'y') flags |= BondType.Flag.Aromatic;
                switch (value_order.value(i)) {
                    case 'delo': flags |= BondType.Flag.Aromatic; break;
                    case 'doub': ord = 2; break;
                    case 'trip': ord = 3; break;
                    case 'quad': ord = 4; break;
                }
                order[j] = ord;
                flag[j] = flags;

                j += 1;
            }

            IndexPairBonds.Provider.set(first, IndexPairBonds.fromData({ pairs: {
                indexA: Column.ofIntArray(indexA),
                indexB: Column.ofIntArray(indexB),
                order: Column.ofIntArray(order),
                flag: Column.ofIntArray(flag)
            }, count: atomCount }));
        }
    }

    return modelsModel;
}
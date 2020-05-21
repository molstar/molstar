/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model } from '../../mol-model/structure/model/model';
import { Task } from '../../mol-task';
import { ModelFormat } from '../format';
import { CifFrame, CIF } from '../../mol-io/reader/cif';
import { mmCIF_Database } from '../../mol-io/reader/cif/schema/mmcif';
import { createModels } from './basic/parser';
import { ModelSymmetry } from './property/symmetry';
import { ModelSecondaryStructure } from './property/secondary-structure';
import { Table } from '../../mol-data/db';
import { AtomSiteAnisotrop } from './property/anisotropic';
import { ComponentBond } from './property/bonds/comp';
import { StructConn } from './property/bonds/struct_conn';

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

//

export { MmcifFormat };

type MmcifFormat = ModelFormat<MmcifFormat.Data>

namespace MmcifFormat {
    export type Data = { db: mmCIF_Database, frame: CifFrame }
    export function is(x: ModelFormat): x is MmcifFormat {
        return x.kind === 'mmCIF';
    }

    export function fromFrame(frame: CifFrame, db?: mmCIF_Database): MmcifFormat {
        if (!db) db = CIF.schema.mmCIF(frame);
        return { kind: 'mmCIF', name: db._name, data: { db, frame } };
    }
}

export function trajectoryFromMmCIF(frame: CifFrame): Task<Model.Trajectory> {
    const format = MmcifFormat.fromFrame(frame);
    return Task.create('Create mmCIF Model', ctx => createModels(format.data.db, format, ctx));
}
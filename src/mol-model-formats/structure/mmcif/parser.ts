/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Table } from '../../../mol-data/db';
import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';
import { RuntimeContext } from '../../../mol-task';
import { Model } from '../../../mol-model/structure/model/model';
import { ModelFormat } from '../format';
import mmCIF_Format = ModelFormat.mmCIF
import { AtomSiteAnisotrop } from '../property/anisotropic';
import { _parse_basic } from '../basic/parser';
import { ModelSymmetry } from '../property/symmetry';
import { ModelSecondaryStructure } from '../property/secondary-structure';
import { ComponentBond } from '../property/bonds/comp';
import { StructConn } from '../property/bonds/struct_conn';
import { ModelCrossLinkRestraint } from '../property/pair-restraints/cross-links';

export async function _parse_mmCif(format: mmCIF_Format, ctx: RuntimeContext) {
    return _parse_basic(format.data, format, ctx)
}

function modelSymmetryFromMmcif(model: Model) {
    if (model.sourceData.kind !== 'mmCIF') return;
    return ModelSymmetry.fromData(model.sourceData.data)
}
ModelSymmetry.Provider.formatRegistry.add('mmCIF', modelSymmetryFromMmcif)

function secondaryStructureFromMmcif(model: Model) {
    if (model.sourceData.kind !== 'mmCIF') return;
    const { struct_conf, struct_sheet_range } = model.sourceData.data
    return ModelSecondaryStructure.fromStruct(struct_conf, struct_sheet_range, model.atomicHierarchy)
}
ModelSecondaryStructure.Provider.formatRegistry.add('mmCIF', secondaryStructureFromMmcif)

function atomSiteAnisotropFromMmcif(model: Model) {
    if (model.sourceData.kind !== 'mmCIF') return;
    const { atom_site_anisotrop } = model.sourceData.data
    const data = Table.ofColumns(mmCIF_Schema['atom_site_anisotrop'], atom_site_anisotrop);
    const elementToAnsiotrop = AtomSiteAnisotrop.getElementToAnsiotrop(model, data)
    return { data, elementToAnsiotrop }
}
AtomSiteAnisotrop.Provider.formatRegistry.add('mmCIF', atomSiteAnisotropFromMmcif)

function componentBondFromMmcif(model: Model) {
    if (model.sourceData.kind !== 'mmCIF') return;
    const { chem_comp_bond } = model.sourceData.data;
    if (chem_comp_bond._rowCount === 0) return;
    return {
        data: chem_comp_bond,
        entries: ComponentBond.getEntriesFromChemCompBond(chem_comp_bond)
    }
}
ComponentBond.Provider.formatRegistry.add('mmCIF', componentBondFromMmcif)

function structConnFromMmcif(model: Model) {
    if (model.sourceData.kind !== 'mmCIF') return;
    const { struct_conn } = model.sourceData.data;
    if (struct_conn._rowCount === 0) return;
    const entries = StructConn.getEntriesFromStructConn(struct_conn, model)
    return {
        data: struct_conn,
        byAtomIndex: StructConn.getAtomIndexFromEntries(entries),
        entries,
    }
}
StructConn.Provider.formatRegistry.add('mmCIF', structConnFromMmcif)

function crossLinkRestraintFromMmcif(model: Model) {
    if (model.sourceData.kind !== 'mmCIF') return;
    const { ihm_cross_link_restraint } = model.sourceData.data;
    if (ihm_cross_link_restraint._rowCount === 0) return;
    return ModelCrossLinkRestraint.fromTable(ihm_cross_link_restraint, model)
}
ModelCrossLinkRestraint.Provider.formatRegistry.add('mmCIF', crossLinkRestraintFromMmcif)
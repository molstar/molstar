/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';
import { Table } from '../../../mol-data/db';

// TODO split into conformation and hierarchy parts

export type Entry = Table<mmCIF_Schema['entry']>
export type Struct = Table<mmCIF_Schema['struct']>
export type StructAsym = Table<mmCIF_Schema['struct_asym']>
export type IhmModelList = Table<mmCIF_Schema['ihm_model_list']>
export type IhmModelGroup = Table<mmCIF_Schema['ihm_model_group']>
export type IhmModelGroupLink = Table<mmCIF_Schema['ihm_model_group_link']>
export type Entity = Table<mmCIF_Schema['entity']>
export type EntityPoly = Table<mmCIF_Schema['entity_poly']>
export type EntityPolySeq = Table<mmCIF_Schema['entity_poly_seq']>
export type EntityBranch = Table<mmCIF_Schema['pdbx_entity_branch']>
export type ChemComp = Table<mmCIF_Schema['chem_comp']>
export type ChemCompIdentifier = Table<mmCIF_Schema['pdbx_chem_comp_identifier']>
export type AtomSite = Table<mmCIF_Schema['atom_site']>
export type IhmSphereObjSite = Table<mmCIF_Schema['ihm_sphere_obj_site']>
export type IhmGaussianObjSite =Table<mmCIF_Schema['ihm_gaussian_obj_site']>
export type UnobsOrZeroOccResidues =Table<mmCIF_Schema['pdbx_unobs_or_zero_occ_residues']>
export type Molecule =Table<mmCIF_Schema['pdbx_molecule']>

export const BasicSchema = {
    entry: mmCIF_Schema.entry,
    struct: mmCIF_Schema.struct,
    struct_asym: mmCIF_Schema.struct_asym,
    ihm_model_list: mmCIF_Schema.ihm_model_list,
    ihm_model_group: mmCIF_Schema.ihm_model_group,
    ihm_model_group_link: mmCIF_Schema.ihm_model_group_link,
    entity: mmCIF_Schema.entity,
    entity_poly: mmCIF_Schema.entity_poly,
    entity_poly_seq: mmCIF_Schema.entity_poly_seq,
    pdbx_entity_branch: mmCIF_Schema.pdbx_entity_branch,
    chem_comp: mmCIF_Schema.chem_comp,
    pdbx_chem_comp_identifier: mmCIF_Schema.pdbx_chem_comp_identifier,
    atom_site: mmCIF_Schema.atom_site,
    ihm_sphere_obj_site: mmCIF_Schema.ihm_sphere_obj_site,
    ihm_gaussian_obj_site: mmCIF_Schema.ihm_gaussian_obj_site,
    pdbx_unobs_or_zero_occ_residues: mmCIF_Schema.pdbx_unobs_or_zero_occ_residues,
    pdbx_molecule: mmCIF_Schema.pdbx_molecule,
};

export interface BasicData {
    entry: Entry
    struct: Struct
    struct_asym: StructAsym
    ihm_model_list: IhmModelList
    ihm_model_group: IhmModelGroup
    ihm_model_group_link: IhmModelGroupLink
    entity: Entity
    entity_poly: EntityPoly
    entity_poly_seq: EntityPolySeq
    pdbx_entity_branch: EntityBranch
    chem_comp: ChemComp
    pdbx_chem_comp_identifier: ChemCompIdentifier
    atom_site: AtomSite
    ihm_sphere_obj_site: IhmSphereObjSite
    ihm_gaussian_obj_site: IhmGaussianObjSite
    pdbx_unobs_or_zero_occ_residues: UnobsOrZeroOccResidues
    pdbx_molecule: Molecule
}

export function createBasic(data: Partial<BasicData>): BasicData {
    const basic = Object.create(null);
    for (const name of Object.keys(BasicSchema)) {
        if (name in data) {
            basic[name] = data[name as keyof typeof BasicSchema];
        } else {
            basic[name] = Table.ofUndefinedColumns(BasicSchema[name as keyof typeof BasicSchema], 0);
        }
    }
    return basic;
}
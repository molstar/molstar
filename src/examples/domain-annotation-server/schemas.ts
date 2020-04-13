/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '../../mol-data/db';

import Type = Column.Schema

export const Sources = {
    id: Type.str,
    count: Type.int
};
export type Sources = typeof Sources

export const Base = {
    id: Type.str,
    identifier: Type.str,
    mapping_group_id: Type.int
};
export type Base = typeof Base

export const mapping = {
    id: Type.int,
    group_id: Type.int,

    label_entity_id: Type.str,
    label_asym_id: Type.str,
    auth_asym_id: Type.str,

    beg_label_seq_id: Type.int,
    beg_auth_seq_id: Type.int,
    pdbx_beg_PDB_ins_code: Type.str,

    end_label_seq_id: Type.int,
    end_auth_seq_id: Type.int,
    pdbx_end_PDB_ins_code: Type.str
};
export type mapping = typeof mapping

export const Pfam = {
    description: Type.str
};
export type Pfam = typeof Pfam

export const InterPro = {
    name: Type.str
};
export type InterPro = typeof InterPro

export const CATH = {
    name: Type.str,
    homology: Type.str,
    architecture: Type.str,
    identifier: Type.str,
    class: Type.str,
    topology: Type.str,
};
export type CATH = typeof CATH

export const EC = {
    accepted_name: Type.str,
    reaction: Type.str,
    systematic_name: Type.str
};
export type EC = typeof EC

export const UniProt = {
    name: Type.str
};
export type UniProt = typeof UniProt

export const SCOP = {
    sccs: Type.str,
    description: Type.str
};
export type SCOP = typeof SCOP

export const GO = {
    category: Type.str,
    definition: Type.str,
    name: Type.str
};
export type GO = typeof GO

export const categories = {
    Pfam,
    InterPro,
    CATH,
    EC,
    UniProt,
    SCOP,
    GO
};
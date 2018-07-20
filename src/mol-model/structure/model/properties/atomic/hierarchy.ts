/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, Table } from 'mol-data/db'
import { Segmentation } from 'mol-data/int'
import { mmCIF_Schema as mmCIF } from 'mol-io/reader/cif/schema/mmcif'
import { ElementSymbol } from '../../types'
import { ChainIndex, EntityIndex, ResidueIndex, ElementIndex } from '../../indexing';
import SortedRanges from 'mol-data/int/sorted-ranges';

export const AtomsSchema = {
    type_symbol: Column.Schema.Aliased<ElementSymbol>(mmCIF.atom_site.type_symbol),
    label_atom_id: mmCIF.atom_site.label_atom_id,
    auth_atom_id: mmCIF.atom_site.auth_atom_id,
    label_alt_id: mmCIF.atom_site.label_alt_id,
    pdbx_formal_charge: mmCIF.atom_site.pdbx_formal_charge
    // id, occupancy and B_iso_or_equiv are part of conformation
};

export type AtomsSchema = typeof AtomsSchema
export interface Atoms extends Table<AtomsSchema> { }

export const ResiduesSchema = {
    group_PDB: mmCIF.atom_site.group_PDB,
    label_comp_id: mmCIF.atom_site.label_comp_id,
    auth_comp_id: mmCIF.atom_site.auth_comp_id,
    label_seq_id: mmCIF.atom_site.label_seq_id,
    auth_seq_id: mmCIF.atom_site.auth_seq_id,
    pdbx_PDB_ins_code: mmCIF.atom_site.pdbx_PDB_ins_code
};
export type ResiduesSchema = typeof ResiduesSchema
export interface Residues extends Table<ResiduesSchema> { }

export const ChainsSchema = {
    label_asym_id: mmCIF.atom_site.label_asym_id,
    auth_asym_id: mmCIF.atom_site.auth_asym_id,
    label_entity_id: mmCIF.atom_site.label_entity_id
}
export type ChainsSchema = typeof ChainsSchema
export interface Chains extends Table<ChainsSchema> { }

export interface AtomicData {
    atoms: Atoms,
    residues: Residues,
    chains: Chains
}

export interface AtomicSegments {
    /** Maps residueIndex to a range of atoms [segments[rI], segments[rI + 1]) */
    residueAtomSegments: Segmentation<ElementIndex, ResidueIndex>,
    /**
     * Maps chainIndex to a range of atoms [segments[cI], segments[cI + 1]),
     *
     * residues of i-th chain are accessed like this:
     * const rI = residueAtomSegments.index, offsets = chainAtomSegments.offsets;
     * const start = rI[offsets[i]], const end = rI[offsets[i + 1] - 1] + 1;
     * for (let j = start; j < end; i++) { }
     */
    chainAtomSegments: Segmentation<ElementIndex, ChainIndex>,

    // TODO: include entity segments?
}

export interface AtomicKeys {
    // TODO: include (lazily computed) "entity/chain/residue" indices?

    /** @returns index or -1 if not present. */
    getEntityKey(cI: ChainIndex): EntityIndex,

    /** @returns index or -1 if not present. */
    findChainKey(entityId: string, label_asym_id: string): ChainIndex,

    /**
     * Unique number for each of the residue. Also the index of the 1st occurence of this residue.
     * @returns index or -1 if not present.
     */
    findResidueKey(entityId: string, label_asym_id: string, label_comp_id: string, auth_seq_id: number, pdbx_PDB_ins_code: string): ResidueIndex
}

export interface AtomicRanges {
    polymerRanges: SortedRanges<ElementIndex>
    gapRanges: SortedRanges<ElementIndex>
    cyclicPolymerMap: Map<ResidueIndex, ResidueIndex>
}

type _Hierarchy = AtomicData & AtomicSegments & AtomicKeys & AtomicRanges
export interface AtomicHierarchy extends _Hierarchy { }

export namespace AtomicHierarchy {
    /** Start residue inclusive */
    export function chainStartResidueIndex(segs: AtomicSegments, cI: ChainIndex) {
        return segs.residueAtomSegments.index[segs.chainAtomSegments.offsets[cI]];
    }

    /** End residue exclusive */
    export function chainEndResidueIndexExcl(segs: AtomicSegments, cI: ChainIndex) {
        return segs.residueAtomSegments.index[segs.chainAtomSegments.offsets[cI + 1] - 1] + 1 as ResidueIndex;
    }
}

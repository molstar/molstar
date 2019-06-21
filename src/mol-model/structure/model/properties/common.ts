/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mmCIF_Database, mmCIF_Schema } from '../../../../mol-io/reader/cif/schema/mmcif'
import { Table } from '../../../../mol-data/db';
import { EntityIndex } from '../indexing';

export interface Entities {
    data: mmCIF_Database['entity'],
    getEntityIndex(id: string): EntityIndex
}

export type ChemicalComponent = Table.Row<mmCIF_Schema['chem_comp']>
export type ChemicalComponentMap = ReadonlyMap<string, ChemicalComponent>

export type MissingResidue = Table.Row<Pick<
    mmCIF_Schema['pdbx_unobs_or_zero_occ_residues'],
    'polymer_flag' | 'occupancy_flag'>
>
export interface MissingResidues {
    has(model_num: number, asym_id: string, seq_id: number): boolean
    get(model_num: number, asym_id: string, seq_id: number): MissingResidue | undefined
    readonly size: number
}
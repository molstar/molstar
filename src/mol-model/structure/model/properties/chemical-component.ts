/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { mmCIF_Schema } from '../../../../mol-io/reader/cif/schema/mmcif';
import { Table } from '../../../../mol-data/db';

export type ChemicalComponent = Table.Row<mmCIF_Schema['chem_comp']>
export type ChemicalComponentMap = ReadonlyMap<string, ChemicalComponent>

// TODO add data for common chemical components
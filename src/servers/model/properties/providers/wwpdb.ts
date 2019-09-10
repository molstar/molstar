/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as fs from 'fs'
import * as util from 'util'
import { AttachModelProperty } from '../../property-provider';
import { CIF } from '../../../../mol-io/reader/cif';
import { getParam } from '../../../common/util';
import { ComponentBond } from '../../../../mol-model-formats/structure/mmcif/bonds';
import { mmCIF_Database, mmCIF_Schema } from '../../../../mol-io/reader/cif/schema/mmcif';

require('util.promisify').shim()
const readFile = util.promisify(fs.readFile)

export const wwPDB_chemCompBond: AttachModelProperty = async ({ model, params }) => {
    const table = await getChemCompBondTable(getTablePath(params))
    return ComponentBond.attachFromExternalData(model, table, true)
}

async function read(path: string) {
    return path.endsWith('.bcif') ? new Uint8Array(await readFile(path)) : readFile(path, 'utf8');
}

let chemCompBondTable: mmCIF_Database['chem_comp_bond']
async function getChemCompBondTable(path: string): Promise<mmCIF_Database['chem_comp_bond']> {
    if (!chemCompBondTable) {
        const parsed = await CIF.parse(await read(path)).run()
        if (parsed.isError) throw new Error(parsed.toString())
        const table = CIF.toDatabase(mmCIF_Schema, parsed.result.blocks[0])
        chemCompBondTable = table.chem_comp_bond
    }
    return chemCompBondTable
}

function getTablePath(params: any) {
    const path = getParam<string>(params, 'wwPDB', 'chemCompBondTablePath');
    if (!path) throw new Error(`wwPDB 'chemCompBondTablePath' not set!`);
    return path;
}
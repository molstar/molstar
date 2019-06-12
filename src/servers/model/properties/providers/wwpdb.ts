/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as fs from 'fs'
import * as util from 'util'
import { AttachModelProperty } from '../../property-provider';
import { ChemCompBond } from '../../../../mol-model-props/wwpdb/chem-comp-bond';
import { Table } from '../../../../mol-data/db';
import { CIF } from '../../../../mol-io/reader/cif';
import { getParam } from '../../../common/util';

require('util.promisify').shim()
const readFile = util.promisify(fs.readFile)

export const wwPDB_chemCompBond: AttachModelProperty = ({ model, params }) => {
    const wwPDB_apiSourceTable = getChemCompBondTableProvider(getTablePath(params))
    return ChemCompBond.attachFromCifOrTable(model, { wwPDB_apiSourceTable });
}

async function read(path: string) {
    return path.endsWith('.bcif') ? new Uint8Array(await readFile(path)) : readFile(path, 'utf8');
}

function getChemCompBondTableProvider(path: string): () => Promise<Table<ChemCompBond.Schema['chem_comp_bond']>> {
    let chemCompBondTable: Table<ChemCompBond.Schema['chem_comp_bond']>
    return async function() {
        if (chemCompBondTable === undefined) {
            const parsed = await CIF.parse(await read(path)).run()
            if (parsed.isError) throw new Error(parsed.toString())
            const table = CIF.toDatabase(ChemCompBond.Schema, parsed.result.blocks[0])
            chemCompBondTable = table.chem_comp_bond
        }
        return chemCompBondTable
    }
}

function getTablePath(params: any) {
    const path = getParam<string>(params, 'wwPDB', 'chemCompBondTablePath');
    if (!path) throw new Error(`wwPDB 'chemCompBondTablePath' not set!`);
    return path;
}
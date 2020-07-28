/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as fs from 'fs';
import * as util from 'util';
import { AttachModelProperty } from '../../property-provider';
import { CIF } from '../../../../mol-io/reader/cif';
import { getParam } from '../../../common/util';
import { mmCIF_Database, mmCIF_Schema } from '../../../../mol-io/reader/cif/schema/mmcif';
import { ComponentBond } from '../../../../mol-model-formats/structure/property/bonds/chem_comp';
import { ComponentAtom } from '../../../../mol-model-formats/structure/property/atoms/chem_comp';
import { CCD_Database, CCD_Schema } from '../../../../mol-io/reader/cif/schema/ccd';

require('util.promisify').shim();
const readFile = util.promisify(fs.readFile);

export const wwPDB_chemCompBond: AttachModelProperty = async ({ model, params }) => {
    const table = await getChemCompBondTable(getBondTablePath(params));
    const data = ComponentBond.chemCompBondFromTable(model, table);
    const entries = ComponentBond.getEntriesFromChemCompBond(data);
    return ComponentBond.Provider.set(model, { entries, data });
};

async function read(path: string) {
    return path.endsWith('.bcif') ? new Uint8Array(await readFile(path)) : readFile(path, 'utf8');
}

let chemCompBondTable: mmCIF_Database['chem_comp_bond'];
async function getChemCompBondTable(path: string): Promise<mmCIF_Database['chem_comp_bond']> {
    if (!chemCompBondTable) {
        const parsed = await CIF.parse(await read(path)).run();
        if (parsed.isError) throw new Error(parsed.toString());
        const table = CIF.toDatabase(mmCIF_Schema, parsed.result.blocks[0]);
        chemCompBondTable = table.chem_comp_bond;
    }
    return chemCompBondTable;
}

function getBondTablePath(params: any) {
    const path = getParam<string>(params, 'wwPDB', 'chemCompBondTablePath');
    if (!path) throw new Error(`wwPDB 'chemCompBondTablePath' not set!`);
    return path;
}

export const wwPDB_chemCompAtom: AttachModelProperty = async ({ model, params }) => {
    const table = await getChemCompAtomTable(getAtomTablePath(params));
    const data = ComponentAtom.chemCompAtomFromTable(model, table);
    const entries = ComponentAtom.getEntriesFromChemCompAtom(data);
    return ComponentAtom.Provider.set(model, { entries, data });
};

let chemCompAtomTable: CCD_Database['chem_comp_atom'];
async function getChemCompAtomTable(path: string): Promise<CCD_Database['chem_comp_atom']> {
    if (!chemCompAtomTable) {
        const parsed = await CIF.parse(await read(path)).run();
        if (parsed.isError) throw new Error(parsed.toString());
        const table = CIF.toDatabase(CCD_Schema, parsed.result.blocks[0]);
        chemCompAtomTable = table.chem_comp_atom;
    }
    return chemCompAtomTable;
}

function getAtomTablePath(params: any) {
    const path = getParam<string>(params, 'wwPDB', 'chemCompAtomTablePath');
    if (!path) throw new Error(`wwPDB 'chemCompAtomTablePath' not set!`);
    return path;
}
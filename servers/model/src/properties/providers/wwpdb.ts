/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import * as fs from 'fs';
import type { AttachModelProperty } from '@molstar/model-server/property-provider';
import { CIF } from '@molstar/io/reader/cif';
import { getParam } from '@molstar/common-server/util';
import { type mmCIF_Database, mmCIF_Schema } from '@molstar/io/reader/cif/schema/mmcif';
import { ComponentBond } from '@molstar/model/formats/structure/property/bonds/chem_comp';
import { ComponentAtom } from '@molstar/model/formats/structure/property/atoms/chem_comp';
import { type CCD_Database, CCD_Schema } from '@molstar/io/reader/cif/schema/ccd';

const readFileAsync = fs.promises.readFile;

export const wwPDB_chemCompBond: AttachModelProperty = async ({ model, params }) => {
  const table = await getChemCompBondTable(getBondTablePath(params));
  const data = ComponentBond.chemCompBondFromTable(model, table);
  const entries = ComponentBond.getEntriesFromChemCompBond(data);
  return ComponentBond.Provider.set(model, { entries, data });
};

async function read(path: string) {
  return path.endsWith('.bcif') ? new Uint8Array(await readFileAsync(path)) : readFileAsync(path, 'utf8');
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

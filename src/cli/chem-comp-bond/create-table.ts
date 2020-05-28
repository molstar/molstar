#!/usr/bin/env node
/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as argparse from 'argparse';
import * as util from 'util';
import * as path from 'path';
import * as fs from 'fs';
import * as zlib from 'zlib';
import fetch from 'node-fetch';
require('util.promisify').shim();
const readFile = util.promisify(fs.readFile);
const writeFile = util.promisify(fs.writeFile);

import { Progress } from '../../mol-task';
import { Database, Table, DatabaseCollection } from '../../mol-data/db';
import { CIF } from '../../mol-io/reader/cif';
import { CifWriter } from '../../mol-io/writer/cif';
import { CCD_Schema } from '../../mol-io/reader/cif/schema/ccd';
import { SetUtils } from '../../mol-util/set';
import { DefaultMap } from '../../mol-util/map';
import { mmCIF_chemCompBond_schema } from '../../mol-io/reader/cif/schema/mmcif-extras';

export async function ensureAvailable(path: string, url: string) {
    if (FORCE_DOWNLOAD || !fs.existsSync(path)) {
        console.log(`downloading ${url}...`);
        const data = await fetch(url);
        if (!fs.existsSync(DATA_DIR)) {
            fs.mkdirSync(DATA_DIR);
        }
        if (url.endsWith('.gz')) {
            await writeFile(path, zlib.gunzipSync(await data.buffer()));
        } else {
            await writeFile(path, await data.text());
        }
        console.log(`done downloading ${url}`);
    }
}

export async function ensureDataAvailable() {
    await ensureAvailable(CCD_PATH, CCD_URL);
    await ensureAvailable(PVCD_PATH, PVCD_URL);
}

export async function readFileAsCollection<S extends Database.Schema>(path: string, schema: S) {
    const parsed = await parseCif(await readFile(path, 'utf8'));
    return CIF.toDatabaseCollection(schema, parsed.result);
}

export async function readCCD() {
    return readFileAsCollection(CCD_PATH, CCD_Schema);
}

export async function readPVCD() {
    return readFileAsCollection(PVCD_PATH, CCD_Schema);
}

async function parseCif(data: string | Uint8Array) {
    const comp = CIF.parse(data);
    console.time('parse cif');
    const parsed = await comp.run(p => console.log(Progress.format(p)), 250);
    console.timeEnd('parse cif');
    if (parsed.isError) throw parsed;
    return parsed;
}

export function getEncodedCif(name: string, database: Database<Database.Schema>, binary = false) {
    const encoder = CifWriter.createEncoder({ binary, encoderName: 'mol*' });
    CifWriter.Encoder.writeDatabase(encoder, name, database);
    return encoder.getData();
}

type CCB = Table<CCD_Schema['chem_comp_bond']>
type CCA = Table<CCD_Schema['chem_comp_atom']>

function ccbKey(compId: string, atomId1: string, atomId2: string) {
    return atomId1 < atomId2 ? `${compId}:${atomId1}-${atomId2}` : `${compId}:${atomId2}-${atomId1}`;
}

function addChemCompBondToSet(set: Set<string>, ccb: CCB) {
    for (let i = 0, il = ccb._rowCount; i < il; ++i) {
        set.add(ccbKey(ccb.comp_id.value(i), ccb.atom_id_1.value(i), ccb.atom_id_2.value(i)));
    }
    return set;
}

function addChemCompAtomToSet(set: Set<string>, cca: CCA) {
    for (let i = 0, il = cca._rowCount; i < il; ++i) {
        set.add(cca.atom_id.value(i));
    }
    return set;
}

function checkAddingBondsFromPVCD(pvcd: DatabaseCollection<CCD_Schema>) {
    const ccbSetByParent = DefaultMap<string, Set<string>>(() => new Set());

    for (const k in pvcd) {
        const { chem_comp, chem_comp_bond } = pvcd[k];
        if (chem_comp_bond._rowCount) {
            const parentIds = chem_comp.mon_nstd_parent_comp_id.value(0);
            if (parentIds.length === 0) {
                const set = ccbSetByParent.getDefault(chem_comp.id.value(0));
                addChemCompBondToSet(set, chem_comp_bond);
            } else {
                for (let i = 0, il = parentIds.length; i < il; ++i) {
                    const parentId = parentIds[i];
                    const set = ccbSetByParent.getDefault(parentId);
                    addChemCompBondToSet(set, chem_comp_bond);
                }
            }
        }
    }

    for (const k in pvcd) {
        const { chem_comp, chem_comp_atom, chem_comp_bond } = pvcd[k];
        if (chem_comp_bond._rowCount) {
            const parentIds = chem_comp.mon_nstd_parent_comp_id.value(0);
            if (parentIds.length > 0) {
                for (let i = 0, il = parentIds.length; i < il; ++i) {
                    const entryBonds = addChemCompBondToSet(new Set<string>(), chem_comp_bond);
                    const entryAtoms = addChemCompAtomToSet(new Set<string>(), chem_comp_atom);
                    const extraBonds = SetUtils.difference(ccbSetByParent.get(parentIds[i])!, entryBonds);
                    extraBonds.forEach(bk => {
                        const [a1, a2] = bk.split('|');
                        if (entryAtoms.has(a1) && entryAtoms.has(a2)) {
                            console.error(`Adding all PVCD bonds would wrongly add bond ${bk} for ${k}`);
                        }
                    });
                }
            }
        }
    }
}

async function createBonds() {
    await ensureDataAvailable();
    const ccd = await readCCD();
    const pvcd = await readPVCD();

    const ccbSet = new Set<string>();

    const comp_id: string[] = [];
    const atom_id_1: string[] = [];
    const atom_id_2: string[] = [];
    const value_order: typeof mmCIF_chemCompBond_schema['value_order']['T'][] = [];
    const pdbx_aromatic_flag: typeof mmCIF_chemCompBond_schema['pdbx_aromatic_flag']['T'][] = [];
    const pdbx_stereo_config: typeof mmCIF_chemCompBond_schema['pdbx_stereo_config']['T'][] = [];
    const molstar_protonation_variant: string[] = [];

    function addBonds(compId: string, ccb: CCB, protonationVariant: boolean) {
        for (let i = 0, il = ccb._rowCount; i < il; ++i) {
            const atomId1 = ccb.atom_id_1.value(i);
            const atomId2 = ccb.atom_id_2.value(i);
            const k = ccbKey(compId, atomId1, atomId2);
            if (!ccbSet.has(k)) {
                atom_id_1.push(atomId1);
                atom_id_2.push(atomId2);
                comp_id.push(compId);
                value_order.push(ccb.value_order.value(i));
                pdbx_aromatic_flag.push(ccb.pdbx_aromatic_flag.value(i));
                pdbx_stereo_config.push(ccb.pdbx_stereo_config.value(i));
                molstar_protonation_variant.push(protonationVariant ? 'Y' : 'N');
                ccbSet.add(k);
            }
        }
    }

    // check adding bonds from PVCD
    checkAddingBondsFromPVCD(pvcd);

    // add bonds from PVCD
    for (const k in pvcd) {
        const { chem_comp, chem_comp_bond } = pvcd[k];
        if (chem_comp_bond._rowCount) {
            const parentIds = chem_comp.mon_nstd_parent_comp_id.value(0);
            if (parentIds.length === 0) {
                addBonds(chem_comp.id.value(0), chem_comp_bond, false);
            } else {
                for (let i = 0, il = parentIds.length; i < il; ++i) {
                    addBonds(parentIds[i], chem_comp_bond, true);
                }
            }
        }
    }

    // add bonds from CCD
    for (const k in ccd) {
        const { chem_comp, chem_comp_bond } = ccd[k];
        if (chem_comp_bond._rowCount) {
            addBonds(chem_comp.id.value(0), chem_comp_bond, false);
        }
    }

    const bondTable = Table.ofArrays(mmCIF_chemCompBond_schema, {
        comp_id, atom_id_1, atom_id_2, value_order,
        pdbx_aromatic_flag, pdbx_stereo_config, molstar_protonation_variant
    });

    const bondDatabase =  Database.ofTables(
        TABLE_NAME,
        { chem_comp_bond: mmCIF_chemCompBond_schema },
        { chem_comp_bond: bondTable }
    );

    return bondDatabase;
}

async function run(out: string, binary = false) {
    const bonds = await createBonds();

    const cif = getEncodedCif(TABLE_NAME, bonds, binary);
    if (!fs.existsSync(path.dirname(out))) {
        fs.mkdirSync(path.dirname(out));
    }
    writeFile(out, cif);
}

const TABLE_NAME = 'CHEM_COMP_BONDS';

const DATA_DIR = path.join(__dirname, '..', '..', '..', '..', 'build/data');
const CCD_PATH = path.join(DATA_DIR, 'components.cif');
const PVCD_PATH = path.join(DATA_DIR, 'aa-variants-v1.cif');
const CCD_URL = 'http://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif';
const PVCD_URL = 'http://ftp.wwpdb.org/pub/pdb/data/monomers/aa-variants-v1.cif';

const parser = new argparse.ArgumentParser({
    addHelp: true,
    description: 'Create a cif file with one big table of all chem_comp_bond entries from the CCD and PVCD.'
});
parser.addArgument('out', {
    help: 'Generated file output path.'
});
parser.addArgument([ '--forceDownload', '-f' ], {
    action: 'storeTrue',
    help: 'Force download of CCD and PVCD.'
});
parser.addArgument([ '--binary', '-b' ], {
    action: 'storeTrue',
    help: 'Output as BinaryCIF.'
});
interface Args {
    out: string
    forceDownload?: boolean
    binary?: boolean
}
const args: Args = parser.parseArgs();

const FORCE_DOWNLOAD = args.forceDownload;

run(args.out, args.binary);

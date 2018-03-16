/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as argparse from 'argparse'
import * as util from 'util'
import * as fs from 'fs'
import fetch from 'node-fetch'
require('util.promisify').shim();
const readFileAsync = util.promisify(fs.readFile);

import { DatabaseCollection, Database, Table } from 'mol-data/db'
import CIF from 'mol-io/reader/cif'
// import { CCD_Schema } from 'mol-io/reader/cif/schema/ccd'
import * as Encoder from 'mol-io/writer/cif'
import Computation from 'mol-util/computation'
import { mmCIF_Schema } from 'mol-io/reader/cif/schema/mmcif';
import { CCD_Schema } from 'mol-io/reader/cif/schema/ccd';
// import { Table } from 'mol-io/reader/csv/data-model';

// import { Model, Structure } from 'mol-model/structure'
// import to_mmCIF from 'mol-model/structure/export/mmcif'

export async function ensureAvailable(path: string, url: string) {
    if (FORCE_DOWNLOAD || !fs.existsSync(path)) {
        console.log(`downloading ${url}...`)
        const data = await fetch(url)
        if (!fs.existsSync(DATA_DIR)) {
            fs.mkdirSync(DATA_DIR);
        }
        fs.writeFileSync(path, await data.text())
        console.log(`done downloading ${url}`)
    }
}

export async function ensureDataAvailable() {
    await ensureAvailable(CCD_PATH, CCD_URL)
    await ensureAvailable(PVCD_PATH, PVCD_URL)
}

function showProgress(tag: string, p: Computation.Progress) {
    console.log(`[${tag}] ${p.message} ${p.isIndeterminate ? '' : (p.current / p.max * 100).toFixed(2) + '% '}(${p.elapsedMs | 0}ms)`)
}

export async function readCCD(bcif = false) {
    const parsed = await parseCif(await readFileAsync(CCD_PATH, 'utf8'))
    const ccd: DatabaseCollection<CCD_Schema> = {}
    for (const data of parsed.result.blocks) {
        // console.log(data.header)
        ccd[data.header] = CIF.schema.CCD(data)
    }
    return ccd;
}

export async function readPVCD(bcif = false) {
    const parsed = await parseCif(await readFileAsync(PVCD_PATH, 'utf8'))
    const ccd: DatabaseCollection<CCD_Schema> = {}
    for (const data of parsed.result.blocks) {
        // console.log(data.header)
        ccd[data.header] = CIF.schema.CCD(data)
    }
    return ccd;
}

export async function getCCD() {
    await ensureDataAvailable()
    return readPVCD()
}

async function parseCif(data: string|Uint8Array) {
    const comp = CIF.parse(data)
    const ctx = Computation.observable({
        updateRateMs: 250,
        observer: p => showProgress(`cif parser ${typeof data === 'string' ? 'string' : 'binary'}`, p)
    });
    console.time('parse cif')
    const parsed = await comp(ctx);
    console.timeEnd('parse cif')
    if (parsed.isError) throw parsed;
    return parsed
}

export function getEncodedCif(name: string, database: Database<Database.Schema>, binary = false) {
    const encoder = Encoder.create({ binary, encoderName: 'mol*' });
    Encoder.writeDatabase(encoder, name, database)
    return encoder.getData();
}

export async function getPdb(pdb: string) {
    console.log(`downloading ${pdb}...`)
    const data = await fetch(`https://files.rcsb.org/download/${pdb}.cif`)
    console.log(`done downloading ${pdb}`)

    const parsed = await parseCif(await data.text())
    return CIF.schema.mmCIF(parsed.result.blocks[0])
}

async function run(pdb: string, out?: string) {
    const ccd = await getCCD()
    // console.log(Object.keys(ccd).length)

    const mmcif = await getPdb(pdb)
    // console.log(mmcif.chem_comp.id.toArray())

    // const chemCompBond_Schema = { chem_comp_bond: mmCIF_Schema.chem_comp_bond }
    // type chemCompBond_Schema = typeof chemCompBond_Schema;
    // type chemCompBond_Database = Database<chemCompBond_Schema>;
    // interface chemCompBond_Database extends Database<chemCompBond_Schema> {}

    const chemCompBondTables: Table<mmCIF_Schema['chem_comp_bond']>[] = []

    for (let i = 0, n = mmcif.chem_comp._rowCount; i < n; ++i) {
        const ccdId = mmcif.chem_comp.id.value(i);
        // console.log(ccdId)
        if (ccdId in ccd) {
            // console.log(`ccdId ${ccdId} has ${ccd[ccdId].chem_comp_atom._rowCount} atoms`)
            // console.log(ccd[ccdId].chem_comp_bond.atom_id_1.toArray())
            chemCompBondTables.push(ccd[ccdId].chem_comp_bond)
        } else {
            console.error(`ccdId ${ccdId} not found`)
        }
    }

    // for (const k of Object.keys(ccd)) {
    //     console.log(k)
    // }

    // const combinedChemCompBonds = Database.ofTables('chemCompBond', chemCompBond_Schema, {
    //     chem_comp_bond: Table.concat(chemCompBondTables, mmCIF_Schema.chem_comp_bond)
    // })
    // console.log('concat done')
    // console.log(getEncodedCif('chemCompBond', combinedChemCompBonds))
    const combinedMmcif = Database.ofTables('chemCompBond', mmCIF_Schema, Object.assign({},
        mmcif,
        { chem_comp_bond: Table.concat(chemCompBondTables, mmCIF_Schema.chem_comp_bond) }
    ))
    console.log(getEncodedCif(pdb, combinedMmcif))
}

const DATA_DIR = './build/data'
const CCD_PATH = `${DATA_DIR}/components.cif`
const PVCD_PATH = `${DATA_DIR}/aa-variants-v1.cif`
const CCD_URL = 'http://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif'
const PVCD_URL = 'http://ftp.wwpdb.org/pub/pdb/data/monomers/aa-variants-v1.cif'

const parser = new argparse.ArgumentParser({
  addHelp: true,
  description: 'Create a mmcif file that includes relevant CCD and BIRD entries'
});
parser.addArgument([ '--pdb', '-p' ], {
    help: 'Pdb entry id'
});
parser.addArgument([ '--out', '-o' ], {
    help: 'Generated file output path, if not given printed to stdout'
});
parser.addArgument([ '--forceDownload', '-f' ], {
    action: 'storeTrue',
    help: 'Force download of CCD and BIRD'
});
interface Args {
    pdb: string
    out: string
    forceDownload: boolean
}
const args: Args = parser.parseArgs();

const FORCE_DOWNLOAD = args.forceDownload

run(args.pdb, args.out)

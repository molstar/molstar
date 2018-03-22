/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as argparse from 'argparse'
import * as util from 'util'
import * as fs from 'fs'
import * as zlib from 'zlib'
import fetch from 'node-fetch'
require('util.promisify').shim();
const readFile = util.promisify(fs.readFile);
const writeFile = util.promisify(fs.writeFile);

import { Database, Table, DatabaseCollection } from 'mol-data/db'
import CIF from 'mol-io/reader/cif'
// import { CCD_Schema } from 'mol-io/reader/cif/schema/ccd'
import * as Encoder from 'mol-io/writer/cif'
import Computation from 'mol-util/computation'
import { mmCIF_Schema, mmCIF_Database } from 'mol-io/reader/cif/schema/mmcif';
import { CCD_Schema } from 'mol-io/reader/cif/schema/ccd';
import { BIRD_Schema } from 'mol-io/reader/cif/schema/bird';
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
        if (url.endsWith('.gz')) {
            await writeFile(path, zlib.gunzipSync(await data.buffer()))
        } else {
            await writeFile(path, await data.text())
        }
        console.log(`done downloading ${url}`)
    }
}

export async function ensureDataAvailable() {
    await ensureAvailable(CCD_PATH, CCD_URL)
    await ensureAvailable(PVCD_PATH, PVCD_URL)
    await ensureAvailable(BIRD_PATH, BIRD_URL)
}

function showProgress(tag: string, p: Computation.Progress) {
    console.log(`[${tag}] ${p.message} ${p.isIndeterminate ? '' : (p.current / p.max * 100).toFixed(2) + '% '}(${p.elapsedMs | 0}ms)`)
}

export async function readFileAsCollection<S extends Database.Schema>(path: string, schema: S) {
    const parsed = await parseCif(await readFile(path, 'utf8'))
    return CIF.toDatabaseCollection(schema, parsed.result)
}

export async function readCCD() {
    return readFileAsCollection(CCD_PATH, CCD_Schema)
}

export async function readPVCD() {
    return readFileAsCollection(PVCD_PATH, CCD_Schema)
}

export async function readBIRD() {
    return readFileAsCollection(BIRD_PATH, BIRD_Schema)
}

export async function getCCD() {
    await ensureDataAvailable()
    return readPVCD()
}

export async function getBIRD() {
    await ensureDataAvailable()
    return readBIRD()
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

type extraTables = {
    chem_comp_bond: Table<mmCIF_Schema['chem_comp_bond']>,
    pdbx_reference_entity_list: Table<mmCIF_Schema['pdbx_reference_entity_list']>,
    pdbx_reference_entity_link: Table<mmCIF_Schema['pdbx_reference_entity_link']>,
    pdbx_reference_entity_poly_link: Table<mmCIF_Schema['pdbx_reference_entity_poly_link']>,
}
type extraTablesLists = {
    [k in keyof extraTables]: extraTables[k][]
}

export function getExtraTables(mmcif: mmCIF_Database, ccd: DatabaseCollection<CCD_Schema>, bird: DatabaseCollection<BIRD_Schema>) {
    const extraTablesLists: extraTablesLists = {
        chem_comp_bond: [],
        pdbx_reference_entity_list: [],
        pdbx_reference_entity_link: [],
        pdbx_reference_entity_poly_link: []
    }

    for (let i = 0, n = mmcif.chem_comp._rowCount; i < n; ++i) {
        const ccdId = mmcif.chem_comp.id.value(i);
        if (ccdId in ccd) {
            extraTablesLists.chem_comp_bond.push(ccd[ccdId].chem_comp_bond)
        } else {
            console.error(`ccdId ${ccdId} not found`)
        }
    }

    for (let i = 0, n = mmcif.pdbx_molecule_features._rowCount; i < n; ++i) {
        const birdId = mmcif.pdbx_molecule_features.prd_id.value(i);
        if (birdId in bird) {
            const e = bird[birdId]
            extraTablesLists.pdbx_reference_entity_list.push(e.pdbx_reference_entity_list)
            extraTablesLists.pdbx_reference_entity_link.push(e.pdbx_reference_entity_link)
            extraTablesLists.pdbx_reference_entity_poly_link.push(e.pdbx_reference_entity_poly_link)
        } else {
            console.error(`birdId ${birdId} not found`)
        }
    }

    const extraTables: extraTables = Object.assign({}, ...Object.keys(extraTablesLists).map(k => {
        // TODO how to avoid type casting?
        return { [k]: Table.concat((extraTablesLists as any)[k], (mmCIF_Schema as any)[k]) }
    }))

    return extraTables
}

type PartialStructConnRow = Partial<Table.Row<mmCIF_Schema['struct_conn']>>

export function getStructConnValueOrder(value: any): mmCIF_Schema['struct_conn']['pdbx_value_order']['T'] {
    return mmCIF_Schema.struct_conn.pdbx_value_order['T'].includes(value) ? value : 'sing'
}

export function getBirdBonds(mmcif: mmCIF_Database) {
    const bonds: PartialStructConnRow[] = []

    const mol = mmcif.pdbx_molecule
    const molFeat = mmcif.pdbx_molecule_features
    for (let i = 0, n = molFeat._rowCount; i < n; ++i) {
        // console.log(Table.getRow(molFeat, i))
        const instancesAsymIdList: { [k: number]: string[] } = {}
        for (let j = 0, m = mol._rowCount; j < m; ++j) {
            if (mol.prd_id.value(j) === molFeat.prd_id.value(i)) {
                // console.log(Table.getRow(mol, j))
                const instanceId = mol.instance_id.value(j)
                if (instancesAsymIdList[instanceId] === undefined) {
                    instancesAsymIdList[instanceId] = []
                }
                instancesAsymIdList[instanceId].push(mol.asym_id.value(j))
            }
        }
        // console.log(instancesAsymIdList)
        const entityLink = mmcif.pdbx_reference_entity_link
        for (const instanceId of Object.keys(instancesAsymIdList)) {
            const asymIdList = instancesAsymIdList[instanceId as any]
            for (let j = 0, m = entityLink._rowCount; j < m; ++j) {
                if (entityLink.prd_id.value(j) === molFeat.prd_id.value(i)) {
                    // console.log(Table.getRow(entityLink, j))
                    const link: PartialStructConnRow = {
                        ptnr1_label_asym_id: asymIdList[ entityLink.component_1.value(j) - 1 ],
                        ptnr1_label_atom_id: entityLink.atom_id_1.value(j),
                        ptnr1_label_comp_id: entityLink.comp_id_1.value(j),
                        ptnr1_label_seq_id: entityLink.entity_seq_num_1.value(j),

                        ptnr2_label_asym_id: asymIdList[ entityLink.component_2.value(j) - 1 ],
                        ptnr2_label_atom_id: entityLink.atom_id_2.value(j),
                        ptnr2_label_comp_id: entityLink.comp_id_2.value(j),
                        ptnr2_label_seq_id: entityLink.entity_seq_num_2.value(j),

                        pdbx_value_order: getStructConnValueOrder(entityLink.value_order.value(j)),
                    }
                    // console.log(link)
                    bonds.push(link)
                }
            }
        }
    }
    return bonds
}

export function getCcdBonds(mmcif: mmCIF_Database) {
    // const bonds: PartialStructConnRow[] = []
}

async function run(pdb: string, out?: string) {
    const ccd = await getCCD()
    const bird = await getBIRD()

    const mmcif = await getPdb(pdb)
    // console.log(mmcif.chem_comp.id.toArray())

    for (const k of Object.keys(bird)) {
        const entity = bird[k].pdbx_reference_entity_list
        for (let i = 0, n = entity._rowCount; i < n; ++i) {
            if (entity.ref_entity_id.value(i) !== entity.component_id.value(i).toString()) {
                console.log(Table.getRow(entity, i))
            }
        }

        const link = bird[k].pdbx_reference_entity_link
        for (let i = 0, n = link._rowCount; i < n; ++i) {
            if (link.value_order.value(i) !== 'sing') {
                console.log(Table.getRow(link, i))
            }
        }

        const polyLink = bird[k].pdbx_reference_entity_poly_link
        for (let i = 0, n = link._rowCount; i < n; ++i) {
            if (polyLink.value_order.value(i) !== 'sing') {
                console.log(Table.getRow(polyLink, i))
            }
        }
    }

    const extraTables = getExtraTables(mmcif, ccd, bird)
    const combinedMmcif = Database.ofTables('mmcif_combined', mmCIF_Schema, Object.assign({}, mmcif, extraTables))
    // console.log(getEncodedCif(pdb, combinedMmcif))
    // console.log(Database.getTablesAsRows(combinedMmcif))

    // console.log(getBirdBonds(combinedMmcif))
    console.log(getCcdBonds(combinedMmcif))
}

const DATA_DIR = './build/data'
const CCD_PATH = `${DATA_DIR}/components.cif`
const PVCD_PATH = `${DATA_DIR}/aa-variants-v1.cif`
const BIRD_PATH = `${DATA_DIR}/prd-all.cif`
const CCD_URL = 'http://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif'
const PVCD_URL = 'http://ftp.wwpdb.org/pub/pdb/data/monomers/aa-variants-v1.cif'
const BIRD_URL = 'http://ftp.wwpdb.org/pub/pdb/data/bird/prd/prd-all.cif.gz'

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

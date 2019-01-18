/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as argparse from 'argparse'
// import * as util from 'util'
import * as fs from 'fs'
import fetch from 'node-fetch'

import Csv from 'mol-io/reader/csv/parser'
import CIF, { CifFrame } from 'mol-io/reader/cif'
import { generateSchema } from './util/cif-dic'
import { generate } from './util/generate'
import { Filter } from './util/schema'

async function runGenerateSchema(name: string, fieldNamesPath?: string, typescript = false, out?: string) {
    await ensureMmcifDicAvailable()
    const mmcifDic = await CIF.parseText(fs.readFileSync(MMCIF_DIC_PATH, 'utf8')).run();
    if (mmcifDic.isError) throw mmcifDic

    await ensureIhmDicAvailable()
    const ihmDic = await CIF.parseText(fs.readFileSync(IHM_DIC_PATH, 'utf8')).run();
    if (ihmDic.isError) throw ihmDic

    await ensureCarbBranchDicAvailable()
    const carbBranchDic = await CIF.parseText(fs.readFileSync(CARB_BRANCH_DIC_PATH, 'utf8')).run();
    if (carbBranchDic.isError) throw carbBranchDic

    await ensureCarbCompDicAvailable()
    const carbCompDic = await CIF.parseText(fs.readFileSync(CARB_COMP_DIC_PATH, 'utf8')).run();
    if (carbCompDic.isError) throw carbCompDic

    const mmcifDicVersion = CIF.schema.dic(mmcifDic.result.blocks[0]).dictionary.version.value(0)
    const ihmDicVersion = CIF.schema.dic(ihmDic.result.blocks[0]).dictionary.version.value(0)
    const carbDicVersion = 'draft'
    const version = `Dictionary versions: mmCIF ${mmcifDicVersion}, IHM ${ihmDicVersion}, CARB ${carbDicVersion}.`

    const frames: CifFrame[] = [...mmcifDic.result.blocks[0].saveFrames, ...ihmDic.result.blocks[0].saveFrames, ...carbBranchDic.result.blocks[0].saveFrames, ...carbCompDic.result.blocks[0].saveFrames]
    const schema = generateSchema(frames)

    const filter = fieldNamesPath ? await getFieldNamesFilter(fieldNamesPath) : undefined
    const output = typescript ? generate(name, version, schema, filter) : JSON.stringify(schema, undefined, 4)

    if (out) {
        fs.writeFileSync(out, output)
    } else {
        console.log(output)
    }
}

async function getFieldNamesFilter(fieldNamesPath: string): Promise<Filter> {
    const fieldNamesStr = fs.readFileSync(fieldNamesPath, 'utf8')
    const parsed = await Csv(fieldNamesStr, { noColumnNames: true }).run();
    if (parsed.isError) throw parser.error
    const csvFile = parsed.result;

    const fieldNamesCol = csvFile.table.getColumn('0')
    if (!fieldNamesCol) throw 'error getting fields columns'
    const fieldNames = fieldNamesCol.toStringArray()

    const filter: Filter = {}
    fieldNames.forEach((name, i) => {
        const [ category, field ] = name.split('.')
        // console.log(category, field)
        if (!filter[ category ]) filter[ category ] = {}
        filter[ category ][ field ] = true
    })
    // console.log(filter)
    return filter
}

async function ensureMmcifDicAvailable() { await ensureDicAvailable(MMCIF_DIC_PATH, MMCIF_DIC_URL) }
async function ensureIhmDicAvailable() { await ensureDicAvailable(IHM_DIC_PATH, IHM_DIC_URL) }
async function ensureCarbBranchDicAvailable() { await ensureDicAvailable(CARB_BRANCH_DIC_PATH, CARB_BRANCH_DIC_URL) }
async function ensureCarbCompDicAvailable() { await ensureDicAvailable(CARB_COMP_DIC_PATH, CARB_COMP_DIC_URL) }

async function ensureDicAvailable(dicPath: string, dicUrl: string) {
    if (FORCE_DIC_DOWNLOAD || !fs.existsSync(dicPath)) {
        const name = dicUrl.substr(dicUrl.lastIndexOf('/') + 1)
        console.log(`downloading ${name}...`)
        const data = await fetch(dicUrl)
        if (!fs.existsSync(DIC_DIR)) {
            fs.mkdirSync(DIC_DIR);
        }
        fs.writeFileSync(dicPath, await data.text())
        console.log(`done downloading ${name}`)
    }
}

const DIC_DIR = './build/dics'
const MMCIF_DIC_PATH = `${DIC_DIR}/mmcif_pdbx_v50.dic`
const MMCIF_DIC_URL = 'http://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic'
const IHM_DIC_PATH = `${DIC_DIR}/ihm-extension.dic`
const IHM_DIC_URL = 'https://raw.githubusercontent.com/ihmwg/IHM-dictionary/master/ihm-extension.dic'
const CARB_BRANCH_DIC_PATH = `${DIC_DIR}/entity_branch-extension.dic`
const CARB_BRANCH_DIC_URL = 'https://raw.githubusercontent.com/pdbxmmcifwg/carbohydrate-extension/master/dict/entity_branch-extension.dic'
const CARB_COMP_DIC_PATH = `${DIC_DIR}/chem_comp-extension.dic`
const CARB_COMP_DIC_URL = 'https://raw.githubusercontent.com/pdbxmmcifwg/carbohydrate-extension/master/dict/chem_comp-extension.dic'

const parser = new argparse.ArgumentParser({
  addHelp: true,
  description: 'Create schema from mmcif dictionary (v50 plus IHM and entity_branch extensions, downloaded from wwPDB)'
});
parser.addArgument([ '--name', '-n' ], {
    defaultValue: 'mmCIF',
    help: 'Schema name'
});
parser.addArgument([ '--out', '-o' ], {
    help: 'Generated schema output path, if not given printed to stdout'
});
parser.addArgument([ '--typescript', '-ts' ], {
    action: 'storeTrue',
    help: 'Output schema as TypeScript instead of as JSON'
});
parser.addArgument([ '--fieldNamesPath', '-fn' ], {
    defaultValue: '',
    help: 'Field names to include'
});
parser.addArgument([ '--forceDicDownload', '-f' ], {
    action: 'storeTrue',
    help: 'Force download of dictionaries'
});
interface Args {
    name: string
    forceDicDownload: boolean
    fieldNamesPath: string
    typescript: boolean
    out: string
}
const args: Args = parser.parseArgs();

const FORCE_DIC_DOWNLOAD = args.forceDicDownload

if (args.name) {
    runGenerateSchema(args.name, args.fieldNamesPath, args.typescript, args.out).catch(e => {
        console.error(e)
    })
}

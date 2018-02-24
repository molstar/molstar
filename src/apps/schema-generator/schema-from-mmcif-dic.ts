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
import CIF from 'mol-io/reader/cif'
import { generateSchema } from './util/cif-dic'
import { generate } from './util/generate'
import { Filter, mergeFilters } from './util/json-schema'

async function runGenerateSchema(name: string, fieldNamesPath?: string, minCount = 0, typescript = false, out?: string) {
    await ensureMmcifDicAvailable()
    const comp = CIF.parseText(fs.readFileSync(MMCIF_DIC_PATH, 'utf8'))
    const parsed = await comp();
    if (parsed.isError) throw parsed

    console.log(fieldNamesPath, minCount)

    let filter: Filter | undefined
    if (minCount && fieldNamesPath) {
        filter = mergeFilters(
            await getUsageCountsFilter(minCount),
            await getFieldNamesFilter(fieldNamesPath)
        )
    } else if (minCount) {
        filter = await getUsageCountsFilter(minCount)
    } else if (fieldNamesPath) {
        console.log('MOIN')
        filter = await getFieldNamesFilter(fieldNamesPath)
    }

    const schema = generateSchema(parsed.result.blocks[0])
    const output = typescript ? generate(name, schema, filter) : JSON.stringify(schema, undefined, 4)

    if (out) {
        fs.writeFileSync(out, output)
    } else {
        console.log(output)
    }
}

async function getFieldNamesFilter(fieldNamesPath: string): Promise<Filter> {
    const fieldNamesStr = fs.readFileSync(fieldNamesPath, 'utf8')
    const parsed = await Csv(fieldNamesStr, { noColumnNames: true })();
    if (parsed.isError) throw parser.error
    const csvFile = parsed.result;

    const fieldNamesCol = csvFile.table.getColumn('0')
    if (!fieldNamesCol) throw 'error getting fields columns'
    const fieldNames = fieldNamesCol.toStringArray()

    const filter: Filter = {}
    fieldNames.forEach((name, i) => {
        const [ category, field ] = name.split('.')
        console.log(category, field)
        if (!filter[ category ]) filter[ category ] = {}
        filter[ category ][ field ] = true
    })
    console.log(filter)
    return filter
}

async function getUsageCountsFilter(minCount: number): Promise<Filter> {
    const usageCountsStr = fs.readFileSync(MMCIF_USAGE_COUNTS_PATH, 'utf8')
    const parsed = await Csv(usageCountsStr, { delimiter: ' ' })();
    if (parsed.isError) throw parser.error
    const csvFile = parsed.result;

    const fieldNamesCol = csvFile.table.getColumn('field_name')
    const usageCountsCol = csvFile.table.getColumn('usage_count')
    if (!fieldNamesCol || !usageCountsCol) throw 'error getting usage columns'
    const fieldNames = fieldNamesCol.toStringArray()
    const usageCounts = usageCountsCol.toIntArray()

    const filter: Filter = {}
    fieldNames.forEach((name, i) => {
        if (usageCounts[i] < minCount) return
        const [ category, field ] = name.substr(1).split('.')
        if (!filter[ category ]) filter[ category ] = {}
        filter[ category ][ field ] = true
    })
    return filter
}

async function ensureMmcifDicAvailable() {
    if (FORCE_MMCIF_DOWNLOAD || !fs.existsSync(MMCIF_DIC_PATH)) {
        console.log('downloading mmcif dic...')
        const data = await fetch(MMCIF_DIC_URL)
        if (!fs.existsSync(MMCIF_DIC_DIR)){
            fs.mkdirSync(MMCIF_DIC_DIR);
        }
        fs.writeFileSync(MMCIF_DIC_PATH, await data.text())
        console.log('done downloading mmcif dic')
    }
}

const MMCIF_USAGE_COUNTS_PATH = './data/mmcif-usage-counts.txt'
const MMCIF_DIC_DIR = './build/dics'
const MMCIF_DIC_PATH = `${MMCIF_DIC_DIR}/mmcif_pdbx_v50.dic`
const MMCIF_DIC_URL = 'http://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic'

const parser = new argparse.ArgumentParser({
  addHelp: true,
  description: 'Create schema from mmcif dictionary (v50, downloaded from wwPDB)'
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
parser.addArgument([ '--minFieldUsageCount', '-mc' ], {
    defaultValue: 0,
    type: parseInt,
    help: 'Minimum mmcif field usage counts'
});
parser.addArgument([ '--fieldNamesPath', '-fn' ], {
    defaultValue: '',
    help: 'Field names to include'
});
parser.addArgument([ '--forceMmcifDicDownload', '-f' ], {
    action: 'storeTrue',
    help: 'Force download of mmcif dictionary'
});
interface Args {
    name: string
    forceMmcifDicDownload: boolean
    fieldNamesPath: string
    minFieldUsageCount: number
    typescript: boolean
    out: string
}
const args: Args = parser.parseArgs();

const FORCE_MMCIF_DOWNLOAD = args.forceMmcifDicDownload

if (args.name) {
    runGenerateSchema(args.name, args.fieldNamesPath, args.minFieldUsageCount, args.typescript, args.out)
}

#!/usr/bin/env node
/**
 * Copyright (c) 2017-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as argparse from 'argparse';
import * as fs from 'fs';
import * as path from 'path';
import fetch from 'node-fetch';

import { parseCsv } from '../../mol-io/reader/csv/parser';
import { CifFrame, CifBlock } from '../../mol-io/reader/cif';
import { parseCifText } from '../../mol-io/reader/cif/text/parser';
import { generateSchema } from './util/cif-dic';
import { generate } from './util/generate';
import { Filter, Database } from './util/schema';
import { parseImportGet } from './util/helper';

function getDicVersion(block: CifBlock) {
    return block.categories.dictionary.getField('version')!.str(0);
}

function getDicNamespace(block: CifBlock) {
    return block.categories.dictionary.getField('namespace')!.str(0);
}

async function runGenerateSchemaMmcif(name: string, fieldNamesPath: string, typescript = false, out: string, moldbImportPath: string, addAliases: boolean) {
    await ensureMmcifDicAvailable();
    const mmcifDic = await parseCifText(fs.readFileSync(MMCIF_DIC_PATH, 'utf8')).run();
    if (mmcifDic.isError) throw mmcifDic;

    await ensureIhmDicAvailable();
    const ihmDic = await parseCifText(fs.readFileSync(IHM_DIC_PATH, 'utf8')).run();
    if (ihmDic.isError) throw ihmDic;

    await ensureMaDicAvailable();
    const maDic = await parseCifText(fs.readFileSync(MA_DIC_PATH, 'utf8')).run();
    if (maDic.isError) throw maDic;

    const mmcifDicVersion = getDicVersion(mmcifDic.result.blocks[0]);
    const ihmDicVersion = getDicVersion(ihmDic.result.blocks[0]);
    const maDicVersion = getDicVersion(maDic.result.blocks[0]);
    const version = `Dictionary versions: mmCIF ${mmcifDicVersion}, IHM ${ihmDicVersion}, MA ${maDicVersion}.`;

    const frames: CifFrame[] = [...mmcifDic.result.blocks[0].saveFrames, ...ihmDic.result.blocks[0].saveFrames, ...maDic.result.blocks[0].saveFrames];
    const schema = generateSchema(frames);

    await runGenerateSchema(name, version, schema, fieldNamesPath, typescript, out, moldbImportPath, addAliases);
}

async function runGenerateSchemaCifCore(name: string, fieldNamesPath: string, typescript = false, out: string, moldbImportPath: string, addAliases: boolean) {
    await ensureCifCoreDicAvailable();
    const cifCoreDic = await parseCifText(fs.readFileSync(CIF_CORE_DIC_PATH, 'utf8')).run();
    if (cifCoreDic.isError) throw cifCoreDic;

    const cifCoreDicVersion = getDicVersion(cifCoreDic.result.blocks[0]);
    const version = `Dictionary versions: CifCore ${cifCoreDicVersion}.`;

    const frames: CifFrame[] = [...cifCoreDic.result.blocks[0].saveFrames];
    const imports = await resolveImports(frames, DIC_DIR);
    const schema = generateSchema(frames, imports);

    await runGenerateSchema(name, version, schema, fieldNamesPath, typescript, out, moldbImportPath, addAliases);
}

async function resolveImports(frames: CifFrame[], baseDir: string): Promise<Map<string, CifFrame[]>> {
    const imports = new Map<string, CifFrame[]>();

    for (const d of frames) {
        if ('import' in d.categories) {
            const importGet = parseImportGet(d.categories['import'].getField('get')!.str(0));
            for (const g of importGet) {
                const { file } = g;
                if (!file) continue;
                if (imports.has(file)) continue;

                const dic = await parseCifText(fs.readFileSync(path.join(baseDir, file), 'utf8')).run();
                if (dic.isError) throw dic;

                imports.set(file, [...dic.result.blocks[0].saveFrames]);
            }
        }
    }

    return imports;
}

async function runGenerateSchemaDic(name: string, dicPath: string, fieldNamesPath: string, typescript = false, out: string, moldbImportPath: string, addAliases: boolean) {
    const dic = await parseCifText(fs.readFileSync(dicPath, 'utf8')).run();
    if (dic.isError) throw dic;

    const dicVersion = getDicVersion(dic.result.blocks[0]);
    const dicName = getDicNamespace(dic.result.blocks[0]);
    const version = `Dictionary versions: ${dicName} ${dicVersion}.`;

    const frames: CifFrame[] = [...dic.result.blocks[0].saveFrames];
    const imports = await resolveImports(frames, path.dirname(dicPath));
    const schema = generateSchema(frames, imports);

    await runGenerateSchema(name, version, schema, fieldNamesPath, typescript, out, moldbImportPath, addAliases);
}

async function runGenerateSchema(name: string, version: string, schema: Database, fieldNamesPath: string, typescript = false, out: string, moldbImportPath: string, addAliases: boolean) {
    const filter = fieldNamesPath ? await getFieldNamesFilter(fieldNamesPath) : undefined;
    const output = typescript ? generate(name, version, schema, filter, moldbImportPath, addAliases) : JSON.stringify(schema, undefined, 4);

    if (out) {
        fs.writeFileSync(out, output);
    } else {
        console.log(output);
    }
}

async function getFieldNamesFilter(fieldNamesPath: string): Promise<Filter> {
    const fieldNamesStr = fs.readFileSync(fieldNamesPath, 'utf8');
    const parsed = await parseCsv(fieldNamesStr, { noColumnNames: true }).run();
    if (parsed.isError) throw parser.error;
    const csvFile = parsed.result;

    const fieldNamesCol = csvFile.table.getColumn('0');
    if (!fieldNamesCol) throw new Error('error getting fields columns');
    const fieldNames = fieldNamesCol.toStringArray();

    const filter: Filter = {};
    fieldNames.forEach((name, i) => {
        const [category, field] = name.split('.');
        // console.log(category, field)
        if (!filter[category]) filter[category] = {};
        filter[category][field] = true;
    });
    return filter;
}

async function ensureMmcifDicAvailable() { await ensureDicAvailable(MMCIF_DIC_PATH, MMCIF_DIC_URL); }
async function ensureIhmDicAvailable() { await ensureDicAvailable(IHM_DIC_PATH, IHM_DIC_URL); }
async function ensureMaDicAvailable() { await ensureDicAvailable(MA_DIC_PATH, MA_DIC_URL); }
async function ensureCifCoreDicAvailable() {
    await ensureDicAvailable(CIF_CORE_DIC_PATH, CIF_CORE_DIC_URL);
    await ensureDicAvailable(CIF_CORE_ENUM_PATH, CIF_CORE_ENUM_URL);
    await ensureDicAvailable(CIF_CORE_ATTR_PATH, CIF_CORE_ATTR_URL);
}

async function ensureDicAvailable(dicPath: string, dicUrl: string) {
    if (FORCE_DIC_DOWNLOAD || !fs.existsSync(dicPath)) {
        const name = dicUrl.substr(dicUrl.lastIndexOf('/') + 1);
        console.log(`downloading ${name}...`);
        const data = await fetch(dicUrl);
        if (!fs.existsSync(DIC_DIR)) {
            fs.mkdirSync(DIC_DIR);
        }
        fs.writeFileSync(dicPath, await data.text());
        console.log(`done downloading ${name}`);
    }
}

const DIC_DIR = path.resolve(__dirname, '../../../../build/dics/');
const MMCIF_DIC_PATH = `${DIC_DIR}/mmcif_pdbx_v50.dic`;
const MMCIF_DIC_URL = 'http://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic';
const IHM_DIC_PATH = `${DIC_DIR}/mmcif_ihm_ext.dic`;
const IHM_DIC_URL = 'https://raw.githubusercontent.com/ihmwg/IHMCIF/master/dist/mmcif_ihm_ext.dic';
const MA_DIC_PATH = `${DIC_DIR}/ma-extension.dic`;
const MA_DIC_URL = 'https://raw.githubusercontent.com/ihmwg/ModelCIF/master/dist/mmcif_ma.dic';

const CIF_CORE_DIC_PATH = `${DIC_DIR}/cif_core.dic`;
const CIF_CORE_DIC_URL = 'https://raw.githubusercontent.com/COMCIFS/cif_core/master/cif_core.dic';
const CIF_CORE_ENUM_PATH = `${DIC_DIR}/templ_enum.cif`;
const CIF_CORE_ENUM_URL = 'https://raw.githubusercontent.com/COMCIFS/cif_core/master/templ_enum.cif';
const CIF_CORE_ATTR_PATH = `${DIC_DIR}/templ_attr.cif`;
const CIF_CORE_ATTR_URL = 'https://raw.githubusercontent.com/COMCIFS/cif_core/master/templ_attr.cif';

const parser = new argparse.ArgumentParser({
    add_help: true,
    description: 'Create schema from mmcif dictionary (v50 plus IHM and entity_branch extensions, downloaded from wwPDB)'
});
parser.add_argument('--preset', '-p', {
    default: '',
    choices: ['', 'mmCIF', 'CCD', 'BIRD', 'CifCore'],
    help: 'Preset name'
});
parser.add_argument('--name', '-n', {
    default: '',
    help: 'Schema name'
});
parser.add_argument('--out', '-o', {
    help: 'Generated schema output path, if not given printed to stdout'
});
parser.add_argument('--targetFormat', '-tf', {
    default: 'typescript-molstar',
    choices: ['typescript-molstar', 'json-internal'],
    help: 'Target format'
});
parser.add_argument('--dicPath', '-d', {
    default: '',
    help: 'Path to dictionary'
});
parser.add_argument('--fieldNamesPath', '-fn', {
    default: '',
    help: 'Field names to include'
});
parser.add_argument('--forceDicDownload', '-f', {
    action: 'store_true',
    help: 'Force download of dictionaries'
});
parser.add_argument('--moldataImportPath', '-mip', {
    default: 'molstar/lib/mol-data',
    help: 'mol-data import path (for typescript target only)'
});
parser.add_argument('--addAliases', '-aa', {
    action: 'store_true',
    help: 'Add field name/path aliases'
});
interface Args {
    name: string
    preset: '' | 'mmCIF' | 'CCD' | 'BIRD' | 'CifCore'
    forceDicDownload: boolean
    dic: '' | 'mmCIF' | 'CifCore'
    dicPath: string,
    fieldNamesPath: string
    targetFormat: 'typescript-molstar' | 'json-internal'
    out: string,
    moldataImportPath: string
    addAliases: boolean
}
const args: Args = parser.parse_args();

const FORCE_DIC_DOWNLOAD = args.forceDicDownload;

switch (args.preset) {
    case 'mmCIF':
        args.name = 'mmCIF';
        args.dic = 'mmCIF';
        args.fieldNamesPath = path.resolve(__dirname, '../../../../data/cif-field-names/mmcif-field-names.csv');
        break;
    case 'CCD':
        args.name = 'CCD';
        args.dic = 'mmCIF';
        args.fieldNamesPath = path.resolve(__dirname, '../../../../data/cif-field-names/ccd-field-names.csv');
        break;
    case 'BIRD':
        args.name = 'BIRD';
        args.dic = 'mmCIF';
        args.fieldNamesPath = path.resolve(__dirname, '../../../../data/cif-field-names/bird-field-names.csv');
        break;
    case 'CifCore':
        args.name = 'CifCore';
        args.dic = 'CifCore';
        args.fieldNamesPath = path.resolve(__dirname, '../../../../data/cif-field-names/cif-core-field-names.csv');
        break;
}

if (args.name) {
    const typescript = args.targetFormat === 'typescript-molstar';
    if (args.dicPath) {
        runGenerateSchemaDic(args.name, args.dicPath, args.fieldNamesPath, typescript, args.out, args.moldataImportPath, args.addAliases).catch(e => {
            console.error(e);
        });
    } else if (args.dic === 'mmCIF') {
        runGenerateSchemaMmcif(args.name, args.fieldNamesPath, typescript, args.out, args.moldataImportPath, args.addAliases).catch(e => {
            console.error(e);
        });
    } else if (args.dic === 'CifCore') {
        runGenerateSchemaCifCore(args.name, args.fieldNamesPath, typescript, args.out, args.moldataImportPath, args.addAliases).catch(e => {
            console.error(e);
        });
    }
}

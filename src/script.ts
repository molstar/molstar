/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as util from 'util'
import * as fs from 'fs'

import Gro from './reader/gro/parser'
import CIF from './reader/cif/index'

import { getSchema } from './reader/cif/schema/utils'

const file = '1crn.gro'
// const file = 'water.gro'
// const file = 'test.gro'
// const file = 'md_1u19_trj.gro'

function showProgress(tag: string, p: Computation.Progress) {
    console.log(`[${tag}] ${p.message} ${p.isIndeterminate ? '' : (p.current / p.max * 100).toFixed(2) + '% '}(${p.elapsedMs | 0}ms)`)
}

async function runGro(input: string) {
    console.time('parseGro');
    const comp = Gro(input);

    const ctx = Computation.observable({ updateRateMs: 150, observer: p => showProgress('GRO', p) });
    const parsed = await comp(ctx);
    console.timeEnd('parseGro');

    if (parsed.isError) {
        console.log(parsed);
        return;
    }

    const groFile = parsed.result

    console.log('structure count: ', groFile.structures.length);

    const data = groFile.structures[0];

    // const header = groFile.blocks[0].getCategory('header')
    const { header, atoms } = data;
    console.log(JSON.stringify(header, null, 2));
    console.log('number of atoms:', atoms.count);

    console.log(`'${atoms.residueNumber.value(1)}'`)
    console.log(`'${atoms.residueName.value(1)}'`)
    console.log(`'${atoms.atomName.value(1)}'`)
    console.log(atoms.z.value(1))
    console.log(`'${atoms.z.value(1)}'`)

    const n = atoms.count;
    console.log('rowCount', n)

    console.time('getFloatArray x')
    const x = atoms.x.toArray({ array: Float32Array })
    console.timeEnd('getFloatArray x')
    console.log(x.length, x[0], x[x.length - 1])

    console.time('getFloatArray y')
    const y = atoms.y.toArray({ array: Float32Array })
    console.timeEnd('getFloatArray y')
    console.log(y.length, y[0], y[y.length - 1])

    console.time('getFloatArray z')
    const z = atoms.z.toArray({ array: Float32Array })
    console.timeEnd('getFloatArray z')
    console.log(z.length, z[0], z[z.length - 1])

    console.time('getIntArray residueNumber')
    const residueNumber = atoms.residueNumber.toArray({ array: Int32Array })
    console.timeEnd('getIntArray residueNumber')
    console.log(residueNumber.length, residueNumber[0], residueNumber[residueNumber.length - 1])
}

export function _gro() {
    fs.readFile(`./examples/${file}`, 'utf8', function (err, input) {
        if (err) {
            return console.log(err);
        }
        runGro(input)
    });
}

// _gro()

async function runCIF(input: string | Uint8Array) {
    console.time('parseCIF');
    const comp = typeof input === 'string' ? CIF.parseText(input) : CIF.parseBinary(input);

    const ctx = Computation.observable({ updateRateMs: 250, observer: p => showProgress('CIF', p) });
    const parsed = await comp(ctx);
    console.timeEnd('parseCIF');
    if (parsed.isError) {
        console.log(parsed);
        return;
    }

    const data = parsed.result.blocks[0];
    const atom_site = data.categories._atom_site;
    console.log(atom_site.getField('Cartn_x')!.float(0));
    //console.log(atom_site.getField('label_atom_id')!.toStringArray());

    const mmcif = CIF.schema.mmCIF(data);
    console.log(mmcif.atom_site.Cartn_x.value(0));
    console.log(mmcif.entity.type.toArray());
    console.log(mmcif.pdbx_struct_oper_list.matrix.value(0));
}

export function _cif() {
    let path = `./examples/1cbs_updated.cif`;
    path = '../test/3j3q.cif'  // lets have a relative path for big test files
    fs.readFile(path, 'utf8', function (err, input) {
        if (err) {
            return console.log(err);
        }
        console.log('------------------');
        console.log('Text CIF:');
        runCIF(input);
    });

    path = `./examples/1cbs_full.bcif`;
    // const path = 'c:/test/quick/3j3q.cif';
    fs.readFile(path, function (err, input) {
        if (err) {
            return console.log(err);
        }
        console.log('------------------');
        console.log('BinaryCIF:');
        const data = new Uint8Array(input.byteLength);
        for (let i = 0; i < input.byteLength; i++) data[i] = input[i];
        runCIF(input);
    });
}

// _cif();

async function runDic(input: string | Uint8Array) {
    console.time('parseDic');
    const comp = typeof input === 'string' ? CIF.parseText(input) : CIF.parseBinary(input);

    const ctx = Computation.observable({ updateRateMs: 250, observer: p => showProgress('DIC', p) });
    const parsed = await comp(ctx);
    console.timeEnd('parseDic');
    if (parsed.isError) {
        console.log(parsed);
        return;
    }

    const schema = getSchema(parsed.result.blocks[0])
    // console.log(util.inspect(schema, {showHidden: false, depth: 1}))
    console.log(util.inspect(Object.keys(schema).length, {showHidden: false, depth: 1}))
}

export function _dic() {
    let path = './build/dics/mmcif_pdbx_v50.dic'
    fs.readFile(path, 'utf8', function (err, input) {
        if (err) {
            return console.log(err);
        }
        console.log('------------------');
        console.log('Text DIC:');
        runDic(input);
    });
}

_dic();

import Computation from './utils/computation'
const comp = Computation.create(async ctx => {
    for (let i = 0; i < 0; i++) {
        await new Promise(res => setTimeout(res, 500));
        if (ctx.requiresUpdate) await ctx.update({ message: 'working', current: i, max: 2 });
    }
    return 42;
});
async function testComp() {
    const ctx = Computation.observable({ observer: p => showProgress('test', p) });
    const ret = await comp(ctx);
    console.log('computation returned', ret);
}
testComp();
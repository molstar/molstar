/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as argparse from 'argparse'
import fetch from 'node-fetch'
require('util.promisify').shim();

import CIF from 'mol-io/reader/cif'
import Computation from 'mol-util/computation'
import { Model, Structure } from 'mol-model/structure'

function showProgress(tag: string, p: Computation.Progress) {
    console.log(`[${tag}] ${p.message} ${p.isIndeterminate ? '' : (p.current / p.max * 100).toFixed(2) + '% '}(${p.elapsedMs | 0}ms)`)
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

export async function getPdb(pdb: string) {
    console.log(`downloading ${pdb}...`)
    const data = await fetch(`https://files.rcsb.org/download/${pdb}.cif`)
    console.log(`done downloading ${pdb}`)

    const parsed = await parseCif(await data.text())
    return CIF.schema.mmCIF(parsed.result.blocks[0])
}

async function run(pdb: string, out?: string) {
    const mmcif = await getPdb(pdb)

    const models = Model.create({ kind: 'mmCIF', data: mmcif });
    const structure = Structure.ofModel(models[0])

    console.log(structure)
    console.log(Model.bonds(models[0]))
}

const parser = new argparse.ArgumentParser({
  addHelp: true,
  description: 'Print info about a structure'
});
parser.addArgument([ '--pdb', '-p' ], {
    help: 'Pdb entry id'
});
interface Args {
    pdb: string
}
const args: Args = parser.parseArgs();

run(args.pdb)

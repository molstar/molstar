/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as argparse from 'argparse'
import fetch from 'node-fetch'
require('util.promisify').shim();

// import { Table } from 'mol-data/db'
import CIF from 'mol-io/reader/cif'
import { Model, Structure, ElementSet, Unit } from 'mol-model/structure'
import { Run, Progress } from 'mol-task'
import { OrderedSet } from 'mol-data/int';

async function parseCif(data: string|Uint8Array) {
    const comp = CIF.parse(data);
    const parsed = await Run(comp, p => console.log(Progress.format(p)), 250);
    if (parsed.isError) throw parsed;
    return parsed
}

async function getPdb(pdb: string) {
    const data = await fetch(`https://files.rcsb.org/download/${pdb}.cif`)
    const parsed = await parseCif(await data.text())
    return CIF.schema.mmCIF(parsed.result.blocks[0])
}

export function atomLabel(model: Model, aI: number) {
    const { atoms, residues, chains, residueSegments, chainSegments } = model.hierarchy
    const { label_atom_id } = atoms
    const { label_comp_id, label_seq_id } = residues
    const { label_asym_id } = chains
    const rI = residueSegments.segmentMap[aI]
    const cI = chainSegments.segmentMap[aI]
    return `${label_asym_id.value(cI)} ${label_comp_id.value(rI)} ${label_seq_id.value(rI)} ${label_atom_id.value(aI)}`
}


function printBonds(structure: Structure) {

    const { units, elements } = structure;
    const unitIds = ElementSet.unitIndices(elements);

    for (let i = 0, _i = OrderedSet.size(unitIds); i < _i; i++) {
        const unit = units[OrderedSet.getAt(unitIds, i)];
        const group = ElementSet.groupFromUnitIndex(elements, OrderedSet.getAt(unitIds, i));

        const { count, offset, neighbor } = Unit.getGroupBonds(unit, group);
        const { model }  = unit;

        for (let j = 0; j < count; ++j) {
            const start = offset[j];
            const end = offset[j + 1];
            for (let bI = start; bI < end; bI++) {
                console.log(`${atomLabel(model, j)} -- ${atomLabel(model, neighbor[bI])}`)
            }
        }
    }
}

async function run(pdb: string) {
    const mmcif = await getPdb(pdb)
    const models = Model.create({ kind: 'mmCIF', data: mmcif });
    const structure = Structure.ofModel(models[0])
    // console.log(structure)
    printBonds(structure)
}

const parser = new argparse.ArgumentParser({
  addHelp: true,
  description: 'Print info about a structure, mainly to test and showcase the mol-model module'
});
parser.addArgument([ '--pdb', '-p' ], {
    help: 'Pdb entry id'
});
interface Args {
    pdb: string
}
const args: Args = parser.parseArgs();

run(args.pdb)

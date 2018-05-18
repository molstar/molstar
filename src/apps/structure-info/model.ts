/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as argparse from 'argparse'
require('util.promisify').shim();

// import { Table } from 'mol-data/db'
import CIF from 'mol-io/reader/cif'
import { Model, Structure, Element, Unit, Queries } from 'mol-model/structure'
// import { Run, Progress } from 'mol-task'
import { OrderedSet } from 'mol-data/int';
import { Table } from 'mol-data/db';
import { mmCIF_Database } from 'mol-io/reader/cif/schema/mmcif';
import { openCif, downloadCif } from './helpers';
import { BitFlags } from 'mol-util';
import { SecondaryStructureType } from 'mol-model/structure/model/types';


async function downloadFromPdb(pdb: string) {
    // `https://files.rcsb.org/download/${pdb}.cif`
    const parsed = await downloadCif(`http://www.ebi.ac.uk/pdbe/static/entry/${pdb}_updated.cif`, false);
    return CIF.schema.mmCIF(parsed.blocks[0]);
}

async function readPdbFile(path: string) {
    const parsed = await openCif(path);
    return CIF.schema.mmCIF(parsed.blocks[0]);
}

export function atomLabel(model: Model, aI: number) {
    const { atoms, residues, chains, residueSegments, chainSegments } = model.atomicHierarchy
    const { label_atom_id } = atoms
    const { label_comp_id, label_seq_id } = residues
    const { label_asym_id } = chains
    const rI = residueSegments.segmentMap[aI]
    const cI = chainSegments.segmentMap[aI]
    return `${label_asym_id.value(cI)} ${label_comp_id.value(rI)} ${label_seq_id.value(rI)} ${label_atom_id.value(aI)}`
}

export function residueLabel(model: Model, rI: number) {
    const { residues, chains, residueSegments, chainSegments } = model.atomicHierarchy
    const { label_comp_id, label_seq_id } = residues
    const { label_asym_id } = chains
    const cI = chainSegments.segmentMap[residueSegments.segments[rI]]
    return `${label_asym_id.value(cI)} ${label_comp_id.value(rI)} ${label_seq_id.value(rI)}`
}

export function printSecStructure(model: Model) {
    console.log('Secondary Structure\n=============');
    const { residues } = model.atomicHierarchy;
    const { type, key } = model.properties.secondaryStructure;

    const count = residues._rowCount;
    let rI = 0;
    while (rI < count) {
        let start = rI;
        while (rI < count && key[start] === key[rI]) rI++;
        rI--;

        if (BitFlags.has(type[start], SecondaryStructureType.Flag.Beta)) {
            console.log(`Sheet: ${residueLabel(model, start)} - ${residueLabel(model, rI)} (key ${key[start]})`);
        } else if (BitFlags.has(type[start], SecondaryStructureType.Flag.Helix)) {
            console.log(`Helix: ${residueLabel(model, start)} - ${residueLabel(model, rI)} (key ${key[start]})`);
        }

        rI++;
    }
}

export function printBonds(structure: Structure) {
    for (const unit of structure.units) {
        if (!Unit.isAtomic(unit)) continue;

        const elements = unit.elements;
        const { count, offset, neighbor } = unit.bonds;
        const { model }  = unit;

        if (!count) continue;

        for (let j = 0; j < offset.length - 1; ++j) {
            const start = offset[j];
            const end = offset[j + 1];

            if (end <= start) continue;

            const aI = elements[j];
            for (let _bI = start; _bI < end; _bI++) {
                const bI = elements[neighbor[_bI]];
                console.log(`${atomLabel(model, aI)} -- ${atomLabel(model, bI)}`);
            }
        }
    }
}

export function printSequence(model: Model) {
    console.log('Sequence\n=============');
    const { byEntityKey } = model.sequence;
    for (const key of Object.keys(byEntityKey)) {
        const seq = byEntityKey[+key];
        console.log(`${seq.entityId} (${seq.num.value(0)}, ${seq.num.value(seq.num.rowCount - 1)}) (${seq.compId.value(0)}, ${seq.compId.value(seq.compId.rowCount - 1)})`);
        // for (let i = 0; i < seq.compId.rowCount; i++) {
        //     console.log(`${seq.entityId} ${seq.num.value(i)} ${seq.compId.value(i)}`);
        // }
    }
    console.log();
}

export function printUnits(structure: Structure) {
    console.log('Units\n=============');
    const l = Element.Location();

    for (const unit of structure.units) {
        l.unit = unit;
        const elements = unit.elements;
        const size = OrderedSet.size(elements);

        if (Unit.isAtomic(l.unit)) {
            console.log(`Atomic unit ${unit.id} ${unit.conformation.operator.name}: ${size} elements`);
        } else if (Unit.isCoarse(l.unit)) {
            console.log(`Coarse unit ${unit.id} ${unit.conformation.operator.name} (${Unit.isSpheres(l.unit) ? 'spheres' : 'gaussians'}): ${size} elements.`);

            const props = Queries.props.coarse;
            const seq = l.unit.model.sequence;

            for (let j = 0, _j = Math.min(size, 3); j < _j; j++) {
                l.element = OrderedSet.getAt(elements, j);

                const residues: string[] = [];
                const start = props.seq_id_begin(l), end = props.seq_id_end(l);
                const compId = seq.byEntityKey[props.entityKey(l)].compId.value;
                for (let e = start; e <= end; e++) residues.push(compId(e));
                console.log(`${props.asym_id(l)}:${start}-${end} (${residues.join('-')}) ${props.asym_id(l)} [${props.x(l).toFixed(2)}, ${props.y(l).toFixed(2)}, ${props.z(l).toFixed(2)}]`);
            }
            if (size > 3) console.log(`...`);
        }
    }
}


export function printIHMModels(model: Model) {
    if (!model.coarseHierarchy.isDefined) return false;
    console.log('IHM Models\n=============');
    console.log(Table.formatToString(model.coarseHierarchy.models));
}

async function run(mmcif: mmCIF_Database) {
    const models = await Model.create({ kind: 'mmCIF', data: mmcif }).run();
    const structure = Structure.ofModel(models[0]);
    printSequence(models[0]);
    printIHMModels(models[0]);
    printUnits(structure);
    // printBonds(structure);
    printSecStructure(models[0]);
}

async function runDL(pdb: string) {
    const mmcif = await downloadFromPdb(pdb)
    run(mmcif);
}

async function runFile(filename: string) {
    const mmcif = await readPdbFile(filename);
    run(mmcif);
}

const parser = new argparse.ArgumentParser({
  addHelp: true,
  description: 'Print info about a structure, mainly to test and showcase the mol-model module'
});
parser.addArgument([ '--download', '-d' ], {
    help: 'Pdb entry id'
});
parser.addArgument([ '--file', '-f' ], {
    help: 'filename'
});
interface Args {
    download?: string,
    file?: string
}
const args: Args = parser.parseArgs();

if (args.download) runDL(args.download)
else if (args.file) runFile(args.file)

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

// import { Table } from 'mol-data/db'
import CIF from 'mol-io/reader/cif'
import { Model, Structure, Element, ElementSet, Unit, ElementGroup, Queries } from 'mol-model/structure'
import { Run, Progress } from 'mol-task'
import { OrderedSet } from 'mol-data/int';
import { Table } from 'mol-data/db';
import { mmCIF_Database } from 'mol-io/reader/cif/schema/mmcif';
import CoarseGrained from 'mol-model/structure/model/properties/coarse-grained';

const readFileAsync = util.promisify(fs.readFile);

async function readFile(path: string) {
    if (path.match(/\.bcif$/)) {
        const input = await readFileAsync(path)
        const data = new Uint8Array(input.byteLength);
        for (let i = 0; i < input.byteLength; i++) data[i] = input[i];
        return data;
    } else {
        return readFileAsync(path, 'utf8');
    }
}

async function readCif(path: string) {
    const data = await readFile(path);
    const parsed = await parseCif(data);
    return CIF.schema.mmCIF(parsed.result.blocks[0])
}

async function parseCif(data: string|Uint8Array) {
    const comp = CIF.parse(data);
    const parsed = await Run(comp, p => console.log(Progress.format(p)), 250);
    if (parsed.isError) throw parsed;
    return parsed
}

async function getPdb(pdb: string) {
    //const data = await fetch(`https://files.rcsb.org/download/${pdb}.cif`)
    const data = await fetch(`http://www.ebi.ac.uk/pdbe/static/entry/${pdb}_updated.cif`);
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


export function printBonds(structure: Structure) {
    const { units, elements } = structure;
    const unitIds = ElementSet.unitIndices(elements);

    for (let i = 0, _i = OrderedSet.size(unitIds); i < _i; i++) {
        const unit = units[OrderedSet.getAt(unitIds, i)];
        const group = ElementSet.groupFromUnitIndex(elements, OrderedSet.getAt(unitIds, i));

        const { count, offset, neighbor } = Unit.getGroupBonds(unit, group);
        const { model }  = unit;

        if (!count) continue;

        for (let j = 0; j < offset.length - 1; ++j) {
            const start = offset[j];
            const end = offset[j + 1];

            if (end <= start) continue;

            const aI = ElementGroup.getAt(group, j);
            for (let _bI = start; _bI < end; _bI++) {
                const bI = ElementGroup.getAt(group, neighbor[_bI])
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
    const { elements, units } = structure;
    const unitIds = ElementSet.unitIndices(elements);
    const l = Element.Location();

    for (let i = 0, _i = unitIds.length; i < _i; i++) {
        const unitId = unitIds[i];
        l.unit = units[unitId];
        const set = ElementSet.groupAt(elements, i).elements;
        const size = OrderedSet.size(set);

        if (Unit.isAtomic(l.unit)) {
            console.log(`Atomic unit ${unitId}: ${size} elements`);
        } else if (Unit.isCoarse(l.unit)) {
            console.log(`Coarse unit ${unitId} (${l.unit.elementType === CoarseGrained.ElementType.Sphere ? 'spheres' : 'gaussians'}): ${size} elements.`);

            const props = Queries.props.coarse_grained;
            const seq = l.unit.model.sequence;

            for (let j = 0, _j = Math.min(size, 10); j < _j; j++) {
                l.element = OrderedSet.getAt(set, j);

                const residues: string[] = [];
                const start = props.seq_id_begin(l), end = props.seq_id_end(l);
                const compId = seq.byEntityKey[props.entityKey(l)].compId.value;
                for (let e = start; e <= end; e++) residues.push(compId(e));
                console.log(`${props.asym_id(l)}:${start}-${end} (${residues.join('-')}) ${props.asym_id(l)} [${props.x(l).toFixed(2)}, ${props.y(l).toFixed(2)}, ${props.z(l).toFixed(2)}]`);
            }
            if (size > 10) console.log(`...`);
        }
    }
}


export function printIHMModels(model: Model) {
    if (!model.coarseGrained.isDefined) return false;
    console.log('IHM Models\n=============');
    console.log(Table.formatToString(model.coarseGrained.modelList));
}

async function run(mmcif: mmCIF_Database) {
    const models = Model.create({ kind: 'mmCIF', data: mmcif });
    const structure = Structure.ofModel(models[0]);
    printSequence(models[0]);
    printIHMModels(models[0]);
    printUnits(structure);
}

async function runDL(pdb: string) {
    const mmcif = await getPdb(pdb)
    run(mmcif);
}

async function runFile(filename: string) {
    const mmcif = await readCif(filename);
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

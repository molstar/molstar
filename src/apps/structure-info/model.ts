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
import { UnitRings } from 'mol-model/structure/structure/unit/rings';


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
    console.log('\nSecondary Structure\n=============');
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

export function printLinks(structure: Structure, showIntra: boolean, showInter: boolean) {
    if (showIntra) {
        console.log('\nIntra Unit Links\n=============');
        for (const unit of structure.units) {
            if (!Unit.isAtomic(unit)) continue;

            const elements = unit.elements;
            const { a, b } = unit.links;
            const { model } = unit;

            if (!a.length) continue;

            for (let bI = 0, _bI = a.length; bI < _bI; bI++) {
                const x = a[bI], y = b[bI];
                if (x >= y) continue;
                console.log(`${atomLabel(model, elements[x])} -- ${atomLabel(model, elements[y])}`);
            }
        }
    }

    if (showInter) {
        console.log('\nInter Unit Links\n=============');
        const links = structure.links;
        for (const unit of structure.units) {
            if (!Unit.isAtomic(unit)) continue;

            for (const pairLinks of links.getLinkedUnits(unit)) {
                if (!pairLinks.areUnitsOrdered) continue;

                const { unitA, unitB } = pairLinks;

                console.log(`${pairLinks.unitA.id} - ${pairLinks.unitB.id}: ${pairLinks.bondCount} bond(s)`);

                for (const aI of pairLinks.linkedElementIndices) {
                    for (const link of pairLinks.getBonds(aI)) {
                        console.log(`${atomLabel(unitA.model, unitA.elements[aI])} -- ${atomLabel(unitB.model, unitB.elements[link.indexB])}`);
                    }
                }
            }
        }
    }
}

export function printSequence(model: Model) {
    console.log('\nSequence\n=============');
    const { byEntityKey } = model.sequence;
    for (const key of Object.keys(byEntityKey)) {
        const seq = byEntityKey[+key];
        console.log(`${seq.entityId} (${seq.sequence.kind} ${seq.num.value(0)} (offset ${seq.sequence.offset}), ${seq.num.value(seq.num.rowCount - 1)}) (${seq.compId.value(0)}, ${seq.compId.value(seq.compId.rowCount - 1)})`);
        console.log(`${seq.sequence.sequence}`);
    }
    console.log();
}

export function printRings(structure: Structure) {
    console.log('\nRings\n=============');
    for (const unit of structure.units) {
        if (!Unit.isAtomic(unit)) continue;
        const { all, byFingerprint } = unit.rings;
        const fps: string[] = [];
        for (let i = 0, _i = Math.min(5, all.length); i < _i; i++) {
            fps[fps.length] = UnitRings.getRingFingerprint(unit, all[i]);
        }
        if (all.length > 5) fps.push('...')
        console.log(`Unit ${unit.id}, ${all.length} ring(s), ${byFingerprint.size} different fingerprint(s).\n  ${fps.join(', ')}`);
    }
    console.log();
}

export function printUnits(structure: Structure) {
    console.log('\nUnits\n=============');
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
    console.log('\nIHM Models\n=============');
    console.log(Table.formatToString(model.coarseHierarchy.models));
}

async function run(mmcif: mmCIF_Database) {
    const models = await Model.create({ kind: 'mmCIF', data: mmcif }).run();
    const structure = Structure.ofModel(models[0]);
    printSequence(models[0]);
    //printIHMModels(models[0]);
    printUnits(structure);
    printRings(structure);
    printLinks(structure, false, true);
    //printSecStructure(models[0]);
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
parser.addArgument(['--download', '-d'], {
    help: 'Pdb entry id'
});
parser.addArgument(['--file', '-f'], {
    help: 'filename'
});
interface Args {
    download?: string,
    file?: string
}
const args: Args = parser.parseArgs();

if (args.download) runDL(args.download)
else if (args.file) runFile(args.file)

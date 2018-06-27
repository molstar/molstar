/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as argparse from 'argparse'
require('util.promisify').shim();

import { CifFrame } from 'mol-io/reader/cif'
import { Model, Structure, Element, Unit, Format, StructureProperties } from 'mol-model/structure'
// import { Run, Progress } from 'mol-task'
import { OrderedSet } from 'mol-data/int';
import { openCif, downloadCif } from './helpers';
import { BitFlags } from 'mol-util';
import { SecondaryStructureType } from 'mol-model/structure/model/types';
import { UnitRings } from 'mol-model/structure/structure/unit/rings';
import { Vec3 } from 'mol-math/linear-algebra';


async function downloadFromPdb(pdb: string) {
    // `https://files.rcsb.org/download/${pdb}.cif`
    const parsed = await downloadCif(`http://www.ebi.ac.uk/pdbe/static/entry/${pdb}_updated.cif`, false);
    return parsed.blocks[0];
}

async function readPdbFile(path: string) {
    const parsed = await openCif(path);
    return parsed.blocks[0];
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
            const { a, b, edgeCount } = unit.links;
            const { model } = unit;

            if (!edgeCount) continue;

            for (let bI = 0, _bI = edgeCount * 2; bI < _bI; bI++) {
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
                if (!pairLinks.areUnitsOrdered || pairLinks.bondCount === 0) continue;

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

export function printModRes(model: Model) {
    console.log('\nModified Residues\n=============');
    const map = model.properties.modifiedResidueNameMap;
    const { label_comp_id, _rowCount } = model.atomicHierarchy.residues;
    for (let i = 0; i < _rowCount; i++) {
        const comp_id = label_comp_id.value(i);
        if (!map.has(comp_id)) continue;
        console.log(`[${i}] ${map.get(comp_id)} -> ${comp_id}`);
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

            const props = StructureProperties.coarse;
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

export function printSymmetryInfo(model: Model) {
    console.log('\nSymmetry Info\n=============');
    const { symmetry } = model;
    const { size, anglesInRadians } = symmetry.spacegroup.cell;
    console.log(`Spacegroup: ${symmetry.spacegroup.name} size: ${Vec3.toString(size)} angles: ${Vec3.toString(anglesInRadians)}`);
    console.log(`Assembly names: ${symmetry.assemblies.map(a => a.id).join(', ')}`);
    // NCS example: 1auy
    console.log(`NCS operators: ${symmetry.ncsOperators && symmetry.ncsOperators.map(a => a.name).join(', ')}`);
}

export function printModelStats(models: ReadonlyArray<Model>) {
    console.log('\nModels\n=============');

    for (const m of models) {
        if (m.coarseHierarchy.isDefined) {
            console.log(`${m.label} ${m.modelNum}: ${m.atomicHierarchy.atoms._rowCount} atom(s), ${m.coarseHierarchy.spheres.count} sphere(s), ${m.coarseHierarchy.gaussians.count} gaussian(s)`);
        } else {
            console.log(`${m.label} ${m.modelNum}: ${m.atomicHierarchy.atoms._rowCount} atom(s)`);
        }
    }
    console.log();
}

async function run(frame: CifFrame, args: Args) {
    const models = await Model.create(Format.mmCIF(frame)).run();
    const structure = Structure.ofModel(models[0]);

    if (args.models) printModelStats(models);
    if (args.seq) printSequence(models[0]);
    if (args.units) printUnits(structure);
    if (args.sym) printSymmetryInfo(models[0]);
    if (args.rings) printRings(structure);
    if (args.intraLinks) printLinks(structure, true, false);
    if (args.interLinks) printLinks(structure, false, true);
    if (args.mod) printModRes(models[0]);
    if (args.sec) printSecStructure(models[0]);
}

async function runDL(pdb: string, args: Args) {
    const mmcif = await downloadFromPdb(pdb)
    run(mmcif, args);
}

async function runFile(filename: string, args: Args) {
    const mmcif = await readPdbFile(filename);
    run(mmcif, args);
}

const parser = new argparse.ArgumentParser({
    addHelp: true,
    description: 'Print info about a structure, mainly to test and showcase the mol-model module'
});
parser.addArgument(['--download', '-d'], { help: 'Pdb entry id' });
parser.addArgument(['--file', '-f'], { help: 'filename' });

parser.addArgument(['--models'], { help: 'print models info', action: 'storeTrue' });
parser.addArgument(['--seq'], { help: 'print sequence', action: 'storeTrue' });
parser.addArgument(['--ihm'], { help: 'print IHM', action: 'storeTrue' });
parser.addArgument(['--units'], { help: 'print units', action: 'storeTrue' });
parser.addArgument(['--sym'], { help: 'print symmetry', action: 'storeTrue' });
parser.addArgument(['--rings'], { help: 'print rings', action: 'storeTrue' });
parser.addArgument(['--intraLinks'], { help: 'print intra unit links', action: 'storeTrue' });
parser.addArgument(['--interLinks'], { help: 'print inter unit links', action: 'storeTrue' });
parser.addArgument(['--mod'], { help: 'print modified residues', action: 'storeTrue' });
parser.addArgument(['--sec'], { help: 'print secoundary structure', action: 'storeTrue' });
interface Args {
    download?: string,
    file?: string,

    models?:boolean,
    seq?: boolean,
    ihm?: boolean,
    units?: boolean,
    sym?: boolean,
    rings?: boolean,
    intraLinks?: boolean,
    interLinks?: boolean,
    mod?: boolean,
    sec?: boolean,
}
const args: Args = parser.parseArgs();

if (args.download) runDL(args.download, args)
else if (args.file) runFile(args.file, args)

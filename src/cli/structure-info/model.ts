#!/usr/bin/env node
/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as argparse from 'argparse';
require('util.promisify').shim();

import { CifFrame } from '../../mol-io/reader/cif';
import { Model, Structure, StructureElement, Unit, StructureProperties, UnitRing } from '../../mol-model/structure';
// import { Run, Progress } from '../../mol-task'
import { OrderedSet } from '../../mol-data/int';
import { openCif, downloadCif } from './helpers';
import { Vec3 } from '../../mol-math/linear-algebra';
import { trajectoryFromMmCIF } from '../../mol-model-formats/structure/mmcif';
import { Sequence } from '../../mol-model/sequence';
import { ModelSecondaryStructure } from '../../mol-model-formats/structure/property/secondary-structure';
import { ModelSymmetry } from '../../mol-model-formats/structure/property/symmetry';


async function downloadFromPdb(pdb: string) {
    // `https://files.rcsb.org/download/${pdb}.cif`
    const parsed = await downloadCif(`http://www.ebi.ac.uk/pdbe/static/entry/${pdb}_updated.cif`, false);
    return parsed.blocks[0];
}

export async function readCifFile(path: string) {
    const parsed = await openCif(path);
    return parsed.blocks[0];
}

export function atomLabel(model: Model, aI: number) {
    const { atoms, residues, chains, residueAtomSegments, chainAtomSegments } = model.atomicHierarchy;
    const { label_atom_id, label_comp_id } = atoms;
    const { label_seq_id } = residues;
    const { label_asym_id } = chains;
    const rI = residueAtomSegments.index[aI];
    const cI = chainAtomSegments.index[aI];
    return `${label_asym_id.value(cI)} ${label_comp_id.value(aI)} ${label_seq_id.value(rI)} ${label_atom_id.value(aI)}`;
}

export function residueLabel(model: Model, rI: number) {
    const { atoms, residues, chains, residueAtomSegments, chainAtomSegments } = model.atomicHierarchy;
    const { label_comp_id } = atoms;
    const { label_seq_id } = residues;
    const { label_asym_id } = chains;
    const aI = residueAtomSegments.offsets[rI];
    const cI = chainAtomSegments.index[aI];
    return `${label_asym_id.value(cI)} ${label_comp_id.value(aI)} ${label_seq_id.value(rI)}`;
}

export function printSecStructure(model: Model) {
    console.log('\nSecondary Structure\n=============');
    const { residues } = model.atomicHierarchy;
    const secondaryStructure = ModelSecondaryStructure.Provider.get(model);
    if (!secondaryStructure) return;

    const { key, elements } = secondaryStructure;
    const count = residues._rowCount;
    let rI = 0;
    while (rI < count) {
        let start = rI;
        while (rI < count && key[start] === key[rI]) rI++;
        rI--;

        const e = elements[key[start]];
        if (e.kind !== 'none') console.log(`${e.kind}: ${residueLabel(model, start)} - ${residueLabel(model, rI)}`);
        rI++;
    }
}

export function printBonds(structure: Structure, showIntra: boolean, showInter: boolean) {
    if (showIntra) {
        console.log('\nIntra Unit Bonds\n=============');
        for (const unit of structure.units) {
            if (!Unit.isAtomic(unit)) continue;

            const elements = unit.elements;
            const { a, b, edgeCount } = unit.bonds;
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
        console.log('\nInter Unit Bonds\n=============');
        const bonds = structure.interUnitBonds;
        for (const unit of structure.units) {
            if (!Unit.isAtomic(unit)) continue;

            for (const pairBonds of bonds.getConnectedUnits(unit)) {
                if (!pairBonds.areUnitsOrdered || pairBonds.edgeCount === 0) continue;

                const { unitA, unitB } = pairBonds;
                console.log(`${pairBonds.unitA.id} - ${pairBonds.unitB.id}: ${pairBonds.edgeCount} bond(s)`);

                for (const aI of pairBonds.connectedIndices) {
                    for (const bond of pairBonds.getEdges(aI)) {
                        console.log(`${atomLabel(unitA.model, unitA.elements[aI])} -- ${atomLabel(unitB.model, unitB.elements[bond.indexB])}`);
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
        const { sequence, entityId } = byEntityKey[+key];
        const { seqId, compId } = sequence;
        console.log(`${entityId} (${sequence.kind} ${seqId.value(0)}, ${seqId.value(seqId.rowCount - 1)}) (${compId.value(0)}, ${compId.value(compId.rowCount - 1)})`);
        console.log(`${Sequence.getSequenceString(sequence)}`);
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
            fps[fps.length] = UnitRing.fingerprint(unit, all[i]);
        }
        if (all.length > 5) fps.push('...');
        console.log(`Unit ${unit.id}, ${all.length} ring(s), ${byFingerprint.size} different fingerprint(s).\n  ${fps.join(', ')}`);
    }
    console.log();
}

export function printUnits(structure: Structure) {
    console.log('\nUnits\n=============');
    const l = StructureElement.Location.create(structure);

    for (const unit of structure.units) {
        l.unit = unit;
        const elements = unit.elements;
        const size = OrderedSet.size(elements);

        if (Unit.isAtomic(l.unit)) {
            console.log(`Atomic unit ${unit.id} ${unit.conformation.operator.name}: ${size} elements`);
        } else if (Unit.isCoarse(l.unit)) {
            console.log(`Coarse unit ${unit.id} ${unit.conformation.operator.name} (${Unit.isSpheres(l.unit) ? 'spheres' : 'gaussians'}): ${size} elements.`);

            const props = StructureProperties.coarse;
            const modelSeq = l.unit.model.sequence;

            for (let j = 0, _j = Math.min(size, 3); j < _j; j++) {
                l.element = OrderedSet.getAt(elements, j);

                const residues: string[] = [];
                const start = props.seq_id_begin(l), end = props.seq_id_end(l);
                const compId = modelSeq.byEntityKey[props.entityKey(l)].sequence.compId.value;
                for (let e = start; e <= end; e++) residues.push(compId(e));
                console.log(`${props.asym_id(l)}:${start}-${end} (${residues.join('-')}) ${props.asym_id(l)} [${props.x(l).toFixed(2)}, ${props.y(l).toFixed(2)}, ${props.z(l).toFixed(2)}]`);
            }
            if (size > 3) console.log(`...`);
        }
    }
}

export function printSymmetryInfo(model: Model) {
    console.log('\nSymmetry Info\n=============');
    const symmetry = ModelSymmetry.Provider.get(model);
    if (!symmetry) return;
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

export async function getModelsAndStructure(frame: CifFrame) {
    const models = await trajectoryFromMmCIF(frame).run();
    const structure = Structure.ofModel(models[0]);
    return { models, structure };
}

async function run(frame: CifFrame, args: Args) {
    const { models, structure } = await getModelsAndStructure(frame);

    if (args.models) printModelStats(models);
    if (args.seq) printSequence(models[0]);
    if (args.units) printUnits(structure);
    if (args.sym) printSymmetryInfo(models[0]);
    if (args.rings) printRings(structure);
    if (args.intraBonds) printBonds(structure, true, false);
    if (args.interBonds) printBonds(structure, false, true);
    if (args.sec) printSecStructure(models[0]);
}

async function runDL(pdb: string, args: Args) {
    const mmcif = await downloadFromPdb(pdb);
    run(mmcif, args);
}

async function runFile(filename: string, args: Args) {
    const mmcif = await readCifFile(filename);
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
parser.addArgument(['--units'], { help: 'print units', action: 'storeTrue' });
parser.addArgument(['--sym'], { help: 'print symmetry', action: 'storeTrue' });
parser.addArgument(['--rings'], { help: 'print rings', action: 'storeTrue' });
parser.addArgument(['--intraBonds'], { help: 'print intra unit bonds', action: 'storeTrue' });
parser.addArgument(['--interBonds'], { help: 'print inter unit bonds', action: 'storeTrue' });
parser.addArgument(['--mod'], { help: 'print modified residues', action: 'storeTrue' });
parser.addArgument(['--sec'], { help: 'print secoundary structure', action: 'storeTrue' });
interface Args {
    download?: string,
    file?: string,

    models?: boolean,
    seq?: boolean,
    ihm?: boolean,
    units?: boolean,
    sym?: boolean,
    rings?: boolean,
    intraBonds?: boolean,
    interBonds?: boolean,
    mod?: boolean,
    sec?: boolean,
}
const args: Args = parser.parseArgs();

if (args.download) runDL(args.download, args);
else if (args.file) runFile(args.file, args);

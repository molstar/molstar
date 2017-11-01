/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as B from 'benchmark'

import * as util from 'util'
import * as fs from 'fs'
import CIF from '../mol-io/reader/cif'

import buildModels from '../mol-data/model/builders/mmcif'
import { ofModel } from '../mol-data/structure/base'
import Model from '../mol-data/Model'
import Structure from '../mol-data/structure'
import OrdSet from '../mol-base/collections/integer/ordered-set'
import Atom from '../mol-data/structure/atom'
import AtomSet from '../mol-data/structure/atom-set'
import Segmentation from '../mol-base/collections/integer/segmentation'

require('util.promisify').shim();
const readFileAsync = util.promisify(fs.readFile);

async function readData(path: string) {
    if (path.match(/\.bcif$/)) {
        const input = await readFileAsync(path)
        const data = new Uint8Array(input.byteLength);
        for (let i = 0; i < input.byteLength; i++) data[i] = input[i];
        return data;
    } else {
        return readFileAsync(path, 'utf8');
    }
}

export async function readCIF(path: string) {
    console.time('readData');
    const input = await readData(path)
    console.timeEnd('readData');

    console.time('parse');
    const comp = typeof input === 'string' ? CIF.parseText(input) : CIF.parseBinary(input);
    const parsed = await comp();
    console.timeEnd('parse');
    if (parsed.isError) {
        throw parsed;
    }

    const data = parsed.result.blocks[0];
    console.time('schema')
    const mmcif = CIF.schema.mmCIF(data);
    console.timeEnd('schema')
    console.time('buildModels')
    const models = buildModels(mmcif);
    console.timeEnd('buildModels')
    const structures = models.map(ofModel);

    return { mmcif, models, structures };
}

export namespace PropertyAccess {
    function baselineRaw(model: Model) {
        const atom_site = model.sourceData.data._frame.categories['_atom_site'];
        const id = atom_site.getField('id')!.int;
        let s = 0;
        for (let i = 0, _i = atom_site.rowCount; i < _i; i++) {
            s += id(i);
        }
        return s;
    }

    function baseline(model: Model) {
        const atom_site = model.sourceData.data.atom_site;
        const id = atom_site.id.value;
        let s = 0;
        for (let i = 0, _i = atom_site._rowCount; i < _i; i++) {
            s += id(i);
        }
        return s;
    }

    function sumProperty(structure: Structure, p: Atom.Property<number>) {
        const { atoms, units } = structure;
        const unitIds = AtomSet.unitIds(atoms);
        const l = Atom.Location();

        let s = 0;

        for (let i = 0, _i = unitIds.length; i < _i; i++) {
            l.unit = units[unitIds[i]];
            const set = AtomSet.unitGetByIndex(atoms, i);

            for (let j = 0, _j = OrdSet.size(set); j < _j; j++) {
                l.atom = OrdSet.getAt(set, j);
                s += p(l);
            }
        }

        return s;
    }

    function sumPropertySegmented(structure: Structure, p: Atom.Property<number>) {
        const { atoms, units } = structure;
        const unitIds = AtomSet.unitIds(atoms);
        const l = Atom.Location();

        let s = 0;

        for (let i = 0, _i = unitIds.length; i < _i; i++) {
            const unit = units[unitIds[i]];
            l.unit = unit;
            const set = AtomSet.unitGetByIndex(atoms, i);

            const chainsIt = Segmentation.transientSegments(unit.hierarchy.chainSegments, set);
            const residues = unit.hierarchy.residueSegments;
            while (chainsIt.hasNext) {
                const chainSegment = chainsIt.move();
                const residuesIt = Segmentation.transientSegments(residues, set, chainSegment);
                while (residuesIt.hasNext) {
                    const residueSegment = residuesIt.move();
                    for (let j = residueSegment.start, _j = residueSegment.end; j < _j; j++) {
                        l.atom = OrdSet.getAt(set, j);
                        s += p(l);
                    }
                }
            }
        }

        return s;
    }

    function sumPropertyResidue(structure: Structure, p: Atom.Property<number>) {
        const { atoms, units } = structure;
        const unitIds = AtomSet.unitIds(atoms);
        const l = Atom.Location();

        let s = 0;

        for (let i = 0, _i = unitIds.length; i < _i; i++) {
            const unit = units[unitIds[i]];
            l.unit = unit;
            const set = AtomSet.unitGetByIndex(atoms, i);
            const residuesIt = Segmentation.transientSegments(unit.hierarchy.residueSegments, set);
            while (residuesIt.hasNext) {
                l.atom = OrdSet.getAt(set, residuesIt.move().start);
                s += p(l);
            }
        }

        return s;
    }

    function sumPropertyAtomSetIt(structure: Structure, p: Atom.Property<number>) {
        const { atoms, units } = structure;

        let s = 0;
        const atomsIt = AtomSet.atoms(atoms);
        const l = Atom.Location();
        while (atomsIt.hasNext) {
            const a = atomsIt.move();
            l.unit = units[Atom.unit(a)];
            l.atom = Atom.index(a);
            s += p(l);
        }
        return s;
    }

    // function sumPropertySegmentedMutable(structure: Structure, p: Property<number>) {
    //     const { atoms, units } = structure;
    //     const unitIds = AtomSet.unitIds(atoms);
    //     const l = Property.createLocation();

    //     let s = 0;

    //     for (let i = 0, _i = unitIds.length; i < _i; i++) {
    //         const unit = units[unitIds[i]];
    //         l.unit = unit;
    //         const set = AtomSet.unitGetByIndex(atoms, i);

    //         const chainsIt = Segmentation.transientSegments(unit.hierarchy.chainSegments, set);
    //         const residuesIt = Segmentation.transientSegments(unit.hierarchy.residueSegments, set);
    //         while (chainsIt.hasNext) {
    //             const chainSegment = chainsIt.move();
    //             residuesIt.updateRange(chainSegment);
    //             while (residuesIt.hasNext) {
    //                 const residueSegment = residuesIt.move();
    //                 for (let j = residueSegment.start, _j = residueSegment.end; j < _j; j++) {
    //                     l.atom = OrdSet.getAt(set, j);
    //                     s += p(l);
    //                 }
    //             }
    //         }
    //     }

    //     return s;
    // }

    function sumDirect(structure: Structure) {
        const { atoms, units } = structure;
        const unitIds = AtomSet.unitIds(atoms);

        let s = 0;

        for (let i = 0, _i = unitIds.length; i < _i; i++) {
            const unitId = unitIds[i];
            const unit = units[unitId];
            const set = AtomSet.unitGetByIndex(atoms, i);
            //const { residueIndex, chainIndex } = unit;
            const p = unit.conformation.atomId.value;
            for (let j = 0, _j = OrdSet.size(set); j < _j; j++) {
                const aI = OrdSet.getAt(set, j);
                s += p(aI);
            }
        }

        return s;
    }

    // function concatProperty(structure: Structure, p: Property<string>) {
    //     const { atoms, units } = structure;
    //     const unitIds = AtomSet.unitIds(atoms);
    //     const l = Property.createLocation(structure);

    //     let s = [];

    //     for (let i = 0, _i = unitIds.length; i < _i; i++) {
    //         const unitId = unitIds[i];
    //         l.unit = units[unitId];
    //         const set = AtomSet.unitGetByIndex(atoms, i);
    //         const { residueIndex, chainIndex } = l.unit;

    //         for (let j = 0, _j = OrdSet.size(set); j < _j; j++) {
    //             const aI = OrdSet.getAt(set, j);
    //             l.atom = aI;
    //             l.residueIndex = residueIndex[aI];
    //             l.chainIndex = chainIndex[aI];
    //             s[s.length] = p(l);
    //         }
    //     }

    //     return s;
    // }

    export async function run() {
        //const { structures, models } = await readCIF('./examples/1cbs_full.bcif');
        const { structures, models } = await readCIF('e:/test/quick/1jj2_full.bcif');
        //const { structures, models } = await readCIF('e:/test/quick/3j3q_updated.cif');

        console.log('parsed');

        console.log(baseline(models[0]));
        console.log(baselineRaw(models[0]));
        console.log(sumProperty(structures[0], l => l.unit.model.conformation.atomId.value(l.atom)));
        console.log(sumPropertySegmented(structures[0], l => l.unit.model.conformation.atomId.value(l.atom)));
        //console.log(sumPropertySegmentedMutable(structures[0], l => l.unit.model.conformation.atomId.value(l.atom)));
        console.log(sumPropertyAtomSetIt(structures[0], l => l.unit.model.conformation.atomId.value(l.atom)));
        //console.log(sumProperty(structures[0], Property.cachedAtomColumn(m => m.conformation.atomId)));
        console.log(sumDirect(structures[0]));
        console.log('r', sumPropertyResidue(structures[0], l => l.unit.hierarchy.residues.auth_seq_id.value(l.unit.residueIndex[l.atom])));

        //const col = models[0].conformation.atomId.value;
        const suite = new B.Suite();
        suite
            //.add('test int', () => sumProperty(structures[0], l => col(l.atom)))
            // .add('baseline raw', () =>  baselineRaw(models[0]))
            .add('sum residue', () => sumPropertyResidue(structures[0], l => l.unit.hierarchy.residues.auth_seq_id.value(l.unit.residueIndex[l.atom])))

            .add('baseline', () =>  baseline(models[0]))
            .add('direct', () =>  sumDirect(structures[0]))
            //.add('normal int', () => sumProperty(structures[0], l => l.unit.model.conformation.atomId.value(l.atom)))
            //.add('atom set it int', () => sumPropertyAtomSetIt(structures[0], l => l.unit.conformation.atomId.value(l.atom)))
            // .add('segmented faster int', () => sumPropertySegmented(structures[0], l => l.unit.conformation.atomId.value(l.atom)))
            // .add('faster int', () => sumProperty(structures[0], l => l.unit.conformation.atomId.value(l.atom)))
            .add('segmented faster _x', () => sumPropertySegmented(structures[0], l => l.unit.conformation.__x[l.atom]))
            .add('faster _x', () => sumProperty(structures[0], l => l.unit.conformation.__x[l.atom] +  l.unit.conformation.__y[l.atom] +  l.unit.conformation.__z[l.atom]))
            //.add('segmented mut faster int', () => sumPropertySegmentedMutable(structures[0], l => l.unit.conformation.atomId.value(l.atom)))
            //.add('normal shortcut int', () => sumProperty(structures[0], l => l.conformation.atomId.value(l.atom)))
            //.add('cached int', () => sumProperty(structures[0], Property.cachedAtomColumn(m => m.conformation.atomId)))
            //.add('concat str', () => concatProperty(structures[0], l => l.unit.model.hierarchy.atoms.auth_atom_id.value(l.atom)))
            //.add('cached concat str', () => concatProperty(structures[0], Property.cachedAtomColumn(m => m.hierarchy.atoms.auth_atom_id)))
            .on('cycle', (e: any) => console.log(String(e.target)))
            .run();
    }
}

PropertyAccess.run();
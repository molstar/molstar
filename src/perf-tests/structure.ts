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
import Property from '../mol-data/structure/property'
import Model from '../mol-data/Model'
import Structure from '../mol-data/structure'
import OrdSet from '../mol-base/collections/integer/ordered-set'
import AtomSet from '../mol-data/structure/atom-set'

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
    const input = await readData(path)
    const comp = typeof input === 'string' ? CIF.parseText(input) : CIF.parseBinary(input);
    const parsed = await comp();
    if (parsed.isError) {
        throw parsed;
    }

    const data = parsed.result.blocks[0];
    const mmcif = CIF.schema.mmCIF(data);
    const models = buildModels(mmcif);
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

    function sumProperty(structure: Structure, p: Property<number>, initial: number) {
        const { atoms, units } = structure;
        const unitIds = AtomSet.unitIds(atoms);
        const l = Property.createLocation(structure);

        let s = initial;

        for (let i = 0, _i = unitIds.length; i < _i; i++) {
            const unitId = unitIds[i];
            l.unit = units[unitId];
            const set = AtomSet.unitGetByIndex(atoms, i);
            const { residueIndex, chainIndex } = l.unit;

            for (let j = 0, _j = OrdSet.size(set); j < _j; j++) {
                const aI = OrdSet.getAt(set, j);
                l.atomIndex = aI;
                l.residueIndex = residueIndex[aI];
                l.chainIndex = chainIndex[aI];
                s += p(l);
            }
        }

        return s;
    }

    function sumDirect(structure: Structure) {
        const { atoms, units } = structure;
        const unitIds = AtomSet.unitIds(atoms);

        let s = 0;

        for (let i = 0, _i = unitIds.length; i < _i; i++) {
            const unitId = unitIds[i];
            const unit = units[unitId];
            const set = AtomSet.unitGetByIndex(atoms, i);
            //const { residueIndex, chainIndex } = unit;
            const p = unit.model.conformation.atomId.value;
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
    //             l.atomIndex = aI;
    //             l.residueIndex = residueIndex[aI];
    //             l.chainIndex = chainIndex[aI];
    //             s[s.length] = p(l);
    //         }
    //     }

    //     return s;
    // }

    export async function run() {
        const { structures, models } = await readCIF('./examples/1cbs_full.bcif');

        console.log(baseline(models[0]));
        console.log(baselineRaw(models[0]));
        console.log(sumProperty(structures[0], l => l.unit.model.conformation.atomId.value(l.atomIndex), 0));
        console.log(sumProperty(structures[0], Property.cachedAtomColumn(m => m.conformation.atomId), 0));
        console.log(sumDirect(structures[0]));

        const col = models[0].conformation.atomId.value;
        const suite = new B.Suite();
        suite
            .add('baseline raw', () =>  baselineRaw(models[0]))
            .add('baseline', () =>  baseline(models[0]))
            .add('direct', () =>  sumDirect(structures[0]))
            .add('normal int', () => sumProperty(structures[0], l => l.unit.model.conformation.atomId.value(l.atomIndex), 0))
            .add('test int', () => sumProperty(structures[0], l => col(l.atomIndex), 0))
            .add('cached int', () => sumProperty(structures[0], Property.cachedAtomColumn(m => m.conformation.atomId), 0))
            //.add('concat str', () => concatProperty(structures[0], l => l.unit.model.hierarchy.atoms.auth_atom_id.value(l.atomIndex)))
            //.add('cached concat str', () => concatProperty(structures[0], Property.cachedAtomColumn(m => m.hierarchy.atoms.auth_atom_id)))
            .on('cycle', (e: any) => console.log(String(e.target)))
            .run();
    }
}

PropertyAccess.run();
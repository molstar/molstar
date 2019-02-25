/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CustomElementProperty } from 'mol-model-props/common/custom-element-property';
import { Model, ElementIndex, ResidueIndex } from 'mol-model/structure';
import { Color } from 'mol-util/color';

// export const StripedResidues = CustomElementProperty.create<number>({
//     isStatic: true,
//     name: 'basic-wrapper-residue-striping',
//     display: 'Residue Stripes',
//     getData(model: Model) {
//         const map = new Map<ElementIndex, number>();
//         const residueIndex = model.atomicHierarchy.residueAtomSegments.index;
//         for (let i = 0, _i = model.atomicHierarchy.atoms._rowCount; i < _i; i++) {
//             map.set(i as ElementIndex, residueIndex[i] % 2);
//         }
//         return map;
//     },
//     coloring: {
//         getColor(e) { return e === 0 ? Color(0xff0000) : Color(0x0000ff) },
//         defaultColor: Color(0x777777)
//     },
//     format(e) {
//         return e === 0 ? 'Odd stripe' : 'Even stripe'
//     }
// });

const EvolutionaryConservationPalette: Color[] = [
    [255, 255, 150], // 9
    [160, 37, 96],
    [240, 125, 171],
    [250, 201, 222],
    [252, 237, 244],
    [255, 255, 255],
    [234, 255, 255],
    [215, 255, 255],
    [140, 255, 255],
    [16, 200, 209] // 1
].reverse().map(([r, g, b]) => Color.fromRgb(r, g, b));

export const EvolutionaryConservation = CustomElementProperty.create<number>({
    isStatic: true,
    name: 'proteopedia-wrapper-evolutionary-conservation',
    display: 'Evolutionary Conservation',
    async getData(model: Model) {
        const id = model.label.toLowerCase();
        const req = await fetch(`https://proteopedia.org/cgi-bin/cnsrf?${id}`);
        const json = await req.json();
        const annotations = (json && json.residueAnnotations) || [];

        const conservationMap = new Map<string, number>();

        for (const e of annotations) {
            for (const r of e.ids) {
                conservationMap.set(r, e.annotation);
            }
        }

        const map = new Map<ElementIndex, number>();

        const { _rowCount: residueCount } = model.atomicHierarchy.residues;
        const { offsets: residueOffsets } = model.atomicHierarchy.residueAtomSegments;
        const chainIndex = model.atomicHierarchy.chainAtomSegments.index;

        for (let rI = 0 as ResidueIndex; rI < residueCount; rI++) {
            const cI = chainIndex[residueOffsets[rI]];
            const key = `${model.atomicHierarchy.chains.auth_asym_id.value(cI)} ${model.atomicHierarchy.residues.auth_seq_id.value(rI)}`;
            if (!conservationMap.has(key)) continue;
            const ann = conservationMap.get(key)!;
            for (let aI = residueOffsets[rI]; aI < residueOffsets[aI + 1]; aI++) {
                map.set(aI, ann);
            }
        }

        return map;
    },
    coloring: {
        getColor(e) { return EvolutionaryConservationPalette[(e - 1) || 0]; },
        defaultColor: Color(0x999999)
    },
    format(e) {
        return e ? `Evolutionary Conservation ${e}` : void 0;
    }
});
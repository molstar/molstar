/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../mol-model/loci';
import { StructureElement, Structure } from '../mol-model/structure';
import { Script } from '../mol-script/script';
import BitFlags from '../mol-util/bit-flags';

export { Clipping };

type Clipping = { readonly layers: ReadonlyArray<Clipping.Layer> }

function Clipping(layers: ReadonlyArray<Clipping.Layer>): Clipping {
    return { layers };
}

namespace Clipping {
    export type Layer = { readonly loci: StructureElement.Loci, readonly groups: Groups }
    export const Empty: Clipping = { layers: [] };

    export type Groups = BitFlags<Groups.Flag>
    export namespace Groups {
        export const is: (g: Groups, f: Flag) => boolean = BitFlags.has;
        export const enum Flag {
            None = 0x0,
            One = 0x1,
            Two = 0x2,
            Three = 0x4,
            Four = 0x8,
            Five = 0x10,
            Six = 0x20,
        }

        export function create(flags: Flag): Groups {
            return BitFlags.create(flags);
        }

        export const Names = {
            'one': Flag.One,
            'two': Flag.Two,
            'three': Flag.Three,
            'four': Flag.Four,
            'five': Flag.Five,
            'six': Flag.Six,
        };
        export type Names = keyof typeof Names

        export function isName(name: string): name is Names {
            return name in Names;
        }

        export function fromName(name: Names): Flag {
            switch (name) {
                case 'one': return Flag.One;
                case 'two': return Flag.Two;
                case 'three': return Flag.Three;
                case 'four': return Flag.Four;
                case 'five': return Flag.Five;
                case 'six': return Flag.Six;
            }
        }

        export function fromNames(names: Names[]): Flag {
            let f = Flag.None;
            for (let i = 0, il = names.length; i < il; ++i) {
                f |= fromName(names[i]);
            }
            return f;
        }

        export function toNames(groups: Groups): Names[] {
            let names: Names[] = [];
            if (is(groups, Flag.One)) names.push('one');
            if (is(groups, Flag.Two)) names.push('two');
            if (is(groups, Flag.Three)) names.push('three');
            if (is(groups, Flag.Four)) names.push('four');
            if (is(groups, Flag.Five)) names.push('five');
            if (is(groups, Flag.Six)) names.push('six');
            return names;
        }
    }

    /** Clip object types */
    export const Type = {
        none: 0, // to switch clipping off
        plane: 1,
        sphere: 2,
        cube: 3,
        cylinder: 4,
        infiniteCone: 5,
    };

    export type Variant = 'instance' | 'pixel'

    export function areEqual(cA: Clipping, cB: Clipping) {
        if (cA.layers.length === 0 && cB.layers.length === 0) return true;
        if (cA.layers.length !== cB.layers.length) return false;
        for (let i = 0, il = cA.layers.length; i < il; ++i) {
            if (cA.layers[i].groups !== cB.layers[i].groups) return false;
            if (!Loci.areEqual(cA.layers[i].loci, cB.layers[i].loci)) return false;
        }
        return true;
    }

    export function isEmpty(clipping: Clipping) {
        return clipping.layers.length === 0;
    }

    export function remap(clipping: Clipping, structure: Structure) {
        const layers: Clipping.Layer[] = [];
        for (const layer of clipping.layers) {
            let { loci, groups } = layer;
            loci = StructureElement.Loci.remap(loci, structure);
            if (!StructureElement.Loci.isEmpty(loci)) {
                layers.push({ loci, groups });
            }
        }
        return { layers };
    }

    export function merge(clipping: Clipping): Clipping {
        if (isEmpty(clipping)) return clipping;
        const { structure } = clipping.layers[0].loci;
        const map = new Map<Groups, StructureElement.Loci>();
        let shadowed = StructureElement.Loci.none(structure);
        for (let i = 0, il = clipping.layers.length; i < il; ++i) {
            let { loci, groups } = clipping.layers[il - i - 1]; // process from end
            loci = StructureElement.Loci.subtract(loci, shadowed);
            shadowed = StructureElement.Loci.union(loci, shadowed);
            if (!StructureElement.Loci.isEmpty(loci)) {
                if (map.has(groups)) {
                    loci = StructureElement.Loci.union(loci, map.get(groups)!);
                }
                map.set(groups, loci);
            }
        }
        const layers: Clipping.Layer[] = [];
        map.forEach((loci, groups) => {
            layers.push({ loci, groups });
        });
        return { layers };
    }

    export function filter(clipping: Clipping, filter: Structure): Clipping {
        if (isEmpty(clipping)) return clipping;
        const { structure } = clipping.layers[0].loci;
        const layers: Clipping.Layer[] = [];
        for (const layer of clipping.layers) {
            let { loci, groups } = layer;
            // filter by first map to the `filter` structure and
            // then map back to the original structure of the clipping loci
            const filtered = StructureElement.Loci.remap(loci, filter);
            loci = StructureElement.Loci.remap(filtered, structure);
            if (!StructureElement.Loci.isEmpty(loci)) {
                layers.push({ loci, groups });
            }
        }
        return { layers };
    }

    export type ScriptLayer = { script: Script, groups: Groups }
    export function ofScript(scriptLayers: ScriptLayer[], structure: Structure): Clipping {
        const layers: Clipping.Layer[] = [];
        for (let i = 0, il = scriptLayers.length; i < il; ++i) {
            const { script, groups } = scriptLayers[i];
            const loci = Script.toLoci(script, structure);
            if (!StructureElement.Loci.isEmpty(loci)) {
                layers.push({ loci, groups });
            }
        }
        return { layers };
    }

    export type BundleLayer = { bundle: StructureElement.Bundle, groups: Groups }
    export function ofBundle(bundleLayers: BundleLayer[], structure: Structure): Clipping {
        const layers: Clipping.Layer[] = [];
        for (let i = 0, il = bundleLayers.length; i < il; ++i) {
            const { bundle, groups } = bundleLayers[i];
            const loci = StructureElement.Bundle.toLoci(bundle, structure.root);
            layers.push({ loci, groups });
        }
        return { layers };
    }

    export function toBundle(clipping: Clipping) {
        const layers: BundleLayer[] = [];
        for (let i = 0, il = clipping.layers.length; i < il; ++i) {
            let { loci, groups } = clipping.layers[i];
            const bundle = StructureElement.Bundle.fromLoci(loci);
            layers.push({ bundle, groups });
        }
        return { layers };
    }
}
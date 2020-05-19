/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../mol-model/loci';
import { StructureElement, Structure } from '../mol-model/structure';
import { Script } from '../mol-script/script';

export { Transparency };

type Transparency = { readonly layers: ReadonlyArray<Transparency.Layer> }

function Transparency(layers: ReadonlyArray<Transparency.Layer>): Transparency {
    return { layers };
}

namespace Transparency {
    export type Layer = { readonly loci: StructureElement.Loci, readonly value: number }
    export const Empty: Transparency = { layers: [] };

    export type Variant = 'single' | 'multi'

    export function areEqual(tA: Transparency, tB: Transparency) {
        if (tA.layers.length === 0 && tB.layers.length === 0) return true;
        if (tA.layers.length !== tB.layers.length) return false;
        for (let i = 0, il = tA.layers.length; i < il; ++i) {
            if (tA.layers[i].value !== tB.layers[i].value) return false;
            if (!Loci.areEqual(tA.layers[i].loci, tB.layers[i].loci)) return false;
        }
        return true;
    }

    export function isEmpty(transparency: Transparency) {
        return transparency.layers.length === 0;
    }

    export function remap(transparency: Transparency, structure: Structure) {
        const layers: Transparency.Layer[] = [];
        for (const layer of transparency.layers) {
            let { loci, value } = layer;
            loci = StructureElement.Loci.remap(loci, structure);
            if (!StructureElement.Loci.isEmpty(loci)) {
                layers.push({ loci, value });
            }
        }
        return { layers };
    }

    export function merge(transparency: Transparency): Transparency {
        if (isEmpty(transparency)) return transparency;
        const { structure } = transparency.layers[0].loci;
        const map = new Map<number, StructureElement.Loci>();
        let shadowed = StructureElement.Loci.none(structure);
        for (let i = 0, il = transparency.layers.length; i < il; ++i) {
            let { loci, value } = transparency.layers[il - i - 1]; // process from end
            loci = StructureElement.Loci.subtract(loci, shadowed);
            shadowed = StructureElement.Loci.union(loci, shadowed);
            if (!StructureElement.Loci.isEmpty(loci)) {
                if (map.has(value)) {
                    loci = StructureElement.Loci.union(loci, map.get(value)!);
                }
                map.set(value, loci);
            }
        }
        const layers: Transparency.Layer[] = [];
        map.forEach((loci, value) => {
            layers.push({ loci, value });
        });
        return { layers };
    }

    export function filter(transparency: Transparency, filter: Structure): Transparency {
        if (isEmpty(transparency)) return transparency;
        const { structure } = transparency.layers[0].loci;
        const layers: Transparency.Layer[] = [];
        for (const layer of transparency.layers) {
            let { loci, value } = layer;
            // filter by first map to the `filter` structure and
            // then map back to the original structure of the transparency loci
            const filtered = StructureElement.Loci.remap(loci, filter);
            loci = StructureElement.Loci.remap(filtered, structure);
            if (!StructureElement.Loci.isEmpty(loci)) {
                layers.push({ loci, value });
            }
        }
        return { layers };
    }

    export type ScriptLayer = { script: Script, value: number }
    export function ofScript(scriptLayers: ScriptLayer[], structure: Structure): Transparency {
        const layers: Transparency.Layer[] = [];
        for (let i = 0, il = scriptLayers.length; i < il; ++i) {
            const { script, value } = scriptLayers[i];
            const loci = Script.toLoci(script, structure);
            if (!StructureElement.Loci.isEmpty(loci)) {
                layers.push({ loci, value });
            }
        }
        return { layers };
    }

    export type BundleLayer = { bundle: StructureElement.Bundle, value: number }
    export function ofBundle(bundleLayers: BundleLayer[], structure: Structure): Transparency {
        const layers: Transparency.Layer[] = [];
        for (let i = 0, il = bundleLayers.length; i < il; ++i) {
            const { bundle, value } = bundleLayers[i];
            const loci = StructureElement.Bundle.toLoci(bundle, structure.root);
            layers.push({ loci, value });
        }
        return { layers };
    }

    export function toBundle(transparency: Transparency) {
        const layers: BundleLayer[] = [];
        for (let i = 0, il = transparency.layers.length; i < il; ++i) {
            let { loci, value } = transparency.layers[i];
            const bundle = StructureElement.Bundle.fromLoci(loci);
            layers.push({ bundle, value });
        }
        return { layers };
    }
}
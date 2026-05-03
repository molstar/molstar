/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../mol-model/loci';
import { Structure, StructureElement } from '../mol-model/structure';
import { Script } from '../mol-script/script';

export { Wiggle };

type Wiggle<T extends Loci = Loci> = {
    readonly kind: T['kind']
    readonly layers: ReadonlyArray<Wiggle.Layer<T>>
}

function Wiggle<T extends Loci>(kind: T['kind'], layers: ReadonlyArray<Wiggle.Layer<T>>): Wiggle {
    return { kind, layers };
}

namespace Wiggle {
    export type Layer<T extends Loci = Loci> = { readonly loci: T, readonly value: number }
    export const Empty: Wiggle = { kind: 'empty-loci', layers: [] };

    export function areEqual(wA: Wiggle, wB: Wiggle) {
        if (wA.layers.length === 0 && wB.layers.length === 0) return true;
        if (wA.layers.length !== wB.layers.length) return false;
        for (let i = 0, il = wA.layers.length; i < il; ++i) {
            if (wA.layers[i].value !== wB.layers[i].value) return false;
            if (!Loci.areEqual(wA.layers[i].loci, wB.layers[i].loci)) return false;
        }
        return true;
    }

    export function isEmpty(wiggle: Wiggle) {
        return wiggle.layers.length === 0;
    }

    export function remap(wiggle: Wiggle, structure: Structure): Wiggle {
        if (wiggle.kind === 'element-loci') {
            const layers: Wiggle.Layer[] = [];
            for (const layer of wiggle.layers) {
                let { loci, value } = layer;
                loci = StructureElement.Loci.remap(loci as StructureElement.Loci, structure);
                if (!StructureElement.Loci.isEmpty(loci)) {
                    layers.push({ loci, value });
                }
            }
            return { kind: 'element-loci', layers };
        } else {
            return wiggle;
        }
    }

    export function merge(wiggle: Wiggle): Wiggle {
        if (isEmpty(wiggle)) return wiggle;
        if (wiggle.kind === 'element-loci') {
            const { structure } = wiggle.layers[0].loci as StructureElement.Loci;
            const map = new Map<number, StructureElement.Loci>();
            let shadowed = StructureElement.Loci.none(structure);
            for (let i = 0, il = wiggle.layers.length; i < il; ++i) {
                let { loci, value } = wiggle.layers[il - i - 1]; // process from end
                loci = StructureElement.Loci.subtract(loci as StructureElement.Loci, shadowed);
                shadowed = StructureElement.Loci.union(loci, shadowed);
                if (!StructureElement.Loci.isEmpty(loci)) {
                    if (map.has(value)) {
                        loci = StructureElement.Loci.union(loci, map.get(value)!);
                    }
                    map.set(value, loci);
                }
            }
            const layers: Wiggle.Layer<StructureElement.Loci>[] = [];
            map.forEach((loci, value) => {
                layers.push({ loci, value });
            });
            return { kind: 'element-loci', layers };
        } else {
            return wiggle;
        }
    }

    export function filter(wiggle: Wiggle, filter: Structure): Wiggle {
        if (isEmpty(wiggle)) return wiggle;
        if (wiggle.kind === 'element-loci') {
            const { structure } = wiggle.layers[0].loci as StructureElement.Loci;
            const layers: Wiggle.Layer[] = [];
            for (const layer of wiggle.layers) {
                let { loci, value } = layer;
                // filter by first map to the `filter` structure and
                // then map back to the original structure of the wiggle loci
                const filtered = StructureElement.Loci.remap(loci as StructureElement.Loci, filter);
                loci = StructureElement.Loci.remap(filtered, structure);
                if (!StructureElement.Loci.isEmpty(loci)) {
                    layers.push({ loci, value });
                }
            }
            return { kind: 'element-loci', layers };
        } else {
            return wiggle;
        }
    }

    export type ScriptLayer = { script: Script, value: number }
    export function ofScript(scriptLayers: ScriptLayer[], structure: Structure): Wiggle {
        const layers: Wiggle.Layer[] = [];
        for (let i = 0, il = scriptLayers.length; i < il; ++i) {
            const { script, value } = scriptLayers[i];
            const loci = Script.toLoci(script, structure);
            if (!StructureElement.Loci.isEmpty(loci)) {
                layers.push({ loci, value });
            }
        }
        return { kind: 'element-loci', layers };
    }

    export type BundleLayer = { bundle: StructureElement.Bundle, value: number }
    export function ofBundle(bundleLayers: BundleLayer[], structure: Structure): Wiggle {
        const layers: Wiggle.Layer[] = [];
        for (let i = 0, il = bundleLayers.length; i < il; ++i) {
            const { bundle, value } = bundleLayers[i];
            const loci = StructureElement.Bundle.toLoci(bundle, structure.root);
            layers.push({ loci, value });
        }
        return { kind: 'element-loci', layers };
    }

    export function toBundle(wiggle: Wiggle<StructureElement.Loci>) {
        const layers: BundleLayer[] = [];
        for (let i = 0, il = wiggle.layers.length; i < il; ++i) {
            const { loci, value } = wiggle.layers[i];
            const bundle = StructureElement.Bundle.fromLoci(loci);
            layers.push({ bundle, value });
        }
        return { kind: 'element-loci', layers };
    }
}

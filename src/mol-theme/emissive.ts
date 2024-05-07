/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../mol-model/loci';
import { Structure, StructureElement } from '../mol-model/structure';
import { Script } from '../mol-script/script';

export { Emissive };

type Emissive<T extends Loci = Loci> = {
    readonly kind: T['kind']
    readonly layers: ReadonlyArray<Emissive.Layer<T>>
}

function Emissive<T extends Loci>(kind: T['kind'], layers: ReadonlyArray<Emissive.Layer<T>>): Emissive {
    return { kind, layers };
}

namespace Emissive {
    export type Layer<T extends Loci = Loci> = { readonly loci: T, readonly value: number }
    export const Empty: Emissive = { kind: 'empty-loci', layers: [] };

    export function areEqual(eA: Emissive, eB: Emissive) {
        if (eA.layers.length === 0 && eB.layers.length === 0) return true;
        if (eA.layers.length !== eB.layers.length) return false;
        for (let i = 0, il = eA.layers.length; i < il; ++i) {
            if (eA.layers[i].value !== eB.layers[i].value) return false;
            if (!Loci.areEqual(eA.layers[i].loci, eB.layers[i].loci)) return false;
        }
        return true;
    }

    export function isEmpty(emissive: Emissive) {
        return emissive.layers.length === 0;
    }

    export function remap(emissive: Emissive, structure: Structure): Emissive {
        if (emissive.kind === 'element-loci') {
            const layers: Emissive.Layer[] = [];
            for (const layer of emissive.layers) {
                let { loci, value } = layer;
                loci = StructureElement.Loci.remap(loci as StructureElement.Loci, structure);
                if (!StructureElement.Loci.isEmpty(loci)) {
                    layers.push({ loci, value });
                }
            }
            return { kind: 'element-loci', layers };
        } else {
            return emissive;
        }
    }

    export function merge(emissive: Emissive): Emissive {
        if (isEmpty(emissive)) return emissive;
        if (emissive.kind === 'element-loci') {
            const { structure } = emissive.layers[0].loci as StructureElement.Loci;
            const map = new Map<number, StructureElement.Loci>();
            let shadowed = StructureElement.Loci.none(structure);
            for (let i = 0, il = emissive.layers.length; i < il; ++i) {
                let { loci, value } = emissive.layers[il - i - 1]; // process from end
                loci = StructureElement.Loci.subtract(loci as StructureElement.Loci, shadowed);
                shadowed = StructureElement.Loci.union(loci, shadowed);
                if (!StructureElement.Loci.isEmpty(loci)) {
                    if (map.has(value)) {
                        loci = StructureElement.Loci.union(loci, map.get(value)!);
                    }
                    map.set(value, loci);
                }
            }
            const layers: Emissive.Layer<StructureElement.Loci>[] = [];
            map.forEach((loci, value) => {
                layers.push({ loci, value });
            });
            return { kind: 'element-loci', layers };
        } else {
            return emissive;
        }
    }

    export function filter(emissive: Emissive, filter: Structure): Emissive {
        if (isEmpty(emissive)) return emissive;
        if (emissive.kind === 'element-loci') {
            const { structure } = emissive.layers[0].loci as StructureElement.Loci;
            const layers: Emissive.Layer[] = [];
            for (const layer of emissive.layers) {
                let { loci, value } = layer;
                // filter by first map to the `filter` structure and
                // then map back to the original structure of the emissive loci
                const filtered = StructureElement.Loci.remap(loci as StructureElement.Loci, filter);
                loci = StructureElement.Loci.remap(filtered, structure);
                if (!StructureElement.Loci.isEmpty(loci)) {
                    layers.push({ loci, value });
                }
            }
            return { kind: 'element-loci', layers };
        } else {
            return emissive;
        }
    }

    export type ScriptLayer = { script: Script, value: number }
    export function ofScript(scriptLayers: ScriptLayer[], structure: Structure): Emissive {
        const layers: Emissive.Layer[] = [];
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
    export function ofBundle(bundleLayers: BundleLayer[], structure: Structure): Emissive {
        const layers: Emissive.Layer[] = [];
        for (let i = 0, il = bundleLayers.length; i < il; ++i) {
            const { bundle, value } = bundleLayers[i];
            const loci = StructureElement.Bundle.toLoci(bundle, structure.root);
            layers.push({ loci, value });
        }
        return { kind: 'element-loci', layers };
    }

    export function toBundle(emissive: Emissive<StructureElement.Loci>) {
        const layers: BundleLayer[] = [];
        for (let i = 0, il = emissive.layers.length; i < il; ++i) {
            const { loci, value } = emissive.layers[i];
            const bundle = StructureElement.Bundle.fromLoci(loci);
            layers.push({ bundle, value });
        }
        return { kind: 'element-loci', layers };
    }
}
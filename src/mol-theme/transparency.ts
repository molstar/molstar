/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../mol-model/loci';
import { StructureElement, Structure } from '../mol-model/structure';
import { Script } from '../mol-script/script';

export { Transparency };

type Transparency<T extends Loci = Loci> = {
    readonly kind: T['kind']
    readonly layers: ReadonlyArray<Transparency.Layer<T>>
}

function Transparency<T extends Loci>(kind: T['kind'], layers: ReadonlyArray<Transparency.Layer<T>>): Transparency<T> {
    return { kind, layers };
}

namespace Transparency {
    export type Layer<T extends Loci = Loci> = { readonly loci: T, readonly value: number }
    export const Empty: Transparency = { kind: 'empty-loci', layers: [] };

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

    export function remap(transparency: Transparency, structure: Structure): Transparency {
        if (transparency.kind === 'element-loci') {
            const layers: Transparency.Layer[] = [];
            for (const layer of transparency.layers) {
                const loci = StructureElement.Loci.remap(layer.loci as StructureElement.Loci, structure);
                if (!StructureElement.Loci.isEmpty(loci)) {
                    layers.push({ loci, value: layer.value });
                }
            }
            return { kind: 'element-loci', layers };
        } else {
            return transparency;
        }
    }

    export function merge(transparency: Transparency): Transparency {
        if (isEmpty(transparency)) return transparency;
        if (transparency.kind === 'element-loci') {
            const { structure } = transparency.layers[0].loci as StructureElement.Loci;
            const map = new Map<number, StructureElement.Loci>();
            let shadowed = StructureElement.Loci.none(structure);
            for (let i = 0, il = transparency.layers.length; i < il; ++i) {
                let { loci, value } = transparency.layers[il - i - 1]; // process from end
                loci = StructureElement.Loci.subtract(loci as StructureElement.Loci, shadowed);
                shadowed = StructureElement.Loci.union(loci, shadowed);
                if (!StructureElement.Loci.isEmpty(loci)) {
                    if (map.has(value)) {
                        loci = StructureElement.Loci.union(loci, map.get(value)!);
                    }
                    map.set(value, loci);
                }
            }
            const layers: Transparency.Layer<StructureElement.Loci>[] = [];
            map.forEach((loci, value) => {
                layers.push({ loci, value });
            });
            return { kind: 'element-loci', layers };
        } else {
            return transparency;
        }
    }

    export function filter(transparency: Transparency, filter: Structure): Transparency {
        if (isEmpty(transparency)) return transparency;
        if (transparency.kind === 'element-loci') {
            const { structure } = transparency.layers[0].loci as StructureElement.Loci;
            const layers: Transparency.Layer<StructureElement.Loci>[] = [];
            for (const layer of transparency.layers) {
                let { loci, value } = layer;
                // filter by first map to the `filter` structure and
                // then map back to the original structure of the transparency loci
                const filtered = StructureElement.Loci.remap(loci as StructureElement.Loci, filter);
                loci = StructureElement.Loci.remap(filtered, structure);
                if (!StructureElement.Loci.isEmpty(loci)) {
                    layers.push({ loci, value });
                }
            }
            return { kind: 'element-loci', layers };
        } else {
            return transparency;
        }
    }

    export type ScriptLayer = { script: Script, value: number }
    export function ofScript(scriptLayers: ScriptLayer[], structure: Structure): Transparency<StructureElement.Loci> {
        const layers: Transparency.Layer<StructureElement.Loci>[] = [];
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
    export function ofBundle(bundleLayers: BundleLayer[], structure: Structure): Transparency<StructureElement.Loci> {
        const layers: Transparency.Layer<StructureElement.Loci>[] = [];
        for (let i = 0, il = bundleLayers.length; i < il; ++i) {
            const { bundle, value } = bundleLayers[i];
            const loci = StructureElement.Bundle.toLoci(bundle, structure.root);
            layers.push({ loci, value });
        }
        return { kind: 'element-loci', layers };
    }

    export function toBundle(transparency: Transparency<StructureElement.Loci>) {
        const layers: BundleLayer[] = [];
        for (let i = 0, il = transparency.layers.length; i < il; ++i) {
            const { loci, value } = transparency.layers[i];
            const bundle = StructureElement.Bundle.fromLoci(loci);
            layers.push({ bundle, value });
        }
        return { kind: 'element-loci', layers };
    }
}
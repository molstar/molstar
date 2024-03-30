/**
 * Copyright (c) 2021-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../mol-model/loci';
import { Structure, StructureElement } from '../mol-model/structure';
import { Script } from '../mol-script/script';
import { Material } from '../mol-util/material';
import { shallowEqual } from '../mol-util/object';

export { Substance };

type Substance<T extends Loci = Loci> = {
    readonly kind: T['kind']
    readonly layers: ReadonlyArray<Substance.Layer<T>>
}

function Substance<T extends Loci>(kind: T['kind'], layers: ReadonlyArray<Substance.Layer<T>>): Substance {
    return { kind, layers };
}

namespace Substance {
    export type Layer<T extends Loci = Loci> = { readonly loci: T, readonly material: Material, readonly clear: boolean }
    export const Empty: Substance = { kind: 'empty-loci', layers: [] };

    export function areEqual(sA: Substance, sB: Substance) {
        if (sA.layers.length === 0 && sB.layers.length === 0) return true;
        if (sA.layers.length !== sB.layers.length) return false;
        for (let i = 0, il = sA.layers.length; i < il; ++i) {
            if (sA.layers[i].clear !== sB.layers[i].clear) return false;
            if (!shallowEqual(sA.layers[i].material, sB.layers[i].material)) return false;
            if (!Loci.areEqual(sA.layers[i].loci, sB.layers[i].loci)) return false;
        }
        return true;
    }

    export function isEmpty(substance: Substance) {
        return substance.layers.length === 0;
    }

    export function remap(substance: Substance, structure: Structure): Substance {
        if (substance.kind === 'element-loci') {
            const layers: Substance.Layer[] = [];
            for (const layer of substance.layers) {
                let { loci, material, clear } = layer;
                loci = StructureElement.Loci.remap(loci as StructureElement.Loci, structure);
                if (!StructureElement.Loci.isEmpty(loci)) {
                    layers.push({ loci, material, clear });
                }
            }
            return { kind: 'element-loci', layers };
        } else {
            return substance;
        }
    }

    export function merge(substance: Substance): Substance {
        if (isEmpty(substance)) return substance;
        if (substance.kind === 'element-loci') {
            const { structure } = substance.layers[0].loci as StructureElement.Loci;
            let clearLoci: StructureElement.Loci | undefined = void 0;
            const map = new Map<Material, StructureElement.Loci>();
            let shadowed = StructureElement.Loci.none(structure);
            for (let i = 0, il = substance.layers.length; i < il; ++i) {
                let { loci, material, clear } = substance.layers[il - i - 1]; // process from end
                loci = StructureElement.Loci.subtract(loci as StructureElement.Loci, shadowed);
                shadowed = StructureElement.Loci.union(loci, shadowed);
                if (!StructureElement.Loci.isEmpty(loci)) {
                    if (clear) {
                        clearLoci = clearLoci
                            ? StructureElement.Loci.union(loci, clearLoci)
                            : loci;
                    } else {
                        if (map.has(material)) {
                            loci = StructureElement.Loci.union(loci, map.get(material)!);
                        }
                        map.set(material, loci);
                    }
                }
            }
            const layers: Substance.Layer[] = [];
            if (clearLoci) {
                layers.push({ loci: clearLoci, material: Material(), clear: true });
            }
            map.forEach((loci, material) => {
                layers.push({ loci, material, clear: false });
            });
            return { kind: 'element-loci', layers };
        } else {
            return substance;
        }
    }

    export function filter(substance: Substance, filter: Structure): Substance {
        if (isEmpty(substance)) return substance;
        if (substance.kind === 'element-loci') {
            const { structure } = substance.layers[0].loci as StructureElement.Loci;
            const layers: Substance.Layer[] = [];
            for (const layer of substance.layers) {
                let { loci, material, clear } = layer;
                // filter by first map to the `filter` structure and
                // then map back to the original structure of the substance loci
                const filtered = StructureElement.Loci.remap(loci as StructureElement.Loci, filter);
                loci = StructureElement.Loci.remap(filtered, structure);
                if (!StructureElement.Loci.isEmpty(loci)) {
                    layers.push({ loci, material, clear });
                }
            }
            return { kind: 'element-loci', layers };
        } else {
            return substance;
        }
    }

    export type ScriptLayer = { script: Script, material: Material, clear: boolean }
    export function ofScript(scriptLayers: ScriptLayer[], structure: Structure): Substance {
        const layers: Substance.Layer[] = [];
        for (let i = 0, il = scriptLayers.length; i < il; ++i) {
            const { script, material, clear } = scriptLayers[i];
            const loci = Script.toLoci(script, structure);
            if (!StructureElement.Loci.isEmpty(loci)) {
                layers.push({ loci, material, clear });
            }
        }
        return { kind: 'element-loci', layers };
    }

    export type BundleLayer = { bundle: StructureElement.Bundle, material: Material, clear: boolean }
    export function ofBundle(bundleLayers: BundleLayer[], structure: Structure): Substance {
        const layers: Substance.Layer[] = [];
        for (let i = 0, il = bundleLayers.length; i < il; ++i) {
            const { bundle, material, clear } = bundleLayers[i];
            const loci = StructureElement.Bundle.toLoci(bundle, structure.root);
            layers.push({ loci, material, clear });
        }
        return { kind: 'element-loci', layers };
    }

    export function toBundle(overpaint: Substance<StructureElement.Loci>) {
        const layers: BundleLayer[] = [];
        for (let i = 0, il = overpaint.layers.length; i < il; ++i) {
            const { loci, material, clear } = overpaint.layers[i];
            const bundle = StructureElement.Bundle.fromLoci(loci);
            layers.push({ bundle, material, clear });
        }
        return { kind: 'element-loci', layers };
    }
}
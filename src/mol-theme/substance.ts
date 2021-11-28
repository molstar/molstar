/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../mol-model/loci';
import { Structure, StructureElement } from '../mol-model/structure';
import { Script } from '../mol-script/script';
import { Material } from '../mol-util/material';

export { Substance };

type Substance = { readonly layers: ReadonlyArray<Substance.Layer> }

function Substance(layers: ReadonlyArray<Substance.Layer>): Substance {
    return { layers };
}

namespace Substance {
    export type Layer = { readonly loci: StructureElement.Loci, readonly material: Material, readonly clear: boolean }
    export const Empty: Substance = { layers: [] };

    export function areEqual(sA: Substance, sB: Substance) {
        if (sA.layers.length === 0 && sB.layers.length === 0) return true;
        if (sA.layers.length !== sB.layers.length) return false;
        for (let i = 0, il = sA.layers.length; i < il; ++i) {
            if (sA.layers[i].clear !== sB.layers[i].clear) return false;
            if (sA.layers[i].material !== sB.layers[i].material) return false;
            if (!Loci.areEqual(sA.layers[i].loci, sB.layers[i].loci)) return false;
        }
        return true;
    }

    export function isEmpty(overpaint: Substance) {
        return overpaint.layers.length === 0;
    }

    export function remap(substance: Substance, structure: Structure) {
        const layers: Substance.Layer[] = [];
        for (const layer of substance.layers) {
            let { loci, material, clear } = layer;
            loci = StructureElement.Loci.remap(loci, structure);
            if (!StructureElement.Loci.isEmpty(loci)) {
                layers.push({ loci, material, clear });
            }
        }
        return { layers };
    }

    export function merge(substance: Substance): Substance {
        if (isEmpty(substance)) return substance;
        const { structure } = substance.layers[0].loci;
        const map = new Map<Material | -1, StructureElement.Loci>();
        let shadowed = StructureElement.Loci.none(structure);
        for (let i = 0, il = substance.layers.length; i < il; ++i) {
            let { loci, material, clear } = substance.layers[il - i - 1]; // process from end
            loci = StructureElement.Loci.subtract(loci, shadowed);
            shadowed = StructureElement.Loci.union(loci, shadowed);
            if (!StructureElement.Loci.isEmpty(loci)) {
                const materialOrClear = clear ? -1 : material;
                if (map.has(materialOrClear)) {
                    loci = StructureElement.Loci.union(loci, map.get(materialOrClear)!);
                }
                map.set(materialOrClear, loci);
            }
        }
        const layers: Substance.Layer[] = [];
        map.forEach((loci, materialOrClear) => {
            const clear = materialOrClear === -1;
            const material = clear ? Material(0) : materialOrClear;
            layers.push({ loci, material, clear });
        });
        return { layers };
    }

    export function filter(substance: Substance, filter: Structure): Substance {
        if (isEmpty(substance)) return substance;
        const { structure } = substance.layers[0].loci;
        const layers: Substance.Layer[] = [];
        for (const layer of substance.layers) {
            let { loci, material, clear } = layer;
            // filter by first map to the `filter` structure and
            // then map back to the original structure of the substance loci
            const filtered = StructureElement.Loci.remap(loci, filter);
            loci = StructureElement.Loci.remap(filtered, structure);
            if (!StructureElement.Loci.isEmpty(loci)) {
                layers.push({ loci, material, clear });
            }
        }
        return { layers };
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
        return { layers };
    }

    export type BundleLayer = { bundle: StructureElement.Bundle, material: Material, clear: boolean }
    export function ofBundle(bundleLayers: BundleLayer[], structure: Structure): Substance {
        const layers: Substance.Layer[] = [];
        for (let i = 0, il = bundleLayers.length; i < il; ++i) {
            const { bundle, material, clear } = bundleLayers[i];
            const loci = StructureElement.Bundle.toLoci(bundle, structure.root);
            layers.push({ loci, material, clear });
        }
        return { layers };
    }

    export function toBundle(overpaint: Substance) {
        const layers: BundleLayer[] = [];
        for (let i = 0, il = overpaint.layers.length; i < il; ++i) {
            const { loci, material, clear } = overpaint.layers[i];
            const bundle = StructureElement.Bundle.fromLoci(loci);
            layers.push({ bundle, material, clear });
        }
        return { layers };
    }
}
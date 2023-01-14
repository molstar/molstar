/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../mol-model/loci';
import { Color } from '../mol-util/color';
import { Structure, StructureElement } from '../mol-model/structure';
import { Script } from '../mol-script/script';

export { Overpaint };

type Overpaint<T extends Loci = Loci> = {
    readonly kind: T['kind']
    readonly layers: ReadonlyArray<Overpaint.Layer<T>>
}

function Overpaint<T extends Loci>(kind: T['kind'], layers: ReadonlyArray<Overpaint.Layer<T>>): Overpaint<T> {
    return { kind, layers };
}

namespace Overpaint {
    export type Layer<T extends Loci = Loci> = { readonly loci: T, readonly color: Color, readonly clear: boolean }
    export const Empty: Overpaint = { kind: 'empty-loci', layers: [] };

    export function areEqual(oA: Overpaint, oB: Overpaint) {
        if (oA.layers.length === 0 && oB.layers.length === 0) return true;
        if (oA.layers.length !== oB.layers.length) return false;
        for (let i = 0, il = oA.layers.length; i < il; ++i) {
            if (oA.layers[i].clear !== oB.layers[i].clear) return false;
            if (oA.layers[i].color !== oB.layers[i].color) return false;
            if (!Loci.areEqual(oA.layers[i].loci, oB.layers[i].loci)) return false;
        }
        return true;
    }

    export function isEmpty(overpaint: Overpaint) {
        return overpaint.layers.length === 0;
    }

    export function remap(overpaint: Overpaint, structure: Structure): Overpaint {
        if (overpaint.kind === 'element-loci') {
            const layers: Overpaint.Layer[] = [];
            for (const layer of overpaint.layers) {
                let { loci, color, clear } = layer;
                loci = StructureElement.Loci.remap(loci as StructureElement.Loci, structure);
                if (!StructureElement.Loci.isEmpty(loci)) {
                    layers.push({ loci, color, clear });
                }
            }
            return { kind: 'element-loci', layers };
        } else {
            return overpaint;
        }
    }

    export function merge(overpaint: Overpaint): Overpaint {
        if (isEmpty(overpaint)) return overpaint;
        if (overpaint.kind === 'element-loci') {
            const { structure } = overpaint.layers[0].loci as StructureElement.Loci;
            const map = new Map<Color | -1, StructureElement.Loci>();
            let shadowed = StructureElement.Loci.none(structure);
            for (let i = 0, il = overpaint.layers.length; i < il; ++i) {
                let { loci, color, clear } = overpaint.layers[il - i - 1]; // process from end
                loci = StructureElement.Loci.subtract(loci as StructureElement.Loci, shadowed);
                shadowed = StructureElement.Loci.union(loci, shadowed);
                if (!StructureElement.Loci.isEmpty(loci)) {
                    const colorOrClear = clear ? -1 : color;
                    if (map.has(colorOrClear)) {
                        loci = StructureElement.Loci.union(loci, map.get(colorOrClear)!);
                    }
                    map.set(colorOrClear, loci);
                }
            }
            const layers: Overpaint.Layer[] = [];
            map.forEach((loci, colorOrClear) => {
                const clear = colorOrClear === -1;
                const color = clear ? Color(0) : colorOrClear;
                layers.push({ loci, color, clear });
            });
            return { kind: 'element-loci', layers };
        } else {
            return overpaint;
        }
    }

    export function filter(overpaint: Overpaint, filter: Structure): Overpaint {
        if (isEmpty(overpaint)) return overpaint;
        if (overpaint.kind === 'element-loci') {
            const { structure } = overpaint.layers[0].loci as StructureElement.Loci;
            const layers: Overpaint.Layer[] = [];
            for (const layer of overpaint.layers) {
                let { loci, color, clear } = layer;
                // filter by first map to the `filter` structure and
                // then map back to the original structure of the overpaint loci
                const filtered = StructureElement.Loci.remap(loci as StructureElement.Loci, filter);
                loci = StructureElement.Loci.remap(filtered, structure);
                if (!StructureElement.Loci.isEmpty(loci)) {
                    layers.push({ loci, color, clear });
                }
            }
            return { kind: 'element-loci', layers };
        } else {
            return overpaint;
        }
    }

    export type ScriptLayer = { script: Script, color: Color, clear: boolean }
    export function ofScript(scriptLayers: ScriptLayer[], structure: Structure): Overpaint {
        const layers: Overpaint.Layer[] = [];
        for (let i = 0, il = scriptLayers.length; i < il; ++i) {
            const { script, color, clear } = scriptLayers[i];
            const loci = Script.toLoci(script, structure);
            if (!StructureElement.Loci.isEmpty(loci)) {
                layers.push({ loci, color, clear });
            }
        }
        return { kind: 'element-loci', layers };
    }

    export type BundleLayer = { bundle: StructureElement.Bundle, color: Color, clear: boolean }
    export function ofBundle(bundleLayers: BundleLayer[], structure: Structure): Overpaint {
        const layers: Overpaint.Layer[] = [];
        for (let i = 0, il = bundleLayers.length; i < il; ++i) {
            const { bundle, color, clear } = bundleLayers[i];
            const loci = StructureElement.Bundle.toLoci(bundle, structure.root);
            layers.push({ loci, color, clear });
        }
        return { kind: 'element-loci', layers };
    }

    export function toBundle(overpaint: Overpaint<StructureElement.Loci>) {
        const layers: BundleLayer[] = [];
        for (let i = 0, il = overpaint.layers.length; i < il; ++i) {
            const { loci, color, clear } = overpaint.layers[i];
            const bundle = StructureElement.Bundle.fromLoci(loci);
            layers.push({ bundle, color, clear });
        }
        return { kind: 'element-loci', layers };
    }
}
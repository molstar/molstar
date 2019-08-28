/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../mol-model/loci';
import { Color } from '../mol-util/color';
import { Structure, StructureElement } from '../mol-model/structure';
import { Script } from '../mol-script/script';

export { Overpaint }

type Overpaint = { layers: ReadonlyArray<Overpaint.Layer>, readonly alpha: number }

namespace Overpaint {
    export type Layer = { readonly loci: Loci, readonly color: Color, readonly clear: boolean }
    export const Empty: Overpaint = { layers: [], alpha: 1 }

    export function areEqual(oA: Overpaint, oB: Overpaint) {
        if (oA.layers.length === 0 && oB.layers.length === 0) return true
        if (oA.layers.length !== oB.layers.length) return false
        if (oA.alpha !== oB.alpha) return false
        for (let i = 0, il = oA.layers.length; i < il; ++i) {
            if (oA.layers[i].clear !== oB.layers[i].clear) return false
            if (oA.layers[i].color !== oB.layers[i].color) return false
            if (!Loci.areEqual(oA.layers[i].loci, oB.layers[i].loci)) return false
        }
        return true
    }

    export function remap(overpaint: Overpaint, structure: Structure) {
        const layers: Overpaint.Layer[] = []
        for (const layer of overpaint.layers) {
            const { loci, color, clear } = layer
            layers.push({ loci: Loci.remap(loci, structure), color, clear })
        }
        return { layers, alpha: overpaint.alpha }
    }

    export function ofScript(scriptLayers: { script: Script, color: Color, clear: boolean }[], alpha: number, structure: Structure): Overpaint {
        const layers: Overpaint.Layer[] = []
        for (let i = 0, il = scriptLayers.length; i < il; ++i) {
            const { script, color, clear } = scriptLayers[i]
            layers.push({ loci: Script.toLoci(script, structure), color, clear })
        }
        return { layers, alpha }
    }

    export function ofBundle(bundleLayers: { bundle: StructureElement.Bundle, color: Color, clear: boolean }[], alpha: number, structure: Structure): Overpaint {
        const layers: Overpaint.Layer[] = []
        for (let i = 0, il = bundleLayers.length; i < il; ++i) {
            const { bundle, color, clear } = bundleLayers[i]
            const loci = StructureElement.Bundle.toLoci(bundle, structure.root)
            layers.push({ loci, color, clear })
        }
        return { layers, alpha }
    }
}
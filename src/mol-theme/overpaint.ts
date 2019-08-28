/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../mol-model/loci';
import { Color } from '../mol-util/color';
import { Structure } from '../mol-model/structure';

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
}
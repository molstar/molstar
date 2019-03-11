/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from 'mol-model/loci';
import { Color } from 'mol-util/color';

export { Overpaint }

type Overpaint = { layers: ReadonlyArray<Overpaint.Layer>, readonly alpha: number }

namespace Overpaint {
    export type Layer = { readonly loci: Loci, readonly color: Color }
    export const Empty: Overpaint = { layers: [], alpha: 1 }

    export function areEqual(oA: Overpaint, oB: Overpaint) {
        if (oA.layers.length === 0 && oB.layers.length === 0) return true
        if (oA.layers.length !== oB.layers.length) return false
        if (oA.alpha !== oB.alpha) return false
        for (let i = 0, il = oA.layers.length; i < il; ++i) {
            if (oA.layers[i].color !== oB.layers[i].color) return false
            if (!Loci.areEqual(oA.layers[i].loci, oB.layers[i].loci)) return false
        }
        return true
    }
}
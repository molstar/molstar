/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from 'mol-model/loci';

export { Transparency }

type Transparency = { layers: ReadonlyArray<Transparency.Layer> }

namespace Transparency {
    export type Layer = { readonly loci: Loci, readonly value: number }
    export const Empty: Transparency = { layers: [] }

    export function areEqual(tA: Transparency, tB: Transparency) {
        if (tA.layers.length === 0 && tB.layers.length === 0) return true
        if (tA.layers.length !== tB.layers.length) return false
        for (let i = 0, il = tA.layers.length; i < il; ++i) {
            if (tA.layers[i].value !== tB.layers[i].value) return false
            if (!Loci.areEqual(tA.layers[i].loci, tB.layers[i].loci)) return false
        }
        return true
    }
}
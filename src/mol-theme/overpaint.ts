/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from 'mol-model/loci';
import { Color } from 'mol-util/color';

export { Overpaint }
namespace Overpaint {
    export type Layer = { readonly loci: Loci, readonly color: Color }
    export type Layers = { list: ReadonlyArray<Layer>, readonly alpha: number }
    export const EmptyLayers: Layers = { list: [], alpha: 1 }

    export function areEqual(layersA: Layers, layersB: Layers) {
        if (layersA.list.length === 0 && layersB.list.length === 0) return true
        if (layersA.list.length !== layersB.list.length) return false
        if (layersA.alpha !== layersB.alpha) return false
        for (let i = 0, il = layersA.list.length; i < il; ++i) {
            if (layersA.list[i].color !== layersB.list[i].color) return false
            if (!Loci.areEqual(layersA.list[i].loci, layersB.list[i].loci)) return false
        }
        return true
    }
}
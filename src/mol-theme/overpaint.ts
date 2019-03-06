/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from 'mol-model/loci';
import { Color } from 'mol-util/color';

export { Overpaint }
namespace Overpaint {
    export interface Layer {
        loci: Loci
        color: Color
        alpha: number
    }
    export type Layers = Layer[]

    export function areEqual(layersA: Layers, layersB: Layers) {
        if (layersA.length === 0 && layersB.length === 0) return true
        if (layersA.length !== layersB.length) return false
        for (let i = 0, il = layersA.length; i < il; ++i) {
            if (layersA[i].alpha !== layersB[i].alpha) return false
            if (layersA[i].color !== layersB[i].color) return false
            if (!Loci.areEqual(layersA[i].loci, layersB[i].loci)) return false
        }
        return true
    }
}
/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Interval } from '../../../mol-data/int';
import { Loci } from '../../../mol-model/loci';
import { MarkerAction, applyMarkerAction } from '../../../mol-util/marker-action';
import { StructureElement } from '../../../mol-model/structure';
import { Color } from '../../../mol-util/color';

export { SequenceWrapper }

abstract class SequenceWrapper<D> {
    abstract seqId(i: number): number
    abstract residueLabel(i: number): string
    abstract residueColor(i: number): Color

    abstract eachResidue(loci: Loci, apply: (interval: Interval) => boolean): boolean
    abstract getLoci(seqId: number): StructureElement.Loci

    markResidue(loci: Loci, action: MarkerAction) {
        return this.eachResidue(loci, (i: Interval) => {
            return applyMarkerAction(this.markerArray, i, action)
        })
    }

    constructor(readonly data: D, readonly markerArray: Uint8Array, readonly length: number) {

    }
}

namespace SequenceWrapper {
    export type Any = SequenceWrapper<any>
}
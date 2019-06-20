/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Interval } from '../../../mol-data/int';
import { Loci } from '../../../mol-model/loci';
import { MarkerAction, applyMarkerAction } from '../../../mol-util/marker-action';
import { StructureElement } from '../../../mol-model/structure';
import { Sequence } from '../../../mol-model/sequence';

export { SequenceWrapper }

abstract class SequenceWrapper<D> {
    label: string
    data: D
    markerArray: Uint8Array
    sequence: Sequence
    abstract eachResidue(loci: Loci, apply: (interval: Interval) => boolean): boolean
    abstract getLoci(seqId: number): StructureElement.Loci

    markResidue(loci: Loci, action: MarkerAction) {
        return this.eachResidue(loci, (i: Interval) => {
            return applyMarkerAction(this.markerArray, i, action)
        })
    }
}

namespace SequenceWrapper {
    export type Any = SequenceWrapper<any>
}
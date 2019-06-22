/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet } from '../../../mol-data/int';
import { Loci } from '../../../mol-model/loci';
import { MarkerAction, applyMarkerAction } from '../../../mol-util/marker-action';
import { StructureElement, Structure, Unit } from '../../../mol-model/structure';
import { Color } from '../../../mol-util/color';

export type StructureUnit = { structure: Structure, unit: Unit }

export { SequenceWrapper }

abstract class SequenceWrapper<D> {
    abstract residueLabel(seqIdx: number): string
    abstract residueColor(seqIdx: number): Color

    abstract eachResidue(loci: Loci, apply: (set: OrderedSet) => boolean): boolean
    abstract getLoci(seqIdx: number): StructureElement.Loci

    markResidue(loci: Loci, action: MarkerAction) {
        return this.eachResidue(loci, (set: OrderedSet) => {
            return applyMarkerAction(this.markerArray, set, action)
        })
    }

    constructor(readonly data: D, readonly markerArray: Uint8Array, readonly length: number) {

    }
}

namespace SequenceWrapper {
    export type Any = SequenceWrapper<any>
}
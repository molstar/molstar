/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Interval, OrderedSet } from '../../mol-data/int';
import { Loci, isEveryLoci } from '../../mol-model/loci';
import { MarkerAction, applyMarkerAction } from '../../mol-util/marker-action';
import { StructureElement, Structure, Unit } from '../../mol-model/structure';
import { Color } from '../../mol-util/color';

export type StructureUnit = { structure: Structure, units: Unit[] }

export { SequenceWrapper };

abstract class SequenceWrapper<D> {
    abstract residueLabel(seqIdx: number): string
    abstract residueColor(seqIdx: number): Color
    abstract residueClass(seqIdx: number): string

    abstract getLoci(seqIdx: number): StructureElement.Loci

    /** Return list of sequence viewer positions that correspond to `loci` */
    abstract getSeqIndices(loci: Loci): OrderedSet;

    mark(loci: Loci, action: MarkerAction): boolean {
        const seqIdxs = this.getSeqIndices(loci);
        if (OrderedSet.size(seqIdxs) === 0) return false;
        return applyMarkerAction(this.markerArray, seqIdxs, action);
    }

    markResidue(loci: Loci, action: MarkerAction | 'focus' | 'unfocus'): boolean {
        if (action === 'focus') return this.markResidueFocus(loci, true);
        if (action === 'unfocus') return this.markResidueFocus(loci, false);

        if (isEveryLoci(loci)) {
            return applyMarkerAction(this.markerArray, Interval.ofLength(this.length), action);
        } else {
            return this.mark(loci, action);
        }
    }

    private markResidueFocus(loci: Loci, focusState: boolean) {
        const value = focusState ? 1 : 0;
        if (isEveryLoci(loci)) {
            this.focusMarkerArray.fill(value, 0, this.length);
            return true;
        } else {
            const seqIdxs = this.getSeqIndices(loci);
            OrderedSet.forEach(seqIdxs, seqIdx => this.focusMarkerArray[seqIdx] = value);
            return OrderedSet.size(seqIdxs) > 0;
        }
    }

    /** Markers for "highlighted" and "selected" (2-bits per position) */
    readonly markerArray: Uint8Array;
    /** Markers for "focused" (1-bit per position) */
    readonly focusMarkerArray: Uint8Array;

    constructor(readonly data: D, readonly length: number) {
        this.markerArray = new Uint8Array(length);
        this.focusMarkerArray = new Uint8Array(length);
    }
}

namespace SequenceWrapper {
    export type Any = SequenceWrapper<any>
}

/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Volume } from '../../mol-model/volume';
import { Loci } from '../../mol-model/loci';
import { Interval, OrderedSet } from '../../mol-data/int';
import { equalEps } from '../../mol-math/linear-algebra/3d/common';

export function eachVolumeLoci(loci: Loci, volume: Volume, isoValue: Volume.IsoValue | undefined, apply: (interval: Interval) => boolean) {
    let changed = false;
    if (Volume.isLoci(loci)) {
        if (!Volume.areEquivalent(loci.volume, volume)) return false;
        if (apply(Interval.ofLength(volume.grid.cells.data.length))) changed = true;
    } else if (Volume.Isosurface.isLoci(loci)) {
        if (!Volume.areEquivalent(loci.volume, volume)) return false;
        if (isoValue) {
            if (!Volume.IsoValue.areSame(loci.isoValue, isoValue, volume.grid.stats)) return false;
            if (apply(Interval.ofLength(volume.grid.cells.data.length))) changed = true;
        } else {
            // TODO find a cheaper way?
            const { stats, cells: { data } } = volume.grid;
            const eps = stats.sigma;
            const v = Volume.IsoValue.toAbsolute(loci.isoValue, stats).absoluteValue;
            for (let i = 0, il = data.length; i < il; ++i) {
                if (equalEps(v, data[i], eps)) {
                    if (apply(Interval.ofSingleton(i))) changed = true;
                }
            }
        }
    } else if (Volume.Cell.isLoci(loci)) {
        if (!Volume.areEquivalent(loci.volume, volume)) return false;
        if (Interval.is(loci.indices)) {
            if (apply(loci.indices)) changed = true;
        } else {
            OrderedSet.forEach(loci.indices, v => {
                if (apply(Interval.ofSingleton(v))) changed = true;
            });
        }
    }
    return changed;
}
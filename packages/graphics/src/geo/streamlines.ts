import { Interval } from '@molstar/core/data/int/interval';
import { OrderedSet } from '@molstar/core/data/int/ordered-set';
import type { PickingId } from './geometry/picking.js';
import { LocationIterator } from './util/location-iterator.js';
import { EmptyLoci, type Loci } from '@molstar/model/model/loci';
import { Volume } from '@molstar/model/model/volume/volume';
import { StreamlinesProvider } from '@molstar/model/props/volume/streamlines';
import { type CommonStreamlinesProps, type StreamlinesIndex, StreamlinesLoci, isStreamlinesLoci, StreamlinesLocation } from '@molstar/model/props/volume/streamlines/shared';

export function getStreamlinesVisualLoci(volume: Volume, _props: CommonStreamlinesProps) {
    const streamlines = StreamlinesProvider.get(volume).value!;
    const indices = Interval.ofLength(streamlines.length as Volume.InstanceIndex);
    const instances = Interval.ofLength(volume.instances.length as Volume.InstanceIndex);
    return StreamlinesLoci(streamlines, volume, [{ indices, instances }]);
}

export function getStreamlinesLoci(pickingId: PickingId, volume: Volume, _key: number, _props: CommonStreamlinesProps, id: number) {
    const { objectId, groupId, instanceId } = pickingId;
    if (id === objectId) {
        const granularity = Volume.PickingGranularity.get(volume);
        const instances = OrderedSet.ofSingleton(instanceId as Volume.InstanceIndex);
        if (granularity === 'volume') return Volume.Loci(volume, instances);

        const streamlines = StreamlinesProvider.get(volume).value!;
        const indices = OrderedSet.ofSingleton(groupId as StreamlinesIndex);
        return StreamlinesLoci(streamlines, volume, [{ indices, instances }]);
    }
    return EmptyLoci;
}

export function eachStreamlines(loci: Loci, volume: Volume, _key: number, _props: CommonStreamlinesProps, apply: (interval: Interval) => boolean) {
    let changed = false;
    const streamlines = StreamlinesProvider.get(volume).value!;
    const count = streamlines.length;
    if (Volume.isLoci(loci)) {
        if (!Volume.areEquivalent(loci.volume, volume)) return false;
        if (Interval.is(loci.instances)) {
            const start = Interval.start(loci.instances) * count;
            const end = Interval.end(loci.instances) * count;
            if (apply(Interval.ofBounds(start, end))) changed = true;
        } else {
            for (let i = 0, il = loci.instances.length; i < il; ++i) {
                const offset = loci.instances[i] * count;
                if (apply(Interval.ofBounds(offset, offset + count))) changed = true;
            }
        }
    } else if (isStreamlinesLoci(loci)) {
        if (!Volume.areEquivalent(loci.data.volume, volume)) return false;
        for (const { indices, instances } of loci.elements) {
            if (Interval.is(indices)) {
                OrderedSet.forEach(instances, j => {
                    const offset = j * count;
                    if (apply(Interval.offset(indices, offset))) changed = true;
                });
            } else {
                OrderedSet.forEach(indices, v => {
                    OrderedSet.forEach(instances, j => {
                        const offset = j * count;
                        if (apply(Interval.ofSingleton(offset + v))) changed = true;
                    });
                });
            }
        }
    }
    return changed;
}

export function createStreamlinesLocationIterator(volume: Volume): LocationIterator {
    const streamlines = StreamlinesProvider.get(volume).value!;
    const groupCount = streamlines.length;
    const instanceCount = volume.instances.length;

    const l = StreamlinesLocation(streamlines, volume);
    const getLocation = (groupIndex: number, instanceIndex: number) => {
        l.element.index = groupIndex as StreamlinesIndex;
        l.element.instance = instanceIndex as Volume.InstanceIndex;
        return l;
    };
    return LocationIterator(groupCount, instanceCount, 1, getLocation);
}

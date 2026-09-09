/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Volume } from '../../../mol-model/volume';
import { BodyId, BodyInfo, LabelStore } from './types';

const PropertyKey = '__body-labels__';

/** Accessor for the `LabelStore` attached to a volume (pattern of `Volume.Segmentation`). */
export const BodyLabels = {
    get(volume: Volume): LabelStore | undefined {
        return volume._propertyData[PropertyKey];
    },
    set(volume: Volume, store: LabelStore) {
        volume._propertyData[PropertyKey] = store;
    },
    /** Returns the existing store or attaches a fresh, all-unassigned one. */
    ensure(volume: Volume): LabelStore {
        let store = BodyLabels.get(volume);
        if (!store) {
            store = {
                labels: new Uint8Array(volume.grid.cells.data.length),
                version: 0,
                bodies: [],
                nextId: 1,
            };
            BodyLabels.set(volume, store);
        }
        return store;
    },
    bump(store: LabelStore) {
        store.version++;
    },
    clear(volume: Volume) {
        delete volume._propertyData[PropertyKey];
    },
    getBody(store: LabelStore, id: BodyId): BodyInfo | undefined {
        return store.bodies.find(b => b.id === id);
    },
};

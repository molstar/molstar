/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Volume } from '../../../mol-model/volume';

/** Which voxels the current view polygons select: `1` selected, `0` not. */
export interface SelectionStore {
    selected: Uint8Array;
    /** Bumped on every change; mirrored into the color theme to force a recolor. */
    version: number;
}

const PropertyKey = '__mask-selection__';

/** Accessor for the `SelectionStore` attached to a volume (pattern of `BodyLabels`). */
export const MaskSelection = {
    get(volume: Volume): SelectionStore | undefined {
        return volume._propertyData[PropertyKey];
    },
    /** Returns the existing store or attaches a fresh, empty one. */
    ensure(volume: Volume): SelectionStore {
        let store = MaskSelection.get(volume);
        if (!store) {
            store = { selected: new Uint8Array(volume.grid.cells.data.length), version: 0 };
            volume._propertyData[PropertyKey] = store;
        }
        return store;
    },
    clear(volume: Volume) {
        delete volume._propertyData[PropertyKey];
    },
};

/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../mol-util/value-cell';
import { Vec2 } from '../../mol-math/linear-algebra';
import { TextureImage, createTextureImage } from '../../mol-gl/renderable/util';

export type MarkerData = {
    tMarker: ValueCell<TextureImage<Uint8Array>>
    uMarkerTexDim: ValueCell<Vec2>
}

export function createMarkers(count: number, markerData?: MarkerData): MarkerData {
    const markers = createTextureImage(Math.max(1, count), 1, Uint8Array, markerData && markerData.tMarker.ref.value.array);
    if (markerData) {
        ValueCell.update(markerData.tMarker, markers);
        ValueCell.update(markerData.uMarkerTexDim, Vec2.create(markers.width, markers.height));
        return markerData;
    } else {
        return {
            tMarker: ValueCell.create(markers),
            uMarkerTexDim: ValueCell.create(Vec2.create(markers.width, markers.height)),
        };
    }
}

const emptyMarkerTexture = { array: new Uint8Array(1), width: 1, height: 1 };
export function createEmptyMarkers(markerData?: MarkerData): MarkerData {
    if (markerData) {
        ValueCell.update(markerData.tMarker, emptyMarkerTexture);
        ValueCell.update(markerData.uMarkerTexDim, Vec2.create(1, 1));
        return markerData;
    } else {
        return {
            tMarker: ValueCell.create(emptyMarkerTexture),
            uMarkerTexDim: ValueCell.create(Vec2.create(1, 1)),
        };
    }
}
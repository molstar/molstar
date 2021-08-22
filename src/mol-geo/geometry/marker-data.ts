/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../mol-util/value-cell';
import { Vec2 } from '../../mol-math/linear-algebra';
import { TextureImage, createTextureImage } from '../../mol-gl/renderable/util';

export type MarkerData = {
    tMarker: ValueCell<TextureImage<Uint8Array>>
    uMarkerTexDim: ValueCell<Vec2>
    markerAverage: ValueCell<number>
    markerStatus: ValueCell<number>
}

export function getMarkersAverage(array: Uint8Array, count: number): number {
    if (count === 0) return 0;
    let sum = 0;
    for (let i = 0; i < count; ++i) {
        if (array[i]) sum += 1;
    }
    return sum / count;
}

export function createMarkers(count: number, markerData?: MarkerData): MarkerData {
    const markers = createTextureImage(Math.max(1, count), 1, Uint8Array, markerData && markerData.tMarker.ref.value.array);
    const average = getMarkersAverage(markers.array, count);
    const status = average === 0 ? 0 : -1;
    if (markerData) {
        ValueCell.update(markerData.tMarker, markers);
        ValueCell.update(markerData.uMarkerTexDim, Vec2.create(markers.width, markers.height));
        ValueCell.updateIfChanged(markerData.markerAverage, average);
        ValueCell.updateIfChanged(markerData.markerStatus, status);
        return markerData;
    } else {
        return {
            tMarker: ValueCell.create(markers),
            uMarkerTexDim: ValueCell.create(Vec2.create(markers.width, markers.height)),
            markerAverage: ValueCell.create(average),
            markerStatus: ValueCell.create(status),
        };
    }
}

const emptyMarkerTexture = { array: new Uint8Array(1), width: 1, height: 1 };
export function createEmptyMarkers(markerData?: MarkerData): MarkerData {
    if (markerData) {
        ValueCell.update(markerData.tMarker, emptyMarkerTexture);
        ValueCell.update(markerData.uMarkerTexDim, Vec2.create(1, 1));
        ValueCell.updateIfChanged(markerData.markerAverage, 0);
        ValueCell.updateIfChanged(markerData.markerStatus, 0);
        return markerData;
    } else {
        return {
            tMarker: ValueCell.create(emptyMarkerTexture),
            uMarkerTexDim: ValueCell.create(Vec2.create(1, 1)),
            markerAverage: ValueCell.create(0),
            markerStatus: ValueCell.create(0),
        };
    }
}
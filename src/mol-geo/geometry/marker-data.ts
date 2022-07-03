/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../mol-util/value-cell';
import { Vec2 } from '../../mol-math/linear-algebra';
import { TextureImage, createTextureImage } from '../../mol-gl/renderable/util';

export type MarkerType = 'instance' | 'groupInstance';

export type MarkerData = {
    uMarker: ValueCell<number>
    tMarker: ValueCell<TextureImage<Uint8Array>>
    uMarkerTexDim: ValueCell<Vec2>
    markerAverage: ValueCell<number>
    markerStatus: ValueCell<number>
    dMarkerType: ValueCell<string>
}

const MarkerCountLut = new Uint8Array(0x0303 + 1);
MarkerCountLut[0x0001] = 1;
MarkerCountLut[0x0002] = 1;
MarkerCountLut[0x0003] = 1;
MarkerCountLut[0x0100] = 1;
MarkerCountLut[0x0200] = 1;
MarkerCountLut[0x0300] = 1;
MarkerCountLut[0x0101] = 2;
MarkerCountLut[0x0201] = 2;
MarkerCountLut[0x0301] = 2;
MarkerCountLut[0x0102] = 2;
MarkerCountLut[0x0202] = 2;
MarkerCountLut[0x0302] = 2;
MarkerCountLut[0x0103] = 2;
MarkerCountLut[0x0203] = 2;
MarkerCountLut[0x0303] = 2;

/**
 * Calculates the average number of entries that have any marker flag set.
 *
 * For alternative implementations and performance tests see
 * `src\perf-tests\markers-average.ts`.
 */
export function getMarkersAverage(array: Uint8Array, count: number): number {
    if (count === 0) return 0;

    const view = new Uint32Array(array.buffer, 0, array.buffer.byteLength >> 2);
    const viewEnd = (count - 4) >> 2;
    const backStart = 4 * viewEnd;

    let sum = 0;
    if (viewEnd < 0) {
        // avoid edge cases with small arrays
        for (let i = 0; i < count; ++i) {
            sum += array[i] && 1;
        }
    } else {
        for (let i = 0; i < viewEnd; ++i) {
            const v = view[i];
            sum += MarkerCountLut[v & 0xFFFF] + MarkerCountLut[v >> 16];
        }
        for (let i = backStart; i < count; ++i) {
            sum += array[i] && 1;
        }
    }
    return sum / count;
}

export function createMarkers(count: number, type: MarkerType, markerData?: MarkerData): MarkerData {
    const markers = createTextureImage(Math.max(1, count), 1, Uint8Array, markerData && markerData.tMarker.ref.value.array);
    const average = getMarkersAverage(markers.array, count);
    const status = average === 0 ? 0 : -1;
    if (markerData) {
        ValueCell.updateIfChanged(markerData.uMarker, 0);
        ValueCell.update(markerData.tMarker, markers);
        ValueCell.update(markerData.uMarkerTexDim, Vec2.create(markers.width, markers.height));
        ValueCell.updateIfChanged(markerData.markerAverage, average);
        ValueCell.updateIfChanged(markerData.markerStatus, status);
        ValueCell.updateIfChanged(markerData.dMarkerType, type);
        return markerData;
    } else {
        return {
            uMarker: ValueCell.create(0),
            tMarker: ValueCell.create(markers),
            uMarkerTexDim: ValueCell.create(Vec2.create(markers.width, markers.height)),
            markerAverage: ValueCell.create(average),
            markerStatus: ValueCell.create(status),
            dMarkerType: ValueCell.create(type),
        };
    }
}

const emptyMarkerTexture = { array: new Uint8Array(1), width: 1, height: 1 };
export function createEmptyMarkers(markerData?: MarkerData): MarkerData {
    if (markerData) {
        ValueCell.updateIfChanged(markerData.uMarker, 0);
        ValueCell.update(markerData.tMarker, emptyMarkerTexture);
        ValueCell.update(markerData.uMarkerTexDim, Vec2.create(1, 1));
        ValueCell.updateIfChanged(markerData.markerAverage, 0);
        ValueCell.updateIfChanged(markerData.markerStatus, 0);
        return markerData;
    } else {
        return {
            uMarker: ValueCell.create(0),
            tMarker: ValueCell.create(emptyMarkerTexture),
            uMarkerTexDim: ValueCell.create(Vec2.create(1, 1)),
            markerAverage: ValueCell.create(0),
            markerStatus: ValueCell.create(0),
            dMarkerType: ValueCell.create('groupInstance'),
        };
    }
}
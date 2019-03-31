/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { Vec2 } from 'mol-math/linear-algebra'
import { TextureImage, createTextureImage } from 'mol-gl/renderable/util';

export type MarkerData = {
    tMarker: ValueCell<TextureImage<Uint8Array>>
    uMarkerTexDim: ValueCell<Vec2>
}

export enum MarkerAction {
    Highlight,
    RemoveHighlight,
    Select,
    Deselect,
    Toggle,
    Clear
}

export function applyMarkerAction(array: Uint8Array, start: number, end: number, action: MarkerAction) {
    let changed = false
    for (let i = start; i < end; ++i) {
        let v = array[i]
        switch (action) {
            case MarkerAction.Highlight:
                if (v % 2 === 0) {
                    v += 1
                }
                break
            case MarkerAction.RemoveHighlight:
                if (v % 2 !== 0) {
                    v -= 1
                }
                break
            case MarkerAction.Select:
                if (v < 2) v += 2
                // v += 2
                break
            case MarkerAction.Deselect:
                // if (v >= 2) {
                //     v -= 2
                // }
                v = v % 2
                break
            case MarkerAction.Toggle:
                if (v >= 2) {
                    v -= 2
                } else {
                    v += 2
                }
                break
            case MarkerAction.Clear:
                v = 0
                break
        }
        changed = array[i] !== v || changed
        array[i] = v
    }
    return changed
}

export function createMarkers(count: number, markerData?: MarkerData): MarkerData {
    const markers = createTextureImage(Math.max(1, count), 1, Uint8Array, markerData && markerData.tMarker.ref.value.array)
    if (markerData) {
        ValueCell.update(markerData.tMarker, markers)
        ValueCell.update(markerData.uMarkerTexDim, Vec2.create(markers.width, markers.height))
        return markerData
    } else {
        return {
            tMarker: ValueCell.create(markers),
            uMarkerTexDim: ValueCell.create(Vec2.create(markers.width, markers.height)),
        }
    }
}

const emptyMarkerTexture = { array: new Uint8Array(1), width: 1, height: 1 }
export function createEmptyMarkers(markerData?: MarkerData): MarkerData {
    if (markerData) {
        ValueCell.update(markerData.tMarker, emptyMarkerTexture)
        ValueCell.update(markerData.uMarkerTexDim, Vec2.create(1, 1))
        return markerData
    } else {
        return {
            tMarker: ValueCell.create(emptyMarkerTexture),
            uMarkerTexDim: ValueCell.create(Vec2.create(1, 1)),
        }
    }
}
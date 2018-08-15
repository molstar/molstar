/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { Vec2 } from 'mol-math/linear-algebra'
import { TextureImage, createTextureImage } from 'mol-gl/renderable/util';

export type MarkerData = {
    tMarker: ValueCell<TextureImage>
    uMarkerTexSize: ValueCell<Vec2>
}

export enum MarkerAction {
    Highlight,
    RemoveHighlight,
    Select,
    Deselect,
    ToggleSelect,
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
                    changed = true
                }
                break
            case MarkerAction.RemoveHighlight:
                if (v % 2 !== 0) {
                    v -= 1
                    changed = true
                }
                break
            case MarkerAction.Select:
                v += 2
                changed = true
                break
            case MarkerAction.Deselect:
                if (v >= 2) {
                    v -= 2
                    changed = true
                }
                break
            case MarkerAction.ToggleSelect:
                if (v === 0) {
                    v = 2
                } else if (v === 1) {
                    v = 3
                } else if (v === 2) {
                    v = 0
                } else {
                    v -= 2
                }
                changed = true
                break
            case MarkerAction.Clear:
                v = 0
                changed = true
                break
        }
        array[i] = v
    }
    return changed
}

export function createMarkers(count: number, markerData?: MarkerData): MarkerData {
    const markers = markerData && markerData.tMarker.ref.value.array.length >= count
        ? markerData.tMarker.ref.value
        : createTextureImage(count, 1)
    if (markerData) {
        ValueCell.update(markerData.tMarker, markers)
        ValueCell.update(markerData.uMarkerTexSize, Vec2.create(markers.width, markers.height))
        return markerData
    } else {
        return {
            tMarker: ValueCell.create(markers),
            uMarkerTexSize: ValueCell.create(Vec2.create(markers.width, markers.height)),
        }
    }
}

const emptyMarkerTexture = { array: new Uint8Array(1), width: 1, height: 1 }
export function createEmptyMarkers(markerData?: MarkerData) {
    if (markerData) {
        ValueCell.update(markerData.tMarker, emptyMarkerTexture)
        ValueCell.update(markerData.uMarkerTexSize, Vec2.create(1, 1))
        return markerData
    } else {
        return {
            tMarker: ValueCell.create(emptyMarkerTexture),
            uMarkerTexSize: ValueCell.create(Vec2.create(1, 1)),
        }
    }
}
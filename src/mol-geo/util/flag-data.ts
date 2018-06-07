/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { Vec2 } from 'mol-math/linear-algebra'
import { TextureImage, createTextureImage } from 'mol-gl/renderable/util';

export type FlagData = {
    tFlag: ValueCell<TextureImage>
    uFlagTexSize: ValueCell<Vec2>
}

export enum FlagAction {
    Highlight,
    RemoveHighlight,
    Select,
    Deselect,
    ToggleSelect,
    Clear
}

export function applyFlagAction(array: Uint8Array, start: number, end: number, action: FlagAction) {
    let changed = false
    for (let i = start; i < end; ++i) {
        let v = array[i]
        switch (action) {
            case FlagAction.Highlight:
                if (v % 2 === 0) {
                    v += 1
                    changed = true
                }
                break
            case FlagAction.RemoveHighlight:
                if (v % 2 !== 0) {
                    v -= 1
                    changed = true
                } 
                break
            case FlagAction.Select:
                v += 2
                changed = true
                break
            case FlagAction.Deselect:
                if (v >= 2) {
                    v -= 2
                    changed = true
                }
                break
            case FlagAction.ToggleSelect:
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
            case FlagAction.Clear:
                v = 0
                changed = true
                break
        }
        array[i] = v
    }
    return changed
}

export function createFlags(count: number, flagData?: FlagData): FlagData {
    const flags = flagData && flagData.tFlag.ref.value.array.length >= count
        ? flagData.tFlag.ref.value
        : createTextureImage(count, 1)
    if (flagData) {
        ValueCell.update(flagData.tFlag, flags)
        ValueCell.update(flagData.uFlagTexSize, Vec2.create(flags.width, flags.height))
        return flagData
    } else {
        return {
            tFlag: ValueCell.create(flags),
            uFlagTexSize: ValueCell.create(Vec2.create(flags.width, flags.height)),
        }
    }
}

const emptyFlagTexture = { array: new Uint8Array(1), width: 1, height: 1 }
export function createEmptyFlags(flagData?: FlagData) {
    if (flagData) {
        ValueCell.update(flagData.tFlag, emptyFlagTexture)
        ValueCell.update(flagData.uFlagTexSize, Vec2.create(1, 1))
        return flagData
    } else {
        return {
            tFlag: ValueCell.create(emptyFlagTexture),
            uFlagTexSize: ValueCell.create(Vec2.create(1, 1)),
        }
    }
}
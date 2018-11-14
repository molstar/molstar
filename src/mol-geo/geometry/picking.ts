/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

function decodeFloatRGBA(r: number, g: number, b: number) {
    r = Math.floor(r)
    g = Math.floor(g)
    b = Math.floor(b)
    return r * 256 * 256 + g * 256 + b
}

export function decodeIdRGB(r: number, g: number, b: number) {
    return decodeFloatRGBA(r, g, b) - 1
}

export interface PickingId {
    objectId: number
    instanceId: number
    groupId: number
}

export namespace PickingId {
    export function areSame(a: PickingId, b: PickingId) {
        return a.objectId === b.objectId && a.instanceId === b.instanceId && a.groupId === b.groupId;
    }
}

export interface PickingInfo {
    label: string
    data?: any
}
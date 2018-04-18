/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export { ElementColor } from './structure/element';

export function hexColorToArray(hexColor: number, array: Helpers.NumberArray, offset: number) {
    array[ offset ] = (hexColor >> 16 & 255) / 255
    array[ offset + 1 ] = (hexColor >> 8 & 255) / 255
    array[ offset + 2 ] = (hexColor & 255) / 255
}
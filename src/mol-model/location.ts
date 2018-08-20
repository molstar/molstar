/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement } from './structure'
import { Link } from './structure/structure/unit/links'

/** A null value Location */
export const NullLocation = { kind: 'null-location' as 'null-location' }
export type NullLocation = typeof NullLocation
export function isNullLocation(x: any): x is NullLocation {
    return !!x && x.kind === 'null-location';
}

/** A custom Location */
export interface CustomLocation<D = any, K = any> {
    readonly kind: 'custom-location'
    data: D
    key: K
}
export function CustomLocation<D, K>(data: D, key: K): CustomLocation<D, K> {
    return { kind: 'custom-location', data, key }
}
export function isCustomLocation(x: any): x is CustomLocation<any, any> {
    return !!x && x.kind === 'custom-location';
}

export type Location = StructureElement | Link.Location | NullLocation | CustomLocation
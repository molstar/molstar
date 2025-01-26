/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement } from './structure';
import { Bond } from './structure/structure/unit/bonds';
import { ShapeGroup } from './shape/shape';
import { PositionLocation } from '../mol-geo/util/location-iterator';
import { Volume } from './volume';

/** A null value Location */
export const NullLocation = { kind: 'null-location' as const };
export type NullLocation = typeof NullLocation
export function isNullLocation(x: any): x is NullLocation {
    return !!x && x.kind === 'null-location';
}

/** A generic data Location */
export interface DataLocation<T = unknown, E = unknown> {
    readonly kind: 'data-location',
    readonly tag: string
    readonly data: T,
    element: E
}
export function DataLocation<T = unknown, E = unknown>(tag: string, data: T, element: E): DataLocation<T, E> {
    return { kind: 'data-location', tag, data, element };
}
export function isDataLocation(x: any): x is DataLocation {
    return !!x && x.kind === 'data-location';
}

/**
 * A direct Location.
 *
 * For it, the location is implicitly clear from context and is not explicitly given.
 * This is used for themes with direct-volume rendering where the location is the volume
 * grid cell itself and coloring is applied in a shader on the GPU.
 */
export const DirectLocation = { kind: 'direct-location' as const };
export type DirectLocation = typeof DirectLocation
export function isDirectLocation(x: any): x is DirectLocation {
    return !!x && x.kind === 'direct-location';
}

export type Location = StructureElement.Location | Bond.Location | ShapeGroup.Location | PositionLocation | DataLocation | NullLocation | DirectLocation | Volume.Cell.Location | Volume.Segment.Location

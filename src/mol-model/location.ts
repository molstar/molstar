/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement } from './structure';
import { Bond } from './structure/structure/unit/bonds';
import { ShapeGroup } from './shape/shape';

/** A null value Location */
export const NullLocation = { kind: 'null-location' as 'null-location' };
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

export type Location = StructureElement.Location | Bond.Location | ShapeGroup.Location | DataLocation | NullLocation
/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mesh } from 'mol-geo/mesh/mesh';
import { Color } from 'mol-util/color';
import { UUID } from 'mol-util';
import { OrderedSet } from 'mol-data/int';

export interface Shape {
    readonly id: UUID
    readonly name: string
    readonly mesh: Mesh
    getColor(group: number): Color
    getLabel(group: number): string
}

export namespace Shape {
    export function create(mesh: Mesh, name: string, getColor: (group: number) => Color, getLabel: (group: number) => string): Shape {
        return { id: UUID.create(), name, mesh, getColor, getLabel }
    }

    export interface Location {
        readonly kind: 'group-location'
        shape: Shape
        group: number
    }

    export function Location(shape?: Shape, group?: number): Location {
        return { kind: 'group-location', shape: shape!, group: group || 0 };
    }

    export function isLocation(x: any): x is Location {
        return !!x && x.kind === 'group-location';
    }

    export interface Loci {
        readonly kind: 'group-loci',
        readonly groups: ReadonlyArray<{
            shape: Shape,
            ids: OrderedSet<number>
        }>
    }

    export function Loci(groups: ArrayLike<{ shape: Shape, ids: OrderedSet<number> }>): Loci {
        return { kind: 'group-loci', groups: groups as Loci['groups'] };
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'group-loci';
    }
}
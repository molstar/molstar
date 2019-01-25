/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from 'mol-util/color';
import { UUID } from 'mol-util';
import { OrderedSet } from 'mol-data/int';
import { Geometry } from 'mol-geo/geometry/geometry';
import { Mat4 } from 'mol-math/linear-algebra';

export interface Shape<G extends Geometry = Geometry> {
    readonly id: UUID
    readonly name: string
    readonly geometry: G
    readonly transforms: Mat4[]
    readonly groupCount: number
    getColor(groupId: number): Color
    getLabel(groupId: number): string
}

export namespace Shape {
    export function create<G extends Geometry>(name: string, geometry: G, getColor: Shape['getColor'], getLabel: Shape['getLabel'], transforms?: Mat4[]): Shape<G> {
        return {
            id: UUID.create22(),
            name,
            geometry,
            transforms: transforms || [Mat4.identity()],
            get groupCount() { return Geometry.getGroupCount(geometry) },
            getColor,
            getLabel
        }
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
        readonly shape: Shape,
        readonly groups: ReadonlyArray<{
            ids: OrderedSet<number>
        }>
    }

    export function Loci(shape: Shape, groups: ArrayLike<{ ids: OrderedSet<number> }>): Loci {
        return { kind: 'group-loci', shape, groups: groups as Loci['groups'] };
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'group-loci';
    }

    export function areLociEqual(a: Loci, b: Loci) {
        if (a.shape !== b.shape) return false
        if (a.groups.length !== b.groups.length) return false
        for (let i = 0, il = a.groups.length; i < il; ++i) {
            const groupA = a.groups[i]
            const groupB = b.groups[i]
            if (!OrderedSet.areEqual(groupA.ids, groupB.ids)) return false
        }
        return true
    }
}
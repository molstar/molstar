/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from 'mol-util/color';
import { UUID } from 'mol-util';
import { OrderedSet } from 'mol-data/int';
import { Geometry } from 'mol-geo/geometry/geometry';
import { Mat4 } from 'mol-math/linear-algebra';

export interface Shape<G extends Geometry = Geometry> {
    /** A uuid to identify a shape object */
    readonly id: UUID
    /** A name to describe the shape */
    readonly name: string
    /** The data used to create the shape */
    readonly sourceData: unknown
    /** The geometry of the shape, e.g. `Mesh` or `Lines` */
    readonly geometry: G
    /** An array of transformation matrices to describe multiple instances of the geometry */
    readonly transforms: Mat4[]
    /** Number of groups in the geometry */
    readonly groupCount: number
    /** Get color for a given group */
    getColor(groupId: number, instanceId: number): Color
    /** Get size for a given group */
    getSize(groupId: number, instanceId: number): number
    /** Get label for a given group */
    getLabel(groupId: number, instanceId: number): string
}

export namespace Shape {
    export function create<G extends Geometry>(name: string, sourceData: unknown, geometry: G, getColor: Shape['getColor'], getSize: Shape['getSize'], getLabel: Shape['getLabel'], transforms?: Mat4[]): Shape<G> {
        return {
            id: UUID.create22(),
            name,
            sourceData,
            geometry,
            transforms: transforms || [Mat4.identity()],
            get groupCount() { return Geometry.getGroupCount(geometry) },
            getColor,
            getSize,
            getLabel
        }
    }

    export interface Loci { readonly kind: 'shape-loci', readonly shape: Shape }
    export function Loci(shape: Shape): Loci { return { kind: 'shape-loci', shape } }
    export function isLoci(x: any): x is Loci { return !!x && x.kind === 'shape-loci' }
    export function areLociEqual(a: Loci, b: Loci) { return a.shape === b.shape }
}

export namespace ShapeGroup {
    export interface Location {
        readonly kind: 'group-location'
        shape: Shape
        group: number
        instance: number
    }

    export function Location(shape?: Shape, group?: number, instance?: number): Location {
        return { kind: 'group-location', shape: shape!, group: group || 0, instance: instance || 0 };
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
        readonly instance: number
    }

    export function Loci(shape: Shape, groups: ArrayLike<{ ids: OrderedSet<number> }>, instance: number): Loci {
        return { kind: 'group-loci', shape, groups: groups as Loci['groups'], instance };
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'group-loci';
    }

    export function areLociEqual(a: Loci, b: Loci) {
        if (a.shape !== b.shape) return false
        if (a.groups.length !== b.groups.length) return false
        if (a.instance !== b.instance) return false
        for (let i = 0, il = a.groups.length; i < il; ++i) {
            const groupA = a.groups[i]
            const groupB = b.groups[i]
            if (!OrderedSet.areEqual(groupA.ids, groupB.ids)) return false
        }
        return true
    }
}
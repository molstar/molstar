/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { Color } from 'mol-util/color';
import { UUID, ValueCell } from 'mol-util';
import { OrderedSet } from 'mol-data/int';
import { arrayMax } from 'mol-util/array';

export interface Shape {
    readonly id: UUID
    readonly name: string
    readonly mesh: Mesh
    readonly colors: ValueCell<Color[]>
    readonly labels: ValueCell<string[]>
    readonly groupCount: number
}

export namespace Shape {
    export function create(name: string, mesh: Mesh, colors: Color[], labels: string[]): Shape {
        let currentGroupBufferVersion = -1
        let currentGroupCount = -1

        return {
            id: UUID.create(),
            name,
            mesh,
            get groupCount() {
                if (mesh.groupBuffer.ref.version !== currentGroupBufferVersion) {
                    currentGroupCount = arrayMax(mesh.groupBuffer.ref.value) + 1
                }
                return currentGroupCount
            },
            colors: ValueCell.create(colors),
            labels: ValueCell.create(labels),
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
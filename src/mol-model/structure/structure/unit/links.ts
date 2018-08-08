/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement } from '../../structure'

export * from './links/data'
export * from './links/intra-compute'
export * from './links/inter-compute'

namespace Link {
    export interface Location {
        readonly kind: 'link-location',
        aUnit: Unit,
        /** Index into aUnit.elements */
        aIndex: StructureElement.UnitIndex,
        bUnit: Unit,
        /** Index into bUnit.elements */
        bIndex: StructureElement.UnitIndex,
    }

    export function Location(aUnit?: Unit, aIndex?: StructureElement.UnitIndex, bUnit?: Unit, bIndex?: StructureElement.UnitIndex): Location {
        return { kind: 'link-location', aUnit: aUnit as any, aIndex: aIndex as any, bUnit: bUnit as any, bIndex: bIndex as any };
    }

    export function isLocation(x: any): x is Location {
        return !!x && x.kind === 'link-location';
    }

    export interface Loci {
        readonly kind: 'link-loci',
        readonly links: ReadonlyArray<Location>
    }

    export function Loci(links: ArrayLike<Location>): Loci {
        return { kind: 'link-loci', links: links as Loci['links'] };
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'link-loci';
    }
}

export { Link }
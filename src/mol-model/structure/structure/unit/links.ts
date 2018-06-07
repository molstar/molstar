/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Unit } from '../../structure'

export * from './links/data'
export * from './links/intra-compute'
export * from './links/inter-compute'

namespace Link {
    export interface Location {
        readonly aUnit: Unit,
        /** Index into aUnit.elements */
        readonly aIndex: number,
        readonly bUnit: Unit,
        /** Index into bUnit.elements */
        readonly bIndex: number,
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
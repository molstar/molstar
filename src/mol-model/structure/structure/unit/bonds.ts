/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Unit } from '../../structure'

export * from './bonds/intra-data'
export * from './bonds/intra-compute'

namespace Bond {
    export interface Location {
        readonly aUnit: Unit,
        /** Index into aUnit.elements */
        readonly aIndex: number,
        readonly bUnit: Unit,
        /** Index into bUnit.elements */
        readonly bIndex: number,
    }

    export interface Loci {
        readonly kind: 'bond-loci',
        readonly bonds: ReadonlyArray<Location>
    }

    export function Loci(bonds: ArrayLike<Location>): Loci {
        return { kind: 'bond-loci', bonds: bonds as Loci['bonds'] };
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'bond-loci';
    }
}

export { Bond }
/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Element } from '../../structure'

export * from './bonds/intra-data'
export * from './bonds/intra-compute'

interface Bond {
    readonly a: Readonly<Element.Location>,
    readonly b: Readonly<Element.Location>
}

namespace Bond {
    export interface Loci {
        readonly kind: 'bond-loci',
        readonly bonds: ReadonlyArray<Bond>
    }

    export function Loci(bonds: ArrayLike<Bond>): Loci {
        return { kind: 'bond-loci', bonds: bonds as Loci['bonds'] };
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'bond-loci';
    }
}

export { Bond }
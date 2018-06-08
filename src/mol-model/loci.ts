/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Element } from './structure'
import { Link } from './structure/structure/unit/links'

/** A Loci that includes every loci */
export const EveryLoci = { kind: 'every-loci' as 'every-loci' }
export type EveryLoci = typeof EveryLoci
export function isEveryLoci(x: any): x is EveryLoci {
    return !!x && x.kind === 'every-loci';
}

export type Loci =  Element.Loci | Link.Loci | EveryLoci
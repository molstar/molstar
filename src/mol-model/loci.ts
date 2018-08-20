/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement } from './structure'
import { Link } from './structure/structure/unit/links'
import { CustomLocation } from './location';

/** A Loci that includes every loci */
export const EveryLoci = { kind: 'every-loci' as 'every-loci' }
export type EveryLoci = typeof EveryLoci
export function isEveryLoci(x: any): x is EveryLoci {
    return !!x && x.kind === 'every-loci';
}

/** A Loci that is empty */
export const EmptyLoci = { kind: 'empty-loci' as 'empty-loci' }
export type EmptyLoci = typeof EmptyLoci
export function isEmptyLoci(x: any): x is EmptyLoci {
    return !!x && x.kind === 'empty-loci';
}

/** A Loci of custom locations */
export interface CustomLoci<D = any, K = any> {
    readonly kind: 'custom-loci'
    readonly locations: ReadonlyArray<CustomLocation>
}
export function CustomLoci<D, K>(locations: CustomLocation<D, K>[]): CustomLoci<D, K> {
    return { kind: 'custom-loci', locations }
}
export function isCustomLoci(x: any): x is CustomLoci {
    return !!x && x.kind === 'custom-loci';
}

export type Loci =  StructureElement.Loci | Link.Loci | EveryLoci | EmptyLoci | CustomLoci
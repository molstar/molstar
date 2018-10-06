/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AttachModelProperties } from '../property-provider';
import { RCSB_assemblySymmetry } from './providers/rcsb';

export const attachModelProperties: AttachModelProperties = (args) => {
    // return a list of promises that start attaching the props in parallel
    // (if there are downloads etc.)
    return [
        RCSB_assemblySymmetry(args)
    ];
}
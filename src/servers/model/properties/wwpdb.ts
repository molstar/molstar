/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AttachModelProperties } from '../property-provider';
import { wwPDB_chemCompBond } from './providers/wwpdb';

export const attachModelProperties: AttachModelProperties = (args) => {
    // return a list of promises that start attaching the props in parallel
    // (if there are downloads etc.)
    return [
        wwPDB_chemCompBond(args)
    ];
};
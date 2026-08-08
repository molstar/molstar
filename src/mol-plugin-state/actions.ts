/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as Structure from './actions/structure';
import * as Volume from './actions/volume';
import * as Particles from './actions/particles';
import * as DataFormat from './actions/file';

export const StateActions = {
    Structure,
    Volume,
    Particles,
    DataFormat
};
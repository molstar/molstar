/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Data from './transforms/data';
import * as Misc from './transforms/misc';
import * as Model from './transforms/model';
import * as Volume from './transforms/volume';
import * as Representation from './transforms/representation';

export const StateTransforms = {
    Data,
    Misc,
    Model,
    Volume,
    Representation
};

export type StateTransforms = typeof StateTransforms
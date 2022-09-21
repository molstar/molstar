

import { loadCheckpoint } from '../mol-util/debug';
loadCheckpoint(`mol-plugin-state/actions.ts::start`);
/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Structure from './actions/structure';
import * as Volume from './actions/volume';
import * as DataFormat from './actions/file';

export const StateActions = {
    Structure,
    Volume,
    DataFormat
};
loadCheckpoint(`mol-plugin-state/actions.ts::end`);

/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from 'mol-state';
import Viewer from 'mol-canvas3d/viewer';

export class PluginContext {
    state = {
        data: State,
        behaviour: State,
        plugin: State
    };

    viewer: Viewer;

    // logger = ;
    // settings = ;
}
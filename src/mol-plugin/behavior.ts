/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export * from './behavior/behavior';

import * as StaticState from './behavior/static/state';
import * as StaticRepresentation from './behavior/static/representation';
import * as StaticCamera from './behavior/static/camera';
import * as StaticMisc from './behavior/static/misc';

import * as DynamicRepresentation from './behavior/dynamic/representation';
import * as DynamicCamera from './behavior/dynamic/camera';
import * as DynamicCustomProps from './behavior/dynamic/custom-props';

export const BuiltInPluginBehaviors = {
    State: StaticState,
    Representation: StaticRepresentation,
    Camera: StaticCamera,
    Misc: StaticMisc
};

export const PluginBehaviors = {
    Representation: DynamicRepresentation,
    Camera: DynamicCamera,
    CustomProps: DynamicCustomProps
};
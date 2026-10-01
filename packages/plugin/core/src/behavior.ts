/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export * from '@molstar/plugin/behavior/behavior';

import * as StaticState from '@molstar/plugin/behavior/static/state';
import * as StaticRepresentation from '@molstar/plugin/behavior/static/representation';
import * as StaticCamera from '@molstar/plugin/behavior/static/camera';
import * as StaticMisc from '@molstar/plugin/behavior/static/misc';

import * as DynamicRepresentation from '@molstar/plugin/behavior/dynamic/representation';
import * as DynamicCamera from '@molstar/plugin/behavior/dynamic/camera';
import * as DynamicState from '@molstar/plugin/behavior/dynamic/state';
import * as DynamicCustomProps from '@molstar/plugin/behavior/dynamic/custom-props';

export const BuiltInPluginBehaviors = {
    State: StaticState,
    Representation: StaticRepresentation,
    Camera: StaticCamera,
    Misc: StaticMisc
};

export const PluginBehaviors = {
    Representation: DynamicRepresentation,
    Camera: DynamicCamera,
    State: DynamicState,
    CustomProps: DynamicCustomProps
};
/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as DynamicRepresentation from '@molstar/plugin/behavior/dynamic/representation';
import * as DynamicCamera from '@molstar/plugin/behavior/dynamic/camera';
import * as DynamicState from '@molstar/plugin/behavior/dynamic/state';
import * as DynamicCustomProps from '@molstar/plugin/behavior/dynamic/custom-props';

export const PluginBehaviors = {
  Representation: DynamicRepresentation,
  Camera: DynamicCamera,
  State: DynamicState,
  CustomProps: DynamicCustomProps,
};

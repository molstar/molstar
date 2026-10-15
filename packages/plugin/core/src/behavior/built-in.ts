/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as StaticState from '@molstar/plugin/behavior/static/state';
import * as StaticRepresentation from '@molstar/plugin/behavior/static/representation';
import * as StaticCamera from '@molstar/plugin/behavior/static/camera';
import * as StaticMisc from '@molstar/plugin/behavior/static/misc';

export const BuiltInPluginBehaviors = {
  State: StaticState,
  Representation: StaticRepresentation,
  Camera: StaticCamera,
  Misc: StaticMisc,
};

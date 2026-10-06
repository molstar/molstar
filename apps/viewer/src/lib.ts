/**
 * Copyright (c) 2025-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Structure from '@molstar/model/model/structure';
import { DataLoci, EveryLoci, Loci } from '@molstar/model/model/loci';
import { Volume } from '@molstar/model/model/volume';
import { Shape, ShapeGroup } from '@molstar/model/model/shape';
import * as LinearAlgebra3D from '@molstar/core/math/linear-algebra/3d';
import { PluginContext } from '@molstar/plugin/context';
import { PluginUIContext } from '@molstar/plugin-ui/context';
import { PluginConfig } from '@molstar/plugin/config';
import { PluginBehavior } from '@molstar/plugin/behavior';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { PluginSpec } from '@molstar/plugin/spec';
import { DefaultPluginUISpec } from '@molstar/plugin-ui/default-spec';
import { PluginStateObject, PluginStateTransform } from '@molstar/plugin/state/objects';
import { StateActions } from '@molstar/plugin/state/actions';
import { StateTransforms } from '@molstar/viewer/state-transforms';
import { PluginExtensions } from '@molstar/viewer/extensions';

export const lib = {
  structure: {
    ...Structure,
  },
  volume: {
    Volume,
  },
  shape: {
    Shape,
    ShapeGroup,
  },
  loci: {
    Loci,
    DataLoci,
    EveryLoci,
  },
  math: {
    LinearAlgebra: {
      ...LinearAlgebra3D,
    },
  },
  plugin: {
    PluginContext,
    PluginUIContext,
    PluginConfig,
    PluginBehavior,
    PluginSpec,
    PluginStateObject,
    PluginStateTransform,
    StateTransforms,
    StateActions,
    DefaultPluginSpec,
    DefaultPluginUISpec,
  },
  extensions: {
    ...PluginExtensions,
  },
};

/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ventura Rivera <venturaxrivera@gmail.com>
 */

import { CreateVolumeStreamingBehavior } from '@molstar/plugin/behavior/dynamic/volume-streaming/transformers';
import { VolumeStreamingCustomControls } from '@molstar/plugin-ui/custom/volume';
import { DefaultStructureTools } from '@molstar/plugin-ui/default-structure-tools';
import type { PluginUISpec } from '@molstar/plugin-ui/spec';

// The UI parts of `DefaultPluginUISpec()` without its registry, for apps that compose their own registry. Importing
// `@molstar/plugin-ui/default-spec` loads `DefaultRegistry` and every built-in catalog.

export const DefaultPluginUIComponents = (): NonNullable<PluginUISpec['components']> => ({
  structureTools: DefaultStructureTools,
});

export const DefaultPluginUICustomParamEditors = (): NonNullable<PluginUISpec['customParamEditors']> => [
  [CreateVolumeStreamingBehavior, VolumeStreamingCustomControls],
];

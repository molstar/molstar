/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ventura Rivera <venturaxrivera@gmail.com>
 */

import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { DefaultPluginUIComponents, DefaultPluginUICustomParamEditors } from '@molstar/plugin-ui/default-ui';
import type { PluginUISpec } from '@molstar/plugin-ui/spec';

export const DefaultPluginUISpec = (): PluginUISpec => ({
  ...DefaultPluginSpec(),
  components: DefaultPluginUIComponents(),
  customParamEditors: DefaultPluginUICustomParamEditors(),
});

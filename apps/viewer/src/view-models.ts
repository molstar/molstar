/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import type { PluginSpec } from '@molstar/plugin/spec';
import { DefaultPluginUISpec } from '@molstar/plugin-ui/default-spec';
import type { PluginUISpec } from '@molstar/plugin-ui/spec';
import { PluginViewModel } from '@molstar/plugin-extension/view-model';
import { PluginUIViewModel } from '@molstar/plugin-extension/ui-view-model';

// The library view models require an explicit spec. The classic Viewer global keeps the 5.x behavior of using the
// default spec when none is given.

export class ViewerPluginViewModel extends PluginViewModel {
  constructor(options?: { spec?: PluginSpec }) {
    super({ spec: options?.spec ?? DefaultPluginSpec() });
  }
}

export class ViewerPluginUIViewModel extends PluginUIViewModel {
  constructor(options?: { spec?: PluginUISpec }) {
    super({ spec: options?.spec ?? DefaultPluginUISpec() });
  }
}

/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { ChemicalComponentPreset, ChemicalCompontentTrajectoryHierarchyPreset } from './representation.js';

export const wwPDBChemicalComponentDictionary = PluginBehavior.create<{}>({
  name: 'wwpdb-chemical-component-dictionary',
  category: 'representation',
  display: {
    name: 'wwPDB Chemical Compontent Dictionary',
    description: 'Custom representation for data loaded from the CCD.',
  },
  ctor: class extends PluginBehavior.Handler<{}> {
    private unregisterEntry: (() => void) | undefined;

    register(): void {
      const entry: PluginRegistryEntry = {
        structure: {
          presets: {
            hierarchy: [ChemicalCompontentTrajectoryHierarchyPreset],
            representation: [ChemicalComponentPreset],
          },
        },
      };
      this.unregisterEntry = this.ctx.register(entry);
    }

    update() {
      return false;
    }

    unregister() {
      this.unregisterEntry?.();
      this.unregisterEntry = undefined;
    }
  },
  params: () => ({}),
});

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 * The slim-plugin acceptance target (.v6/designs/plugin-composition.md §12): a plugin that lists only what it needs
 * and renders an SDF ligand as ball-and-stick. The bundle excludes every module in `scripts/workspace/slim-exclusions.json`
 * (mmCIF, CCP4, cartoon, the volume representations, the other presets, the query catalog, script transpilers, MP4 export);
 * `pnpm check:workspace` verifies this.
 */

import { HighlightLoci, SelectLoci } from '@molstar/plugin/behavior/dynamic/representation';
import { FocusLoci as CameraFocusLoci } from '@molstar/plugin/behavior/dynamic/camera';
import { PluginConfig } from '@molstar/plugin/config';
import { PluginSpec } from '@molstar/plugin/spec';
import { DefaultHierarchyPresetEntry } from '@molstar/plugin/state/builder/structure/hierarchy-presets/default';
import { BallAndStickPresetEntry } from '@molstar/plugin/state/builder/structure/representation-presets/ball-and-stick';
import { Sdf } from '@molstar/plugin/state/formats/trajectory/sdf';
import { Asset } from '@molstar/core/util/assets';
import { createPluginUI } from '@molstar/plugin-ui';
import type { PluginUIContext } from '@molstar/plugin-ui/context';
import { renderReact18 } from '@molstar/plugin-ui/react18';
import type { PluginUISpec } from '@molstar/plugin-ui/spec';
import '@molstar/plugin-ui/skin/light.scss';
import './index.html';
import './ligand.sdf';

const spec: PluginUISpec = {
  registry: [Sdf, DefaultHierarchyPresetEntry, BallAndStickPresetEntry],
  behaviors: [
    PluginSpec.Behavior(HighlightLoci),
    PluginSpec.Behavior(SelectLoci),
    PluginSpec.Behavior(CameraFocusLoci),
  ],
  config: [[PluginConfig.Structure.DefaultRepresentationPreset, 'preset-structure-representation-ball-and-stick']],
};

async function init(target: string | HTMLElement) {
  const plugin = await createPluginUI({
    target: typeof target === 'string' ? document.getElementById(target)! : target,
    spec,
    render: renderReact18,
  });

  const data = await plugin.builders.data.download({ url: Asset.Url('ligand.sdf') }, { state: { isGhost: true } });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'sdf');
  await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'preset-trajectory-default');

  window.slimPlugin = plugin;
  return plugin;
}

declare global {
  interface Window {
    SlimPlugin: { init(target: string | HTMLElement): Promise<PluginUIContext> };
    /** The plugin, set once the ligand is loaded. For smoke checks. */
    slimPlugin?: PluginUIContext;
    /** Resolves with the plugin once the ligand is loaded, rejects when loading fails. For smoke checks. */
    slimPluginReady?: Promise<PluginUIContext>;
  }
}

window.SlimPlugin = {
  init(target) {
    return (window.slimPluginReady = init(target));
  },
};

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultHierarchyPresetEntry } from '@molstar/plugin/state/builder/structure/hierarchy-presets/default';
import { BallAndStickPresetEntry } from '@molstar/plugin/state/builder/structure/representation-presets/ball-and-stick';
import { PolymerCartoonPresetEntry } from '@molstar/plugin/state/builder/structure/representation-presets/polymer-cartoon';
import { Mmcif } from '@molstar/plugin/state/formats/trajectory/mmcif';
import type { PrismViewerSpec } from '@molstar/prism-viewer/spec';

/**
 * The smallest Prism viewer composition that loads an mmCIF structure and shows it as cartoon or ball-and-stick: not
 * `DefaultPrismViewerSpec()` minus parts, but the spec written out from nothing (viewer-spec.md §1).
 */
export function MinimalPrismViewerSpec(): PrismViewerSpec {
  return {
    plugin: {
      // The engine's registry (viewer-spec.md §2.1): the mmCIF format with its parsers, the default hierarchy preset
      // (which makes the structure), and the two representation presets; each preset entry brings the representation
      // and the themes it builds with. No behaviors: no highlight, selection or focus.
      registry: [Mmcif, DefaultHierarchyPresetEntry, PolymerCartoonPresetEntry, BallAndStickPresetEntry],
      behaviors: [],
    },
    // None of the viewer's features (loaders, snapshots, animation, export, …) and no extensions (§3).
    extensions: [],
    // No sequence viewer (Seq* is not bundled) and no tools or assistant (§2.3, §2.4).
    sequence: false,
    agent: false,
    // The canvas's two remaining buttons (the pick mode switch, the trajectory bar), off by option (§2.2).
    options: { ui: { buttons: { selectionMode: false, trajectory: false } } },
    // Nothing is stored in the browser: no `persist` (§2.6).
    // The shell's parts, removed by id (§4): every tab and rail button, so only the canvas is left.
    overrides: {
      'core.tab.home': null,
      'core.tab.structure': null,
      'core.tab.selection': null,
      'core.tab.log': null,
      'core.tab.history': null,
      'core.tab.settings': null,
      'core.tab.help': null,
      'core.tab.sequence': null,
      'core.toolbar.theme': null,
      'core.toolbar.palette': null,
      'core.toolbar.reset': null,
      'core.toolbar.expand': null,
      'core.toolbar.export': null,
      'core.toolbar.sequence': null,
      'core.overlay.export': null,
      'core.overlay.sequence-legend': null,
    },
  };
}

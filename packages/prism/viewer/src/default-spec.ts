/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { AnnotationsFeature } from './features/annotations/index.js';
import { HelpFeature } from './features/help/index.js';
import { HomeFeature } from './features/home/index.js';
import { SelectionPanelFeature } from './features/selection-panel/index.js';
import { StructurePanelFeature } from './features/structure-panel/index.js';
import { StructureFeature } from './features/structure/index.js';
import { TrajectoryFeature } from './features/trajectory/index.js';
import { VolumeFeature } from './features/volume/index.js';
import { AllQueryLanguages } from './targets/languages/all.js';
import { AnimationFeature } from './features/animation/index.js';
import { CanvasSettingsFeature } from './features/canvas-settings/index.js';
import { ExportImageFeature } from './features/export-image/index.js';
import { ExportSessionFeature } from './features/export-session/index.js';
import { LazyParticles, LazyRemoteStates, LazyStateTree, LazySuperposition } from './features/lazy/index.js';
import { LoadersFeature } from './features/loaders/index.js';
import { ClippingFeature } from './features/clipping/index.js';
import { MeasurementsFeature } from './features/measurements/index.js';
import { SnapshotsFeature } from './features/snapshots/index.js';
import { prismViewerPluginParts } from './sources/base.js';
import { DefaultSequenceTracks } from './sequence/tracks.js';
import type { PrismViewerSpec } from './spec.js';

// The viewer's default spec (viewer-spec.md §5): its base, the default sequence and agent (no assistant), its own
// built-in features as extensions, nothing kept. Mol*'s extensions are
// the app's: they depend on the viewer.

/** A new object each call. */
export function DefaultPrismViewerSpec(): PrismViewerSpec {
  return {
    plugin: prismViewerPluginParts(),
    extensions: [
      // The scene's API (structure, components, representations, styles), maps and particles, trajectories,
      // annotations: first, as the rest use them.
      StructureFeature,
      VolumeFeature,
      LazyParticles,
      TrajectoryFeature,
      AnnotationsFeature,
      // The panels: Home, Structure, Selection, Help (with the core's tab ids).
      HomeFeature,
      StructurePanelFeature,
      SelectionPanelFeature,
      HelpFeature,
      LoadersFeature,
      CanvasSettingsFeature,
      MeasurementsFeature,
      ClippingFeature,
      LazySuperposition,
      LazyStateTree,
      LazyRemoteStates,
      SnapshotsFeature,
      AnimationFeature,
      ExportImageFeature,
      ExportSessionFeature,
    ],
    // PyMOL, VMD, Jmol and MolScript for target text.
    languages: AllQueryLanguages,
    // Seq*, its own chunk (loaded when a structure is first shown).
    sequence: { runtime: () => import('./sequence/runtime.js'), behaviors: [], tracks: DefaultSequenceTracks },
    agent: {},
  };
}

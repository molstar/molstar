/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolViewSpecBehavior } from '@molstar/mvs/behavior';
import { PluginBehaviors } from '@molstar/plugin/behavior';
import { StructureFocusRepresentation } from '@molstar/plugin/behavior/dynamic/selection/structure-focus-representation';
import { DefaultRegistry } from '@molstar/plugin/default-registry';
import { PluginSpec } from '@molstar/plugin/spec';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import type { PrismViewerPluginParts } from '../spec.js';

// The viewer's base (viewer-spec.md §2.1): the registry entries and behaviors its own
// features need, exported as named parts so an app narrows or swaps them by identity. Behaviors that only drive the
// classic UI are left out; nothing of `PluginUISpec`.

/**
 * `DefaultRegistry`, for parity. It holds everything `MVSRuntimeRegistry` lists; the `MolViewSpec` behavior registers
 * MVS's own themes and label representations.
 */
export const PrismViewerRegistry: readonly PluginRegistryEntry[] = DefaultRegistry;

/** Highlight, select, focus, camera, snapshot keys, the focus representation, the default custom properties, MVS. */
export const PrismViewerBehaviors = {
  highlight: PluginSpec.Behavior(PluginBehaviors.Representation.HighlightLoci),
  select: PluginSpec.Behavior(PluginBehaviors.Representation.SelectLoci),
  lociLabels: PluginSpec.Behavior(PluginBehaviors.Representation.DefaultLociLabelProvider),
  focusRepresentation: PluginSpec.Behavior(PluginBehaviors.Representation.FocusLoci),
  cameraFocus: PluginSpec.Behavior(PluginBehaviors.Camera.FocusLoci),
  cameraAxes: PluginSpec.Behavior(PluginBehaviors.Camera.CameraAxisHelper),
  cameraControls: PluginSpec.Behavior(PluginBehaviors.Camera.CameraControls),
  snapshotControls: PluginSpec.Behavior(PluginBehaviors.State.SnapshotControls),
  structureFocus: PluginSpec.Behavior(StructureFocusRepresentation),
  structureInfo: PluginSpec.Behavior(PluginBehaviors.CustomProps.StructureInfo),
  accessibleSurfaceArea: PluginSpec.Behavior(PluginBehaviors.CustomProps.AccessibleSurfaceArea),
  sequenceMapping: PluginSpec.Behavior(PluginBehaviors.CustomProps.BestDatabaseSequenceMapping),
  interactions: PluginSpec.Behavior(PluginBehaviors.CustomProps.Interactions),
  secondaryStructure: PluginSpec.Behavior(PluginBehaviors.CustomProps.SecondaryStructure),
  valenceModel: PluginSpec.Behavior(PluginBehaviors.CustomProps.ValenceModel),
  crossLinks: PluginSpec.Behavior(PluginBehaviors.CustomProps.CrossLinkRestraint),
  streamlines: PluginSpec.Behavior(PluginBehaviors.CustomProps.Streamlines),
  molViewSpec: MolViewSpecBehavior,
} as const;

/** The base parts, a new object each call. */
export function prismViewerPluginParts(): PrismViewerPluginParts {
  return {
    registry: PrismViewerRegistry,
    behaviors: Object.values(PrismViewerBehaviors),
  };
}

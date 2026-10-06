/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginBehaviors } from '@molstar/plugin/behavior';
import { StructureFocusRepresentation } from '@molstar/plugin/behavior/dynamic/selection/structure-focus-representation';
import { DefaultActions, DefaultAnimations, DefaultRegistry } from '@molstar/plugin/default-registry';
import { PluginSpec } from '@molstar/plugin/spec';

export const DefaultPluginSpec = (): PluginSpec => ({
  registry: DefaultRegistry,
  // Kept until the spec fields are removed; the entries registered through `registry` already list the same
  // actions and animations, so these only count again.
  actions: DefaultActions.actions!.map((a) => PluginSpec.Action(a as Parameters<typeof PluginSpec.Action>[0])),
  animations: [...DefaultAnimations.animations!],
  behaviors: [
    PluginSpec.Behavior(PluginBehaviors.Representation.HighlightLoci),
    PluginSpec.Behavior(PluginBehaviors.Representation.SelectLoci),
    PluginSpec.Behavior(PluginBehaviors.Representation.DefaultLociLabelProvider),
    PluginSpec.Behavior(PluginBehaviors.Representation.FocusLoci),
    PluginSpec.Behavior(PluginBehaviors.Camera.FocusLoci),
    PluginSpec.Behavior(PluginBehaviors.Camera.CameraAxisHelper),
    PluginSpec.Behavior(PluginBehaviors.Camera.CameraControls),
    PluginSpec.Behavior(PluginBehaviors.State.SnapshotControls),
    PluginSpec.Behavior(StructureFocusRepresentation),

    PluginSpec.Behavior(PluginBehaviors.CustomProps.StructureInfo),
    PluginSpec.Behavior(PluginBehaviors.CustomProps.AccessibleSurfaceArea),
    PluginSpec.Behavior(PluginBehaviors.CustomProps.BestDatabaseSequenceMapping),
    PluginSpec.Behavior(PluginBehaviors.CustomProps.Interactions),
    PluginSpec.Behavior(PluginBehaviors.CustomProps.SecondaryStructure),
    PluginSpec.Behavior(PluginBehaviors.CustomProps.ValenceModel),
    PluginSpec.Behavior(PluginBehaviors.CustomProps.CrossLinkRestraint),
    PluginSpec.Behavior(PluginBehaviors.CustomProps.Streamlines),
  ],
});

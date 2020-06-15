/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import * as ReactDOM from 'react-dom';
import { StateActions } from '../mol-plugin-state/actions';
import { AnimateAssemblyUnwind, AnimateModelIndex, AnimateStateInterpolation, AnimateUnitsExplode } from '../mol-plugin-state/animation/built-in';
import { StateTransforms } from '../mol-plugin-state/transforms';
import { VolumeStreamingCustomControls } from '../mol-plugin-ui/custom/volume';
import { Plugin } from '../mol-plugin-ui/plugin';
import { PluginBehaviors } from './behavior';
import { StructureFocusRepresentation } from './behavior/dynamic/selection/structure-focus-representation';
import { BoxifyVolumeStreaming, CreateVolumeStreamingBehavior, InitVolumeStreaming } from './behavior/dynamic/volume-streaming/transformers';
import { PluginContext } from './context';
import { PluginSpec } from './spec';
import { AssignColorVolume } from '../mol-plugin-state/actions/volume';

export const DefaultPluginSpec: PluginSpec = {
    actions: [
        PluginSpec.Action(StateActions.Structure.DownloadStructure),
        PluginSpec.Action(StateActions.Structure.AddTrajectory),
        PluginSpec.Action(StateActions.Volume.DownloadDensity),
        PluginSpec.Action(StateActions.DataFormat.DownloadFile),
        PluginSpec.Action(StateActions.DataFormat.OpenFiles),
        PluginSpec.Action(StateActions.Structure.EnableModelCustomProps),
        PluginSpec.Action(StateActions.Structure.EnableStructureCustomProps),

        // Volume streaming
        PluginSpec.Action(InitVolumeStreaming),
        PluginSpec.Action(BoxifyVolumeStreaming),
        PluginSpec.Action(CreateVolumeStreamingBehavior),

        PluginSpec.Action(StateTransforms.Data.Download),
        PluginSpec.Action(StateTransforms.Data.ParseCif),
        PluginSpec.Action(StateTransforms.Data.ParseCcp4),
        PluginSpec.Action(StateTransforms.Data.ParseDsn6),

        PluginSpec.Action(StateTransforms.Model.TrajectoryFromMmCif),
        PluginSpec.Action(StateTransforms.Model.TrajectoryFromCifCore),
        PluginSpec.Action(StateTransforms.Model.TrajectoryFromPDB),
        PluginSpec.Action(StateTransforms.Model.TransformStructureConformation),
        PluginSpec.Action(StateTransforms.Model.StructureFromModel),
        PluginSpec.Action(StateTransforms.Model.StructureFromTrajectory),
        PluginSpec.Action(StateTransforms.Model.ModelFromTrajectory),
        PluginSpec.Action(StateTransforms.Model.StructureSelectionFromScript),
        PluginSpec.Action(StateTransforms.Representation.StructureRepresentation3D),
        PluginSpec.Action(StateTransforms.Representation.StructureSelectionsDistance3D),
        PluginSpec.Action(StateTransforms.Representation.StructureSelectionsAngle3D),
        PluginSpec.Action(StateTransforms.Representation.StructureSelectionsDihedral3D),
        PluginSpec.Action(StateTransforms.Representation.StructureSelectionsLabel3D),
        PluginSpec.Action(StateTransforms.Representation.StructureSelectionsOrientation3D),
        PluginSpec.Action(StateTransforms.Representation.ModelUnitcell3D),
        PluginSpec.Action(StateTransforms.Representation.ExplodeStructureRepresentation3D),
        PluginSpec.Action(StateTransforms.Representation.UnwindStructureAssemblyRepresentation3D),
        PluginSpec.Action(StateTransforms.Representation.OverpaintStructureRepresentation3DFromScript),
        PluginSpec.Action(StateTransforms.Representation.TransparencyStructureRepresentation3DFromScript),

        PluginSpec.Action(AssignColorVolume),
        PluginSpec.Action(StateTransforms.Volume.VolumeFromCcp4),
        PluginSpec.Action(StateTransforms.Volume.VolumeFromDsn6),
        PluginSpec.Action(StateTransforms.Volume.VolumeFromCube),
        PluginSpec.Action(StateTransforms.Volume.VolumeFromDx),
        PluginSpec.Action(StateTransforms.Representation.VolumeRepresentation3D),
    ],
    behaviors: [
        PluginSpec.Behavior(PluginBehaviors.Representation.HighlightLoci),
        PluginSpec.Behavior(PluginBehaviors.Representation.SelectLoci),
        PluginSpec.Behavior(PluginBehaviors.Representation.DefaultLociLabelProvider),
        PluginSpec.Behavior(PluginBehaviors.Representation.FocusLoci),
        PluginSpec.Behavior(PluginBehaviors.Camera.FocusLoci),
        PluginSpec.Behavior(StructureFocusRepresentation),

        PluginSpec.Behavior(PluginBehaviors.CustomProps.StructureInfo),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.AccessibleSurfaceArea),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.Interactions),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.SecondaryStructure),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.ValenceModel),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.CrossLinkRestraint),
    ],
    customParamEditors: [
        [CreateVolumeStreamingBehavior, VolumeStreamingCustomControls]
    ],
    animations: [
        AnimateModelIndex,
        AnimateAssemblyUnwind,
        AnimateUnitsExplode,
        AnimateStateInterpolation
    ]
};

export function createPlugin(target: HTMLElement, spec?: PluginSpec): PluginContext {
    const ctx = new PluginContext(spec || DefaultPluginSpec);
    ReactDOM.render(React.createElement(Plugin, { plugin: ctx }), target);
    return ctx;
}
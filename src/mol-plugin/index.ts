/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from './context';
import { Plugin } from '../mol-plugin-ui/plugin'
import * as React from 'react';
import * as ReactDOM from 'react-dom';
import { PluginSpec } from './spec';
import { StateTransforms } from './state/transforms';
import { PluginBehaviors } from './behavior';
import { AnimateModelIndex, AnimateAssemblyUnwind, AnimateUnitsExplode, AnimateStateInterpolation } from './state/animation/built-in';
import { StateActions } from './state/actions';
import { InitVolumeStreaming, BoxifyVolumeStreaming, CreateVolumeStreamingBehavior } from './behavior/dynamic/volume-streaming/transformers';
import { StructureRepresentationInteraction } from './behavior/dynamic/selection/structure-representation-interaction';
import { TransformStructureConformation } from './state/actions/structure';
import { VolumeStreamingCustomControls } from '../mol-plugin-ui/custom/volume';
import { PluginConfig } from './config';

export const DefaultPluginSpec: PluginSpec = {
    actions: [
        PluginSpec.Action(StateActions.Structure.DownloadStructure),
        PluginSpec.Action(StateActions.Structure.AddTrajectory),
        PluginSpec.Action(StateActions.Volume.DownloadDensity),
        PluginSpec.Action(StateActions.DataFormat.OpenFiles),
        PluginSpec.Action(StateActions.Structure.Create3DRepresentationPreset),
        PluginSpec.Action(StateActions.Structure.Remove3DRepresentationPreset),
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
        PluginSpec.Action(StateTransforms.Model.TrajectoryFromPDB),
        PluginSpec.Action(TransformStructureConformation),
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

        PluginSpec.Action(StateTransforms.Volume.VolumeFromCcp4),
        PluginSpec.Action(StateTransforms.Representation.VolumeRepresentation3D),

        PluginSpec.Action(StateActions.Structure.StructureFromSelection),
    ],
    behaviors: [
        PluginSpec.Behavior(PluginBehaviors.Representation.HighlightLoci),
        PluginSpec.Behavior(PluginBehaviors.Representation.SelectLoci),
        PluginSpec.Behavior(PluginBehaviors.Representation.DefaultLociLabelProvider),
        PluginSpec.Behavior(PluginBehaviors.Camera.FocusLoci),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.AccessibleSurfaceArea),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.Interactions),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.SecondaryStructure),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.ValenceModel),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.PDBeStructureQualityReport, { autoAttach: true, showTooltip: true }),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.RCSBAssemblySymmetry),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.RCSBValidationReport),
        PluginSpec.Behavior(StructureRepresentationInteraction)
    ],
    customParamEditors: [
        [CreateVolumeStreamingBehavior, VolumeStreamingCustomControls]
    ],
    animations: [
        AnimateModelIndex,
        AnimateAssemblyUnwind,
        AnimateUnitsExplode,
        AnimateStateInterpolation
    ],
    config: new Map([
        [PluginConfig.State.DefaultServer, 'https://webchem.ncbr.muni.cz/molstar-state']
    ])
}

export function createPlugin(target: HTMLElement, spec?: PluginSpec): PluginContext {
    const ctx = new PluginContext(spec || DefaultPluginSpec);
    ReactDOM.render(React.createElement(Plugin, { plugin: ctx }), target);
    return ctx;
}
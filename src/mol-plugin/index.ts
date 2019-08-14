/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from './context';
import { Plugin } from './ui/plugin'
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
import { VolumeStreamingCustomControls } from './ui/custom/volume';

export const DefaultPluginSpec: PluginSpec = {
    actions: [
        PluginSpec.Action(StateActions.Structure.DownloadStructure),
        PluginSpec.Action(StateActions.Volume.DownloadDensity),
        PluginSpec.Action(StateActions.DataFormat.OpenFile),
        PluginSpec.Action(StateActions.Structure.CreateComplexRepresentation),
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
        PluginSpec.Action(StateTransforms.Model.StructureAssemblyFromModel),
        PluginSpec.Action(StateTransforms.Model.StructureSymmetryFromModel),
        PluginSpec.Action(TransformStructureConformation),
        PluginSpec.Action(StateTransforms.Model.StructureFromModel),
        PluginSpec.Action(StateTransforms.Model.StructureFromTrajectory),
        PluginSpec.Action(StateTransforms.Model.ModelFromTrajectory),
        PluginSpec.Action(StateTransforms.Model.UserStructureSelection),
        PluginSpec.Action(StateTransforms.Representation.StructureRepresentation3D),
        PluginSpec.Action(StateTransforms.Representation.StructureLabels3D),
        PluginSpec.Action(StateTransforms.Representation.ExplodeStructureRepresentation3D),
        PluginSpec.Action(StateTransforms.Representation.UnwindStructureAssemblyRepresentation3D),
        PluginSpec.Action(StateTransforms.Representation.OverpaintStructureRepresentation3D),
        PluginSpec.Action(StateTransforms.Representation.TransparencyStructureRepresentation3D),

        PluginSpec.Action(StateTransforms.Volume.VolumeFromCcp4),
        PluginSpec.Action(StateTransforms.Representation.VolumeRepresentation3D),

        PluginSpec.Action(StateActions.Structure.StructureFromSelection),
    ],
    behaviors: [
        PluginSpec.Behavior(PluginBehaviors.Representation.HighlightLoci),
        PluginSpec.Behavior(PluginBehaviors.Representation.SelectLoci),
        PluginSpec.Behavior(PluginBehaviors.Representation.DefaultLociLabelProvider),
        PluginSpec.Behavior(PluginBehaviors.Camera.FocusLociOnSelect, { minRadius: 8, extraRadius: 4 }),
        // PluginSpec.Behavior(PluginBehaviors.Labels.SceneLabels),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.MolstarSecondaryStructure, { autoAttach: true }),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.PDBeStructureQualityReport, { autoAttach: true }),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.RCSBAssemblySymmetry, { autoAttach: true }),
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
    ]
}

export function createPlugin(target: HTMLElement, spec?: PluginSpec): PluginContext {
    const ctx = new PluginContext(spec || DefaultPluginSpec);
    ReactDOM.render(React.createElement(Plugin, { plugin: ctx }), target);
    return ctx;
}
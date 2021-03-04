/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateActions } from '../mol-plugin-state/actions';
import { AssignColorVolume } from '../mol-plugin-state/actions/volume';
import { StateTransforms } from '../mol-plugin-state/transforms';
import { StateTransformParameters } from '../mol-plugin-ui/state/common';
import { DefaultPluginSpec, PluginSpec } from '../mol-plugin/spec';
import { StateAction, StateTransformer } from '../mol-state';
import { BoxifyVolumeStreaming, CreateVolumeStreamingBehavior, InitVolumeStreaming } from '../mol-plugin/behavior/dynamic/volume-streaming/transformers';

export { PluginUISpec };

interface PluginUISpec extends PluginSpec {
    actions: PluginUISpec.Action[],
    customParamEditors?: [StateAction | StateTransformer, StateTransformParameters.Class][],
    components?: {
        controls?: PluginUISpec.LayoutControls
        remoteState?: 'none' | 'default',
        structureTools?: React.ComponentClass,
        viewport?: {
            view?: React.ComponentClass,
            controls?: React.ComponentClass
        },
        hideTaskOverlay?: boolean
    },
}

namespace PluginUISpec {
    export interface Action {
        action: StateAction | StateTransformer,
        customControl?: StateTransformParameters.Class,
        autoUpdate?: boolean
    }

    export function Action(action: StateAction | StateTransformer, params?: { customControl?: StateTransformParameters.Class, autoUpdate?: boolean }): Action {
        return { action, customControl: params && params.customControl, autoUpdate: params && params.autoUpdate };
    }

    export interface LayoutControls {
        top?: React.ComponentClass | 'none',
        left?: React.ComponentClass | 'none',
        right?: React.ComponentClass | 'none',
        bottom?: React.ComponentClass | 'none'
    }
}

export const DefaultPluginUISpec = (): PluginUISpec => ({
    ...DefaultPluginSpec(),
    actions: [
        PluginUISpec.Action(StateActions.Structure.DownloadStructure),
        PluginUISpec.Action(StateActions.Structure.AddTrajectory),
        PluginUISpec.Action(StateActions.Volume.DownloadDensity),
        PluginUISpec.Action(StateActions.DataFormat.DownloadFile),
        PluginUISpec.Action(StateActions.DataFormat.OpenFiles),
        PluginUISpec.Action(StateActions.Structure.EnableModelCustomProps),
        PluginUISpec.Action(StateActions.Structure.EnableStructureCustomProps),

        // Volume streaming
        PluginUISpec.Action(InitVolumeStreaming),
        PluginUISpec.Action(BoxifyVolumeStreaming),
        PluginUISpec.Action(CreateVolumeStreamingBehavior),

        PluginUISpec.Action(StateTransforms.Data.Download),
        PluginUISpec.Action(StateTransforms.Data.ParseCif),
        PluginUISpec.Action(StateTransforms.Data.ParseCcp4),
        PluginUISpec.Action(StateTransforms.Data.ParseDsn6),

        PluginUISpec.Action(StateTransforms.Model.TrajectoryFromMmCif),
        PluginUISpec.Action(StateTransforms.Model.TrajectoryFromCifCore),
        PluginUISpec.Action(StateTransforms.Model.TrajectoryFromPDB),
        PluginUISpec.Action(StateTransforms.Model.TransformStructureConformation),
        PluginUISpec.Action(StateTransforms.Model.StructureFromModel),
        PluginUISpec.Action(StateTransforms.Model.StructureFromTrajectory),
        PluginUISpec.Action(StateTransforms.Model.ModelFromTrajectory),
        PluginUISpec.Action(StateTransforms.Model.StructureSelectionFromScript),
        PluginUISpec.Action(StateTransforms.Representation.StructureRepresentation3D),
        PluginUISpec.Action(StateTransforms.Representation.StructureSelectionsDistance3D),
        PluginUISpec.Action(StateTransforms.Representation.StructureSelectionsAngle3D),
        PluginUISpec.Action(StateTransforms.Representation.StructureSelectionsDihedral3D),
        PluginUISpec.Action(StateTransforms.Representation.StructureSelectionsLabel3D),
        PluginUISpec.Action(StateTransforms.Representation.StructureSelectionsOrientation3D),
        PluginUISpec.Action(StateTransforms.Representation.ModelUnitcell3D),
        PluginUISpec.Action(StateTransforms.Representation.ExplodeStructureRepresentation3D),
        PluginUISpec.Action(StateTransforms.Representation.UnwindStructureAssemblyRepresentation3D),
        PluginUISpec.Action(StateTransforms.Representation.OverpaintStructureRepresentation3DFromScript),
        PluginUISpec.Action(StateTransforms.Representation.TransparencyStructureRepresentation3DFromScript),

        PluginUISpec.Action(AssignColorVolume),
        PluginUISpec.Action(StateTransforms.Volume.VolumeFromCcp4),
        PluginUISpec.Action(StateTransforms.Volume.VolumeFromDsn6),
        PluginUISpec.Action(StateTransforms.Volume.VolumeFromCube),
        PluginUISpec.Action(StateTransforms.Volume.VolumeFromDx),
        PluginUISpec.Action(StateTransforms.Representation.VolumeRepresentation3D),
    ]
});
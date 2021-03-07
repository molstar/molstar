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
    ]
});
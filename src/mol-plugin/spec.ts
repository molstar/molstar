/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PartialCanvas3DProps } from '../mol-canvas3d/canvas3d';
import { AnimateAssemblyUnwind } from '../mol-plugin-state/animation/built-in/assembly-unwind';
import { AnimateCameraSpin } from '../mol-plugin-state/animation/built-in/camera-spin';
import { AnimateModelIndex } from '../mol-plugin-state/animation/built-in/model-index';
import { AnimateStateSnapshots } from '../mol-plugin-state/animation/built-in/state-snapshots';
import { PluginStateAnimation } from '../mol-plugin-state/animation/model';
import { DataFormatProvider } from '../mol-plugin-state/formats/provider';
import { StateAction, StateTransformer } from '../mol-state';
import { PluginBehaviors } from './behavior';
import { StructureFocusRepresentation } from './behavior/dynamic/selection/structure-focus-representation';
import { PluginConfigItem } from './config';
import { PluginLayoutStateProps } from './layout';
import { StateActions } from '../mol-plugin-state/actions';
import { AssignColorVolume } from '../mol-plugin-state/actions/volume';
import { StateTransforms } from '../mol-plugin-state/transforms';
import { BoxifyVolumeStreaming, CreateVolumeStreamingBehavior, InitVolumeStreaming } from '../mol-plugin/behavior/dynamic/volume-streaming/transformers';
import { AnimateStateInterpolation } from '../mol-plugin-state/animation/built-in/state-interpolation';
import { AnimateStructureSpin } from '../mol-plugin-state/animation/built-in/spin-structure';
import { AnimateCameraRock } from '../mol-plugin-state/animation/built-in/camera-rock';

export { PluginSpec };

interface PluginSpec {
    actions?: PluginSpec.Action[],
    behaviors: PluginSpec.Behavior[],
    animations?: PluginStateAnimation[],
    customFormats?: [string, DataFormatProvider][],
    canvas3d?: PartialCanvas3DProps,
    layout?: {
        initial?: Partial<PluginLayoutStateProps>,
    },
    config?: [PluginConfigItem, unknown][]
}

namespace PluginSpec {
    export interface Action {
        action: StateAction | StateTransformer,
        /* constructible react component with <action.customControl /> */
        customControl?: any,
        autoUpdate?: boolean
    }

    export function Action(action: StateAction | StateTransformer, params?: { customControl?: any /* constructible react component with <action.customControl /> */, autoUpdate?: boolean }): Action {
        return { action, customControl: params && params.customControl, autoUpdate: params && params.autoUpdate };
    }

    export interface Behavior {
        transformer: StateTransformer,
        defaultParams?: any
    }

    export function Behavior<T extends StateTransformer>(transformer: T, defaultParams: Partial<StateTransformer.Params<T>> = {}): Behavior {
        return { transformer, defaultParams };
    }
}

export const DefaultPluginSpec = (): PluginSpec => ({
    actions: [
        PluginSpec.Action(StateActions.Structure.DownloadStructure),
        PluginSpec.Action(StateActions.Volume.DownloadDensity),
        PluginSpec.Action(StateActions.DataFormat.DownloadFile),
        PluginSpec.Action(StateActions.DataFormat.OpenFiles),
        PluginSpec.Action(StateActions.Structure.LoadTrajectory),
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
        PluginSpec.Action(StateTransforms.Representation.StructureBoundingBox3D),
        PluginSpec.Action(StateTransforms.Representation.ExplodeStructureRepresentation3D),
        PluginSpec.Action(StateTransforms.Representation.SpinStructureRepresentation3D),
        PluginSpec.Action(StateTransforms.Representation.UnwindStructureAssemblyRepresentation3D),
        PluginSpec.Action(StateTransforms.Representation.OverpaintStructureRepresentation3DFromScript),
        PluginSpec.Action(StateTransforms.Representation.TransparencyStructureRepresentation3DFromScript),
        PluginSpec.Action(StateTransforms.Representation.ClippingStructureRepresentation3DFromScript),
        PluginSpec.Action(StateTransforms.Representation.SubstanceStructureRepresentation3DFromScript),
        PluginSpec.Action(StateTransforms.Representation.ThemeStrengthRepresentation3D),

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
        PluginSpec.Behavior(PluginBehaviors.Camera.CameraAxisHelper),
        PluginSpec.Behavior(PluginBehaviors.Camera.CameraControls),
        PluginSpec.Behavior(StructureFocusRepresentation),

        PluginSpec.Behavior(PluginBehaviors.CustomProps.StructureInfo),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.AccessibleSurfaceArea),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.BestDatabaseSequenceMapping),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.Interactions),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.SecondaryStructure),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.ValenceModel),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.CrossLinkRestraint),
    ],
    animations: [
        AnimateModelIndex,
        AnimateCameraSpin,
        AnimateCameraRock,
        AnimateStateSnapshots,
        AnimateAssemblyUnwind,
        AnimateStructureSpin,
        AnimateStateInterpolation
    ]
});
/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateTransformer, StateAction } from '../mol-state';
import { StateTransformParameters } from '../mol-plugin-ui/state/common';
import { PluginLayoutStateProps } from './layout';
import { PluginStateAnimation } from '../mol-plugin-state/animation/model';
import { PluginConfigItem } from './config';
import { PartialCanvas3DProps } from '../mol-canvas3d/canvas3d';
import { DataFormatProvider } from '../mol-plugin-state/formats/provider';
import { StateActions } from '../mol-plugin-state/actions';
import { StateTransforms } from '../mol-plugin-state/transforms';
import { VolumeStreamingCustomControls } from '../mol-plugin-ui/custom/volume';
import { PluginBehaviors } from './behavior';
import { StructureFocusRepresentation } from './behavior/dynamic/selection/structure-focus-representation';
import { BoxifyVolumeStreaming, CreateVolumeStreamingBehavior, InitVolumeStreaming } from './behavior/dynamic/volume-streaming/transformers';
import { AssignColorVolume } from '../mol-plugin-state/actions/volume';
import { AnimateModelIndex } from '../mol-plugin-state/animation/built-in/model-index';
import { AnimateAssemblyUnwind } from '../mol-plugin-state/animation/built-in/assembly-unwind';
import { AnimateCameraSpin } from '../mol-plugin-state/animation/built-in/camera-spin';
import { AnimateStateSnapshots } from '../mol-plugin-state/animation/built-in/state-snapshots';

export { PluginSpec };

interface PluginSpec {
    actions: PluginSpec.Action[],
    behaviors: PluginSpec.Behavior[],
    animations?: PluginStateAnimation[],
    customParamEditors?: [StateAction | StateTransformer, StateTransformParameters.Class][],
    customFormats?: [string, DataFormatProvider][]
    layout?: {
        initial?: Partial<PluginLayoutStateProps>,
        controls?: PluginSpec.LayoutControls
    },
    components?: {
        remoteState?: 'none' | 'default',
        structureTools?: React.ComponentClass,
        viewport?: {
            view?: React.ComponentClass,
            controls?: React.ComponentClass,
            canvas3d?: PartialCanvas3DProps
        },
        hideTaskOverlay?: boolean
    },
    config?: [PluginConfigItem, unknown][]
}

namespace PluginSpec {
    export interface Action {
        action: StateAction | StateTransformer,
        customControl?: StateTransformParameters.Class,
        autoUpdate?: boolean
    }

    export function Action(action: StateAction | StateTransformer, params?: { customControl?: StateTransformParameters.Class, autoUpdate?: boolean }): Action {
        return { action, customControl: params && params.customControl, autoUpdate: params && params.autoUpdate };
    }

    export interface Behavior {
        transformer: StateTransformer,
        defaultParams?: any
    }

    export function Behavior<T extends StateTransformer>(transformer: T, defaultParams: Partial<StateTransformer.Params<T>> = {}): Behavior {
        return { transformer, defaultParams };
    }

    export interface LayoutControls {
        top?: React.ComponentClass | 'none',
        left?: React.ComponentClass | 'none',
        right?: React.ComponentClass | 'none',
        bottom?: React.ComponentClass | 'none'
    }
}

export const DefaultPluginSpec = (): PluginSpec => ({
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
        AnimateCameraSpin,
        AnimateStateSnapshots,
        AnimateAssemblyUnwind
    ]
});
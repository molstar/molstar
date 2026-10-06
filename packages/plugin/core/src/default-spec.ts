/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AnimateAssemblyUnwind } from '@molstar/plugin/state/animation/built-in/assembly-unwind';
import { AnimateCameraSpin } from '@molstar/plugin/state/animation/built-in/camera-spin';
import { AnimateModelIndex } from '@molstar/plugin/state/animation/built-in/model-index';
import { AnimateParticleTrajectory } from '@molstar/plugin/state/animation/built-in/particles';
import {
  AnimateStateSnapshotTransition,
  AnimateStateSnapshots,
} from '@molstar/plugin/state/animation/built-in/state-snapshots';
import { PluginBehaviors } from '@molstar/plugin/behavior';
import { StructureFocusRepresentation } from '@molstar/plugin/behavior/dynamic/selection/structure-focus-representation';
import { StateActions } from '@molstar/plugin/state/actions';
import { AssignColorVolume } from '@molstar/plugin/state/actions/volume';
import { StateTransforms } from '@molstar/plugin/state/transforms';
import {
  BoxifyVolumeStreaming,
  CreateVolumeStreamingBehavior,
  InitVolumeStreaming,
} from '@molstar/plugin/behavior/dynamic/volume-streaming/transformers';
import { AnimateStateInterpolation } from '@molstar/plugin/state/animation/built-in/state-interpolation';
import { AnimateStructureSpin } from '@molstar/plugin/state/animation/built-in/spin-structure';
import { AnimateCameraRock } from '@molstar/plugin/state/animation/built-in/camera-rock';
import { AnimateTime } from '@molstar/plugin/state/animation/built-in/time';
import { PluginSpec } from '@molstar/plugin/spec';

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
    PluginSpec.Action(StateTransforms.Model.StructureInstances),
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
    PluginSpec.Action(StateTransforms.Representation.WiggleStructureRepresentation3DFromScript),
    PluginSpec.Action(StateTransforms.Representation.ThemeStrengthRepresentation3D),

    PluginSpec.Action(AssignColorVolume),
    PluginSpec.Action(StateTransforms.Volume.VolumeFromCcp4),
    PluginSpec.Action(StateTransforms.Volume.VolumeFromDsn6),
    PluginSpec.Action(StateTransforms.Volume.VolumeFromCube),
    PluginSpec.Action(StateTransforms.Volume.VolumeFromDx),
    PluginSpec.Action(StateTransforms.Representation.VolumeRepresentation3D),
    PluginSpec.Action(StateTransforms.Volume.VolumeTransform),
    PluginSpec.Action(StateTransforms.Volume.VolumeInstances),

    PluginSpec.Action(StateTransforms.Particles.ParticleListFromRelionStar),
    PluginSpec.Action(StateTransforms.Particles.ParticleListFromDynamoTbl),
    PluginSpec.Action(StateTransforms.Particles.ParticleListFromCryoEtDataPortalNdjson),
    PluginSpec.Action(StateTransforms.Particles.ParticleListFromArtiatomiEm),
    PluginSpec.Action(StateTransforms.Particles.ParticleListFromMmcifAssembly),
    PluginSpec.Action(StateTransforms.Particles.ParticleTrajectoryFromSimularium),
    PluginSpec.Action(StateTransforms.Particles.ParticleListFromTrajectory),
    PluginSpec.Action(StateTransforms.Particles.ParticleListWithTargets),
    PluginSpec.Action(StateTransforms.Particles.ParticleListUnitcell3D),
    PluginSpec.Action(StateTransforms.Particles.ParticlesRepresentation3D),
  ],
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
  animations: [
    AnimateModelIndex,
    AnimateParticleTrajectory,
    AnimateCameraSpin,
    AnimateCameraRock,
    AnimateStateSnapshots,
    AnimateStateSnapshotTransition,
    AnimateAssemblyUnwind,
    AnimateStructureSpin,
    AnimateStateInterpolation,
    AnimateTime,
  ],
});

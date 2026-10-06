/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import '@molstar/model/script/transpilers/all';
import '@molstar/plugin/state/transforms/catalog';
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
import { Download } from '@molstar/plugin/state/transforms/data/fetch';
import { ParseCif } from '@molstar/plugin/state/formats/cif';
import { ParseCcp4, VolumeFromCcp4 } from '@molstar/plugin/state/formats/volume/ccp4';
import { ParseDsn6, VolumeFromDsn6 } from '@molstar/plugin/state/formats/volume/dsn6';
import { TrajectoryFromMmCif } from '@molstar/plugin/state/formats/trajectory/mmcif';
import { TrajectoryFromCifCore } from '@molstar/plugin/state/formats/trajectory/cif-core';
import { TrajectoryFromPDB } from '@molstar/plugin/state/formats/trajectory/pdb';
import {
  ModelFromTrajectory,
  StructureFromModel,
  StructureFromTrajectory,
  StructureInstances,
  TransformStructureConformation,
} from '@molstar/plugin/state/transforms/structure/hierarchy';
import { StructureSelectionFromScript } from '@molstar/plugin/state/transforms/structure/selection';
import { StructureRepresentation3D } from '@molstar/plugin/state/transforms/structure/representation';
import {
  StructureSelectionsAngle3D,
  StructureSelectionsDihedral3D,
  StructureSelectionsDistance3D,
  StructureSelectionsLabel3D,
  StructureSelectionsOrientation3D,
} from '@molstar/plugin/state/transforms/structure/measurement';
import { ModelUnitcell3D } from '@molstar/plugin/state/transforms/structure/unitcell';
import { StructureBoundingBox3D } from '@molstar/plugin/state/transforms/structure/bounding-box';
import {
  ExplodeStructureRepresentation3D,
  SpinStructureRepresentation3D,
  UnwindStructureAssemblyRepresentation3D,
} from '@molstar/plugin/state/transforms/structure/animation';
import { OverpaintStructureRepresentation3DFromScript } from '@molstar/plugin/state/transforms/structure/effects/overpaint';
import { TransparencyStructureRepresentation3DFromScript } from '@molstar/plugin/state/transforms/structure/effects/transparency';
import { ClippingStructureRepresentation3DFromScript } from '@molstar/plugin/state/transforms/structure/effects/clipping';
import { SubstanceStructureRepresentation3DFromScript } from '@molstar/plugin/state/transforms/structure/effects/substance';
import { WiggleStructureRepresentation3DFromScript } from '@molstar/plugin/state/transforms/structure/effects/wiggle';
import { ThemeStrengthRepresentation3D } from '@molstar/plugin/state/transforms/structure/effects/theme-strength';
import { VolumeFromCube } from '@molstar/plugin/state/formats/volume/cube';
import { VolumeFromDx } from '@molstar/plugin/state/formats/volume/dx';
import { VolumeRepresentation3D } from '@molstar/plugin/state/transforms/volume/representation';
import { VolumeInstances, VolumeTransform } from '@molstar/plugin/state/transforms/volume/ops';
import { ParticleListFromRelionStar } from '@molstar/plugin/state/formats/particles/star';
import { ParticleListFromDynamoTbl } from '@molstar/plugin/state/formats/particles/tbl';
import { ParticleListFromCryoEtDataPortalNdjson } from '@molstar/plugin/state/formats/particles/ndjson';
import { ParticleListFromArtiatomiEm } from '@molstar/plugin/state/formats/particles/em';
import { ParticleListFromMmcifAssembly } from '@molstar/plugin/state/formats/particles/mmcif-assembly';
import { ParticleTrajectoryFromSimularium } from '@molstar/plugin/state/formats/particles/simularium';
import { ParticleListFromTrajectory, ParticleListWithTargets } from '@molstar/plugin/state/transforms/particles/ops';
import { ParticleListUnitcell3D } from '@molstar/plugin/state/transforms/particles/unitcell';
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
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

    PluginSpec.Action(Download),
    PluginSpec.Action(ParseCif),
    PluginSpec.Action(ParseCcp4),
    PluginSpec.Action(ParseDsn6),

    PluginSpec.Action(TrajectoryFromMmCif),
    PluginSpec.Action(TrajectoryFromCifCore),
    PluginSpec.Action(TrajectoryFromPDB),
    PluginSpec.Action(TransformStructureConformation),
    PluginSpec.Action(StructureInstances),
    PluginSpec.Action(StructureFromModel),
    PluginSpec.Action(StructureFromTrajectory),
    PluginSpec.Action(ModelFromTrajectory),
    PluginSpec.Action(StructureSelectionFromScript),
    PluginSpec.Action(StructureRepresentation3D),
    PluginSpec.Action(StructureSelectionsDistance3D),
    PluginSpec.Action(StructureSelectionsAngle3D),
    PluginSpec.Action(StructureSelectionsDihedral3D),
    PluginSpec.Action(StructureSelectionsLabel3D),
    PluginSpec.Action(StructureSelectionsOrientation3D),
    PluginSpec.Action(ModelUnitcell3D),
    PluginSpec.Action(StructureBoundingBox3D),
    PluginSpec.Action(ExplodeStructureRepresentation3D),
    PluginSpec.Action(SpinStructureRepresentation3D),
    PluginSpec.Action(UnwindStructureAssemblyRepresentation3D),
    PluginSpec.Action(OverpaintStructureRepresentation3DFromScript),
    PluginSpec.Action(TransparencyStructureRepresentation3DFromScript),
    PluginSpec.Action(ClippingStructureRepresentation3DFromScript),
    PluginSpec.Action(SubstanceStructureRepresentation3DFromScript),
    PluginSpec.Action(WiggleStructureRepresentation3DFromScript),
    PluginSpec.Action(ThemeStrengthRepresentation3D),

    PluginSpec.Action(AssignColorVolume),
    PluginSpec.Action(VolumeFromCcp4),
    PluginSpec.Action(VolumeFromDsn6),
    PluginSpec.Action(VolumeFromCube),
    PluginSpec.Action(VolumeFromDx),
    PluginSpec.Action(VolumeRepresentation3D),
    PluginSpec.Action(VolumeTransform),
    PluginSpec.Action(VolumeInstances),

    PluginSpec.Action(ParticleListFromRelionStar),
    PluginSpec.Action(ParticleListFromDynamoTbl),
    PluginSpec.Action(ParticleListFromCryoEtDataPortalNdjson),
    PluginSpec.Action(ParticleListFromArtiatomiEm),
    PluginSpec.Action(ParticleListFromMmcifAssembly),
    PluginSpec.Action(ParticleTrajectoryFromSimularium),
    PluginSpec.Action(ParticleListFromTrajectory),
    PluginSpec.Action(ParticleListWithTargets),
    PluginSpec.Action(ParticleListUnitcell3D),
    PluginSpec.Action(ParticlesRepresentation3D),
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

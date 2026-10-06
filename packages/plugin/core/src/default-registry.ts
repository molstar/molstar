/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import '@molstar/model/script/transpilers/all';
import '@molstar/plugin/state/transforms/catalog';
import type { ColorTheme } from '@molstar/graphics/theme/color';
import type { SizeTheme } from '@molstar/graphics/theme/size';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { BuiltInColorThemes } from '@molstar/graphics/theme/color/catalog';
import { BuiltInSizeThemes } from '@molstar/graphics/theme/size/catalog';
import { BuiltInStructureRepresentations } from '@molstar/graphics/repr/structure/catalog';
import { BuiltInVolumeRepresentations } from '@molstar/graphics/repr/volume/catalog';
import { BuiltInParticleRepresentations } from '@molstar/graphics/repr/particles/catalog';
import { BuiltInVolumeFormats } from '@molstar/plugin/state/formats/volume/catalog';
import { BuiltInTopologyFormats } from '@molstar/plugin/state/formats/topology/catalog';
import { BuiltInCoordinatesFormats } from '@molstar/plugin/state/formats/coordinates/catalog';
import { BuiltInShapeFormats } from '@molstar/plugin/state/formats/shape/catalog';
import { BuiltInParticlesFormats } from '@molstar/plugin/state/formats/particles/catalog';
import { BuiltInTrajectoryFormats } from '@molstar/plugin/state/formats/trajectory/catalog';
import { PresetTrajectoryHierarchy } from '@molstar/plugin/state/builder/structure/hierarchy-presets/catalog';
import { PresetStructureRepresentations } from '@molstar/plugin/state/builder/structure/representation-presets/catalog';
import { StructureSelectionQueries } from '@molstar/plugin/state/queries/structure/catalog';
import {
  AminoAcidSelectionQueries,
  NucleicBaseSelectionQueries,
} from '@molstar/plugin/state/queries/structure/residue';
import { BuiltInMarkdownExtension } from '@molstar/plugin/state/markdown/catalog';
import { ExternalColorThemes } from '@molstar/plugin/themes/external';
import { AnimateAssemblyUnwind } from '@molstar/plugin/state/animation/built-in/assembly-unwind';
import { AnimateCameraSpin } from '@molstar/plugin/state/animation/built-in/camera-spin';
import { AnimateModelIndex } from '@molstar/plugin/state/animation/built-in/model-index';
import { AnimateParticleTrajectory } from '@molstar/plugin/state/animation/built-in/particles';
import {
  AnimateStateSnapshotTransition,
  AnimateStateSnapshots,
} from '@molstar/plugin/state/animation/built-in/state-snapshots';
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

/*
 * The default plugin as named entries, in the order of `DefaultRegistry`. Applications compose subsets by identity,
 * for example `DefaultRegistry.filter((e) => e !== DefaultAnimations)`.
 */

const ColorThemes = [
  ...Object.values(BuiltInColorThemes),
  ...(ExternalColorThemes.structure?.themes?.color ?? []),
] as readonly ColorTheme.Provider<any, any>[];
const SizeThemes = Object.values(BuiltInSizeThemes) as readonly SizeTheme.Provider<any, any>[];

/** Every built-in color and size theme in all three scopes, plus the external structure and volume color themes. */
export const DefaultThemes: PluginRegistryEntry = {
  structure: { themes: { color: ColorThemes, size: SizeThemes } },
  volume: { themes: { color: ColorThemes, size: SizeThemes } },
  particles: { themes: { color: ColorThemes, size: SizeThemes } },
};

export const DefaultStructureRepresentations: PluginRegistryEntry = {
  structure: { representations: Object.values(BuiltInStructureRepresentations) },
};

export const DefaultVolumeRepresentations: PluginRegistryEntry = {
  volume: { representations: Object.values(BuiltInVolumeRepresentations) },
};

export const DefaultParticleRepresentations: PluginRegistryEntry = {
  particles: { representations: Object.values(BuiltInParticleRepresentations) },
};

export const DefaultActions: PluginRegistryEntry = {
  actions: [
    StateActions.Structure.DownloadStructure,
    StateActions.Volume.DownloadDensity,
    StateActions.DataFormat.DownloadFile,
    StateActions.DataFormat.OpenFiles,
    StateActions.Structure.LoadTrajectory,
    StateActions.Structure.EnableModelCustomProps,
    StateActions.Structure.EnableStructureCustomProps,

    // Volume streaming
    InitVolumeStreaming,
    BoxifyVolumeStreaming,
    CreateVolumeStreamingBehavior,

    Download,
    ParseCif,
    ParseCcp4,
    ParseDsn6,

    TrajectoryFromMmCif,
    TrajectoryFromCifCore,
    TrajectoryFromPDB,
    TransformStructureConformation,
    StructureInstances,
    StructureFromModel,
    StructureFromTrajectory,
    ModelFromTrajectory,
    StructureSelectionFromScript,
    StructureRepresentation3D,
    StructureSelectionsDistance3D,
    StructureSelectionsAngle3D,
    StructureSelectionsDihedral3D,
    StructureSelectionsLabel3D,
    StructureSelectionsOrientation3D,
    ModelUnitcell3D,
    StructureBoundingBox3D,
    ExplodeStructureRepresentation3D,
    SpinStructureRepresentation3D,
    UnwindStructureAssemblyRepresentation3D,
    OverpaintStructureRepresentation3DFromScript,
    TransparencyStructureRepresentation3DFromScript,
    ClippingStructureRepresentation3DFromScript,
    SubstanceStructureRepresentation3DFromScript,
    WiggleStructureRepresentation3DFromScript,
    ThemeStrengthRepresentation3D,

    AssignColorVolume,
    VolumeFromCcp4,
    VolumeFromDsn6,
    VolumeFromCube,
    VolumeFromDx,
    VolumeRepresentation3D,
    VolumeTransform,
    VolumeInstances,

    ParticleListFromRelionStar,
    ParticleListFromDynamoTbl,
    ParticleListFromCryoEtDataPortalNdjson,
    ParticleListFromArtiatomiEm,
    ParticleListFromMmcifAssembly,
    ParticleTrajectoryFromSimularium,
    ParticleListFromTrajectory,
    ParticleListWithTargets,
    ParticleListUnitcell3D,
    ParticlesRepresentation3D,
  ],
};

export const DefaultFormats: PluginRegistryEntry = {
  formats: [
    ...BuiltInVolumeFormats,
    ...BuiltInTopologyFormats,
    ...BuiltInCoordinatesFormats,
    ...BuiltInShapeFormats,
    ...BuiltInParticlesFormats,
    ...BuiltInTrajectoryFormats,
  ],
};

export const DefaultPresets: PluginRegistryEntry = {
  structure: {
    presets: {
      hierarchy: Object.values(PresetTrajectoryHierarchy),
      representation: Object.values(PresetStructureRepresentations),
    },
  },
};

export const DefaultSelectionQueries: PluginRegistryEntry = {
  structure: {
    selectionQueries: [
      ...Object.values(StructureSelectionQueries),
      ...AminoAcidSelectionQueries,
      ...NucleicBaseSelectionQueries,
    ],
  },
};

export const DefaultMarkdownExtensions: PluginRegistryEntry = {
  markdownExtensions: BuiltInMarkdownExtension,
};

/**
 * Empty for now: the open-anything fallback handler moves here when `DragAndDropManager` stops registering it itself
 * (plugin composition step 3).
 */
export const DefaultDragAndDrop: PluginRegistryEntry = {
  dragAndDrop: [],
};

export const DefaultAnimations: PluginRegistryEntry = {
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
};

/** The default plugin: every entry above, in registration order. */
export const DefaultRegistry: readonly PluginRegistryEntry[] = [
  DefaultThemes,
  DefaultStructureRepresentations,
  DefaultVolumeRepresentations,
  DefaultParticleRepresentations,
  DefaultActions, // before DefaultFormats: format entries carry actions, which keep the 5.x registration order
  DefaultFormats,
  DefaultPresets,
  DefaultSelectionQueries,
  DefaultMarkdownExtensions,
  DefaultDragAndDrop,
  DefaultAnimations,
];

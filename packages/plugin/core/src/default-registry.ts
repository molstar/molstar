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
import { Cartoon } from '@molstar/plugin/registry/structure/cartoon';
import { Backbone } from '@molstar/plugin/registry/structure/backbone';
import { BallAndStick } from '@molstar/plugin/registry/structure/ball-and-stick';
import { BlobSurface } from '@molstar/plugin/registry/structure/blob-surface';
import { Carbohydrate } from '@molstar/plugin/registry/structure/carbohydrate';
import { Ellipsoid } from '@molstar/plugin/registry/structure/ellipsoid';
import { GaussianSurface } from '@molstar/plugin/registry/structure/gaussian-surface';
import { GaussianVolume } from '@molstar/plugin/registry/structure/gaussian-volume';
import { Label } from '@molstar/plugin/registry/structure/label';
import { Line } from '@molstar/plugin/registry/structure/line';
import { MolecularSurface } from '@molstar/plugin/registry/structure/molecular-surface';
import { Orientation } from '@molstar/plugin/registry/structure/orientation';
import { Plane } from '@molstar/plugin/registry/structure/plane';
import { Point } from '@molstar/plugin/registry/structure/point';
import { Putty } from '@molstar/plugin/registry/structure/putty';
import { Spacefill } from '@molstar/plugin/registry/structure/spacefill';
import { Polyhedron } from '@molstar/plugin/registry/structure/polyhedron';
import { DirectVolume } from '@molstar/plugin/registry/volume/direct-volume';
import { Dot } from '@molstar/plugin/registry/volume/dot';
import { Isosurface } from '@molstar/plugin/registry/volume/isosurface';
import { Segment } from '@molstar/plugin/registry/volume/segment';
import { Slice } from '@molstar/plugin/registry/volume/slice';
import { ParticleSpacefill } from '@molstar/plugin/registry/particles/spacefill';
import { ParticleOrientation } from '@molstar/plugin/registry/particles/orientation';
import { ParticleFibers } from '@molstar/plugin/registry/particles/fibers';
import { ParticleTarget } from '@molstar/plugin/registry/particles/target';
import { BuiltInVolumeFormats } from '@molstar/plugin/state/formats/volume/catalog';
import { BuiltInTopologyFormats } from '@molstar/plugin/state/formats/topology/catalog';
import { BuiltInCoordinatesFormats } from '@molstar/plugin/state/formats/coordinates/catalog';
import { BuiltInShapeFormats } from '@molstar/plugin/state/formats/shape/catalog';
import { BuiltInParticlesFormats } from '@molstar/plugin/state/formats/particles/catalog';
import { BuiltInTrajectoryFormats } from '@molstar/plugin/state/formats/trajectory/catalog';
import { DefaultHierarchyPresetEntry } from '@molstar/plugin/state/builder/structure/hierarchy-presets/default';
import { AllModelsHierarchyPresetEntry } from '@molstar/plugin/state/builder/structure/hierarchy-presets/all-models';
import { UnitcellHierarchyPresetEntry } from '@molstar/plugin/state/builder/structure/hierarchy-presets/unitcell';
import { SupercellHierarchyPresetEntry } from '@molstar/plugin/state/builder/structure/hierarchy-presets/supercell';
import { CrystalContactsHierarchyPresetEntry } from '@molstar/plugin/state/builder/structure/hierarchy-presets/crystal-contacts';
import { EmptyPresetEntry } from '@molstar/plugin/state/builder/structure/representation-presets/empty';
import { AutoPresetEntry } from '@molstar/plugin/state/builder/structure/representation-presets/auto';
import { AtomicDetailPresetEntry } from '@molstar/plugin/state/builder/structure/representation-presets/atomic-detail';
import { PolymerCartoonPresetEntry } from '@molstar/plugin/state/builder/structure/representation-presets/polymer-cartoon';
import { PolymerAndLigandPresetEntry } from '@molstar/plugin/state/builder/structure/representation-presets/polymer-and-ligand';
import { ProteinAndNucleicPresetEntry } from '@molstar/plugin/state/builder/structure/representation-presets/protein-and-nucleic';
import { CoarseSurfacePresetEntry } from '@molstar/plugin/state/builder/structure/representation-presets/coarse-surface';
import { IllustrativePresetEntry } from '@molstar/plugin/state/builder/structure/representation-presets/illustrative';
import { MolecularSurfacePresetEntry } from '@molstar/plugin/state/builder/structure/representation-presets/molecular-surface';
import { AutoLodPresetEntry } from '@molstar/plugin/state/builder/structure/representation-presets/auto-lod';
import { MesoscalePresetEntry } from '@molstar/plugin/state/builder/structure/representation-presets/mesoscale';
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
import { openDroppedFiles } from '@molstar/plugin/state/actions/file';
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

/*
 * The representation entries each carry the representation and its default themes; the default entries below list the
 * representation providers alone, in the order of the 5.x built-in catalogs. The themes come in with `DefaultThemes`.
 */

export const DefaultStructureRepresentations: PluginRegistryEntry = {
  structure: {
    representations: [
      Cartoon,
      Backbone,
      BallAndStick,
      BlobSurface,
      Carbohydrate,
      Ellipsoid,
      GaussianSurface,
      GaussianVolume,
      Label,
      Line,
      MolecularSurface,
      Orientation,
      Plane,
      Point,
      Putty,
      Spacefill,
      Polyhedron,
    ].flatMap((e) => e.structure!.representations!),
  },
};

export const DefaultVolumeRepresentations: PluginRegistryEntry = {
  volume: {
    representations: [DirectVolume, Dot, Isosurface, Segment, Slice].flatMap((e) => e.volume!.representations!),
  },
};

export const DefaultParticleRepresentations: PluginRegistryEntry = {
  particles: {
    representations: [ParticleSpacefill, ParticleOrientation, ParticleFibers, ParticleTarget].flatMap(
      (e) => e.particles!.representations!,
    ),
  },
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

// A preset entry lists its own preset first, followed by the presets it delegates to.
const ownPreset = <P>(presets: readonly P[] | undefined) => presets![0];

/*
 * The hierarchy and representation presets of 5.x in their 5.x order. Each preset entry also carries the
 * representations and themes the preset builds; the default entries above register those.
 */
export const DefaultPresets: PluginRegistryEntry = {
  structure: {
    presets: {
      hierarchy: [
        DefaultHierarchyPresetEntry,
        AllModelsHierarchyPresetEntry,
        UnitcellHierarchyPresetEntry,
        SupercellHierarchyPresetEntry,
        CrystalContactsHierarchyPresetEntry,
      ].map((e) => ownPreset(e.structure!.presets!.hierarchy)),
      representation: [
        EmptyPresetEntry,
        AutoPresetEntry,
        AtomicDetailPresetEntry,
        PolymerCartoonPresetEntry,
        PolymerAndLigandPresetEntry,
        ProteinAndNucleicPresetEntry,
        CoarseSurfacePresetEntry,
        IllustrativePresetEntry,
        MolecularSurfacePresetEntry,
        AutoLodPresetEntry,
        MesoscalePresetEntry,
      ].map((e) => ownPreset(e.structure!.presets!.representation)),
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

/** The open-anything handler: runs after session handling and every other handler, whatever the entry order. */
export const DefaultDragAndDrop: PluginRegistryEntry = {
  dragAndDrop: [{ name: 'open-files', handle: openDroppedFiles, fallback: true }],
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

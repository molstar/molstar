/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParseBlob, ParseCif } from '@molstar/plugin/state/formats/cif';
import { CoordinatesFromDcd } from '@molstar/plugin/state/formats/coordinates/dcd';
import { CoordinatesFromLammpstraj } from '@molstar/plugin/state/formats/coordinates/lammps';
import { CoordinatesFromNctraj } from '@molstar/plugin/state/formats/coordinates/nctraj';
import { CoordinatesFromTrr } from '@molstar/plugin/state/formats/coordinates/trr';
import { CoordinatesFromXtc } from '@molstar/plugin/state/formats/coordinates/xtc';
import { ParseArtiatomiEm, ParticleListFromArtiatomiEm } from '@molstar/plugin/state/formats/particles/em';
import { ParticleListFromMmcifAssembly } from '@molstar/plugin/state/formats/particles/mmcif-assembly';
import {
  ParseCryoEtDataPortalNdjson,
  ParticleListFromCryoEtDataPortalNdjson,
} from '@molstar/plugin/state/formats/particles/ndjson';
import { ParseSimularium, ParticleTrajectoryFromSimularium } from '@molstar/plugin/state/formats/particles/simularium';
import { ParticleListFromRelionStar } from '@molstar/plugin/state/formats/particles/star';
import { ParseDynamoTbl, ParticleListFromDynamoTbl } from '@molstar/plugin/state/formats/particles/tbl';
import { ParseObj, ShapeFromObj } from '@molstar/plugin/state/formats/shape/obj';
import { ParsePly, ShapeFromPly } from '@molstar/plugin/state/formats/shape/ply';
import { ParseVtp, ShapeFromVtp } from '@molstar/plugin/state/formats/shape/vtp';
import { ParsePrmtop, TopologyFromPrmtop } from '@molstar/plugin/state/formats/topology/prmtop';
import { ParsePsf, TopologyFromPsf } from '@molstar/plugin/state/formats/topology/psf';
import { ParseTop, TopologyFromTop } from '@molstar/plugin/state/formats/topology/top';
import { TrajectoryFromCifCore } from '@molstar/plugin/state/formats/trajectory/cif-core';
import { TrajectoryFromGRO } from '@molstar/plugin/state/formats/trajectory/gro';
import {
  TrajectoryFromLammpsData,
  TrajectoryFromLammpsTrajData,
} from '@molstar/plugin/state/formats/trajectory/lammps';
import { TrajectoryFromBlob, TrajectoryFromMmCif } from '@molstar/plugin/state/formats/trajectory/mmcif';
import { TrajectoryFromMOL } from '@molstar/plugin/state/formats/trajectory/mol';
import { TrajectoryFromMOL2 } from '@molstar/plugin/state/formats/trajectory/mol2';
import { TrajectoryFromPDB } from '@molstar/plugin/state/formats/trajectory/pdb';
import { TrajectoryFromSDF } from '@molstar/plugin/state/formats/trajectory/sdf';
import { TrajectoryFromXYZ } from '@molstar/plugin/state/formats/trajectory/xyz';
import { ParseCcp4, VolumeFromCcp4 } from '@molstar/plugin/state/formats/volume/ccp4';
import { ParseCube, TrajectoryFromCube, VolumeFromCube } from '@molstar/plugin/state/formats/volume/cube';
import { VolumeFromDensityServerCif } from '@molstar/plugin/state/formats/volume/density-server';
import { ParseDsn6, VolumeFromDsn6 } from '@molstar/plugin/state/formats/volume/dsn6';
import { ParseDx, VolumeFromDx } from '@molstar/plugin/state/formats/volume/dx';
import { ParseMtz, VolumeFromMtz } from '@molstar/plugin/state/formats/volume/mtz';
import { VolumeFromSegmentationCif } from '@molstar/plugin/state/formats/volume/segmentation';
import { VolumeFromStructureFactorsCif } from '@molstar/plugin/state/formats/volume/structure-factors';
import {
  DeflateData,
  Download,
  DownloadBlob,
  ImportString,
  LazyVolume,
  RawData,
  ReadFile,
} from '@molstar/plugin/state/transforms/data/fetch';
import { ImportJson, ParseJson } from '@molstar/plugin/state/transforms/data/json';
import { CreateGroup } from '@molstar/plugin/state/transforms/misc/group';
import { ParticleListFromTrajectory, ParticleListWithTargets } from '@molstar/plugin/state/transforms/particles/ops';
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
import { ParticleListUnitcell3D } from '@molstar/plugin/state/transforms/particles/unitcell';
import { BoxShape3D, getBoxMesh } from '@molstar/plugin/state/transforms/shape/box';
import { ShapeRepresentation3D } from '@molstar/plugin/state/transforms/shape/representation';
import {
  ExplodeStructureRepresentation3D,
  SpinStructureRepresentation3D,
  UnwindStructureAssemblyRepresentation3D,
} from '@molstar/plugin/state/transforms/structure/animation';
import { StructureBoundingBox3D } from '@molstar/plugin/state/transforms/structure/bounding-box';
import {
  ClippingStructureRepresentation3DFromBundle,
  ClippingStructureRepresentation3DFromScript,
} from '@molstar/plugin/state/transforms/structure/effects/clipping';
import {
  EmissiveStructureRepresentation3DFromBundle,
  EmissiveStructureRepresentation3DFromScript,
} from '@molstar/plugin/state/transforms/structure/effects/emissive';
import {
  OverpaintStructureRepresentation3DFromBundle,
  OverpaintStructureRepresentation3DFromScript,
} from '@molstar/plugin/state/transforms/structure/effects/overpaint';
import {
  SubstanceStructureRepresentation3DFromBundle,
  SubstanceStructureRepresentation3DFromScript,
} from '@molstar/plugin/state/transforms/structure/effects/substance';
import { ThemeStrengthRepresentation3D } from '@molstar/plugin/state/transforms/structure/effects/theme-strength';
import {
  TransparencyStructureRepresentation3DFromBundle,
  TransparencyStructureRepresentation3DFromScript,
} from '@molstar/plugin/state/transforms/structure/effects/transparency';
import {
  WiggleStructureRepresentation3DFromBundle,
  WiggleStructureRepresentation3DFromScript,
} from '@molstar/plugin/state/transforms/structure/effects/wiggle';
import {
  CustomModelProperties,
  CustomStructureProperties,
  ModelFromTrajectory,
  ModelWithCoordinates,
  StructureFromModel,
  StructureFromTrajectory,
  StructureInstances,
  TrajectoryFromModelAndCoordinates,
  TransformStructureConformation,
} from '@molstar/plugin/state/transforms/structure/hierarchy';
import {
  StructureSelectionsAngle3D,
  StructureSelectionsDihedral3D,
  StructureSelectionsDistance3D,
  StructureSelectionsLabel3D,
  StructureSelectionsOrientation3D,
  StructureSelectionsPlane3D,
} from '@molstar/plugin/state/transforms/structure/measurement';
import { StructureRepresentation3D } from '@molstar/plugin/state/transforms/structure/representation';
import {
  MultiStructureSelectionFromBundle,
  MultiStructureSelectionFromExpression,
  StructureComplexElement,
  StructureComplexElementTypes,
  StructureComponent,
  StructureSelectionFromBundle,
  StructureSelectionFromExpression,
  StructureSelectionFromScript,
} from '@molstar/plugin/state/transforms/structure/selection';
import { getTrajectory } from '@molstar/plugin/state/transforms/structure/trajectory-helpers';
import { ModelUnitcell3D } from '@molstar/plugin/state/transforms/structure/unitcell';
import {
  AssignColorVolume,
  CustomVolumeProperties,
  VolumeInstances,
  VolumeTransform,
} from '@molstar/plugin/state/transforms/volume/ops';
import { VolumeRepresentation3D } from '@molstar/plugin/state/transforms/volume/representation';
import { VolumeRepresentation3DHelpers } from '@molstar/plugin/state/transforms/volume/representation-helpers';

// The `molstar.lib.plugin.StateTransforms` global exists for classic-script compatibility: pages that load the
// Viewer bundle with a plain `<script>` tag cannot import leaf modules. Library code imports each transformer from
// its defining module instead of an aggregate object like this one.
export const StateTransforms = {
  Data: {
    Download,
    DownloadBlob,
    DeflateData,
    RawData,
    ReadFile,
    ParseBlob,
    ParseCif,
    ParseCube,
    ParsePsf,
    ParsePrmtop,
    ParseTop,
    ParsePly,
    ParseObj,
    ParseVtp,
    ParseCcp4,
    ParseDsn6,
    ParseMtz,
    ParseDx,
    ParseDynamoTbl,
    ParseCryoEtDataPortalNdjson,
    ParseArtiatomiEm,
    ParseSimularium,
    ImportString,
    ImportJson,
    ParseJson,
    LazyVolume,
  },
  Misc: {
    CreateGroup,
  },
  Model: {
    CoordinatesFromDcd,
    CoordinatesFromXtc,
    CoordinatesFromTrr,
    CoordinatesFromNctraj,
    CoordinatesFromLammpstraj,
    TopologyFromPsf,
    TopologyFromPrmtop,
    TopologyFromTop,
    TrajectoryFromModelAndCoordinates,
    TrajectoryFromBlob,
    TrajectoryFromMmCif,
    TrajectoryFromPDB,
    TrajectoryFromGRO,
    TrajectoryFromXYZ,
    TrajectoryFromLammpsData,
    TrajectoryFromLammpsTrajData,
    TrajectoryFromMOL,
    TrajectoryFromSDF,
    TrajectoryFromMOL2,
    TrajectoryFromCube,
    TrajectoryFromCifCore,
    ModelFromTrajectory,
    ModelWithCoordinates,
    StructureFromTrajectory,
    StructureFromModel,
    TransformStructureConformation,
    StructureInstances,
    StructureSelectionFromExpression,
    MultiStructureSelectionFromExpression,
    MultiStructureSelectionFromBundle,
    StructureSelectionFromScript,
    StructureSelectionFromBundle,
    StructureComplexElement,
    StructureComponent,
    CustomModelProperties,
    CustomStructureProperties,
    getTrajectory,
    StructureComplexElementTypes,
  },
  Particles: {
    ParticleListFromRelionStar,
    ParticleListFromDynamoTbl,
    ParticleListFromCryoEtDataPortalNdjson,
    ParticleListFromArtiatomiEm,
    ParticleListFromMmcifAssembly,
    ParticleTrajectoryFromSimularium,
    ParticleListFromTrajectory,
    ParticlesRepresentation3D,
    ParticleListWithTargets,
    ParticleListUnitcell3D,
  },
  Volume: {
    VolumeFromCcp4,
    VolumeFromDsn6,
    VolumeFromCube,
    VolumeFromDx,
    AssignColorVolume,
    VolumeFromDensityServerCif,
    VolumeFromSegmentationCif,
    VolumeFromStructureFactorsCif,
    VolumeFromMtz,
    VolumeTransform,
    VolumeInstances,
    CustomVolumeProperties,
  },
  Representation: {
    StructureRepresentation3D,
    ExplodeStructureRepresentation3D,
    SpinStructureRepresentation3D,
    UnwindStructureAssemblyRepresentation3D,
    OverpaintStructureRepresentation3DFromScript,
    OverpaintStructureRepresentation3DFromBundle,
    TransparencyStructureRepresentation3DFromScript,
    TransparencyStructureRepresentation3DFromBundle,
    EmissiveStructureRepresentation3DFromScript,
    EmissiveStructureRepresentation3DFromBundle,
    SubstanceStructureRepresentation3DFromScript,
    SubstanceStructureRepresentation3DFromBundle,
    ClippingStructureRepresentation3DFromScript,
    ClippingStructureRepresentation3DFromBundle,
    WiggleStructureRepresentation3DFromScript,
    WiggleStructureRepresentation3DFromBundle,
    ThemeStrengthRepresentation3D,
    VolumeRepresentation3D,
    ShapeRepresentation3D,
    ModelUnitcell3D,
    StructureBoundingBox3D,
    StructureSelectionsDistance3D,
    StructureSelectionsAngle3D,
    StructureSelectionsDihedral3D,
    StructureSelectionsLabel3D,
    StructureSelectionsOrientation3D,
    StructureSelectionsPlane3D,
    VolumeRepresentation3DHelpers,
  },
  Shape: {
    BoxShape3D,
    ShapeFromPly,
    ShapeFromObj,
    ShapeFromVtp,
    getBoxMesh,
  },
};

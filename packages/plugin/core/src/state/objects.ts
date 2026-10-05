/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Ccp4File } from '@molstar/io/reader/ccp4/schema';
import type { CifFile } from '@molstar/io/reader/cif';
import type { ArtiatomiEmFile } from '@molstar/io/reader/artiatomi/em';
import type { CryoEtDataPortalNdjsonFile } from '@molstar/io/reader/cryoet/ndjson';
import type { SimulariumFile } from '@molstar/io/reader/simularium/schema';
import type { MtzFile } from '@molstar/io/reader/mtz/schema';
import type { DcdFile } from '@molstar/io/reader/dcd/parser';
import type { DynamoTblFile } from '@molstar/io/reader/dynamo/tbl';
import type { Dsn6File } from '@molstar/io/reader/dsn6/schema';
import type { PlyFile } from '@molstar/io/reader/ply/schema';
import type { ObjFile } from '@molstar/io/reader/obj/schema';
import type { VtpFile } from '@molstar/io/reader/vtp/schema';
import type { PsfFile } from '@molstar/io/reader/psf/parser';
import type { ParticleList } from '@molstar/model/model/particles/particle-list';
import type { ParticleTrajectory } from '@molstar/model/model/particles/particle-trajectory';
import type { ShapeProvider } from '@molstar/graphics/geo/shape/provider';
import type {
  Coordinates as _Coordinates,
  Model as _Model,
  Structure as _Structure,
  Trajectory as _Trajectory,
  StructureElement,
  Topology as _Topology,
} from '@molstar/model/model/structure';
import type { Volume as _Volume } from '@molstar/model/model/volume';
import type { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import type { Representation } from '@molstar/graphics/repr/representation';
import type { ShapeRepresentation } from '@molstar/graphics/repr/shape/representation';
import type {
  StructureRepresentation,
  StructureRepresentationState,
} from '@molstar/graphics/repr/structure/representation';
import type { VolumeRepresentation } from '@molstar/graphics/repr/volume/representation';
import type { ParticleRepresentation } from '@molstar/graphics/repr/particles/representation';
import { StateObject, StateTransformer } from '@molstar/core/state';
import type { CubeFile } from '@molstar/io/reader/cube/parser';
import type { DxFile } from '@molstar/io/reader/dx/parser';
import type { Color } from '@molstar/core/util/color/color';
import type { Asset } from '@molstar/core/util/assets';
import type { PrmtopFile } from '@molstar/io/reader/prmtop/parser';
import type { TopFile } from '@molstar/io/reader/top/parser';
import type { StringLike } from '@molstar/core/util/string-like';

export type TypeClass = 'root' | 'data' | 'prop';

export namespace PluginStateObject {
  export type Any = StateObject<any, TypeInfo>;

  export type TypeClass = 'Root' | 'Group' | 'Data' | 'Object' | 'Representation3D' | 'Behavior';
  export interface TypeInfo {
    name: string;
    typeClass: TypeClass;
  }

  export const Create = StateObject.factory<TypeInfo>();

  export function isRepresentation3D(o?: Any): o is StateObject<Representation3DData<Representation.Any>, TypeInfo> {
    return !!o && o.type.typeClass === 'Representation3D';
  }

  export function isBehavior(o?: Any): o is StateObject<PluginBehavior, TypeInfo> {
    return !!o && o.type.typeClass === 'Behavior';
  }

  export interface Representation3DData<T extends Representation.Any, S = any> {
    repr: T;
    sourceData: S;
  }
  export function CreateRepresentation3D<T extends Representation.Any, S = any>(type: { name: string }) {
    return Create<Representation3DData<T, S>>({ ...type, typeClass: 'Representation3D' });
  }

  export function CreateBehavior<T extends PluginBehavior>(type: { name: string }) {
    return Create<T>({ ...type, typeClass: 'Behavior' });
  }

  export class Root extends Create({ name: 'Root', typeClass: 'Root' }) {}
  export class Group extends Create({ name: 'Group', typeClass: 'Group' }) {}

  export namespace Data {
    export class String extends Create<StringLike>({ name: 'String Data', typeClass: 'Data' }) {}
    export class Binary extends Create<Uint8Array<ArrayBuffer>>({ name: 'Binary Data', typeClass: 'Data' }) {}

    export type BlobEntry = { id: string } & (
      | { kind: 'string'; data: string }
      | { kind: 'binary'; data: Uint8Array<ArrayBuffer> }
    );
    export type BlobData = BlobEntry[];
    export class Blob extends Create<BlobData>({ name: 'Data Blob', typeClass: 'Data' }) {}
  }

  export namespace Format {
    export class Json extends Create<any>({ name: 'JSON Data', typeClass: 'Data' }) {}
    export class Cif extends Create<CifFile>({ name: 'CIF File', typeClass: 'Data' }) {}
    export class ArtiatomiEm extends Create<ArtiatomiEmFile>({ name: 'Artiatomi EM File', typeClass: 'Data' }) {}
    export class DynamoTbl extends Create<DynamoTblFile>({ name: 'Dynamo TBL File', typeClass: 'Data' }) {}
    export class CryoEtDataPortalNdjson extends Create<CryoEtDataPortalNdjsonFile>({
      name: 'CryoET Data Portal NDJSON File',
      typeClass: 'Data',
    }) {}
    export class Simularium extends Create<SimulariumFile>({ name: 'Simularium File', typeClass: 'Data' }) {}
    export class Cube extends Create<CubeFile>({ name: 'Cube File', typeClass: 'Data' }) {}
    export class Psf extends Create<PsfFile>({ name: 'PSF File', typeClass: 'Data' }) {}
    export class Prmtop extends Create<PrmtopFile>({ name: 'PRMTOP File', typeClass: 'Data' }) {}
    export class Top extends Create<TopFile>({ name: 'TOP File', typeClass: 'Data' }) {}
    export class Ply extends Create<PlyFile>({ name: 'PLY File', typeClass: 'Data' }) {}
    export class Obj extends Create<ObjFile>({ name: 'OBJ File', typeClass: 'Data' }) {}
    export class Vtp extends Create<VtpFile>({ name: 'VTP File', typeClass: 'Data' }) {}
    export class Ccp4 extends Create<Ccp4File>({ name: 'CCP4/MRC/MAP File', typeClass: 'Data' }) {}
    export class Dsn6 extends Create<Dsn6File>({ name: 'DSN6/BRIX File', typeClass: 'Data' }) {}
    export class Dx extends Create<DxFile>({ name: 'DX File', typeClass: 'Data' }) {}
    export class Mtz extends Create<MtzFile>({ name: 'MTZ File', typeClass: 'Data' }) {}

    export type BlobEntry = { id: string } & (
      | { kind: 'json'; data: unknown }
      | { kind: 'string'; data: string }
      | { kind: 'binary'; data: Uint8Array<ArrayBuffer> }
      | { kind: 'cif'; data: CifFile }
      | { kind: 'pdb'; data: CifFile }
      | { kind: 'gro'; data: CifFile }
      | { kind: 'dcd'; data: DcdFile }
      | { kind: 'ccp4'; data: Ccp4File }
      | { kind: 'dsn6'; data: Dsn6File }
      | { kind: 'dx'; data: DxFile }
      | { kind: 'ply'; data: PlyFile }
      // For non-build in extensions
      | { kind: 'custom'; data: unknown; tag: string }
    );
    export type BlobData = BlobEntry[];
    export class Blob extends Create<BlobData>({ name: 'Format Blob', typeClass: 'Data' }) {}
  }

  export namespace Molecule {
    export class Coordinates extends Create<_Coordinates>({ name: 'Coordinates', typeClass: 'Object' }) {}
    export class Topology extends Create<_Topology>({ name: 'Topology', typeClass: 'Object' }) {}
    export class Model extends Create<_Model>({ name: 'Model', typeClass: 'Object' }) {}
    export class Trajectory extends Create<_Trajectory>({ name: 'Trajectory', typeClass: 'Object' }) {}
    export class Structure extends Create<_Structure>({ name: 'Structure', typeClass: 'Object' }) {}

    export namespace Structure {
      export class Representation3D extends CreateRepresentation3D<StructureRepresentation<any>, _Structure>({
        name: 'Structure 3D',
      }) {}

      export interface Representation3DStateData {
        repr: StructureRepresentation<any>;
        /** used to restore state when the obj is removed */
        initialState: Partial<StructureRepresentationState>;
        state: Partial<StructureRepresentationState>;
        info?: unknown;
      }
      export class Representation3DState extends Create<Representation3DStateData>({
        name: 'Structure 3D State',
        typeClass: 'Object',
      }) {}

      export interface SelectionEntry {
        key: string;
        structureRef: string;
        groupId?: string;
        loci: StructureElement.Loci;
      }
      export class Selections extends Create<ReadonlyArray<SelectionEntry>>({
        name: 'Selections',
        typeClass: 'Object',
      }) {}
    }
  }

  export namespace Volume {
    export interface LazyInfo {
      url: string | Asset.Url;
      isBinary: boolean;
      format: string;
      entryId?: string | string[];
      isovalues: {
        type: 'absolute' | 'relative';
        value: number;
        color: Color;
        alpha?: number;
        volumeIndex?: number;
      }[];
    }

    export class Data extends Create<_Volume>({ name: 'Volume', typeClass: 'Object' }) {}
    export class Lazy extends Create<LazyInfo>({ name: 'Lazy Volume', typeClass: 'Object' }) {}
    export class Representation3D extends CreateRepresentation3D<VolumeRepresentation<any>, _Volume>({
      name: 'Volume 3D',
    }) {}
  }

  export namespace Shape {
    export class Provider extends Create<ShapeProvider<any, any, any>>({
      name: 'Shape Provider',
      typeClass: 'Object',
    }) {}
    export class Representation3D extends CreateRepresentation3D<ShapeRepresentation<any, any, any>, unknown>({
      name: 'Shape 3D',
    }) {}
  }

  export namespace Particle {
    export class Trajectory extends Create<ParticleTrajectory>({ name: 'Particle Trajectory', typeClass: 'Object' }) {}
    export class List extends Create<ParticleList>({ name: 'Particle List', typeClass: 'Object' }) {}
    export class Representation3D extends CreateRepresentation3D<ParticleRepresentation<any>, ParticleList>({
      name: 'Particle 3D',
    }) {}
  }
}

export namespace PluginStateTransform {
  export const CreateBuiltIn = StateTransformer.factory('ms-plugin');
  export const BuiltIn = StateTransformer.builderFactory('ms-plugin');
}

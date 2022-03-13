/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Ccp4File } from '../mol-io/reader/ccp4/schema';
import { CifFile } from '../mol-io/reader/cif';
import { DcdFile } from '../mol-io/reader/dcd/parser';
import { Dsn6File } from '../mol-io/reader/dsn6/schema';
import { PlyFile } from '../mol-io/reader/ply/schema';
import { PsfFile } from '../mol-io/reader/psf/parser';
import { ShapeProvider } from '../mol-model/shape/provider';
import { Coordinates as _Coordinates, Model as _Model, Structure as _Structure, Trajectory as _Trajectory, StructureElement, Topology as _Topology } from '../mol-model/structure';
import { Volume as _Volume } from '../mol-model/volume';
import { PluginBehavior } from '../mol-plugin/behavior/behavior';
import { Representation } from '../mol-repr/representation';
import { ShapeRepresentation } from '../mol-repr/shape/representation';
import { StructureRepresentation, StructureRepresentationState } from '../mol-repr/structure/representation';
import { VolumeRepresentation } from '../mol-repr/volume/representation';
import { StateObject, StateTransformer } from '../mol-state';
import { CubeFile } from '../mol-io/reader/cube/parser';
import { DxFile } from '../mol-io/reader/dx/parser';
import { Color } from '../mol-util/color/color';
import { Asset } from '../mol-util/assets';
import { PrmtopFile } from '../mol-io/reader/prmtop/parser';
import { TopFile } from '../mol-io/reader/top/parser';

export type TypeClass = 'root' | 'data' | 'prop'

export namespace PluginStateObject {
    export type Any = StateObject<any, TypeInfo>

    export type TypeClass = 'Root' | 'Group' | 'Data' | 'Object' | 'Representation3D' | 'Behavior'
    export interface TypeInfo { name: string, typeClass: TypeClass }

    export const Create = StateObject.factory<TypeInfo>();

    export function isRepresentation3D(o?: Any): o is StateObject<Representation3DData<Representation.Any>, TypeInfo> {
        return !!o && o.type.typeClass === 'Representation3D';
    }

    export function isBehavior(o?: Any): o is StateObject<PluginBehavior, TypeInfo> {
        return !!o && o.type.typeClass === 'Behavior';
    }

    export interface Representation3DData<T extends Representation.Any, S = any> { repr: T, sourceData: S }
    export function CreateRepresentation3D<T extends Representation.Any, S = any>(type: { name: string }) {
        return Create<Representation3DData<T, S>>({ ...type, typeClass: 'Representation3D' });
    }

    export function CreateBehavior<T extends PluginBehavior>(type: { name: string }) {
        return Create<T>({ ...type, typeClass: 'Behavior' });
    }

    export class Root extends Create({ name: 'Root', typeClass: 'Root' }) { }
    export class Group extends Create({ name: 'Group', typeClass: 'Group' }) { }

    export namespace Data {
        export class String extends Create<string>({ name: 'String Data', typeClass: 'Data', }) { }
        export class Binary extends Create<Uint8Array>({ name: 'Binary Data', typeClass: 'Data' }) { }

        export type BlobEntry = { id: string } & (
            { kind: 'string', data: string } |
            { kind: 'binary', data: Uint8Array }
        )
        export type BlobData = BlobEntry[]
        export class Blob extends Create<BlobData>({ name: 'Data Blob', typeClass: 'Data' }) { }
    }

    export namespace Format {
        export class Json extends Create<any>({ name: 'JSON Data', typeClass: 'Data' }) { }
        export class Cif extends Create<CifFile>({ name: 'CIF File', typeClass: 'Data' }) { }
        export class Cube extends Create<CubeFile>({ name: 'Cube File', typeClass: 'Data' }) { }
        export class Psf extends Create<PsfFile>({ name: 'PSF File', typeClass: 'Data' }) { }
        export class Prmtop extends Create<PrmtopFile>({ name: 'PRMTOP File', typeClass: 'Data' }) { }
        export class Top extends Create<TopFile>({ name: 'TOP File', typeClass: 'Data' }) { }
        export class Ply extends Create<PlyFile>({ name: 'PLY File', typeClass: 'Data' }) { }
        export class Ccp4 extends Create<Ccp4File>({ name: 'CCP4/MRC/MAP File', typeClass: 'Data' }) { }
        export class Dsn6 extends Create<Dsn6File>({ name: 'DSN6/BRIX File', typeClass: 'Data' }) { }
        export class Dx extends Create<DxFile>({ name: 'DX File', typeClass: 'Data' }) { }

        export type BlobEntry = { id: string } & (
            { kind: 'json', data: unknown } |
            { kind: 'string', data: string } |
            { kind: 'binary', data: Uint8Array } |
            { kind: 'cif', data: CifFile } |
            { kind: 'pdb', data: CifFile } |
            { kind: 'gro', data: CifFile } |
            { kind: 'dcd', data: DcdFile } |
            { kind: 'ccp4', data: Ccp4File } |
            { kind: 'dsn6', data: Dsn6File } |
            { kind: 'dx', data: DxFile } |
            { kind: 'ply', data: PlyFile } |
            // For non-build in extensions
            { kind: 'custom', data: unknown, tag: string }
        )
        export type BlobData = BlobEntry[]
        export class Blob extends Create<BlobData>({ name: 'Format Blob', typeClass: 'Data' }) { }
    }

    export namespace Molecule {
        export class Coordinates extends Create<_Coordinates>({ name: 'Coordinates', typeClass: 'Object' }) { }
        export class Topology extends Create<_Topology>({ name: 'Topology', typeClass: 'Object' }) { }
        export class Model extends Create<_Model>({ name: 'Model', typeClass: 'Object' }) { }
        export class Trajectory extends Create<_Trajectory>({ name: 'Trajectory', typeClass: 'Object' }) { }
        export class Structure extends Create<_Structure>({ name: 'Structure', typeClass: 'Object' }) { }

        export namespace Structure {
            export class Representation3D extends CreateRepresentation3D<StructureRepresentation<any>, _Structure>({ name: 'Structure 3D' }) { }

            export interface Representation3DStateData {
                repr: StructureRepresentation<any>,
                /** used to restore state when the obj is removed */
                initialState: Partial<StructureRepresentationState>,
                state: Partial<StructureRepresentationState>,
                info?: unknown
            }
            export class Representation3DState extends Create<Representation3DStateData>({ name: 'Structure 3D State', typeClass: 'Object' }) { }

            export interface SelectionEntry { key: string, groupId?: string, loci: StructureElement.Loci }
            export class Selections extends Create<ReadonlyArray<SelectionEntry>>({ name: 'Selections', typeClass: 'Object' }) {}
        }
    }

    export namespace Volume {
        export interface LazyInfo {
            url: string | Asset.Url,
            isBinary: boolean,
            format: string,
            entryId?: string | string[],
            isovalues: {
                type: 'absolute' | 'relative',
                value: number,
                color: Color,
                alpha?: number,
                volumeIndex?: number,
            }[]
        }

        export class Data extends Create<_Volume>({ name: 'Volume', typeClass: 'Object' }) { }
        export class Lazy extends Create<LazyInfo>({ name: 'Lazy Volume', typeClass: 'Object' }) { }
        export class Representation3D extends CreateRepresentation3D<VolumeRepresentation<any>, _Volume>({ name: 'Volume 3D' }) { }
    }

    export namespace Shape {
        export class Provider extends Create<ShapeProvider<any, any, any>>({ name: 'Shape Provider', typeClass: 'Object' }) { }
        export class Representation3D extends CreateRepresentation3D<ShapeRepresentation<any, any, any>, unknown>({ name: 'Shape 3D' }) { }
    }
}

export namespace PluginStateTransform {
    export const CreateBuiltIn = StateTransformer.factory('ms-plugin');
    export const BuiltIn = StateTransformer.builderFactory('ms-plugin');
}
/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CifFile } from 'mol-io/reader/cif';
import { Model as _Model, Structure as _Structure } from 'mol-model/structure';
import { VolumeData } from 'mol-model/volume';
import { PluginBehavior } from 'mol-plugin/behavior/behavior';
import { Representation } from 'mol-repr/representation';
import { StructureRepresentation, StructureRepresentationState } from 'mol-repr/structure/representation';
import { VolumeRepresentation } from 'mol-repr/volume/representation';
import { StateObject, StateTransformer } from 'mol-state';
import { Ccp4File } from 'mol-io/reader/ccp4/schema';
import { Dsn6File } from 'mol-io/reader/dsn6/schema';
import { ShapeRepresentation } from 'mol-repr/shape/representation';

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

    export interface Representation3DData<T extends Representation.Any, S extends StateObject = StateObject> { repr: T, source: S }
    export function CreateRepresentation3D<T extends Representation.Any, S extends StateObject = StateObject>(type: { name: string }) {
        return Create<Representation3DData<T, S>>({ ...type, typeClass: 'Representation3D' });
    }

    export function CreateBehavior<T extends PluginBehavior>(type: { name: string }) {
        return Create<T>({ ...type, typeClass: 'Behavior' })
    }

    export class Root extends Create({ name: 'Root', typeClass: 'Root' }) { }
    export class Group extends Create({ name: 'Group', typeClass: 'Group' }) { }

    export namespace Data {
        export class String extends Create<string>({ name: 'String Data', typeClass: 'Data', }) { }
        export class Binary extends Create<Uint8Array>({ name: 'Binary Data', typeClass: 'Data' }) { }

        export type BlobEntry = { id: string } &
            ( { kind: 'string', data: string }
            | { kind: 'binary', data: Uint8Array })
        export type BlobData = BlobEntry[]
        export class Blob extends Create<BlobData>({ name: 'Data Blob', typeClass: 'Data' }) { }
    }

    export namespace Format {
        export class Json extends Create<any>({ name: 'JSON Data', typeClass: 'Data' }) { }
        export class Cif extends Create<CifFile>({ name: 'CIF File', typeClass: 'Data' }) { }
        export class Ccp4 extends Create<Ccp4File>({ name: 'CCP4/MRC/MAP File', typeClass: 'Data' }) { }
        export class Dsn6 extends Create<Dsn6File>({ name: 'DSN6/BRIX File', typeClass: 'Data' }) { }

        export type BlobEntry = { id: string } &
            ( { kind: 'json', data: unknown }
            | { kind: 'string', data: string }
            | { kind: 'binary', data: Uint8Array }
            | { kind: 'cif', data: CifFile }
            | { kind: 'ccp4', data: Ccp4File }
            | { kind: 'dsn6', data: Dsn6File }
            // For non-build in extensions
            | { kind: 'custom', data: unknown, tag: string })
        export type BlobData = BlobEntry[]
        export class Blob extends Create<BlobData>({ name: 'Format Blob', typeClass: 'Data' }) { }
    }

    export namespace Molecule {
        export class Trajectory extends Create<ReadonlyArray<_Model>>({ name: 'Trajectory', typeClass: 'Object' }) { }
        export class Model extends Create<_Model>({ name: 'Model', typeClass: 'Object' }) { }
        export class Structure extends Create<_Structure>({ name: 'Structure', typeClass: 'Object' }) { }

        export namespace Structure {
            export class Representation3D extends CreateRepresentation3D<StructureRepresentation<any> | ShapeRepresentation<any, any, any>, Structure>({ name: 'Structure 3D' }) { }

            export interface Representation3DStateData {
                source: Representation3D,
                /** used to restore state when the obj is removed */
                initialState: Partial<StructureRepresentationState>,
                state: Partial<StructureRepresentationState>,
                info?: unknown
            }
            export class Representation3DState extends Create<Representation3DStateData>({ name: 'Structure 3D State', typeClass: 'Object' }) { }
        }
    }

    export namespace Volume {
        export class Data extends Create<VolumeData>({ name: 'Volume Data', typeClass: 'Object' }) { }
        export class Representation3D extends CreateRepresentation3D<VolumeRepresentation<any>>({ name: 'Volume 3D' }) { }
    }
}

export namespace PluginStateTransform {
    export const CreateBuiltIn = StateTransformer.factory('ms-plugin');
    export const BuiltIn = StateTransformer.builderFactory('ms-plugin');
}
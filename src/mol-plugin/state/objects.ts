/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateObject } from './base';
import { CifFile } from 'mol-io/reader/cif';
import { Model as _Model, Structure as _Structure } from 'mol-model/structure'
import { StructureRepresentation } from 'mol-repr/structure/index';
import { VolumeRepresentation } from 'mol-repr/volume';
import { VolumeData } from 'mol-model/volume';

const _create = PluginStateObject.Create, _createRepr3D = PluginStateObject.CreateRepresentation3D

namespace PluginStateObjects {
    export class Root extends _create({ name: 'Root', typeClass: 'Root' }) { }

    export class Group extends _create({ name: 'Group', typeClass: 'Group' }) { }

    export namespace Data {
        export class String extends _create<string>({ name: 'String Data', typeClass: 'Data', }) { }
        export class Binary extends _create<Uint8Array>({ name: 'Binary Data', typeClass: 'Data' }) { }
        export class Json extends _create<any>({ name: 'JSON Data', typeClass: 'Data' }) { }
        export class Cif extends _create<CifFile>({ name: 'CIF File', typeClass: 'Data' }) { }

        // TODO
        // export class MultipleRaw extends _create<{
        //     [key: string]: { type: 'String' | 'Binary', data: string | Uint8Array }
        // }>({ name: 'Data', typeClass: 'Data', shortName: 'MD', description: 'Multiple Keyed Data.' }) { }
    }

    export namespace Molecule {
        export class Trajectory extends _create<ReadonlyArray<_Model>>({ name: 'Trajectory', typeClass: 'Object' }) { }
        export class Model extends _create<_Model>({ name: 'Model', typeClass: 'Object' }) { }
        export class Structure extends _create<_Structure>({ name: 'Structure', typeClass: 'Object' }) { }
        export class Representation3D extends _createRepr3D<StructureRepresentation<any>>({ name: 'Structure 3D' }) { }
    }

    export namespace Volume {
        export class Data extends _create<VolumeData>({ name: 'Volume Data', typeClass: 'Object' }) { }
        export class Representation3D extends _createRepr3D<VolumeRepresentation<any>>({ name: 'Volume 3D' }) { }
    }

}

export { PluginStateObjects }
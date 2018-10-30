/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateObject } from './base';
import { CifFile } from 'mol-io/reader/cif';
import { Model as _Model, Structure as _Structure } from 'mol-model/structure'

const _create = PluginStateObject.Create

namespace PluginStateObjects {
    export class Root extends _create({ name: 'Root', shortName: 'R', typeClass: 'Root', description: 'Where everything begins.' }) { }
    export class Group extends _create({ name: 'Group', shortName: 'G', typeClass: 'Group', description: 'A group on entities.' }) { }

    export namespace Data {
        export class String extends _create<string>({ name: 'String Data', typeClass: 'Data', shortName: 'S_D', description: 'A string.' }) { }
        export class Binary extends _create<Uint8Array>({ name: 'Binary Data', typeClass: 'Data', shortName: 'B_D', description: 'A binary blob.' }) { }
        export class Json extends _create<any>({ name: 'JSON Data', typeClass: 'Data', shortName: 'JS_D', description: 'Represents JSON data.' }) { }
        export class Cif extends _create<CifFile>({ name: 'Cif File', typeClass: 'Data', shortName: 'CF', description: 'Represents parsed CIF data.' }) { }

        // TODO
        // export class MultipleRaw extends _create<{
        //     [key: string]: { type: 'String' | 'Binary', data: string | Uint8Array }
        // }>({ name: 'Data', typeClass: 'Data', shortName: 'MD', description: 'Multiple Keyed Data.' }) { }
    }

    export class Model extends _create<_Model>({ name: 'Molecule Model', typeClass: 'Object', shortName: 'M_M', description: 'A model of a molecule.' }) { }
    export class Structure extends _create<_Structure>({ name: 'Molecule Structure', typeClass: 'Object', shortName: 'M_S', description: 'A structure of a molecule.' }) { }


    export class StructureRepresentation extends _create<{
        // TODO
    }>({ name: 'Molecule Structure Representation', typeClass: 'Representation', shortName: 'S_R', description: 'A representation of a molecular structure.' }) { }
}

export { PluginStateObjects }
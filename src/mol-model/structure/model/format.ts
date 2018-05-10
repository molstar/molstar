/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// import { File as GroFile } from 'mol-io/reader/gro/schema'
import { mmCIF_Database } from 'mol-io/reader/cif/schema/mmcif'

type Format =
    // | Format.gro
    | Format.mmCIF

namespace Format {
    // export interface gro { kind: 'gro', data: GroFile }
    export interface mmCIF { kind: 'mmCIF', data: mmCIF_Database }
}

export default Format
/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mmCIF_Database } from 'mol-io/reader/cif/schema/mmcif'

type Format =
    | Format.mmCIF

namespace Format {
    export interface mmCIF { kind: 'mmCIF', data: mmCIF_Database }
}

export default Format
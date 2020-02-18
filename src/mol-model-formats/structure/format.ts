/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { mmCIF_Database } from '../../mol-io/reader/cif/schema/mmcif';
import { CIF, CifFrame } from '../../mol-io/reader/cif';

interface Format { readonly kind: string, name: string }

type ModelFormat =
    | ModelFormat.mmCIF

namespace ModelFormat {
    export interface mmCIF extends Format {
        readonly kind: 'mmCIF', data: mmCIF_Database, frame: CifFrame
    }
    export function mmCIF(frame: CifFrame, data?: mmCIF_Database): mmCIF {
        if (!data) data = CIF.schema.mmCIF(frame)
        return { kind: 'mmCIF', name: data._name, data, frame };
    }
}

export { ModelFormat }
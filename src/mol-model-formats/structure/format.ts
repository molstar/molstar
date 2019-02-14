/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mmCIF_Database } from 'mol-io/reader/cif/schema/mmcif';
import CIF, { CifFrame } from 'mol-io/reader/cif';

type ModelFormat =
    | ModelFormat.mmCIF

namespace ModelFormat {
    export interface mmCIF { kind: 'mmCIF', data: mmCIF_Database, frame: CifFrame }
    export function mmCIF(frame: CifFrame, data?: mmCIF_Database): mmCIF { return { kind: 'mmCIF', data: data || CIF.schema.mmCIF(frame), frame }; }
}

export { ModelFormat }
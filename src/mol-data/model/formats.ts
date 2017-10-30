/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Frame as mmCIF_Frame } from '../../mol-io/reader/cif/schema/mmcif'

export type RawData =
    | { source: 'mmCIF', data: mmCIF_Frame }


import { loadCheckpoint } from '../mol-util/debug';
loadCheckpoint(`mol-model-formats/format.ts::start`);
/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export interface ModelFormat<T = unknown> { readonly kind: string, name: string, data: T }
loadCheckpoint(`mol-model-formats/format.ts::end`);

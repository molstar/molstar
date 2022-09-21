

import { loadCheckpoint } from '../mol-util/debug';
loadCheckpoint(`mol-util/sleep.ts::start`);
/**
 * Copyright (c) 2019 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export function sleep(milliseconds: number) {
    return new Promise(resolve => setTimeout(resolve, milliseconds));
}
loadCheckpoint(`mol-util/sleep.ts::end`);

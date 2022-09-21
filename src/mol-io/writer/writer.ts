

import { loadCheckpoint } from '../../mol-util/debug';
loadCheckpoint(`mol-io/writer/writer.ts::start`);
/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

interface Writer {
    writeString(data: string): boolean,
    writeBinary(data: Uint8Array): boolean
}

namespace Writer {

}

export { Writer };
loadCheckpoint(`mol-io/writer/writer.ts::end`);

/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Writer } from './writer.js';

export interface Encoder {
    encode(): void,
    writeTo(writer: Writer): void,
    getSize(): number
}

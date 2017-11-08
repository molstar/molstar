/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import TextCIFEncoder from './cif/encoder/text'
import BinaryCIFEncoder from './cif/encoder/binary'

export * from './cif/encoder'

export function create(params?: { binary?: boolean, encoderName?: string }) {
    const { binary = false, encoderName = 'mol*' } = params || {};
    return binary ? new BinaryCIFEncoder(encoderName) : new TextCIFEncoder();
}
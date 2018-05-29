/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import TextCIFEncoder from './cif/encoder/text'
import BinaryCIFEncoder from './cif/encoder/binary'
import * as _Encoder from './cif/encoder'

// export * from './cif/encoder'

// export function createCIFEncoder(params?: { binary?: boolean, encoderName?: string }): CIFEncoder {
//     const { binary = false, encoderName = 'mol*' } = params || {};
//     return binary ? new BinaryCIFEncoder(encoderName) : new TextCIFEncoder();
//}

export namespace CifWriter {
    export import Encoder = _Encoder.CIFEncoder
    export import Category = _Encoder.CIFCategory
    export import Field = _Encoder.CIFField

    export function createEncoder(params?: { binary?: boolean, encoderName?: string }): Encoder {
        const { binary = false, encoderName = 'mol*' } = params || {};
        return binary ? new BinaryCIFEncoder(encoderName) : new TextCIFEncoder();
    }
}
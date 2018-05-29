/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import TextEncoder from './cif/encoder/text'
import BinaryEncoder from './cif/encoder/binary'
import * as _Encoder from './cif/encoder'
import { ArrayEncoding } from '../common/binary-cif';

export namespace CifWriter {
    export import Encoder = _Encoder.Encoder
    export import Category = _Encoder.Category
    export import Field = _Encoder.Field
    export import Encoding = ArrayEncoding

    export function createEncoder(params?: { binary?: boolean, encoderName?: string }): Encoder {
        const { binary = false, encoderName = 'mol*' } = params || {};
        return binary ? new BinaryEncoder(encoderName) : new TextEncoder();
    }

    import E = Encoding
    export const Encodings = {
        deltaRLE: E.by(E.delta).and(E.runLength).and(E.integerPacking),
        fixedPoint2: E.by(E.fixedPoint(100)).and(E.delta).and(E.integerPacking),
        fixedPoint3: E.by(E.fixedPoint(1000)).and(E.delta).and(E.integerPacking),
    };
}
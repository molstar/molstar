/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import TextEncoder from './cif/encoder/text';
import BinaryEncoder, { BinaryEncodingProvider } from './cif/encoder/binary';
import * as _Encoder from './cif/encoder';
import { ArrayEncoding, ArrayEncoder } from '../common/binary-cif';
import { CifFrame } from '../reader/cif';

export namespace CifWriter {
    export import Encoder = _Encoder.Encoder
    export import Category = _Encoder.Category
    export import Field = _Encoder.Field
    export import Encoding = ArrayEncoding

    export interface EncoderParams {
        binary?: boolean,
        encoderName?: string,
        binaryEncodingPovider?: BinaryEncodingProvider,
        binaryAutoClassifyEncoding?: boolean
    }

    export function createEncoder(params?: EncoderParams): Encoder {
        const { binary = false, encoderName = 'mol*' } = params || {};
        return binary ? new BinaryEncoder(encoderName, params ? params.binaryEncodingPovider : void 0, params ? !!params.binaryAutoClassifyEncoding : false) : new TextEncoder();
    }

    export function fields<K = number, D = any, N extends string = string>() {
        return Field.build<K, D, N>();
    }

    import E = Encoding
    export const Encodings = {
        deltaRLE: E.by(E.delta).and(E.runLength).and(E.integerPacking),
        fixedPoint2: E.by(E.fixedPoint(100)).and(E.delta).and(E.integerPacking),
        fixedPoint3: E.by(E.fixedPoint(1000)).and(E.delta).and(E.integerPacking),
    };

    export function categoryInstance<Key, Data>(fields: Field<Key, Data>[], source: Category.DataSource): Category.Instance {
        return { fields, source: [source] };
    }

    export function createEncodingProviderFromCifFrame(frame: CifFrame): BinaryEncodingProvider {
        return {
            get(c, f) {
                const cat = frame.categories[c];
                if (!cat) return void 0;
                const ff = cat.getField(f);
                return ff && ff.binaryEncoding ? ArrayEncoder.fromEncoding(ff.binaryEncoding) : void 0;
            }
        };
    };

    export function createEncodingProviderFromJsonConfig(hints: EncodingStrategyHint[]): BinaryEncodingProvider {
        return {
            get(c, f) {
                for (let i = 0; i < hints.length; i++) {
                    const hint = hints[i];
                    if (hint.categoryName === c && hint.columnName === f) {
                        return resolveEncoding(hint);
                    }
                }
            }
        };
    }

    function resolveEncoding(hint: EncodingStrategyHint): ArrayEncoder {
        const precision: number | undefined = hint.precision;
        if (precision !== void 0) {
            const multiplier = Math.pow(10, precision);
            const fixedPoint = E.by(E.fixedPoint(multiplier));
            switch (hint.encoding) {
                case 'pack':
                    return fixedPoint.and(E.integerPacking);
                case 'rle':
                    return fixedPoint.and(E.runLength).and(E.integerPacking);
                case 'delta':
                    return fixedPoint.and(E.delta).and(E.integerPacking);
                case 'delta-rle':
                    return fixedPoint.and(E.delta).and(E.runLength).and(E.integerPacking);
            };
        } else {
            switch (hint.encoding) {
                case 'pack':
                    return E.by(E.integerPacking);
                case 'rle':
                    return E.by(E.runLength).and(E.integerPacking);
                case 'delta':
                    return E.by(E.delta).and(E.integerPacking);
                case 'delta-rle':
                    return E.by(E.delta).and(E.runLength).and(E.integerPacking);
            }
        }
        throw new Error('cannot be reached');
    }
}

/**
 * Defines the information needed to encode certain fields: category and column name as well as encoding tag, precision is optional and identifies float columns.
 */
export interface EncodingStrategyHint {
    categoryName: string,
    columnName: string,
    // TODO would be nice to infer strategy and precision if needed
    encoding: EncodingType,
    /**
     * number of decimal places to keep - must be specified to float columns
     */
    precision?: number
}

type EncodingType = 'pack' | 'rle' | 'delta' | 'delta-rle'
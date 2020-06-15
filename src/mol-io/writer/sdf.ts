/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { SdfEncoder } from './sdf/encoder';
import { Encoder } from './cif/encoder';

export namespace SdfWriter {
    export interface EncoderParams {
        encoderName?: string,
        hideMetaInformation?: boolean
    }

    export function createEncoder(params?: EncoderParams): Encoder {
        const { encoderName = 'mol*', hideMetaInformation = false } = params || {};
        return new SdfEncoder(encoderName, hideMetaInformation);
    }
}
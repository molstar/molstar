/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { MolEncoder } from './mol/encoder';
import { Encoder } from './cif/encoder';

export namespace MolWriter {
    export interface EncoderParams {
        encoderName?: string,
        // whether to write hydrogen atoms
        hydrogens?: boolean
    }

    export function createEncoder(params?: EncoderParams): Encoder {
        const { encoderName = 'mol*', hydrogens = true } = params || {};
        return new MolEncoder(encoderName, false, hydrogens);
    }
}
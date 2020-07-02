/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { Mol2Encoder } from './mol2/encoder';
import { Encoder } from './cif/encoder';

export namespace Mol2Writer {
    export interface EncoderParams {
        encoderName?: string,
        // whether to write ModelServer meta-information (query & params)
        metaInformation?: boolean,
        // whether to write hydrogen atoms
        hydrogens?: boolean
    }

    export function createEncoder(params?: EncoderParams): Encoder {
        const { encoderName = 'mol*', metaInformation = true, hydrogens = true } = params || {};
        return new Mol2Encoder(encoderName, metaInformation, hydrogens);
    }
}
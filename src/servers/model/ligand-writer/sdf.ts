/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { MolEncoder } from './mol/encoder';
import { Encoder } from '../../../mol-io/writer/cif/encoder';

export namespace SdfWriter {
    export interface EncoderParams {
        encoderName?: string,
        // whether to write ModelServer meta-information (query & params)
        metaInformation?: boolean,
        // whether to write hydrogen atoms
        hydrogens?: boolean
    }

    export function createEncoder(params?: EncoderParams): Encoder {
        const { encoderName = 'mol*', metaInformation = false, hydrogens = false } = params || {};
        return new MolEncoder(encoderName, metaInformation, hydrogens, '$$$$');
    }
}
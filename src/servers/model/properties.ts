/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model } from 'mol-model/structure';
import { PDBe_structureQualityReport } from './properties/pdbe';
import { RCSB_assemblySymmetry } from './properties/rcsb';

export function attachModelProperties(model: Model): Promise<any>[] {
    // return a list of promises that start attaching the props in parallel
    // (if there are downloads etc.)
    return [
        PDBe_structureQualityReport(model),
        RCSB_assemblySymmetry(model)
    ];
}
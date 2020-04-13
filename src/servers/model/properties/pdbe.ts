/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PDBe_structureQualityReport, PDBe_preferredAssembly, PDBe_structRefDomain } from './providers/pdbe';
import { AttachModelProperties } from '../property-provider';

export const attachModelProperties: AttachModelProperties = (args) => {
    // return a list of promises that start attaching the props in parallel
    // (if there are downloads etc.)
    return [
        PDBe_structureQualityReport(args),
        PDBe_preferredAssembly(args),
        PDBe_structRefDomain(args)
    ];
};
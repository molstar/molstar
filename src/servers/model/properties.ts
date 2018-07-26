/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Model } from 'mol-model/structure';
import { StructureQualityReport } from './properties/structure-quality-report';

export function attachModelProperties(model: Model): Promise<any>[] {
    // return a list of promises that start attaching the props in parallel
    // (if there are downloads etc.)
    return [
        StructureQualityReport.attachFromPDBeApi(model)
    ];
}
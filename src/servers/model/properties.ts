/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Model } from 'mol-model/structure';
import { StructureQualityReport } from './properties/structure-quality-report';

export async function attachModelProperties(model: Model) {
    await StructureQualityReport.attachFromPDBeApi(model);
}
/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

 import { Model } from 'mol-model/structure';
import { StructureQualityReport } from 'mol-model-props/pdbe/structure-quality-report';
import { fetchRetry } from '../../utils/fetch-retry';
import { UUID } from 'mol-util';

const cacheKey = UUID.create();
export function PDBe_structureQualityReport(model: Model, cache: any) {
    return StructureQualityReport.attachFromCifOrApi(model, {
        PDBe_apiSourceJson: async model => {
            if (cache[cacheKey]) return cache[cacheKey];
            const rawData = await fetchRetry(`https://www.ebi.ac.uk/pdbe/api/validation/residuewise_outlier_summary/entry/${model.label.toLowerCase()}`, 1500, 5);
            const json = await rawData.json();
            cache[cacheKey] = json;
            return json;
        }
    });
}
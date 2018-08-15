/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

 import { Model } from 'mol-model/structure';
import { StructureQualityReport } from 'mol-model-props/pdbe/structure-quality-report';
import { fetchRetry } from '../utils/fetch-retry';

export function PDBe_structureQualityReport(model: Model) {
    return StructureQualityReport.attachFromCifOrApi(model, {
        PDBe_apiSourceJson: async model => {
            const rawData = await fetchRetry(`https://www.ebi.ac.uk/pdbe/api/validation/residuewise_outlier_summary/entry/${model.label.toLowerCase()}`, 1500, 5);
            return await rawData.json();
        }
    });
}
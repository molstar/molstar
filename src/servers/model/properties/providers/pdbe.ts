/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'fs'
import { Model } from 'mol-model/structure';
import { StructureQualityReport } from 'mol-model-props/pdbe/structure-quality-report';
import { fetchRetry } from '../../utils/fetch-retry';
import { UUID } from 'mol-util';
import { PDBePreferredAssembly } from 'mol-model-props/pdbe/preferred-assembly';
import { PDBeStructRefDomain } from 'mol-model-props/pdbe/struct-ref-domain';

const USE_FILE_SOURCE = false;

export function PDBe_structureQualityReport(model: Model, cache: any) {
    return StructureQualityReport.attachFromCifOrApi(model, {
        PDBe_apiSourceJson: USE_FILE_SOURCE
            ? residuewise_outlier_summary.getDataFromAggregateFile
            : residuewise_outlier_summary.getDataFromApiProvider(cache)
    });
}

export function PDBe_preferredAssembly(model: Model, cache: any) {
    return PDBePreferredAssembly.attachFromCifOrApi(model, {
        PDBe_apiSourceJson: USE_FILE_SOURCE
            ? void 0
            : preferred_assembly.getDataFromApiProvider(cache)
    });
}

export function PDBe_structRefDomain(model: Model, cache: any) {
    return PDBeStructRefDomain.attachFromCifOrApi(model, {
        PDBe_apiSourceJson: USE_FILE_SOURCE
            ? void 0
            : struct_ref_domain.getDataFromApiProvider(cache)
    });
}

namespace preferred_assembly {
    export const getDataFromApiProvider = apiQueryProvider('https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary');
}

namespace struct_ref_domain {
    export const getDataFromApiProvider = apiQueryProvider('https://www.ebi.ac.uk/pdbe/api/mappings/sequence_domains');
}

namespace residuewise_outlier_summary {
    const json = new Map<string, any>();

    export async function getDataFromAggregateFile(model: Model) {
        const key = `${model.label[1]}${model.label[2]}`;
        if (!json.has(key)) {
            const fn = `e:/test/mol-star/model/props/${key}.json`;
            if (!fs.existsSync(fn)) json.set(key, { });
            // TODO: use async readFile?
            else json.set(key, JSON.parse(fs.readFileSync(fn, 'utf8')));
        }
        return json.get(key)![model.label.toLowerCase()] || { };
    }

    export const getDataFromApiProvider = apiQueryProvider('https://www.ebi.ac.uk/pdbe/api/validation/residuewise_outlier_summary/entry');
}

function apiQueryProvider(urlPrefix: string) {
    return (cache: any) => {
        const cacheKey = UUID.create();
        return async (model: Model) => {
            if (cache[cacheKey]) return cache[cacheKey];
            const rawData = await fetchRetry(`${urlPrefix}/${model.label.toLowerCase()}`, 1500, 5);
            const json = (await rawData.json())[model.label.toLowerCase()] || { };
            cache[cacheKey] = json;
            return json;
        }
    }
}
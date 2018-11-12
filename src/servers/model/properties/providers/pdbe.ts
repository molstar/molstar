/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'fs'
import * as path from 'path'
import { Model } from 'mol-model/structure';
import { StructureQualityReport } from 'mol-model-props/pdbe/structure-quality-report';
import { fetchRetry } from '../../utils/fetch-retry';
import { UUID } from 'mol-util';
import { PDBePreferredAssembly } from 'mol-model-props/pdbe/preferred-assembly';
import { PDBeStructRefDomain } from 'mol-model-props/pdbe/struct-ref-domain';
import { AttachModelProperty } from '../../property-provider';
import { ConsoleLogger } from 'mol-util/console-logger';

export const PDBe_structureQualityReport: AttachModelProperty = ({ model, params, cache }) => {
    const PDBe_apiSourceJson = useFileSource(params)
        ? residuewise_outlier_summary.getDataFromAggregateFile(getFilePrefix(params, 'residuewise_outlier_summary'))
        : apiQueryProvider(getApiUrl(params, 'residuewise_outlier_summary', 'https://www.ebi.ac.uk/pdbe/api/validation/residuewise_outlier_summary/entry'), cache);

    return StructureQualityReport.attachFromCifOrApi(model, { PDBe_apiSourceJson });
}

export const PDBe_preferredAssembly: AttachModelProperty = ({ model, params, cache }) => {
    const PDBe_apiSourceJson = apiQueryProvider(getApiUrl(params, 'preferred_assembly', 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary'), cache);
    return PDBePreferredAssembly.attachFromCifOrApi(model, { PDBe_apiSourceJson });
}

export const PDBe_structRefDomain: AttachModelProperty = ({ model, params, cache }) => {
    const PDBe_apiSourceJson = apiQueryProvider(getApiUrl(params, 'struct_ref_domain', 'https://www.ebi.ac.uk/pdbe/api/mappings/sequence_domains'), cache);
    return PDBeStructRefDomain.attachFromCifOrApi(model, { PDBe_apiSourceJson });
}

namespace residuewise_outlier_summary {
    const json = new Map<string, any>();
    export function getDataFromAggregateFile(pathPrefix: string) {
        // This is for "testing" purposes and should probably only read
        // a single file with the appropriate prop in the "production" version.
        return async (model: Model) => {
            const key = `${model.label[1]}${model.label[2]}`;
            if (!json.has(key)) {
                const fn = path.join(pathPrefix, `${key}.json`);
                if (!fs.existsSync(fn)) json.set(key, { });
                // TODO: use async readFile?
                else json.set(key, JSON.parse(fs.readFileSync(fn, 'utf8')));
            }
            return json.get(key)![model.label.toLowerCase()] || { };
        }
    }
}

function getApiUrl(params: any, name: string, fallback: string) {
    const url = getParam<string>(params, 'PDBe', 'API', name);
    if (!url) return fallback;
    if (url[url.length - 1] === '/') return url.substring(0, url.length - 1);
    return url;
}

function getFilePrefix(params: any, name: string) {
    const ret = getParam<string>(params, 'PDBe', 'File', name);
    if (!ret) throw new Error(`PDBe file prefix '${name}' not set!`);
    return ret;
}

function useFileSource(params: any) {
    return !!getParam<boolean>(params, 'PDBe', 'UseFileSource')
}

function getParam<T>(params: any, ...path: string[]): T | undefined {
    try {
        let current = params;
        for (const p of path) {
            if (typeof current === 'undefined') return;
            current = current[p];
        }
        return current;
    } catch (e) {
        ConsoleLogger.error('Config', `Unable to retrieve property ${path.join('.')} from ${JSON.stringify(params)}`);
    }
}


function apiQueryProvider(urlPrefix: string, cache: any) {
    const cacheKey = UUID.create22();
    return async (model: Model) => {
        try {
            if (cache[cacheKey]) return cache[cacheKey];
            const rawData = await fetchRetry(`${urlPrefix}/${model.label.toLowerCase()}`, 1500, 5);
            // TODO: is this ok?
            if (rawData.status !== 200) return { };
            const json = (await rawData.json())[model.label.toLowerCase()] || { };
            cache[cacheKey] = json;
            return json;
        } catch (e) {
            // TODO: handle better
            ConsoleLogger.warn('Props', `Count not retrieve prop @${`${urlPrefix}/${model.label.toLowerCase()}`}`);
            return { };
        }
    }
}
/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Model } from 'mol-model/structure';
import Config from './config';

export type ModelPropertiesProvider = (model: Model, cache: object) => Promise<any>[]

export function createModelPropertiesProviderFromConfig(): ModelPropertiesProvider {
    return createModelPropertiesProviderFromSources(Config.customPropertyProviders);
}

export function createModelPropertiesProviderFromSources(sources: string[]): ModelPropertiesProvider {
    if (!sources || sources.length === 0) return () => [];

    const ps: ModelPropertiesProvider[] = [];
    for (const p of sources) {
        ps.push(require(p).attachModelProperties);
    }

    return (model, cache) => {
        const ret: Promise<any>[] = [];
        for (const p of ps) {
            for (const e of p(model, cache)) ret.push(e);
        }
        return ret;
    }
}

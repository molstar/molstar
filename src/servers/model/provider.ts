/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Model } from 'mol-model/structure';
import Config from './config';

export type ModelPropertiesProvider = (model: Model) => Promise<any>[]

export function createModelPropertiesProviderFromConfig(): ModelPropertiesProvider {
    if (!Config.customPropertyProviders || Config.customPropertyProviders.length === 0) return () => [];

    const ps: ModelPropertiesProvider[] = [];
    for (const p of Config.customPropertyProviders) {
        ps.push(require(p).attachModelProperties);
    }

    return model => {
        const ret: Promise<any>[] = [];
        for (const p of ps) {
            for (const e of p(model)) ret.push(e);
        }
        return ret;
    }
}


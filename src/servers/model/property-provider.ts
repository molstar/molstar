/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'fs'
import { Model } from 'mol-model/structure';
import Config from './config';
import { ConsoleLogger } from 'mol-util/console-logger';

export interface ModelPropertyProviderConfig {
    sources: string[],
    params?: { [name: string]: any }
}

export type AttachModelProperty = (args: { model: Model, params: any, cache: any }) => Promise<any>
export type AttachModelProperties = (args: { model: Model, params: any, cache: any }) => Promise<any>[]
export type ModelPropertiesProvider = (model: Model, cache: any) => Promise<any>[]

export function createModelPropertiesProviderFromConfig(): ModelPropertiesProvider {
    return createModelPropertiesProvider(Config.customProperties);
}

export function createModelPropertiesProvider(configOrPath: ModelPropertyProviderConfig | string | undefined): ModelPropertiesProvider {
    let config: ModelPropertyProviderConfig;
    if (typeof configOrPath === 'string') {
        try {
            config = JSON.parse(fs.readFileSync(configOrPath, 'utf8'));
        } catch {
            ConsoleLogger.error('Config', `Could not read property provider config file '${configOrPath}', ignoring.`);
            return () => [];
        }
    } else {
        config = configOrPath!;
    }

    if (!config || !config.sources || config.sources.length === 0) return () => [];

    const ps: AttachModelProperties[] = [];
    for (const p of config.sources) {
        ps.push(require(p).attachModelProperties);
    }

    return (model, cache) => {
        const ret: Promise<any>[] = [];
        for (const p of ps) {
            for (const e of p({ model, cache, params: config.params })) ret.push(e);
        }
        return ret;
    }
}

/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'fs';
import { Model } from '../../mol-model/structure';
import { ModelServerConfig as Config } from './config';
import { ConsoleLogger } from '../../mol-util/console-logger';

// TODO enable dynamic imports again
import * as pdbeProps from './properties/pdbe';
import * as wwpdbProps from './properties/wwpdb';

const attachModelProperties: { [k: string]: AttachModelProperties } = {
    pdbe: pdbeProps.attachModelProperties,
    wwpdb: wwpdbProps.attachModelProperties
};

export interface ModelPropertyProviderConfig {
    sources: string[],
    params?: { [name: string]: any }
}

export type AttachModelProperty = (args: { model: Model, params: any, cache: any }) => Promise<any>
export type AttachModelProperties = (args: { model: Model, params: any, cache: any }) => Promise<any>[]
export type ModelPropertiesProvider = (model: Model, cache: any) => Promise<any>[]

export function createModelPropertiesProviderFromConfig() {
    return createModelPropertiesProvider(Config.customProperties);
}

export function createModelPropertiesProvider(configOrPath: ModelPropertyProviderConfig | string | undefined): ModelPropertiesProvider | undefined {
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

    if (!config || !config.sources || config.sources.length === 0) return void 0;

    const ps: AttachModelProperties[] = [];
    for (const p of config.sources) {
        if (p in attachModelProperties) {
            ps.push(attachModelProperties[p]);
        } else {
            ConsoleLogger.error('Config', `Could not find property provider '${p}', ignoring.`);
        }
    }

    return (model, cache) => {
        const ret: Promise<any>[] = [];
        for (const p of ps) {
            for (const e of p({ model, cache, params: config.params })) ret.push(e);
        }
        return ret;
    };
}

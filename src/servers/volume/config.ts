/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import * as argparse from 'argparse';
import { ObjectKeys } from '../../mol-util/type-helpers';
import { VOLUME_SERVER_HEADER } from './server/version';
import * as fs from 'fs';

const DefaultServerConfig = {
    apiPrefix: '/VolumeServer',
    defaultPort: 1337,
    shutdownTimeoutMinutes: 24 * 60, /* a day */
    shutdownTimeoutVarianceMinutes: 60,
    idMap: [] as [string, string][],
    healthCheckPath: [] as string[],
};

function addLimitsArgs(parser: argparse.ArgumentParser) {
    parser.add_argument('--maxRequestBlockCount', {
        default: DefaultLimitsConfig.maxRequestBlockCount,
        metavar: 'COUNT',
        help: `Maximum number of blocks that could be read in 1 query.
This is somewhat tied to the maxOutputSizeInVoxelCountByPrecisionLevel
in that the <maximum number of voxel> = maxRequestBlockCount * <block size>^3.
The default block size is 96 which corresponds to 28,311,552 voxels with 32 max blocks.`
    });
    parser.add_argument('--maxFractionalBoxVolume', {
        default: DefaultLimitsConfig.maxFractionalBoxVolume,
        metavar: 'VOLUME',
        help: `The maximum fractional volume of the query box (to prevent queries that are too big).`
    });
    parser.add_argument('--maxOutputSizeInVoxelCountByPrecisionLevel', {
        nargs: '+',
        default: DefaultLimitsConfig.maxOutputSizeInVoxelCountByPrecisionLevel,
        metavar: 'LEVEL',
        help: `What is the (approximate) maximum desired size in voxel count by precision level
Rule of thumb: <response gzipped size> in [<voxel count> / 8, <voxel count> / 4].
The maximum number of voxels is tied to maxRequestBlockCount.`
    });
}

function addServerArgs(parser: argparse.ArgumentParser) {
    parser.add_argument('--apiPrefix', {
        default: DefaultServerConfig.apiPrefix,
        metavar: 'PREFIX',
        help: `Specify the prefix of the API, i.e. <host>/<apiPrefix>/<API queries>`
    });
    parser.add_argument('--defaultPort', {
        default: DefaultServerConfig.defaultPort,
        metavar: 'PORT',
        type: 'int',
        help: `Specify the port the server is running on`
    });

    parser.add_argument('--shutdownTimeoutMinutes', {
        default: DefaultServerConfig.shutdownTimeoutMinutes,
        metavar: 'TIME',
        type: 'int',
        help: `0 for off, server will shut down after this amount of minutes.`
    });
    parser.add_argument('--shutdownTimeoutVarianceMinutes', {
        default: DefaultServerConfig.shutdownTimeoutVarianceMinutes,
        metavar: 'VARIANCE',
        type: 'int',
        help: `modifies the shutdown timer by +/- timeoutVarianceMinutes (to avoid multiple instances shutting at the same time)`
    });
    parser.add_argument('--idMap', {
        nargs: 2,
        action: 'append',
        metavar: ['TYPE', 'PATH'] as any,
        help: [
            'Map `id`s for a `type` to a file path.',
            'Example: x-ray \'../../data/mdb/xray/${id}-ccp4.mdb\'',
            '',
            '  - JS expressions can be used inside ${}, e.g. \'${id.substr(1, 2)}/${id}.mdb\'',
            '  - Can be specified multiple times.',
            '  - The `TYPE` variable (e.g. `x-ray`) is arbitrary and depends on how you plan to use the server.',
            '    By default, Mol* Viewer uses `x-ray` and `em`, but any particular use case may vary. '
        ].join('\n'),
    });
    parser.add_argument('--healthCheckPath', {
        default: DefaultServerConfig.healthCheckPath,
        action: 'append',
        metavar: 'PATH',
        help: `File path(s) to use for health-checks. Will test if all files are accessible and report a failed health-check if that's not the case.`,
    });
}

function addJsonConfigArgs(parser: argparse.ArgumentParser) {
    parser.add_argument('--cfg', {
        help: [
            'JSON config file path',
            'If a property is not specified, cmd line param/OS variable/default value are used.'
        ].join('\n'),
        required: false,
    });
    parser.add_argument('--printCfg', { help: 'Print current config for validation and exit.', required: false, action: 'store_true' });
    parser.add_argument('--cfgTemplate', { help: 'Prints default JSON config template to be modified and exits.', required: false, action: 'store_true' });
}

export interface ServerJsonConfig {
    cfg?: string,
    printCfg?: any,
    cfgTemplate?: any
}

export type ServerConfig = typeof DefaultServerConfig
export const ServerConfig = { ...DefaultServerConfig };

function setServerConfig(config: ServerConfig) {
    for (const k of ObjectKeys(ServerConfig)) {
        if (config[k] !== void 0) (ServerConfig as any)[k] = config[k];
    }
}

function validateServerConfig() {
    if (!ServerConfig.idMap || ServerConfig.idMap.length === 0) {
        throw new Error(`Please provide 'idMap' configuration. See [-h] for available options.`);
    }
}

const DefaultLimitsConfig = {
    maxRequestBlockCount: 32,
    maxFractionalBoxVolume: 1024,
    maxOutputSizeInVoxelCountByPrecisionLevel: [
        0.5 * 1024 * 1024, // ~ 80*80*80
        1 * 1024 * 1024,
        2 * 1024 * 1024,
        4 * 1024 * 1024,
        8 * 1024 * 1024,
        16 * 1024 * 1024, // ~ 256*256*256
        24 * 1024 * 1024
    ]
};
export type LimitsConfig = typeof DefaultLimitsConfig
export const LimitsConfig = { ...DefaultLimitsConfig };

function setLimitsConfig(config: LimitsConfig) {
    for (const k of ObjectKeys(LimitsConfig)) {
        if (config[k] !== void 0) (LimitsConfig as any)[k] = config[k];
    }
}

type FullServerConfig = ServerConfig & LimitsConfig

function setConfig(config: FullServerConfig) {
    setServerConfig(config);
    setLimitsConfig(config);
}

const ServerConfigTemplate: FullServerConfig = {
    ...DefaultServerConfig,
    idMap: [
        ['x-ray', './path-to-xray-data/${id.substr(1, 2)}/${id}.mdb'],
        ['em', './path-to-em-data/emd-${id}.mdb']
    ] as [string, string][],
    ...DefaultLimitsConfig
};

export function configureServer() {
    const parser = new argparse.ArgumentParser({
        // version: VOLUME_SERVER_VERSION,
        add_help: true,
        description: VOLUME_SERVER_HEADER
    });
    addJsonConfigArgs(parser);
    addServerArgs(parser);
    addLimitsArgs(parser);
    const config = parser.parse_args() as FullServerConfig & ServerJsonConfig;

    if (!!config.cfgTemplate) {
        console.log(JSON.stringify(ServerConfigTemplate, null, 2));
        process.exit(0);
    }

    try {
        setConfig(config); // sets the config for global use

        if (config.cfg) {
            const cfg = JSON.parse(fs.readFileSync(config.cfg, 'utf8')) as FullServerConfig;
            setConfig(cfg);
        }

        if (config.printCfg) {
            console.log(JSON.stringify({ ...ServerConfig, ...LimitsConfig }, null, 2));
            process.exit(0);
        }

        validateServerConfig();
    } catch (e) {
        console.error('' + e);
        process.exit(1);
    }
}

export function configureLocal() {
    const parser = new argparse.ArgumentParser({
        // version: VOLUME_SERVER_VERSION,
        add_help: true,
        description: VOLUME_SERVER_HEADER
    });
    parser.add_argument('--jobs', { help: `Path to a JSON file with job specification.`, required: false });
    parser.add_argument('--jobsTemplate', { help: 'Print example template for jobs.json and exit.', required: false, action: 'store_true' });
    addJsonConfigArgs(parser);
    addLimitsArgs(parser);

    const config = parser.parse_args() as LimitsConfig & ServerJsonConfig;

    if (config.cfgTemplate) {
        console.log(JSON.stringify(DefaultLimitsConfig, null, 2));
        process.exit(0);
    }

    try {
        setLimitsConfig(config); // sets the config for global use

        if (config.cfg) {
            const cfg = JSON.parse(fs.readFileSync(config.cfg, 'utf8')) as FullServerConfig;
            setLimitsConfig(cfg);
        }

        if (config.printCfg) {
            console.log(JSON.stringify(LimitsConfig, null, 2));
            process.exit(0);
        }

        return config as LimitsConfig & { jobs: string, jobsTemplate: any };
    } catch (e) {
        console.error('' + e);
        process.exit(1);
    }
}
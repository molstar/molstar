/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as argparse from 'argparse'

export function addLimitsArgs(parser: argparse.ArgumentParser) {
    parser.addArgument([ '--maxRequestBlockCount' ], {
        defaultValue: DefaultLimitsConfig.maxRequestBlockCount,
        metavar: 'COUNT',
        help: `Maximum number of blocks that could be read in 1 query.
This is somewhat tied to the maxOutputSizeInVoxelCountByPrecisionLevel
in that the <maximum number of voxel> = maxRequestBlockCount * <block size>^3.
The default block size is 96 which corresponds to 28,311,552 voxels with 32 max blocks.`
    });
    parser.addArgument([ '--maxFractionalBoxVolume' ], {
        defaultValue: DefaultLimitsConfig.maxFractionalBoxVolume,
        metavar: 'VOLUME',
        help: `The maximum fractional volume of the query box (to prevent queries that are too big).`
    });
    parser.addArgument([ '--maxOutputSizeInVoxelCountByPrecisionLevel' ], {
        nargs: '+',
        defaultValue: DefaultLimitsConfig.maxOutputSizeInVoxelCountByPrecisionLevel,
        metavar: 'LEVEL',
        help: `What is the (approximate) maximum desired size in voxel count by precision level
Rule of thumb: <response gzipped size> in [<voxel count> / 8, <voxel count> / 4].
The maximum number of voxels is tied to maxRequestBlockCount.`
    });
}

export function addServerArgs(parser: argparse.ArgumentParser) {
    parser.addArgument([ '--apiPrefix' ], {
        defaultValue: DefaultServerConfig.apiPrefix,
        metavar: 'PREFIX',
        help: `Specify the prefix of the API, i.e. <host>/<apiPrefix>/<API queries>`
    });
    parser.addArgument([ '--defaultPort' ], {
        defaultValue: DefaultServerConfig.defaultPort,
        metavar: 'PORT',
        help: `Specify the prefix of the API, i.e. <host>/<apiPrefix>/<API queries>`
    });

    parser.addArgument([ '--shutdownTimeoutMinutes' ], {
        defaultValue: DefaultServerConfig.shutdownTimeoutMinutes,
        metavar: 'TIME',
        help: `0 for off, server will shut down after this amount of minutes.`
    });
    parser.addArgument([ '--shutdownTimeoutVarianceMinutes' ], {
        defaultValue: DefaultServerConfig.shutdownTimeoutVarianceMinutes,
        metavar: 'VARIANCE',
        help: `modifies the shutdown timer by +/- timeoutVarianceMinutes (to avoid multiple instances shutting at the same time)`
    });
    parser.addArgument([ '--idMap' ], {
        nargs: 2,
        action: 'append',
        metavar: ['TYPE', 'PATH'] as any,
        help: [
            'Map `id`s for a `type` to a file path.',
            'Example: x-ray \'../../data/mdb/xray/${id}-ccp4.mdb\'',
            'Note: Can be specified multiple times.'
        ].join('\n'),
    });
}

const DefaultServerConfig = {
    apiPrefix: '/VolumeServer',
    defaultPort: 1337,
    shutdownTimeoutMinutes: 24 * 60, /* a day */
    shutdownTimeoutVarianceMinutes: 60,
    idMap: [] as [string, string][]
}
export type ServerConfig = typeof DefaultServerConfig
export const ServerConfig = { ...DefaultServerConfig }
export function setServerConfig(config: ServerConfig) {
    for (const name in DefaultServerConfig) {
        ServerConfig[name as keyof ServerConfig] = config[name as keyof ServerConfig]
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
}
export type LimitsConfig = typeof DefaultLimitsConfig
export const LimitsConfig = { ...DefaultLimitsConfig }
export function setLimitsConfig(config: LimitsConfig) {
    for (const name in DefaultLimitsConfig) {
        LimitsConfig[name as keyof LimitsConfig] = config[name as keyof LimitsConfig]
    }
}

export function setConfig(config: ServerConfig & LimitsConfig) {
    setServerConfig(config)
    setLimitsConfig(config)
}
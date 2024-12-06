/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import * as argparse from 'argparse';
import { VERSION } from './version';
import { ObjectKeys } from '../../mol-util/type-helpers';

const DefaultMembraneServerConfig = {
    /**
     * Node (V8) sometimes exhibits GC related issues  that significantly slow down the execution
     * (https://github.com/nodejs/node/issues/8670).
     *
     * Therefore an option is provided that automatically shuts down the server.
     * For this to work, the server must be run using a deamon (i.e. forever.js on Linux
     * or IISnode on Windows) so that the server is automatically restarted when the shutdown happens.
     */

    // 0 for off, server will shut down after this amount of minutes.
    shutdownTimeoutMinutes: 24 * 60, /* a day */
    // modifies the shutdown timer by +/- timeoutVarianceMinutes (to avoid multiple instances shutting at the same time)
    shutdownTimeoutVarianceMinutes: 60,

    defaultPort: 1340,

    /**
     * Specify the prefix of the API, i.e.
     * <host>/<apiPrefix>/<API queries>
     */
    apiPrefix: '/MembraneServer',

    /**
     * The maximum number of ms the server spends on a request
     */
    requestTimeoutMs: 60 * 1000,

    /**
     * Default URL from which BinaryCIF data will be pulled.
     */
    bcifSource: (id: string) => `https://models.rcsb.org/${id}.bcif`,
};

function addServerArgs(parser: argparse.ArgumentParser) {
    parser.add_argument('--apiPrefix', {
        default: DefaultMembraneServerConfig.apiPrefix,
        metavar: 'PREFIX',
        help: `Specify the prefix of the API, i.e. <host>/<apiPrefix>/<API queries>`
    });
    parser.add_argument('--defaultPort', {
        default: DefaultMembraneServerConfig.defaultPort,
        metavar: 'PORT',
        type: 'int',
        help: `Specify the port the server is running on`
    });
    parser.add_argument('--shutdownTimeoutMinutes', {
        default: DefaultMembraneServerConfig.shutdownTimeoutMinutes,
        metavar: 'TIME',
        type: 'int',
        help: `0 for off, server will shut down after this amount of minutes.`
    });
    parser.add_argument('--shutdownTimeoutVarianceMinutes', {
        default: DefaultMembraneServerConfig.shutdownTimeoutVarianceMinutes,
        metavar: 'VARIANCE',
        type: 'int',
        help: `modifies the shutdown timer by +/- timeoutVarianceMinutes (to avoid multiple instances shutting at the same time)`
    });
    parser.add_argument('--bcifSource', {
        default: DefaultMembraneServerConfig.bcifSource,
        metavar: 'DEFAULT_SOURCE',
        help: `Where 3D structure data is loaded from.`
    });
}

export type MembraneServerConfig = typeof DefaultMembraneServerConfig
export const MembraneServerConfig = { ...DefaultMembraneServerConfig };

function setConfig(config: MembraneServerConfig) {
    for (const k of ObjectKeys(MembraneServerConfig)) {
        if (config[k] !== void 0) (MembraneServerConfig as any)[k] = config[k];
    }
}

function parseConfigArguments() {
    const parser = new argparse.ArgumentParser({
        add_help: true,
        description: `Mol* MembraneServer ${VERSION}, (c) 2024, Mol* contributors`
    });
    addServerArgs(parser);
    return parser.parse_args() as MembraneServerConfig;
}

export function configureServer() {
    const config = parseConfigArguments();
    setConfig(config); // sets the config for global use
}
﻿/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import * as argparse from 'argparse';
import { ObjectKeys } from '../../mol-util/type-helpers';
import { VERSION } from './version';
import * as fs from 'fs';
import { ModelPropertyProviderConfig } from './property-provider';

const DefaultModelServerConfig = {
    // 0 for off
    cacheMaxSizeInBytes: 2 * 1014 * 1024 * 1024, // 2 GB
    cacheEntryTimeoutMs: 10 * 60 * 1000, // 10 minutes

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

    defaultPort: 1337,

    /**
     * Specify the prefix of the API, i.e.
     * <host>/<apiPrefix>/<API queries>
     */
    apiPrefix: '/ModelServer',

    /**
     * The maximum time the server dedicates to executing a query.
     * Does not include the time it takes to read and export the data.
     */
    queryTimeoutMs: 5 * 1000,

    /**
     * The maximum number of ms the server spends on a request
     */
    requestTimeoutMs: 60 * 1000,

    /** Maximum number of requests before "server busy" */
    maxQueueLength: 30,

    /** The maximum number of queries allowed by the query-many at a time */
    maxQueryManyQueries: 50,

    /**
     * Provide a property config or a path a JSON file with the config.
     */
    // TODO: finish customProperty support
    customProperties: <ModelPropertyProviderConfig | string>{
        sources: [
            // 'pdbe',
            // 'rcsb',
            // 'wwpdb'
        ],
        params: {
            // PDBe: {
            //     UseFileSource: false,
            //     API: {
            //         residuewise_outlier_summary: 'https://www.ebi.ac.uk/pdbe/api/validation/residuewise_outlier_summary/entry',
            //         preferred_assembly: 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary',
            //         struct_ref_domain: 'https://www.ebi.ac.uk/pdbe/api/mappings/sequence_domains'
            //     },
            //     File: {
            //         residuewise_outlier_summary: 'e:/test/mol-star/model/props/'
            //     }
            // },
            // RCSB: {
            //     API: {
            //         assembly_symmetry: 'https://rest-staging.rcsb.org/graphql'
            //     }
            // },
            // wwPDB: {
            //     chemCompBondTablePath: '',
            //     chemCompAtomTablePath: ''
            // }
        }
    },


    /**
     * Default source for fileMapping.
     */
    defaultSource: 'pdb-cif' as string,

    /**
     * Maps a request identifier to either:
     * - filename [source, mapping]
     * - URI [source, mapping, format]
     *
     * Mapping is provided 'source' and 'id' variables to interpolate.
     *
     * /static query uses 'pdb-cif' and 'pdb-bcif' source names.
     */
    sourceMap: [
        ['pdb-cif', 'e:/test/quick/${id}_updated.cif'],
        // ['pdb-bcif', 'e:/test/quick/${id}.bcif'],
    ] as ([string, string] | [string, string, ModelServerFetchFormats])[],

    /**
     * Optionally point to files. The service health-check will assert that all are readable and fail otherwise.
     */
    healthCheckPath: [] as string[],
};

export const ModelServerFetchFormats = ['cif', 'bcif', 'cif.gz', 'bcif.gz'] as const;
export type ModelServerFetchFormats = (typeof ModelServerFetchFormats)[number]

export let mapSourceAndIdToFilename: (source: string, id: string) => [string, ModelServerFetchFormats] = () => {
    throw new Error('call setupConfig & validateConfigAndSetupSourceMap to initialize this function');
};

function addServerArgs(parser: argparse.ArgumentParser) {
    parser.add_argument('--apiPrefix', {
        default: DefaultModelServerConfig.apiPrefix,
        metavar: 'PREFIX',
        help: `Specify the prefix of the API, i.e. <host>/<apiPrefix>/<API queries>`
    });
    parser.add_argument('--defaultPort', {
        default: DefaultModelServerConfig.defaultPort,
        metavar: 'PORT',
        type: 'int',
        help: `Specify the port the server is running on`
    });
    parser.add_argument('--cacheMaxSizeInBytes', {
        default: DefaultModelServerConfig.cacheMaxSizeInBytes,
        metavar: 'CACHE_SIZE',
        type: 'int',
        help: `0 for off.`
    });
    parser.add_argument('--cacheEntryTimeoutMs', {
        default: DefaultModelServerConfig.cacheEntryTimeoutMs,
        metavar: 'CACHE_TIMEOUT',
        type: 'int',
        help: `Specify in ms how long to keep entries in cache.`
    });
    parser.add_argument('--requestTimeoutMs', {
        default: DefaultModelServerConfig.requestTimeoutMs,
        metavar: 'REQUEST_TIMEOUT',
        type: 'int',
        help: `The maximum number of ms the server spends on a request.`
    });
    parser.add_argument('--queryTimeoutMs', {
        default: DefaultModelServerConfig.queryTimeoutMs,
        metavar: 'QUERY_TIMEOUT',
        type: 'int',
        help: `The maximum time the server dedicates to executing a query in ms.\nDoes not include the time it takes to read and export the data.`
    });
    parser.add_argument('--shutdownTimeoutMinutes', {
        default: DefaultModelServerConfig.shutdownTimeoutMinutes,
        metavar: 'TIME',
        type: 'int',
        help: `0 for off, server will shut down after this amount of minutes.`
    });
    parser.add_argument('--shutdownTimeoutVarianceMinutes', {
        default: DefaultModelServerConfig.shutdownTimeoutVarianceMinutes,
        metavar: 'VARIANCE',
        type: 'int',
        help: `modifies the shutdown timer by +/- timeoutVarianceMinutes (to avoid multiple instances shutting at the same time)`
    });
    parser.add_argument('--maxQueryManyQueries', {
        default: DefaultModelServerConfig.maxQueryManyQueries,
        metavar: 'QUERY_MANY_LIMIT',
        type: 'int',
        help: `maximum number of queries allowed by the query-many at a time`
    });
    parser.add_argument('--defaultSource', {
        default: DefaultModelServerConfig.defaultSource,
        metavar: 'DEFAULT_SOURCE',
        help: `modifies which 'sourceMap' source to use by default`
    });
    parser.add_argument('--sourceMap', {
        nargs: 2,
        action: 'append',
        metavar: ['SOURCE', 'PATH'] as any,
        help: [
            'Map `id`s for a `source` to a file path.',
            'Example: pdb-bcif \'../../data/bcif/${id}.bcif\' ',
            'JS expressions can be used inside ${}, e.g. \'${id.substr(1, 2)}/${id}.mdb\'',
            'Can be specified multiple times.',
            'The `SOURCE` variable (e.g. `pdb-bcif`) is arbitrary and depends on how you plan to use the server.',
            `Supported formats: ${ModelServerFetchFormats.join(', ')}`
        ].join('\n'),
    });
    parser.add_argument('--sourceMapUrl', {
        nargs: 3,
        action: 'append',
        metavar: ['SOURCE', 'PATH', 'SOURCE_MAP_FORMAT'] as any,
        help: [
            'Same as --sourceMap but for URL. \'--sourceMapUrl src url format\'',
            'Example: \'pdb-cif "https://www.ebi.ac.uk/pdbe/entry-files/download/${id}_updated.cif" cif\'',
            `Supported formats: ${ModelServerFetchFormats.join(', ')}`
        ].join('\n'),
    });
    parser.add_argument('--healthCheckPath', {
        default: DefaultModelServerConfig.healthCheckPath,
        action: 'append',
        metavar: 'PATH',
        help: `File path(s) to use for health-checks. Will test if all files are accessible and report a failed health-check if that's not the case.`,
    });
}

export type ModelServerConfig = typeof DefaultModelServerConfig
export const ModelServerConfig = { ...DefaultModelServerConfig };

export const ModelServerConfigTemplate: ModelServerConfig = {
    ...DefaultModelServerConfig,
    defaultSource: 'pdb-bcif',
    sourceMap: [
        ['pdb-bcif', './path-to-binary-cif/${id.substr(1, 2)}/${id}.bcif'],
        ['pdb-cif', './path-to-text-cif/${id.substr(1, 2)}/${id}.cif'],
        ['pdb-updated', 'https://www.ebi.ac.uk/pdbe/entry-files/download/${id}_updated.cif', 'cif']
    ] as [string, string][]
};

// TODO: include once properly supported
// cast to any to prevent "The operand of a 'delete' operator must be optional"
delete (ModelServerConfigTemplate as any).customProperties;

interface ServerJsonConfig {
    cfg?: string,
    printCfg?: any,
    cfgTemplate?: any
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

function setConfig(config: ModelServerConfig) {
    for (const k of ObjectKeys(ModelServerConfig)) {
        if (config[k] !== void 0) (ModelServerConfig as any)[k] = config[k];
    }

    if ((config as any).sourceMapUrl) {
        if (!ModelServerConfig.sourceMap) ModelServerConfig.sourceMap = [];
        ModelServerConfig.sourceMap.push(...(config as any).sourceMapUrl);
    }
}

function validateConfigAndSetupSourceMap() {
    if (!ModelServerConfig.sourceMap || ModelServerConfig.sourceMap.length === 0) {
        throw new Error(`Please provide 'sourceMap' configuration. See [-h] for available options.`);
    }

    mapSourceAndIdToFilename = new Function('source', 'id', [
        'switch (source.toLowerCase()) {',
        ...ModelServerConfig.sourceMap.map(([source, path, format]) => `case '${source.toLowerCase()}': return [\`${path}\`, '${format}'];`),
        '}',
    ].join('\n')) as any;
}

function parseConfigArguments() {
    const parser = new argparse.ArgumentParser({
        add_help: true,
        description: `ModelServer ${VERSION}, (c) 2018-2020, Mol* contributors`
    });
    addJsonConfigArgs(parser);
    addServerArgs(parser);
    return parser.parse_args() as ModelServerConfig & ServerJsonConfig;
}

export function configureServer() {
    const config = parseConfigArguments();

    if (!!config.cfgTemplate) {
        console.log(JSON.stringify(ModelServerConfigTemplate, null, 2));
        process.exit(0);
    }

    try {
        setConfig(config); // sets the config for global use

        if (config.cfg) {
            const cfg = JSON.parse(fs.readFileSync(config.cfg, 'utf8')) as ModelServerConfig;
            setConfig(cfg);
        }

        if (!!config.printCfg) {
            console.log(JSON.stringify(ModelServerConfig, null, 2));
            process.exit(0);
        }

        validateConfigAndSetupSourceMap();
    } catch (e) {
        console.error('' + e);
        process.exit(1);
    }
}
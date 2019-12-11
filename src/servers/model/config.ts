/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

const DefaultModelServerConfig = {
    /**
     * Determine if and how long to cache entries after a request.
     */
    cacheParams: {
        useCache: true,
        maxApproximateSizeInBytes: 2 * 1014 * 1024 * 1024, // 2 GB
        entryTimeoutInMs: 10 * 60 * 1000 // 10 minutes
    },

    /**
     * Node (V8) sometimes exhibits GC related issues  that significantly slow down the execution
     * (https://github.com/nodejs/node/issues/8670).
     *
     * Therefore an option is provided that automatically shuts down the server.
     * For this to work, the server must be run using a deamon (i.e. forever.js on Linux
     * or IISnode on Windows) so that the server is automatically restarted when the shutdown happens.
     */
    shutdownParams: {
        // 0 for off, server will shut down after this amount of minutes.
        timeoutMinutes: 24 * 60 /* a day */,

        // modifies the shutdown timer by +/- timeoutVarianceMinutes (to avoid multiple instances shutting at the same time)
        timeoutVarianceMinutes: 60
    },

    defaultPort: 1337,

    /**
     * Specify the prefix of the API, i.e.
     * <host>/<apiPrefix>/<API queries>
     */
    appPrefix: '/ModelServer',

    /**
     * The maximum time the server dedicates to executing a query.
     * Does not include the time it takes to read and export the data.
     */
    maxQueryTimeInMs: 5 * 1000,

    /** Maximum number of requests before "server busy" */
    maxQueueLength: 30,

    /**
     * Provide a property config or a path a JSON file with the config.
     */
    // TODO: finish customProperty support
    customProperties: <import('./property-provider').ModelPropertyProviderConfig | string>{
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
            //     chemCompBondTablePath: ''
            // }
        }
    },


    /**
     * Default source for fileMapping.
     */
    defaultSource: 'pdb-cif' as string,

    /**
     * Maps a request identifier to a filename given a 'source' and 'id' variables.
     *
     * /static query uses 'pdb-cif' and 'pdb-bcif' source names.
     */
    fileMapping: [
        ['pdb-cif', 'e:/test/quick/${id}_updated.cif'],
        // ['pdb-bcif', 'e:/test/quick/${id}.bcif'],
    ] as [string, string][]
};

export type ModelServerConfig = typeof DefaultModelServerConfig
export const ModelServerConfig = DefaultModelServerConfig

export let mapSourceAndIdToFilename: (source: string, id: string) => string = () => {
    throw new Error('call setupConfig to initialize this function');
}

export function setupConfig(cfg?: ModelServerConfig) {
    if (cfg) Object.assign(ModelServerConfig, cfg);
    if (!ModelServerConfig.fileMapping) return;

    mapSourceAndIdToFilename = new Function('source', 'id', [
        'switch (source.toLowerCase()) {',
        ...ModelServerConfig.fileMapping.map(([source, path]) => `case '${source.toLowerCase()}': return \`${path}\`;`),
        '}',
    ].join('\n')) as any;
}
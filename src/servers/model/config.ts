/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

const config = {
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
     * Paths (relative to the root directory of the model server) to JavaScript files that specify custom properties
     */
    customPropertyProviders: [
        './properties/pdbe',
        './properties/rcsb'
    ],

    /**
     * Maps a request identifier to a filename.
     *
     * @param source
     *   Source of the data.
     * @param id
     *   Id provided in the request.
     */
    mapFile(source: string, id: string) {
        switch (source.toLowerCase()) {
            // case 'pdb': return `e:/test/quick/${id}_updated.cif`;
            case 'pdb': return `e:/test/mol-star/model/out/${id}_updated.bcif`;
            default: return void 0;
        }
    }
};

export default config;
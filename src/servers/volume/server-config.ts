
const Config = {
    limits: {
        /**
         * Maximum number of blocks that could be read in 1 query.
         * This is somewhat tied to the maxOutputSizeInVoxelCountByPrecisionLevel
         * in that the <maximum number of voxel> = maxRequestBlockCount * <block size>^3.
         * The default block size is 96 which corresponds to 28,311,552 voxels with 32 max blocks.
         */
        maxRequestBlockCount: 32,

        /**
         * The maximum fractional volume of the query box (to prevent queries that are too big).
         */
        maxFractionalBoxVolume: 1024,

        /**
         * What is the (approximate) maximum desired size in voxel count by precision level
         * Rule of thumb: <response gzipped size> \in [<voxel count> / 8, <voxel count> / 4];
         *
         * The maximum number of voxels is tied to maxRequestBlockCount.
         */
        maxOutputSizeInVoxelCountByPrecisionLevel: [
            0.5 * 1024 * 1024, // ~ 80*80*80
            1 * 1024 * 1024,
            2 * 1024 * 1024,
            4 * 1024 * 1024,
            8 * 1024 * 1024,
            16 * 1024 * 1024, // ~ 256*256*256
            24 * 1024 * 1024
        ]
    },

    /**
     * Specify the prefix of the API, i.e.
     * <host>/<apiPrefix>/<API queries>
     */
    apiPrefix: '/VolumeServer',

    /**
     * If not specified otherwise by the 'port' environment variable, use this port.
     */
    defaultPort: 1337,

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

    /**
     * Maps a request identifier to a filename.
     * 
     * @param source 
     *   Source of the data.
     * @param id
     *   Id provided in the request. For xray, PDB id, for emd, EMDB id number. 
     */
    mapFile(source: string, id: string) {
        switch (source.toLowerCase()) {
            case 'x-ray': return `g:/test/mdb/xray-${id.toLowerCase()}.mdb`;
            case 'emd': return `g:/test/mdb/${id.toLowerCase()}.mdb`;
            default: return void 0;
        }
    }
}

export default Config;
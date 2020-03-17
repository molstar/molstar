/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as compression from 'compression'
import * as express from 'express'
import { ConsoleLogger } from '../../mol-util/console-logger'
import { PerformanceMonitor } from '../../mol-util/performance-monitor'
import { configureServer, ModelServerConfig as ServerConfig } from './config'
import { initWebApi } from './server/api-web'
import Version from './version'

function setupShutdown() {
    if (ServerConfig.shutdownTimeoutVarianceMinutes > ServerConfig.shutdownTimeoutMinutes) {
        ConsoleLogger.log('Server', 'Shutdown timeout variance is greater than the timer itself, ignoring.');
    } else {
        let tVar = 0;
        if (ServerConfig.shutdownTimeoutVarianceMinutes > 0) {
            tVar = 2 * (Math.random() - 0.5) * ServerConfig.shutdownTimeoutVarianceMinutes;
        }
        let tMs = (ServerConfig.shutdownTimeoutMinutes + tVar) * 60 * 1000;

        console.log(`----------------------------------------------------------------------------`);
        console.log(`  The server will shut down in ${PerformanceMonitor.format(tMs)} to prevent slow performance.`);
        console.log(`  Please make sure a daemon is running that will automatically restart it.`);
        console.log(`----------------------------------------------------------------------------`);
        console.log();

        setTimeout(() => {
            // if (WebApi.ApiState.pendingQueries > 0) {
            //     WebApi.ApiState.shutdownOnZeroPending = true;
            // } else {
            ConsoleLogger.log('Server', `Shut down due to timeout.`);
            process.exit(0);
            // }
        }, tMs);
    }
}

configureServer();

function startServer() {
    let app = express();
    app.use(compression(<any>{ level: 6, memLevel: 9, chunkSize: 16 * 16384, filter: () => true }));

    initWebApi(app);

    const port = process.env.port || ServerConfig.defaultPort;
    app.listen(port);

    console.log(`Mol* ModelServer ${Version}`);
    console.log(``);
    console.log(`The server is running on port ${port}.`);
    console.log(``);
}

startServer();

if (ServerConfig.shutdownTimeoutMinutes > 0) {
    setupShutdown();
}
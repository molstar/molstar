#!/usr/bin/env node
/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import compression from 'compression';
import express from 'express';
import { ConsoleLogger } from '../../mol-util/console-logger';
import { configureServer, ServerConfig } from './config';
import { State } from './server/state';
import { VOLUME_SERVER_HEADER } from './server/version';
import init from './server/web-api';


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
        console.log(`  The server will shut down in ${ConsoleLogger.formatTime(tMs)} to prevent slow performance.`);
        console.log(`  Please make sure a daemon is running that will automatically restart it.`);
        console.log(`----------------------------------------------------------------------------`);
        console.log();

        setTimeout(() => {
            if (State.pendingQueries > 0) {
                State.shutdownOnZeroPending = true;
            } else {
                ConsoleLogger.log('Server', `Shut down due to timeout.`);
                process.exit(0);
            }
        }, tMs);
    }
}

configureServer();

const port = process.env.port || ServerConfig.defaultPort;

const app = express();
app.use(compression({ level: 6, memLevel: 9, chunkSize: 16 * 16384, filter: () => true }));
init(app);

app.listen(port);

console.log(VOLUME_SERVER_HEADER);
console.log(``);
console.log(`The server is running on port ${port}.`);
console.log(``);

if (ServerConfig.shutdownTimeoutMinutes > 0) {
    setupShutdown();
}
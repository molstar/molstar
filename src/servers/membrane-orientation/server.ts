#!/usr/bin/env node
/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import compression from 'compression';
import express from 'express';
import { VERSION } from './version';
import { configureServer, MembraneServerConfig } from './config';
import { ConsoleLogger } from '../../mol-util/console-logger';
import { PerformanceMonitor } from '../../mol-util/performance-monitor';
import { initWebApi } from './web-api';

function setupShutdown() {
    if (MembraneServerConfig.shutdownTimeoutVarianceMinutes > MembraneServerConfig.shutdownTimeoutMinutes) {
        ConsoleLogger.log('Server', 'Shutdown timeout variance is greater than the timer itself, ignoring.');
    } else {
        let tVar = 0;
        if (MembraneServerConfig.shutdownTimeoutVarianceMinutes > 0) {
            tVar = 2 * (Math.random() - 0.5) * MembraneServerConfig.shutdownTimeoutVarianceMinutes;
        }
        const tMs = (MembraneServerConfig.shutdownTimeoutMinutes + tVar) * 60 * 1000;

        console.log(`----------------------------------------------------------------------------`);
        console.log(`  The server will shut down in ${PerformanceMonitor.format(tMs)} to prevent slow performance.`);
        console.log(`  Please make sure a daemon is running that will automatically restart it.`);
        console.log(`----------------------------------------------------------------------------`);
        console.log();

        setTimeout(() => {
            ConsoleLogger.log('Server', `Shut down due to timeout.`);
            process.exit(0);
        }, tMs);
    }
}

configureServer();

function startServer() {
    const app = express();
    app.use(compression({
        level: 6, memLevel: 9, chunkSize: 16 * 16384,
        filter: (req, res) => {
            const ct = res.getHeader('Content-Type');
            if (typeof ct === 'string' && ct.indexOf('tar+gzip') > 0) return false;
            return true;
        }
    }));

    initWebApi(app);

    const port = process.env.port || MembraneServerConfig.defaultPort;
    app.listen(port).setTimeout(MembraneServerConfig.requestTimeoutMs);

    console.log(`Mol* Membrane Server ${VERSION}`);
    console.log(``);
    console.log(`The server is running on port ${port}.`);
    console.log(``);
}

startServer();

if (MembraneServerConfig.shutdownTimeoutMinutes > 0) {
    setupShutdown();
}
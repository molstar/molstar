/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as express from 'express'
import * as compression from 'compression'
import * as fs from 'fs'
import * as argparse from 'argparse'
import { ModelServerConfig as ServerConfig, setupConfig, ModelServerConfig } from './config'
import { ConsoleLogger } from '../../mol-util/console-logger';
import { PerformanceMonitor } from '../../mol-util/performance-monitor';
import { initWebApi } from './server/api-web';
import Version from './version'

function setupShutdown() {
    if (ServerConfig.shutdownParams.timeoutVarianceMinutes > ServerConfig.shutdownParams.timeoutMinutes) {
        ConsoleLogger.log('Server', 'Shutdown timeout variance is greater than the timer itself, ignoring.');
    } else {
        let tVar = 0;
        if (ServerConfig.shutdownParams.timeoutVarianceMinutes > 0) {
            tVar = 2 * (Math.random() - 0.5) * ServerConfig.shutdownParams.timeoutVarianceMinutes;
        }
        let tMs = (ServerConfig.shutdownParams.timeoutMinutes + tVar) * 60 * 1000;

        console.log(`----------------------------------------------------------------------------`);
        console.log(`  The server will shut down in ${PerformanceMonitor.format(tMs)} to prevent slow performance.`);
        console.log(`  Please make sure a daemon is running that will automatically restart it.`);
        console.log(`----------------------------------------------------------------------------`);
        console.log();

        setTimeout(() => {
            /*if (WebApi.ApiState.pendingQueries > 0) {
                WebApi.ApiState.shutdownOnZeroPending = true;
            } else*/ {
                ConsoleLogger.log('Server', `Shut down due to timeout.`);
                process.exit(0);
            }
        }, tMs);
    }
}

const cmdParser = new argparse.ArgumentParser({
    addHelp: true
});

cmdParser.addArgument(['--cfg'], { help: 'Config file path.', required: false });
cmdParser.addArgument(['--printcfg'], { help: 'Prints out current config and exits.', required: false, nargs: 0 });

interface CmdArgs {
    cfg?: string,
    printcfg?: boolean
}

const cmdArgs = cmdParser.parseArgs() as CmdArgs;

const cfg = cmdArgs.cfg ? JSON.parse(fs.readFileSync(cmdArgs.cfg, 'utf8')) : void 0;
setupConfig(cfg);

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

if (cmdArgs.printcfg !== null) {
    console.log(JSON.stringify(ModelServerConfig, null, 2));
} else {
    startServer();

    if (ServerConfig.shutdownParams && ServerConfig.shutdownParams.timeoutMinutes > 0) {
        setupShutdown();
    }
}
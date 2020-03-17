/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as express from 'express'
import * as compression from 'compression'

import init from './server/web-api'
import VERSION from './server/version'
import { ConsoleLogger } from '../../mol-util/console-logger'
import { State } from './server/state'
import { addServerArgs, addLimitsArgs, LimitsConfig, setConfig, ServerConfig, addJsonConfigArgs, ServerJsonConfig, ServerConfigTemplate, validateServerConfig } from './config';
import * as argparse from 'argparse'
import * as fs from 'fs'

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

const parser = new argparse.ArgumentParser({
    addHelp: true,
    description: `VolumeServer ${VERSION}, (c) 2018-2019, Mol* contributors`
});
addJsonConfigArgs(parser);
addServerArgs(parser)
addLimitsArgs(parser)

const config: ServerConfig & LimitsConfig & ServerJsonConfig = parser.parseArgs()

if (config.printCfgTemplate !== null) {
    console.log(JSON.stringify(ServerConfigTemplate, null, 2));
    process.exit(0);
}

try {
    setConfig(config) // sets the config for global use

    if (config.cfg) {
        const cfg = JSON.parse(fs.readFileSync(config.cfg, 'utf8')) as ServerConfig & LimitsConfig;
        setConfig(cfg);
    }

    if (config.printCfg !== null) {
        console.log(JSON.stringify({ ...ServerConfig, ...LimitsConfig }, null, 2));
        process.exit(0);
    }

    validateServerConfig();
} catch (e) {
    console.error('' + e);
    process.exit(1);
}

const port = process.env.port || ServerConfig.defaultPort;

const app = express();
app.use(compression({ level: 6, memLevel: 9, chunkSize: 16 * 16384, filter: () => true }));
init(app);

app.listen(port);

console.log(`VolumeServer ${VERSION}, (c) 2018-2019, Mol* contributors`);
console.log(``);
console.log(`The server is running on port ${port}.`);
console.log(``);

if (config.shutdownTimeoutMinutes > 0) {
    setupShutdown();
}
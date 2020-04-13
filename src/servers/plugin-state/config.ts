/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as argparse from 'argparse';

export interface Config {
    working_folder: string,
    port?: string | number,
    api_prefix: string,
    max_states: number
}

export function getConfig() {
    const cmdParser = new argparse.ArgumentParser({
        addHelp: true
    });
    cmdParser.addArgument(['--working-folder'], { help: 'Working forlder path.', required: true });
    cmdParser.addArgument(['--port'], { help: 'Server port. Altenatively use ENV variable PORT.', type: 'int', required: false });
    cmdParser.addArgument(['--api-prefix'], { help: 'Server API prefix.', defaultValue: '', required: false });
    cmdParser.addArgument(['--max-states'], { help: 'Maxinum number of states that could be saved.', defaultValue: 40, type: 'int', required: false });

    const config = cmdParser.parseArgs() as Config;
    if (!config.port) config.port = process.env.port || 1339;
    return config;
}

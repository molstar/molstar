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
        add_help: true
    });
    cmdParser.add_argument('--working-folder', { help: 'Working forlder path.', required: true });
    cmdParser.add_argument('--port', { help: 'Server port. Altenatively use ENV variable PORT.', type: 'int', required: false });
    cmdParser.add_argument('--api-prefix', { help: 'Server API prefix.', default: '', required: false });
    cmdParser.add_argument('--max-states', { help: 'Maxinum number of states that could be saved.', default: 40, type: 'int', required: false });

    const config = cmdParser.parse_args() as Config;
    if (!config.port) config.port = process.env.port || 1339;
    return config;
}

/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'fs';

const VERSION = JSON.parse(fs.readFileSync('./package.json', 'utf8')).version;
const TIMESTAMP = Date.now();

const fileContents = {
    './lib/mol-plugin/version.js': [
        `// This file was replaced by write-version.mjs`,
        `export var PLUGIN_VERSION = '${VERSION}';`,
        `export var PLUGIN_VERSION_DATE = new Date(${TIMESTAMP});`,
    ],
    './lib/commonjs/mol-plugin/version.js': [
        `// This file was replaced by write-version.mjs`,
        `"use strict";`,
        `Object.defineProperty(exports, "__esModule", { value: true });`,
        `exports.PLUGIN_VERSION = '${VERSION}';`,
        `exports.PLUGIN_VERSION_DATE = new Date(${TIMESTAMP});`,
    ],
};

for (const filename in fileContents) {
    if (fs.existsSync(filename)) {
        fs.writeFileSync(filename, fileContents[filename].join('\n'));
    }
}

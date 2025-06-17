/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'fs';

const VERSION = JSON.parse(fs.readFileSync('./package.json', 'utf8')).version;
const TIMESTAMP = Date.now();
const file = `export var PLUGIN_VERSION = '${VERSION}';\nexport var PLUGIN_VERSION_DATE = new Date(${TIMESTAMP})`;
const files = ['./lib/mol-plugin/version.js', './lib/commonjs/mol-plugin/version.js'];
for (const f of files) {
    if (!fs.existsSync(f)) continue;
    fs.writeFileSync(f, file);
}
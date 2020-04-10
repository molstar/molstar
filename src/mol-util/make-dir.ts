/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'fs';

export function makeDir(path: string, root?: string): boolean {
    let dirs = path.split(/\/|\\/g),
        dir = dirs.shift();

    root = (root || '') + dir + '/';

    try {
        fs.mkdirSync(root);
    } catch (e) {
        if (!fs.statSync(root).isDirectory()) throw new Error(e);
    }

    return !dirs.length || makeDir(dirs.join('/'), root);
}
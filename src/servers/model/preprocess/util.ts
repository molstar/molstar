/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Progress } from '../../../mol-task';

export function showProgress(p: Progress) {
    process.stdout.write(`\r${new Array(80).join(' ')}`);
    process.stdout.write(`\r${Progress.format(p)}`);
}

export function clearLine() {
    process.stdout.write(`\r${new Array(80).join(' ')}`);
    process.stdout.write(`\r`);
}
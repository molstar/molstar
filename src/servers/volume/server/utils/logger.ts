/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export function formatTime(t: number) {
    if (isNaN(t)) return 'n/a';

    let h = Math.floor(t / (60 * 60 * 1000)),
        m = Math.floor(t / (60 * 1000) % 60),
        s = Math.floor(t / 1000 % 60),
        ms = Math.floor(t % 1000).toString();

    while (ms.length < 3) ms = '0' + ms;

    if (h > 0) return `${h}h${m}m${s}.${ms}s`;
    if (m > 0) return `${m}m${s}.${ms}s`;
    if (s > 0) return `${s}.${ms}s`;
    return `${t.toFixed(0)}ms`;
}

export function logPlain(tag: string, msg: string) {
    console.log(`[${tag}] ${msg}`);
}

export function log(guid: string, tag: string, msg: string) {
    console.log(`[${guid}][${tag}] ${msg}`);
}

export function errorPlain(ctx: string, e: any) {
    console.error(`[Error] (${ctx}) ${e}`);
    if (e.stack) console.error(e.stack);
}


export function error(guid: string, e: any) {
    console.error(`[${guid}][Error] ${e}`);
    if (e.stack) console.error(e.stack);
}

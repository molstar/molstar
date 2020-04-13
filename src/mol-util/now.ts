/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

declare const process: any;
declare const window: any;

const now: () => now.Timestamp = (function () {
    if (typeof window !== 'undefined' && window.performance) {
        const perf = window.performance;
        return () => perf.now();
    } else if (typeof process !== 'undefined' && process.hrtime !== 'undefined' && typeof process.hrtime === 'function') {
        return () => {
            const t = process.hrtime();
            return t[0] * 1000 + t[1] / 1000000;
        };
    } else if (Date.now) {
        return () => Date.now();
    } else {
        return () => +new Date();
    }
}());

namespace now {
    export type Timestamp = number & { '@type': 'now-timestamp' }
}


function formatTimespan(t: number, includeMsZeroes = true) {
    if (isNaN(t)) return 'n/a';

    let h = Math.floor(t / (60 * 60 * 1000)),
        m = Math.floor(t / (60 * 1000) % 60),
        s = Math.floor(t / 1000 % 60),
        ms = Math.floor(t % 1000).toString();

    while (ms.length < 3) ms = '0' + ms;
    while (!includeMsZeroes && ms.length > 1 && ms[ms.length - 1] === '0') ms = ms.substr(0, ms.length - 1);

    if (h > 0) return `${h}h${m}m${s}.${ms}s`;
    if (m > 0) return `${m}m${s}.${ms}s`;
    if (s > 0) return `${s}.${ms}s`;
    return `${t.toFixed(0)}ms`;
}

export { now, formatTimespan };
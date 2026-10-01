/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * on node `process.env.NODE_ENV` is available, in webpack build it is automatically set
 * by the DefinePlugin to the webpack `mode` value
 */

let isProductionMode = function () {
    try {
        return process.env.NODE_ENV === 'production';
    } catch {
        return false;
    }
}();

/**
 * set to true to enable more comprehensive checks and assertions,
 * mostly used in `mol-gl` and in valence-model calculation
 */
let isDebugMode = function getIsDebug() {
    try {
        const val = process.env.DEBUG;
        return val === '*' || val === 'molstar';
    } catch {
        return false;
    }
}();

/**
 * set to true to gather timings, mostly used in `mol-gl`
 */
let isTimingMode = false;

export { isProductionMode, isDebugMode, isTimingMode };

export function setProductionMode(value?: boolean) {
    if (typeof value !== 'undefined') isProductionMode = value;
}

export function setDebugMode(value?: boolean) {
    if (typeof value !== 'undefined') isDebugMode = value;
}

export function setTimingMode(value?: boolean) {
    if (typeof value !== 'undefined') isTimingMode = value;
}

//

type ConsoleStatsProvider = () => void
const consoleStatsProviders: ConsoleStatsProvider[] = [];

export function addConsoleStatsProvider(p: ConsoleStatsProvider) {
    if (!consoleStatsProviders.includes(p)) consoleStatsProviders.push(p);
}

export function removeConsoleStatsProvider(p: ConsoleStatsProvider) {
    const idx = consoleStatsProviders.indexOf(p);
    if (idx !== -1) consoleStatsProviders.splice(idx, 1);
}

export function consoleStats() {
    for (const p of consoleStatsProviders) {
        p();
    }
}

/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export function dateToUtcString(date: Date) {
    return date.toISOString().replace(/T/, ' ').replace(/\..+/, '');
}
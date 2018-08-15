/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export default function toUpperCase(value: any): string {
    if (!value) return '';
    return typeof value === 'string' ? value.toUpperCase() : `${value}`.toUpperCase();
}
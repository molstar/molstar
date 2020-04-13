/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export function urlQueryParameter (id: string) {
    if (typeof window === 'undefined') return undefined;
    const a = new RegExp(`${id}=([^&#=]*)`);
    const m = a.exec(window.location.search);
    return m ? decodeURIComponent(m[1]) : undefined;
}

export function urlCombine(base: string, query: string) {
    return `${base}${base[base.length - 1] === '/' || query[0] === '/' ? '' : '/'}${query}`;
}
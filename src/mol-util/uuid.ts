/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { now } from '../mol-util/now';

/** A UUID, either standard 36 characters or 22 characters base64 encoded. */
type UUID = string & { '@type': 'uuid' }

namespace UUID {
    const _btoa = typeof btoa !== 'undefined' ? btoa : (s: string) => Buffer.from(s).toString('base64');

    const chars: string[] = [];
    /** Creates a 22 characters 'base64' encoded UUID */
    export function create22(): UUID {
        let d = (+new Date()) + now();
        for (let i = 0; i < 16; i++) {
            chars[i] = String.fromCharCode((d + Math.random() * 0xff) % 0xff | 0);
            d = Math.floor(d / 0xff);
        }
        return _btoa(chars.join('')).replace(/\+/g, '-').replace(/\//g, '_').substr(0, 22) as UUID;
    }

    export function createv4(): UUID {
        let d = (+new Date()) + now();
        const uuid = 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
            const r = (d + Math.random() * 16) % 16 | 0;
            d = Math.floor(d / 16);
            return (c === 'x' ? r : (r & 0x3 | 0x8)).toString(16);
        });
        return uuid as any;
    }

    export function is(x: any): x is UUID {
        return typeof x === 'string';
    }
}

export default UUID;
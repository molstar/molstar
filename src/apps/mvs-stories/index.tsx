/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { getMVSStoriesContext } from './context';
import './index.html';
import './elements/snapshot-markdown';
import './elements/viewer';
import { MVSData } from '../../extensions/mvs/mvs-data';
import './elements/snapshot-markdown';
import './elements/viewer';

import './favicon.ico';
import '../../mol-plugin-ui/skin/light.scss';
import './styles.scss';
import './index.html';

export function getContext(name?: string) {
    return getMVSStoriesContext({ name });
}

export function loadFromURL(url: string, options?: { format: 'mvsx' | 'mvsj', contextName?: string }) {
    setTimeout(() => {
        getContext(options?.contextName).dispatch({
            kind: 'load-mvs',
            format: options?.format ?? 'mvsj',
            url,
        });
    }, 0);
}

export function loadFromData(data: MVSData | string | Uint8Array, options?: { format: 'mvsx' | 'mvsj', contextName?: string }) {
    setTimeout(() => {
        getContext(options?.contextName).dispatch({
            kind: 'load-mvs',
            format: options?.format ?? 'mvsj',
            data,
        });
    }, 0);
}

export { MVSData };
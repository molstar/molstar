/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import fetch, { AbortError, Response } from 'node-fetch';
import { retryIf } from '../../../mol-util/retry-if';
import { downloadGs } from '../../common/google-cloud-storage';


const RETRIABLE_NETWORK_ERRORS = [
    'ECONNRESET', 'ENOTFOUND', 'ESOCKETTIMEDOUT', 'ETIMEDOUT',
    'ECONNREFUSED', 'EHOSTUNREACH', 'EPIPE', 'EAI_AGAIN'
];

function isRetriableNetworkError(error: any) {
    return error && RETRIABLE_NETWORK_ERRORS.includes(error.code);
}

export async function fetchRetry(url: string, timeout: number, retryCount: number, onRetry?: () => void): Promise<Response> {
    const controller = new AbortController();
    const id = setTimeout(() => controller.abort(), timeout);
    const signal = controller.signal as any; // TODO: fix type
    const result = await retryIf(() => wrapFetch(url, { signal }), {
        retryThenIf: r => r.status === 408 /** timeout */ || r.status === 429 /** too many requests */ || (r.status >= 500 && r.status < 600),
        // TODO test retryCatchIf
        retryCatchIf: e => isRetriableNetworkError(e),
        onRetry,
        retryCount
    });
    clearTimeout(id);

    return result;
}

/** Like `fetch` but supports Google Cloud Storage (gs://) protocol. */
export function wrapFetch(url: string, init?: fetch.RequestInit): Promise<Response> {
    if (url.startsWith('gs://')) return fetchGS(url, init);
    else return fetch(url, init);
}

async function fetchGS(url: string, init?: fetch.RequestInit): Promise<Response> {
    if (init?.signal?.aborted) throw new AbortError('The user aborted a request.');
    const data = await downloadGs(url);
    return new Response(data, init);
}

/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fetch, { Response } from 'node-fetch';
import { retryIf } from '../../../mol-util/retry-if';

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
    const result = await retryIf(() => fetch(url, { signal: controller.signal }), {
        retryThenIf: r => r.status === 408 /** timeout */ || r.status === 429 /** too many requests */ || (r.status >= 500 && r.status < 600),
        // TODO test retryCatchIf
        retryCatchIf: e => isRetriableNetworkError(e),
        onRetry,
        retryCount
    });
    clearTimeout(id);

    return result;
}

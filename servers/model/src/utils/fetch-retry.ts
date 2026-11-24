/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Paul Pillot <paul.pillot@tandemai.com>
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import { retryIf } from '@molstar/core/util/retry-if';
import { downloadGs } from '@molstar/common-server/google-cloud-storage';

const RETRIABLE_NETWORK_ERRORS = [
  'ECONNRESET',
  'ENOTFOUND',
  'ESOCKETTIMEDOUT',
  'ETIMEDOUT',
  'ECONNREFUSED',
  'EHOSTUNREACH',
  'EPIPE',
  'EAI_AGAIN',
  // Node's fetch (undici): the server closed a reused keep-alive connection, or connecting timed out
  'UND_ERR_SOCKET',
  'UND_ERR_CONNECT_TIMEOUT',
];

export function isRetriableNetworkError(error: any) {
  // Node's fetch rejects with `TypeError: fetch failed`, the system error is its `cause`
  const code = error?.code ?? error?.cause?.code;
  return !!code && RETRIABLE_NETWORK_ERRORS.includes(code);
}

/** The error with its cause, which says what went wrong with a `fetch failed`. */
export function describeFetchError(error: any) {
  const cause = error?.cause;
  if (!cause) return '' + error;
  return `${error} (${cause.code ?? cause.message ?? cause})`;
}

export function isRetriableStatus(status: number) {
  return status === 408 /** timeout */ || status === 429 /** too many requests */ || (status >= 500 && status < 600);
}

/**
 * `fetch` that retries on network errors and retriable statuses. Each attempt gets its own `timeout` (ms), and the
 * signal stays attached to the response, so the timeout bounds reading the body as well: pick it for the whole
 * download, not the first byte. A timed-out attempt is not retried.
 */
export async function fetchRetry(
  url: string,
  timeout: number,
  retryCount: number,
  onRetry?: () => void,
): Promise<Response> {
  return retryIf(() => wrapFetch(url, { signal: AbortSignal.timeout(timeout) }), {
    retryThenIf: (r) => {
      const retry = isRetriableStatus(r.status);
      // an unread body holds on to its connection until it is collected
      if (retry) void r.body?.cancel().catch(() => {});
      return retry;
    },
    retryCatchIf: (e) => isRetriableNetworkError(e),
    onRetry,
    retryCount,
  });
}

/** Like `fetch` but supports Google Cloud Storage (gs://) protocol. */
export function wrapFetch(url: string, init?: RequestInit): Promise<Response> {
  if (url.startsWith('gs://')) return fetchGS(url, init);
  else return fetch(url, init);
}

async function fetchGS(url: string, init?: RequestInit): Promise<Response> {
  if (init?.signal?.aborted) throw new Error('The user aborted a request.');
  const data = await downloadGs(url);
  return new Response(data, init);
}

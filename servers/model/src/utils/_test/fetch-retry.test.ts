/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { fetchRetry, isRetriableNetworkError } from '../fetch-retry.js';

function response(status: number) {
  const cancel = vi.fn(() => Promise.resolve());
  const r = new Response(status === 204 ? null : 'body', { status });
  Object.defineProperty(r, 'body', { value: { cancel } });
  return { response: r, cancel };
}

describe('model server fetchRetry', () => {
  it('retries retriable statuses with a fresh signal and drops the unread bodies', async () => {
    const failed = response(503);
    const ok = response(200);
    const signals: AbortSignal[] = [];
    const fetch = vi.fn(async (_url: string, init?: RequestInit) => {
      signals.push(init!.signal!);
      return fetch.mock.calls.length === 1 ? failed.response : ok.response;
    });
    vi.stubGlobal('fetch', fetch);

    expect(await fetchRetry('https://example.org/1cbs.bcif', 1000, 3)).toBe(ok.response);
    expect(fetch).toHaveBeenCalledTimes(2);
    expect(failed.cancel).toHaveBeenCalled();
    expect(ok.cancel).not.toHaveBeenCalled();
    // One signal per attempt: an attempt that timed out must not abort the next one before it starts.
    expect(signals[0]).not.toBe(signals[1]);
  });

  it('does not retry a not found', async () => {
    const notFound = response(404);
    const fetch = vi.fn(async () => notFound.response);
    vi.stubGlobal('fetch', fetch);

    expect((await fetchRetry('https://example.org/xxxx.bcif', 1000, 3)).status).toBe(404);
    expect(fetch).toHaveBeenCalledTimes(1);
  });

  it('retries network errors, which Node reports as the cause of a TypeError', async () => {
    const fetch = vi.fn(async () => {
      if (fetch.mock.calls.length === 1) throw new TypeError('fetch failed', { cause: { code: 'ECONNRESET' } });
      return response(200).response;
    });
    vi.stubGlobal('fetch', fetch);

    expect((await fetchRetry('https://example.org/1cbs.bcif', 1000, 3)).status).toBe(200);
    expect(fetch).toHaveBeenCalledTimes(2);
    expect(isRetriableNetworkError({ code: 'ETIMEDOUT' })).toBe(true);
    expect(isRetriableNetworkError(new TypeError('fetch failed', { cause: { code: 'UND_ERR_SOCKET' } }))).toBe(true);
    expect(isRetriableNetworkError(new DOMException('timed out', 'TimeoutError'))).toBe(false);
  });
});

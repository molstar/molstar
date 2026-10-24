/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { fetchData } from '../data-source.js';

// `fetchData` as in a browser: `RUNNING_IN_NODEJS` is fixed at import, so this file mocks the module and runs in a
// registry of its own.
vi.mock('../nodejs-shims.js', async (importOriginal) => ({
  ...(await importOriginal<typeof import('../nodejs-shims.js')>()),
  RUNNING_IN_NODEJS: false,
}));

describe('fetchData in a browser', () => {
  it('fails on a network error with the generic message, the cause attached', async () => {
    const failed = new TypeError('Failed to fetch');
    vi.stubGlobal(
      'fetch',
      vi.fn(() => Promise.reject(failed)),
    );
    const run = fetchData('https://example.org/x').run();
    await expect(run).rejects.toThrow(
      'Failed to download data. Possible reasons: Resource is not available, or CORS is not allowed on the server.',
    );
    await expect(run).rejects.toHaveProperty('cause', failed);
  });
});

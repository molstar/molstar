/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';

describe('pipeline stage runner', () => {
  it('records a throwing stage as failed, without numbers, and exits non-zero (perf-audit-3 12/B4)', async () => {
    // stage.ts needs --expose-gc; the collections only settle the memory columns, not checked here
    vi.stubGlobal('gc', () => {});
    vi.spyOn(console, 'error').mockImplementation(() => {});
    const { createStageRunner, exitOnFailures, summarize } = await import('../stage.js');
    const runner = createStageRunner({ only: [], skip: [], alloc: false });
    const results: Record<string, any> = {};
    await runner.stage('ok', results, () => 'fine');
    await runner.stage('throws', results, () => {
      throw new Error('boom');
    });

    expect(results.ok.error).toBeUndefined();
    expect(results.throws.error).toBe('boom');
    expect(runner.failureCount()).toBe(1);

    const rows = summarize({ ok: [results.ok], throws: [results.throws] });
    expect(rows[0]).toMatchObject({ stage: 'ok', info: 'fine' });
    expect(typeof rows[0].medMs).toBe('number');
    expect(rows[1]).toMatchObject({ stage: 'throws', error: 'boom', coldMs: null, medMs: null, allocMB: null });

    const exitCode = process.exitCode;
    try {
      exitOnFailures(0);
      expect(process.exitCode).toBe(exitCode);
      exitOnFailures(runner.failureCount());
      expect(process.exitCode).toBe(1);
    } finally {
      process.exitCode = exitCode;
    }
  });
});

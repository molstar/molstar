/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { PluginContext } from '@molstar/plugin/context';

// Every evaluated transform logs its timing: a long session (an agent's) grew the log by tens of MB (perf-audit-2 12/B3)
const LogSize = 1000;

describe('the plugin log', () => {
  it('keeps the last 1000 entries (LogSize in context.ts)', async () => {
    const plugin = new PluginContext({ behaviors: [] });
    await plugin.init();
    try {
      const n = LogSize + 50;
      for (let i = 0; i < n; i++) plugin.log.info(`entry ${i}`);
      const { entries } = plugin.log;
      expect(entries.size).toBe(LogSize);
      expect(entries.last()!.message).toBe(`entry ${n - 1}`);
      expect(entries.first()!.message).toBe(`entry ${n - LogSize}`);
    } finally {
      plugin.dispose();
    }
  });
});

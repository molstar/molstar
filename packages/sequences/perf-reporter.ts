/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick tier's report (`test:perf`; `perfTests` in vitest.shared.ts): a table of every case (bytes a step, the
// baseline, the delta, the measurement's note), the moves under the failing threshold listed, and under `--update` the
// package's baseline rewritten from the run (cases that did not run keep theirs; positive controls are not baselined).
// A Vitest reporter by shape: this folder has no vitest of its own to import its types from.

import { readFileSync, writeFileSync } from 'node:fs';

/** What `quickCase` (core/_test/model/quick.ts) leaves on a test's meta. */
interface AllocRecord {
  id: string;
  bytesPerStep: number;
  baseline?: number;
  control?: number;
  note: string;
}

/** A move this large (bytes a step, or 2 % of a larger baseline) is listed though it does not fail. */
const WARN = 1;

export class PerfReporter {
  private readonly records: AllocRecord[] = [];
  private readonly baselinePath: string;
  private readonly update: boolean;
  private readonly threshold: number;
  constructor(baselinePath: string, update: boolean, threshold: number) {
    this.baselinePath = baselinePath;
    this.update = update;
    this.threshold = threshold;
  }

  onTestCaseResult(test: { meta(): object }): void {
    const alloc = (test.meta() as { alloc?: AllocRecord[] }).alloc;
    if (alloc) this.records.push(...alloc);
  }

  onTestRunEnd(): void {
    if (this.records.length === 0) return;
    const rows = [...this.records].sort((a, b) => a.id.localeCompare(b.id));
    const w = Math.max(...rows.map((r) => r.id.length));
    const num = (x: number | undefined) => (x === undefined ? '-' : x.toFixed(2)).padStart(9);
    const lines = [
      `${'case'.padEnd(w)}  ${'B/step'.padStart(9)} ${'baseline'.padStart(9)} ${'delta'.padStart(9)}  note`,
    ];
    const moved: string[] = [];
    for (const r of rows) {
      const delta = r.baseline === undefined ? undefined : r.bytesPerStep - r.baseline;
      const flag =
        r.control !== undefined
          ? `control +${r.control.toFixed(1)} B`
          : delta === undefined
            ? 'new'
            : delta > Math.max(this.threshold, 0.1 * (r.baseline ?? 0))
              ? 'FAIL'
              : Math.abs(delta) >= Math.max(WARN, 0.02 * (r.baseline ?? 0))
                ? 'moved'
                : '';
      if (flag === 'moved' || flag === 'FAIL') moved.push(`${r.id}: ${num(delta).trim()} B/step`);
      lines.push(`${r.id.padEnd(w)}  ${num(r.bytesPerStep)} ${num(r.baseline)} ${num(delta)}  ${flag} ${r.note}`);
    }
    console.log(
      `\nperf:quick (bytes allocated a step; fails over the baseline + ${this.threshold}, or + 10 % of a larger one)\n${lines.join('\n')}`,
    );
    if (moved.length)
      console.log(`moved: ${moved.join('; ')}${this.update ? '' : ' (record an intended move with --update)'}`);
    if (!this.update) return;
    let baseline: Record<string, number> = {};
    try {
      baseline = JSON.parse(readFileSync(this.baselinePath, 'utf8'));
    } catch {}
    for (const r of rows)
      if (r.control === undefined) baseline[r.id] = Math.max(0, Math.round(r.bytesPerStep * 10) / 10);
    const sorted = Object.fromEntries(Object.entries(baseline).sort(([a], [b]) => a.localeCompare(b)));
    writeFileSync(this.baselinePath, `${JSON.stringify(sorted, null, 2)}\n`);
    console.log(`perf:quick: baseline written (${this.baselinePath})`);
  }
}

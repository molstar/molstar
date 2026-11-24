/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The search kernels' allocation (_spec/README.md "Search", Performance): over 100 kb of seeded letters, `degenerate`
// GANTC, `exact` GAATTC and `regex` GG.{2,8}CC run warm with a visitor that counts, 0 bytes a hit (nothing allocated a
// hit: hundreds of hits a step). And the step budget's calibration: the worst patterns (`.{1,126}.*`, `(.|A)*`,
// `[ACGT]*A`, Horspool's `C` and 999 `A`) over 2^20 `A`s end at the default `maxSteps` within about a second, the
// design's measure of 2^28 steps. By hand, never in CI: `pnpm -C packages/seq/core test:perf`.

import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import {
  classKernel,
  classMasks,
  DefaultSearchBounds,
  ExactMatcher,
  exactKernel,
  glushkov,
  parseRegex,
  type RegexParsed,
  RegexMatcher,
  regexKernel,
  regexTables,
} from '../index.js';
import { runSearch, textTarget } from './run.js';

const N = 100_000;
const letters = (() => {
  let s = 42;
  const out = new Uint8Array(N);
  for (let i = 0; i < N; i++) {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    out[i] = 'ACGT'.charCodeAt((s >>> 24) & 3);
  }
  return out;
})();
const codes = (s: string) => Uint8Array.from(s, (c) => c.charCodeAt(0));
const Sets: Record<string, number> = { A: 1, C: 2, G: 4, T: 8, N: 15 };

describe('the search kernels over 100 kb', () => {
  const budget = { steps: 0, limit: 2 ** 40 };
  let hits = 0;
  const visit = () => {
    hits++;
  };
  const cases = {
    'exact GAATTC': exactKernel(codes('GAATTC')),
    'degenerate GANTC': classKernel(
      classMasks(5, (i, c) => {
        const s = Sets[String.fromCharCode(c)] ?? 0;
        return s !== 0 && (s & ~Sets['GANTC'[i]]) === 0;
      }),
      5,
    ),
    'regex GG.{2,8}CC': (() => {
      const p = parseRegex('GG.{2,8}CC', DefaultSearchBounds) as RegexParsed;
      return regexKernel(regexTables(glushkov(p)), p, 12, (set, c) => set[c] === 1);
    })(),
  };
  for (const [name, kernel] of Object.entries(cases)) {
    it(`${name}: 0 bytes a step warm, hits found`, ({ task }) => {
      hits = 0;
      kernel.find(letters, N, 0, N, 3, budget, visit);
      expect(hits).toBeGreaterThan(10);
      const { r } = quickCase(
        task,
        `search ${name} over 100 kb`,
        () => kernel.find(letters, N, 0, N, 3, budget, visit),
        {
          steps: 20,
        },
      );
      console.info(`search ${name} over 100 kb: ${r.bytesPerStep.toFixed(0)} bytes a step (${r.note})`);
      expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    });
  }
});

describe('the step budget over 2^20 As: the worst patterns end at it within about a second', () => {
  const polyA = textTarget('A'.repeat(2 ** 20));
  const cases: [typeof RegexMatcher, string][] = [
    [RegexMatcher, '.{1,126}.*'],
    [RegexMatcher, '(.|A)*'],
    [RegexMatcher, '[ACGT]*A'],
    [ExactMatcher as unknown as typeof RegexMatcher, `C${'A'.repeat(999)}`],
  ];
  for (const [matcher, pattern] of cases) {
    it(`${pattern.slice(0, 12)}: ends at the budget within 1.5 s`, async () => {
      const t0 = performance.now();
      const r = await runSearch(matcher, pattern, [polyA]);
      const ms = performance.now() - t0;
      console.info(`search ${pattern.slice(0, 12)} over 2^20 As: ${ms.toFixed(0)} ms to the step budget`);
      expect(r.issues.map((i) => i.code)).toContain('seq.search.steps.limit');
      expect(ms).toBeLessThan(1500);
    });
  }
});

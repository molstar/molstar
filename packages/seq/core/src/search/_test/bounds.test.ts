/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The regular expression's safety and every bound with its code (_spec/README.md "Search", Regular expressions, Bounds):
// the refusals, a catastrophic pattern costing what a plain one does in counted steps, the bounds each refused before the
// work it sizes, counted past the cap, cancellation mid-run.

import { describe, expect, it } from 'vitest';
import {
  DegenerateMatcher,
  DefaultSearchBounds,
  ExactMatcher,
  exactKernel,
  glushkov,
  parseRegex,
  type RegexParsed,
  RegexMatcher,
  RegexStepOverhead,
  regexKernel,
  regexTables,
  type SearchBatch,
  searchBounds,
} from '../index.js';
import { runSearch, textTarget } from './run.js';

const compileCodes = (m: typeof RegexMatcher | typeof ExactMatcher, pattern: string, over = {}) => {
  const r = m.compile(pattern, {}, searchBounds(over));
  return Array.isArray(r) ? r.map((i) => i.code) : [];
};

describe('a regular expression’s refusals', () => {
  it.each([
    ['(A)\\1', 'a backreference'],
    ['A(?=C)', 'a lookahead'],
    ['A(?!C)', 'a negative lookahead'],
    ['(?<=A)C', 'a lookbehind'],
    ['A*?', 'a lazy repeat'],
    ['A++', 'a possessive repeat'],
    ['\\bAC', '\\b'],
    ['\\d', '\\d'],
    ['\\w+', '\\w'],
    ['(?i)ac', 'a flag'],
    ['A^C', 'an anchor inside'],
  ])('%s: unsupported (%s), its position said', (pattern) => {
    const r = parseRegex(pattern, DefaultSearchBounds);
    expect(Array.isArray(r)).toBe(true);
    const [issue] = r as readonly { code: string; offset?: number; hint?: string }[];
    expect(issue.code).toBe('seq.search.pattern.unsupported');
    expect(issue.offset).toBeTypeOf('number');
  });
  it('a pattern that does not read is a syntax error', () => {
    for (const p of ['(AC', '[AC', 'A)', '*A', 'A{3,2}'])
      expect(compileCodes(RegexMatcher, p)).toEqual(['seq.search.pattern.syntax']);
  });
  it('a pattern that matches only the empty string', () => {
    for (const p of ['()', 'A{0}', '(?:)']) expect(compileCodes(RegexMatcher, p)).toEqual(['seq.search.pattern.empty']);
  });
  it('the language: classes, groups, alternation, counted repeats, escapes', () => {
    for (const p of ['[ACG]T', '[^T]+', '(?:GA|TC){2,3}', 'A{2,}', 'A\\.C', 'N.{0,4}[ST]', '^M'])
      expect(compileCodes(RegexMatcher, p)).toEqual([]);
  });
});

describe('the bounds, each with its code', () => {
  it('a pattern past maxPattern letters (exact, degenerate)', () => {
    expect(compileCodes(ExactMatcher, 'A'.repeat(1001))).toEqual(['seq.search.pattern.size']);
    expect(compileCodes(DegenerateMatcher as typeof ExactMatcher, 'N'.repeat(5), { maxPattern: 4 })).toEqual([
      'seq.search.pattern.size',
    ]);
  });
  it('a regular expression past its characters, positions after expansion, nesting (EC-76)', () => {
    expect(compileCodes(RegexMatcher, 'A'.repeat(257))).toEqual(['seq.search.pattern.size']);
    // 129 positions once the count is copied out, refused before any table is allocated
    expect(compileCodes(RegexMatcher, '(ACG){43}')).toEqual(['seq.search.pattern.size']);
    expect(compileCodes(RegexMatcher, 'A{0,999999}')).toEqual(['seq.search.pattern.size']);
    expect(compileCodes(RegexMatcher, `${'('.repeat(65)}A${')'.repeat(65)}`)).toEqual(['seq.search.pattern.size']);
    expect(compileCodes(RegexMatcher, `${'('.repeat(64)}A${')'.repeat(64)}`)).toEqual([]);
  });
  it('a compile never hangs: a counted repeat of an empty body is nothing to copy, the automaton bounded in nodes', () => {
    const t0 = performance.now();
    // 21 characters that copied out are about 10^12 empty groups
    expect(compileCodes(RegexMatcher, 'A((){999999}){999999}')).toEqual([]);
    expect(compileCodes(RegexMatcher, '((){999}){999}')).toEqual(['seq.search.pattern.empty']);
    // copies of an optional group: within the positions, past the automaton's nodes
    expect(compileCodes(RegexMatcher, `${'('.repeat(30)}A${')?'.repeat(30)}{100}`)).toEqual([
      'seq.search.pattern.size',
    ]);
    expect(performance.now() - t0).toBeLessThan(2000);
  });
  it('a regex hit’s length: the longest match from a start stops at maxHitLength', async () => {
    const r = await runSearch(
      RegexMatcher,
      'A+',
      [textTarget('A'.repeat(50), { alphabet: 'protein' })],
      {},
      { maxHitLength: 10 },
    );
    expect(r.hits[0]).toEqual([0, 0, 10, 0]);
    expect(r.hits.every((h) => h[2] - h[1] <= 10)).toBe(true);
  });
  it('letters a run reads: the target that would pass maxLetters is skipped, the others searched', async () => {
    const r = await runSearch(
      ExactMatcher,
      'AC',
      [
        textTarget('ACAC', { system: '/a' }),
        textTarget('ACACACAC', { system: '/b' }),
        textTarget('AC', { system: '/c' }),
      ],
      { strands: 'forward' },
      { maxLetters: 7 },
    );
    expect(r.hits.map((h) => h[0])).toEqual([0, 0, 2]);
    expect(r.issues.map((i) => [i.code, i.address, i.hint])).toEqual([
      ['seq.search.letters.limit', '/b', 'Search within a region (`within`)'],
    ]);
  });
  it('word steps a run takes: the run ends at maxSteps, incomplete, said', async () => {
    const r = await runSearch(
      ExactMatcher,
      'AC',
      [textTarget('AC'.repeat(100))],
      { strands: 'forward' },
      { maxSteps: 150, chunkLetters: 64 },
    );
    expect(r.issues.map((i) => i.code)).toEqual(['seq.search.steps.limit']);
    expect(r.hits.length).toBeGreaterThan(0);
    expect(r.hits.length).toBeLessThan(100);
  });
  it('hits held: counting goes on past maxHits, the first held, said once', async () => {
    const r = await runSearch(
      ExactMatcher,
      'A',
      [textTarget('A'.repeat(100))],
      { strands: 'forward' },
      { maxHits: 30, chunkLetters: 16 },
    );
    expect(r.hits.length).toBe(30);
    expect(r.hits.length + r.counted).toBe(100);
    expect(r.hits.map((h) => h[1])).toEqual(Array.from({ length: 30 }, (_, i) => i));
    expect(r.issues.map((i) => i.code)).toEqual(['seq.search.hits.limit']);
  });
  it('a failed read skips its target, said', async () => {
    const failing = {
      ...textTarget('ACGT', { system: '/bad' }),
      letters: { read: () => Promise.reject(new Error('gone')) },
    };
    const r = await runSearch(ExactMatcher, 'AC', [failing, textTarget('AC')], { strands: 'forward' });
    expect(r.hits).toEqual([[1, 0, 2, 1]]);
    expect(r.issues.map((i) => [i.code, i.address])).toEqual([['seq.search.letters.failed', '/bad']]);
  });
});

describe('never backtracking', () => {
  it('(A|AA)*B over 100,000 As costs what AB costs, in counted steps', async () => {
    const letters = 'A'.repeat(100_000);
    const steps = async (pattern: string) => {
      // the step budget is the measure: the largest budget a run fits in, found by halving
      let lo = 0;
      let hi = 2 ** 22;
      while (hi - lo > 1) {
        const mid = (lo + hi) >> 1;
        const r = await runSearch(
          RegexMatcher,
          pattern,
          [textTarget(letters, { alphabet: 'protein' })],
          {},
          { maxSteps: mid },
        );
        if (r.issues.some((i) => i.code === 'seq.search.steps.limit')) lo = mid;
        else hi = mid;
      }
      return hi;
    };
    const plain = await steps('AB');
    const catastrophic = await steps('(A|AA)*B');
    // a letter a step of 2 positions: one word, its mask and test, the step's own (`RegexStepOverhead`)
    expect(plain).toBe(100_000 * (1 + 2 + RegexStepOverhead));
    expect(catastrophic).toBe(plain);
    // (A+)+B too
    expect(await steps('(A+)+B')).toBe(plain);
  });
  it('a step is charged the work it does: a regex step chunks × w, a Horspool check its comparisons', () => {
    const n = 2000;
    const polyA = new Uint8Array(n).fill(65);
    const maxLength = 300;
    // every start of poly-A runs on to its longest: min(maxLength, n − s) steps
    let visits = 0;
    for (let s = 0; s < n; s++) visits += Math.min(maxLength, n - s);
    for (const pattern of ['.{1,126}.*', '(.|A)*', '[ACGT]*A']) {
      const p = parseRegex(pattern, DefaultSearchBounds) as RegexParsed;
      const g = glushkov(p);
      const kernel = regexKernel(regexTables(g), p, maxLength, (set, c) => set[c] === 1);
      const budget = { steps: 0, limit: 2 ** 40 };
      expect(kernel.find(polyA, n, 0, n, 3, budget, () => {})).toBe(true);
      // a step reads a table a byte of the state: `.{1,126}.*` 16 bytes of 4 words, 64 words a step
      const stepWords = Math.ceil(g.positions / 8) * g.words + RegexStepOverhead;
      expect(budget.steps, pattern).toBeGreaterThanOrEqual(visits * stepWords);
    }
    // C and 999 A over poly-A: every alignment compares 999 letters (O(n·m)), charged as compared
    const kernel = exactKernel(Uint8Array.from(`C${'A'.repeat(999)}`, (c) => c.charCodeAt(0)));
    const budget = { steps: 0, limit: 2 ** 40 };
    kernel.find(polyA, n, 0, n, 3, budget, () => {});
    expect(budget.steps).toBeGreaterThanOrEqual((n - 999) * 999);
    // and a check past the budget ends the find
    expect(kernel.find(polyA, n, 0, n, 3, { steps: 0, limit: 10_000 }, () => {})).toBe(false);
  });
  it('a pattern long at every start (A* over poly-A) ends at the budget, said, never a hang', async () => {
    const r = await runSearch(
      RegexMatcher,
      'A*',
      [textTarget('A'.repeat(20_000), { alphabet: 'protein' })],
      {},
      {
        maxSteps: 1_000_000,
      },
    );
    expect(r.issues.map((i) => i.code)).toEqual(['seq.search.steps.limit']);
  });
});

describe('cancellation', () => {
  it('a run cancelled mid-way ends with no error and reads no more', async () => {
    const reads: number[][] = [];
    const ac = new AbortController();
    const target = textTarget('ACGT'.repeat(1000), { reads });
    const batches: SearchBatch[] = [];
    const compiled = ExactMatcher.compile('ACG', {}, DefaultSearchBounds);
    if (Array.isArray(compiled)) throw new Error('compile');
    await new Promise<void>((resolve, reject) => {
      ExactMatcher.run(compiled as Exclude<typeof compiled, readonly unknown[]>, [target], {
        signal: ac.signal,
        bounds: searchBounds({ chunkLetters: 100 }),
        yield: async () => {
          if (batches.length === 3) ac.abort();
        },
      }).subscribe({ next: (b) => batches.push(b), error: reject, complete: resolve });
    });
    expect(batches.length).toBe(3);
    expect(reads.length).toBe(3);
  });
  it('a run cancelled while a read is pending ends with no error, emits nothing, and its reader sees the signal', async () => {
    for (const reader of ['resolves late', 'rejects on abort'] as const) {
      const ac = new AbortController();
      let seen: AbortSignal | undefined;
      let late: (() => void) | undefined;
      const codes = Uint8Array.from('ACGT'.repeat(10), (c) => c.charCodeAt(0));
      const target = {
        ...textTarget('ACGT'.repeat(10)),
        letters: {
          read: (start: number, end: number, signal?: AbortSignal) =>
            new Promise<Uint8Array>((resolve, reject) => {
              seen = signal;
              late = () => resolve(codes.slice(start, end));
              if (reader === 'rejects on abort')
                signal?.addEventListener('abort', () => reject(new DOMException('Aborted', 'AbortError')));
            }),
        },
      };
      const batches: SearchBatch[] = [];
      const compiled = ExactMatcher.compile('ACG', {}, DefaultSearchBounds);
      if (Array.isArray(compiled)) throw new Error('compile');
      const done = new Promise<void>((resolve, reject) => {
        ExactMatcher.run(compiled as Exclude<typeof compiled, readonly unknown[]>, [target], {
          signal: ac.signal,
          bounds: DefaultSearchBounds,
        }).subscribe({ next: (b) => batches.push(b), error: reject, complete: resolve });
      });
      await Promise.resolve();
      expect(seen, reader).toBeDefined();
      ac.abort();
      expect(seen!.aborted, reader).toBe(true);
      late?.();
      await done;
      // no hits from the late letters, no `seq.search.letters.failed` for the cancel
      expect(batches, reader).toEqual([]);
    }
  });
  it('unsubscribing cancels', async () => {
    const reads: number[][] = [];
    const target = textTarget('ACGT'.repeat(1000), { reads });
    const compiled = ExactMatcher.compile('ACG', {}, DefaultSearchBounds) as Exclude<
      ReturnType<typeof ExactMatcher.compile>,
      readonly unknown[]
    >;
    const sub = ExactMatcher.run(compiled, [target], {
      signal: new AbortController().signal,
      bounds: searchBounds({ chunkLetters: 100 }),
    }).subscribe(() => sub.unsubscribe());
    for (let i = 0; i < 50; i++) await Promise.resolve();
    expect(reads.length).toBeLessThanOrEqual(2);
  });
});

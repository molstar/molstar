/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @module-tag heavy

// The spec targets (performance.md, Time): validate and resolve a spec of about 1,000 nodes (inline data aside) under
// 20 ms (the model, then the views' resolution: the part a viewer waits for). The spec is generated: 100 short sequences, each with ranges and values over it, and a section of three tracks
// for each. Measured warm (the median of 200 runs after 100); the first, cold run is printed too.

import { parseJson, printJson } from '@molstar/core/repr/codecs/json';
import type { Node } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { checkSeqSpec, createSeqSpec, SeqLanguage } from '../../lang/index.js';
import { resolveViews } from '../../view/index.js';
import { resolveSpec } from '../index.js';

const count = (n: Node): number => 1 + (n.children?.reduce((a, c) => a + count(c), 0) ?? 0);

/** 100 sequences of 60 residues, each with a ranges and a values annotation, one section of three tracks each. */
function bigSpec(sequences = 100): Node {
  const b = createSeqSpec({ title: 'generated' });
  const ids: {
    s: ReturnType<typeof b.data.sequence>;
    r: ReturnType<typeof b.data.ranges>;
    v: ReturnType<typeof b.data.values>;
  }[] = [];
  for (let i = 0; i < sequences; i++) {
    const s = b.data.sequence(`s${i}`, { alphabet: 'protein', letters: 'ACDEFGHIKLMNPQRSTVWY'.repeat(3) });
    const r = b.data.ranges(`r${i}`, {
      target: s,
      rows: { start: new Int32Array([0, 20, 40]), end: new Int32Array([15, 35, 60]), id: ['a', 'b', 'c'] },
    });
    const v = b.data.values(`v${i}`, { target: s, values: new Float32Array(60).fill(0.5) });
    ids.push({ s, r, v });
  }
  b.views.view('main', (view) => {
    ids.forEach(({ s, r, v }, i) => {
      const sec = view.section(`sec${i}`, { of: s });
      sec.track('letters').layer('letters', { type: 'letters' });
      sec.track('ranges').layer('ranges', { type: 'ranges', data: r });
      sec.track('values').layer('values', { type: 'values', data: v });
    });
  });
  return b.build();
}

const median = (t: number[]) => [...t].sort((a, b) => a - b)[t.length >> 1];

function timeRuns(f: () => unknown, runs = 200, warm = 100): number[] {
  for (let i = 0; i < warm; i++) f();
  const out: number[] = [];
  for (let i = 0; i < runs; i++) {
    const t0 = performance.now();
    f();
    out.push(performance.now() - t0);
  }
  return out;
}

describe('a spec of about 1,000 nodes', () => {
  const spec = bigSpec();

  it('has about 1,000 nodes', { tags: ['heavy'] }, () => {
    const n = count(spec);
    console.info(`generated spec: ${n} nodes, ${printJson(spec, { language: SeqLanguage }).length} bytes of JSON`);
    expect(n).toBeGreaterThan(900);
    expect(n).toBeLessThan(1100);
  });

  it('validate, resolve and resolve the views: under 20 ms (median, warm)', { tags: ['heavy'] }, () => {
    const t0 = performance.now();
    checkSeqSpec(spec);
    resolveSpec(spec);
    console.info(`validate + resolve, first call in this process: ${(performance.now() - t0).toFixed(1)} ms`);
    const validate = median(timeRuns(() => checkSeqSpec(spec)));
    const resolve = median(timeRuns(() => resolveSpec(spec)));
    const views = median(timeRuns(() => resolveViews(resolveSpec(spec))));
    const both = median(
      timeRuns(() => {
        const issues = checkSeqSpec(spec);
        resolveViews(resolveSpec(spec));
        return issues;
      }),
    );
    console.info(
      `1,000 nodes, medians: validate ${validate.toFixed(2)} ms, resolve ${resolve.toFixed(2)} ms, resolve + resolve views ${views.toFixed(2)} ms, validate + resolve + views ${both.toFixed(2)} ms`,
    );
    expect(both).toBeLessThan(20);
  });

  it('decode JSON, validate, resolve: reported', { tags: ['heavy'] }, () => {
    const json = printJson(spec);
    const decode = median(timeRuns(() => parseJson(json)));
    const all = median(
      timeRuns(() => {
        const d = parseJson(json);
        checkSeqSpec(d.node as Node);
        return resolveSpec(d.node as Node);
      }),
    );
    console.info(
      `from ${json.length} bytes of JSON, medians: decode ${decode.toFixed(2)} ms, decode + validate + resolve ${all.toFixed(2)} ms`,
    );
    expect(all).toBeLessThan(200);
  });
});

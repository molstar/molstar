/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { validateAndNormalize } from '@molstar/core/repr/schema';
import type { Node } from '@molstar/core/repr/tree';
import { SeqLanguage, SeqLanguageVersion, SeqValidateOptions } from '@molstar/seqx-core/lang';
import { expect } from 'vitest';
import {
  type AssetReader,
  DefaultOpenLimits,
  type FormatRuntime,
  isAbortError,
  type OpenedSource,
  type RegionRequest,
  type SourceEntry,
} from '../protocol.js';
import { checkCoverage } from '../seq-format.js';
import { checkAnswer, sliceTable, tableRows } from '../tables.js';

// The source protocol's conformance (_spec/README.md "Data", Protocol): what every adapter, a default one or a host's,
// is held to. Its coverage is well formed and its nodes are nodes the language takes; every listed shape read at level
// 1 over a whole system is an answer in the protocol's shape (a failure too), with a values channel asked or not; an
// aborted read rejects with an `AbortError`; `changes`, when it has it, tells a change made since the opening, and stops
// telling once unsubscribed; dispose ends it.

export interface AdapterCase {
  readonly runtime: FormatRuntime;
  readonly entry: SourceEntry;
  readonly assets: AssetReader;
  /** Edits the data the adapter serves (a host's adapter): `changes` must tell. */
  readonly edit?: () => void;
}

/** Runs the protocol's checks on an adapter; resolves with what it opened (disposed). */
export async function checkAdapter(c: AdapterCase): Promise<OpenedSource> {
  const signal = new AbortController().signal;
  const opened = await c.runtime.open(c.entry, { assets: c.assets, signal, limits: DefaultOpenLimits });
  const cov = opened.coverage;
  expect(checkCoverage(cov)).toBeUndefined();
  for (const { node } of cov.nodes ?? []) {
    const spec: Node = {
      kind: 'spec',
      params: { version: SeqLanguageVersion },
      children: [{ kind: 'data', id: 'data', children: [node] }],
    };
    expect(validateAndNormalize(SeqLanguage, spec, SeqValidateOptions).issues).toEqual([]);
  }
  for (const s of cov.shapes) {
    const length = cov.systems.find((x) => x.name === s.system)?.length ?? 1;
    for (const field of [undefined, 'b_factor', 'no-such-field']) {
      const request = {
        system: s.system,
        start: 0,
        end: length,
        level: 1,
        shapes: [s.shape],
        ...(field ? { field } : {}),
      };
      expect(checkAnswer(await opened.read(request, signal), request)).toBeUndefined();
    }
  }
  if (cov.shapes.length) {
    const s = cov.shapes[0];
    const aborted = new AbortController();
    aborted.abort();
    const request = { system: s.system, start: 0, end: 1, level: 1, shapes: [s.shape] };
    const rejected = await opened.read(request, aborted.signal).then(
      () => undefined,
      (e: unknown) => e,
    );
    expect(isAbortError(rejected)).toBe(true);
  }
  if (opened.changes) {
    let told = 0;
    // a change between the opening and the subscription is told too (an edit while it was being attached)
    c.edit?.();
    const stop = opened.changes(() => told++);
    expect(typeof stop).toBe('function');
    await Promise.resolve();
    if (c.edit) expect(told).toBe(1);
    c.edit?.();
    if (c.edit) expect(told).toBeGreaterThan(0);
    stop();
    const before = told;
    c.edit?.();
    expect(told).toBe(before);
  }
  opened.dispose();
  return opened;
}

/**
 * The region shapes' conformance (_spec/README.md "Data", Protocol: region reads): each listed shape read by region
 * equals the whole read clipped, at every level the coverage lists (a grid by its x and y regions, a sided system on
 * its side); a pair is found from either end. `tileBins` bins a region (small, to cut many).
 */
export async function checkRegionReads(opened: OpenedSource, tileBins = 8): Promise<void> {
  const signal = new AbortController().signal;
  const cov = opened.coverage;
  const lengthOf = (name: string, side?: 'a' | 'b') =>
    cov.systems.find((s) => s.name === name && s.side === side)?.length ?? 0;
  const read = async (request: RegionRequest) => {
    const a = await opened.read(request, signal);
    expect(checkAnswer(a, request)).toBeUndefined();
    return a;
  };
  for (const info of cov.shapes) {
    const sides = cov.systems.filter((s) => s.name === info.system).map((s) => s.side);
    for (const side of sides)
      for (const level of cov.levels) {
        const length = lengthOf(info.system, side);
        const ys =
          info.shape === 'grid'
            ? info.y === '*'
              ? cov.systems.map((s) => s.name)
              : (info.y ?? [info.system])
            : [undefined];
        for (const ySystem of ys) {
          const yLength = ySystem ? lengthOf(ySystem) : 0;
          const y = ySystem ? { system: ySystem, start: 0, end: yLength } : undefined;
          const base = { system: info.system, level, shapes: [info.shape], ...(side && { side }) };
          const whole = await read({ ...base, start: 0, end: length, ...(y && { y }) });
          expect(whole.ok).toBe(true);
          if (!whole.ok) continue;
          const table = whole.tables[info.shape];
          const w = level * tileBins;
          for (let start = 0; start < length; start += w) {
            const end = Math.min(length, start + w);
            // a grid by both axes: every y part, off the origin and across the diagonal too
            for (let ys = 0; ys < (y ? yLength : 1); ys += w) {
              const yr = y && { system: y.system, start: ys, end: Math.min(yLength, ys + w) };
              const part = await read({ ...base, start, end, ...(yr && { y: yr }) });
              expect(part.ok).toBe(true);
              if (!part.ok || !table) continue;
              const cut = { start, end, system: info.system, ...(side && { side }), ...(yr && { y: yr }) };
              const expected = sliceTable(table, level, cut);
              const got = part.tables[info.shape];
              if (tableRows(expected) === 0) expect(got === undefined || tableRows(got) === 0).toBe(true);
              else expect(got).toEqual(expected);
            }
          }
        }
      }
  }
  for (const level of cov.levels) await checkEitherEnd(opened, level, tileBins);
}

/**
 * Every pair read on one end's system is read on the other end's too (EC-62: found from either end), at `level`: a
 * list's pair by its ends and id, a bin pair by its bins; a sided source's `a` end read on side `a`, found on side `b`.
 * At least one pair checked spans two parts of `tileBins` bins (its ends in different tiles), when the data has one.
 */
async function checkEitherEnd(opened: OpenedSource, level: number, tileBins: number): Promise<void> {
  const signal = new AbortController().signal;
  const cov = opened.coverage;
  const sided = cov.systems.some((s) => s.side);
  const w = level * tileBins;
  for (const info of cov.shapes.filter((s) => s.shape === 'pairs')) {
    if (sided && !cov.systems.some((s) => s.name === info.system && s.side === 'a')) continue;
    const length = cov.systems.find((s) => s.name === info.system && (!sided || s.side === 'a'))?.length ?? 0;
    const asked = { system: info.system, start: 0, end: length, level, shapes: ['pairs' as const] };
    const a = await opened.read(sided ? { ...asked, side: 'a' } : asked, signal);
    if (!a.ok || !a.tables.pairs) continue;
    const t = a.tables.pairs;
    const rows = t.form === 'list' ? t.aStart.length : t.count.length;
    let spanning = 0;
    for (let i = 0; i < Math.min(rows, 40); i++) {
      const other = t.bSystem?.get(i) ?? info.system;
      const [as, bs, be] =
        t.form === 'list'
          ? [t.aStart[i], t.bStart[i], t.bEnd[i]]
          : [t.aBin[i] * level, t.bBin[i] * level, (t.bBin[i] + 1) * level];
      if (Math.floor(as / w) !== Math.floor(bs / w)) spanning++;
      const request: RegionRequest = {
        system: other,
        start: bs,
        end: Math.max(be, bs + 1),
        level,
        shapes: ['pairs'],
        ...(sided && { side: 'b' as const }),
      };
      const b = await opened.read(request, signal);
      const found = b.ok ? b.tables.pairs : undefined;
      let hit = false;
      if (found?.form === 'list' && t.form === 'list')
        for (let j = 0; j < found.aStart.length && !hit; j++)
          hit = found.aStart[j] === t.aStart[i] && found.bStart[j] === t.bStart[i] && found.id?.get(j) === t.id?.get(i);
      else if (found?.form === 'bins' && t.form === 'bins')
        for (let j = 0; j < found.count.length && !hit; j++)
          hit = found.aBin[j] === t.aBin[i] && found.bBin[j] === t.bBin[i];
      expect(hit, `pair ${i} at level ${level} found from its other end`).toBe(true);
    }
    if (rows > 40 && length > 2 * w) expect(spanning, `a pair across two parts at level ${level}`).toBeGreaterThan(0);
  }
}

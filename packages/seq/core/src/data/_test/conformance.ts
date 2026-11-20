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
  type SourceEntry,
} from '../protocol.js';
import { checkCoverage } from '../seq-format.js';
import { checkAnswer } from '../tables.js';

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

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `search-matcher` registry entry (_spec/README.md "Context", Registry entries; "Search", Matchers): the built-ins
// unnamespaced and in the default registry (loaded on first use), an application's composed beside them by its
// `<id>:<name>`, a built-in never shadowed, an extension's in its namespace.

import { describe, expect, it } from 'vitest';
import { composeSeqContext } from '../../context/index.js';
import { createSeqContext } from '../../context/default.js';
import { defineExtension } from '../../extension.js';
import {
  BuiltInSearchRegistry,
  compileFailed,
  exactKernel,
  kernelMatcher,
  loadMatcher,
  type SearchCompiled,
} from '../index.js';
import { runSearch, textTarget } from './run.js';

/** An application's matcher: a kernel and nothing else (a fixed motif, any pattern ignored). */
const Motif = kernelMatcher<SearchCompiled>({
  name: 'app:ecori',
  title: 'EcoRI site',
  compile: (pattern, options) => ({ pattern, strands: options.strands ?? 'both' }),
  kernelOf: () => exactKernel(Uint8Array.from('GAATTC', (c) => c.charCodeAt(0))),
});

describe('the search-matcher entry', () => {
  it('the built-ins are in the default registry, lazy, and load on first use', async () => {
    const ctx = createSeqContext();
    expect(Object.keys(ctx.matchers)).toEqual(['degenerate', 'exact', 'regex']);
    const m = await loadMatcher(ctx.matchers.degenerate);
    expect(m.name).toBe('degenerate');
    const c = m.compile('GANTC', {}, { ...(await import('../index.js')).DefaultSearchBounds });
    expect(compileFailed(c)).toBe(false);
    ctx.retire();
  });
  it('composes an application’s matcher beside the built-ins, run by its name', async () => {
    const ctx = composeSeqContext({
      registry: [...BuiltInSearchRegistry, { kind: 'search-matcher', matcher: Motif }],
    });
    expect(Object.keys(ctx.matchers)).toEqual(['degenerate', 'exact', 'regex', 'app:ecori']);
    const m = await loadMatcher(ctx.matchers['app:ecori']);
    const r = await runSearch(m, 'anything', [textTarget('AAGAATTCAA')]);
    expect(r.hits).toEqual([[0, 2, 8, 0]]);
  });
  it('a built-in is never shadowed: an application’s matcher is namespaced', () => {
    const shadow = { ...Motif, name: 'exact' };
    expect(() =>
      composeSeqContext({ registry: [...BuiltInSearchRegistry, { kind: 'search-matcher', matcher: shadow }] }),
    ).toThrow(/different search matchers under one name: exact/);
    expect(() =>
      composeSeqContext({ registry: [{ kind: 'search-matcher', matcher: { ...Motif, name: 'plain' } }] }),
    ).toThrow(/namespaced.*plain/);
  });
  it('an extension adds one as a registry entry, in its own namespace', () => {
    const ok = defineExtension({
      id: 'app',
      version: '1.0.0',
      requires: '0.10',
      registry: [{ kind: 'search-matcher', matcher: Motif }],
    });
    expect(Object.keys(composeSeqContext({ extensions: [ok] }).matchers)).toEqual(['app:ecori']);
    expect(() =>
      defineExtension({
        id: 'other',
        version: '1.0.0',
        requires: '0.10',
        registry: [{ kind: 'search-matcher', matcher: Motif }],
      }),
    ).toThrow(/other: app:ecori/);
  });
});

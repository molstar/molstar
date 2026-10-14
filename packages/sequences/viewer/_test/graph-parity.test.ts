/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer's graph against core's `compile` (seqspec:viewer#updates-and-status; VS-24, VS-25): the graph stages what
// `compile` does node by node, its edges from core's `nodeRefs`. Over every case of core's `cases.ts`, the set it
// assembles holds what `compile`'s does (the same keys in every table), and its errors are `compile`'s. A query that
// names a computed node of another provider reads its answer, as `compile` runs every other computed node first.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import {
  compile,
  RequiredField,
  ref,
  type SeqViewSpecDoc,
  SimpleParamsSchema,
  type StageContext,
} from '@molstar/seq-core/format';
import type { AnnotationInput, SeqSet } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { CASES } from '../../core/_test/format/cases.js';
import { fixtureLoaders, MolstarFixture } from '../../core/_test/format/load-fixtures.js';
import { Q } from '../../core/src/query/index.js';
import { createSeqViewer, DefaultRegistry, DefaultSeqViewerSpec, manualScheduler } from '../src/index.js';
import { InteractionsViewerRegistry } from '../src/interactions/viewer/index.js';
import { fixtureMeasurer } from '../src/testing/index.js';

/** A set's keys, table by table, sorted. */
function keys(set: SeqSet) {
  const of = (list: readonly ({ key: number } | null | undefined)[]) =>
    list.flatMap((x) => (x ? [set.strings.get(x.key)] : [])).sort();
  const units: string[] = [];
  for (let u = 0; u < set.units.count; u++) if (set.units.alive[u]) units.push(set.strings.get(set.units.key[u]));
  return {
    molecules: of(set.molecules),
    units: units.sort(),
    complexes: of(set.complexes),
    alignments: of(set.alignments),
    mappings: of(set.mappings),
    annotations: of(set.annotations),
  };
}
const errors = (issues: readonly { severity: string; code: string }[]) =>
  [...new Set(issues.filter((i) => i.severity === 'error').map((i) => i.code))].sort();

describe("the graph's set against `compile`'s, over core's cases", () => {
  for (const [name, make] of Object.entries(CASES))
    it(name, { timeout: 60_000 }, async () => {
      const ctx = createSeqViewer(
        { ...DefaultSeqViewerSpec(), registry: [DefaultRegistry, { svs: MolstarFixture }, InteractionsViewerRegistry] },
        { scheduler: manualScheduler(), loaders: fixtureLoaders(), text: fixtureMeasurer() },
      );
      try {
        // The same format registries (the interactions extension's provider among them).
        const compiled = await compile(make(), ctx.registry.svs, fixtureLoaders());
        const r = await ctx.setDocument(make());
        await ctx.graph.update();
        const set = ctx.graph.assembled.value?.set;
        expect(set, name).toBeTruthy();
        expect(keys(set as SeqSet), name).toEqual(keys(compiled.set));
        expect(errors(r.issues), name).toEqual(errors(compiled.issues));
      } finally {
        ctx.dispose();
      }
    });
});

describe('a query naming a computed node of another provider', () => {
  // A provider whose options name a sequence: one ranges item [start, stop) on it.
  const span = {
    name: 'test:span',
    options: SimpleParamsSchema({
      of: RequiredField(ref('sequence'), ''),
      start: RequiredField({ name: 'int', check: () => {} } as never, ''),
      stop: RequiredField({ name: 'int', check: () => {} } as never, ''),
    }),
    geometry: 'ranges' as const,
    async compute(_ctx: StageContext, o: unknown) {
      const { of, start, stop } = o as { of: string; start: number; stop: number };
      return {
        key: '',
        geometry: 'ranges',
        columns: { target: [of], start: [start], stop: [stop] },
      } as AnnotationInput;
    },
  };
  const doc = (): SeqViewSpecDoc => {
    const b = createSeqViewSpec();
    b.data.sequence({ ref: 'p', value: 'MACAACAACAAC' });
    b.data.computed({ ref: 'sp', provider: 'test:span', options: { of: 'p', start: 0, stop: 6 } } as never);
    b.data.computed({ ref: 'q', provider: 'query', query: Q.read('(and sp (letters C))') } as never);
    b.view({ ref: 'v', label: 'V' }).track({ label: 'L' }).letters();
    return b.getDocument();
  };
  const answer = (set: SeqSet | undefined) => {
    const a = set?.annotations.find((x) => x && set.strings.get(x.key) === 'q');
    return a ? { count: a.count, start: Array.from(a.start), stop: Array.from(a.stop) } : null;
  };
  it('the graph answers it as `compile` does, with no issue', async () => {
    const ctx = createSeqViewer(
      { ...DefaultSeqViewerSpec(), registry: [DefaultRegistry, { svs: { providers: [span] } } as never] },
      { scheduler: manualScheduler(), loaders: fixtureLoaders(), text: fixtureMeasurer() },
    );
    try {
      const compiled = await compile(doc(), ctx.registry.svs, fixtureLoaders());
      const r = await ctx.setDocument(doc());
      await ctx.graph.update();
      expect(answer(compiled.set)).toEqual({ count: 2, start: [2, 5], stop: [3, 6] });
      expect({ q: answer(ctx.graph.assembled.value?.set), issues: r.issues.map((i) => i.code) }).toEqual({
        q: answer(compiled.set),
        issues: compiled.issues.map((i) => i.code),
      });
    } finally {
      ctx.dispose();
    }
  });
});

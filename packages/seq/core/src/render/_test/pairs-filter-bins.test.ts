/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A pairs layer's data filter on bins (_spec/README.md "Render", The pairs layer, Themes and filter): bins read by region
// carry no fields, so the filter is not applied to them: said once, `seq.render.pairs.filter` (info). Zoomed in to the
// items the filter applies and nothing is said; a layer with no filter says nothing either.

import { SeqRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import type { PairsLay } from '../pairs.js';
import { settle } from './helpers.js';
import { LinksLength as Length, linksSource } from './region-render-fixtures.js';
import { seqToValue } from '../../data/seq-format.js';

const Size = { width: 720, height: 300 };
const Code = 'seq.render.pairs.filter';

function spec(filter: boolean) {
  const b = createSeqSpec({ title: 'links by region, filtered' });
  const src = b.data.source('src', { format: 'seq', value: seqToValue(linksSource(6000).data) as never });
  const chr = b.data.sequence('chr', { alphabet: 'dna', source: src, name: 'chr' });
  const links = b.data.pairs('links', { a: chr, source: src } as never);
  b.views.view('main', (v) => {
    const s = v.section('chr', { of: chr });
    s.track('links').layer('links', {
      type: 'pairs',
      data: filter ? `${links}/(except :of (item :id l1))` : links,
    } as never);
  });
  return b.build();
}

const layOf = (f: { layout: { sections: readonly { sectionLays: readonly { lay: unknown }[] }[] } }) =>
  f.layout.sections[0].sectionLays[0].lay as PairsLay;

describe('a pairs layer’s data filter on bins', () => {
  it('is said once when the frame draws bins', async () => {
    const s = await settle(spec(true), [SeqRuntime]);
    const f = await s.frame({ size: Size, viewport: { start: 0, end: Length } });
    expect(layOf(f).prep!.level).toBeDefined();
    const issues = f.layout.issues.filter((i) => i.code === Code);
    expect(issues).toHaveLength(1);
    expect(issues[0].severity).toBe('info');
  });

  it('is not said on the items (the filter applies) nor with no filter', async () => {
    const s = await settle(spec(true), [SeqRuntime]);
    const f = await s.frame({ size: Size, viewport: { start: 0, end: 400_000 } });
    expect(layOf(f).prep!.level).toBeUndefined();
    expect(f.layout.issues.map((i) => i.code)).not.toContain(Code);
    const t = await settle(spec(false), [SeqRuntime]);
    const g = await t.frame({ size: Size, viewport: { start: 0, end: Length } });
    expect(layOf(g).prep!.level).toBeDefined();
    expect(g.layout.issues.map((i) => i.code)).not.toContain(Code);
  });
});

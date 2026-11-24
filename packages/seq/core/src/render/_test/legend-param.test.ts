/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { viewLegends } from '../index.js';
import { renderSvg } from '../render-svg.js';
import { settle } from './helpers.js';

// A layer's `legend` param (_spec/README.md "Themes", Legends): `none` keys nothing in the legends and the export's band,
// `auto` (and absent) as before; a legend two layers share stays while one of them keys it.

function spec(legend: Record<string, string | undefined>, second: Record<string, string | undefined> = legend) {
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQ' });
  const d = b.data.ranges('d', {
    target: p,
    rows: { start: [0, 6], end: [5, 12], id: ['A', 'B'], category: ['kinase', 'linker'] },
  });
  b.views.view('main', (v) => {
    const s = v.section('p', { of: p });
    s.track('letters').layer('l', { type: 'letters' });
    s.track('domains', { label: 'Domains' }).layer('d', { type: 'ranges', data: d, ...legend });
    s.track('again', { label: 'Domains' }).layer('e', { type: 'ranges', data: d, ...second });
  });
  return b.build();
}

async function legends(params: Record<string, string | undefined>, second = params) {
  const st = await settle(spec(params, second));
  const f = await st.frame({ size: { width: 640, height: 200 }, dpr: 1, settled: true });
  return viewLegends(f.layout).legends;
}

describe('the layer’s legend param', () => {
  it('`auto` and absent key the layer; `none` does not', async () => {
    const absent = await legends({});
    expect((await legends({ legend: 'auto' })).map((l) => l.key)).toEqual(absent.map((l) => l.key));
    expect(absent.some((l) => l.layers.some((k) => k.endsWith('/d')))).toBe(true);
    const none = await legends({ legend: 'none' });
    expect(none.some((l) => l.layers.some((k) => k.endsWith('/d') || k.endsWith('/e')))).toBe(false);
  });

  it('a legend two layers share stays while one of them keys it, keyed by that layer alone', async () => {
    const shared = await legends({ legend: 'none' }, { legend: 'auto' });
    const l = shared.find((x) => x.layers.some((k) => k.endsWith('/e')));
    expect(l).toBeDefined();
    expect(l!.layers.some((k) => k.endsWith('/d'))).toBe(false);
  });

  it('the export’s band omits what the layer says none of', async () => {
    const texts = async (legend: Record<string, string | undefined>) => {
      const r = await renderSvg(spec(legend), { size: { width: 640, height: 200 }, legend: 'below' });
      return [...r.svg.matchAll(/<text[^>]*>([^<]*)<\/text>/g)].map((m) => m[1]);
    };
    expect(await texts({})).toContain('kinase');
    expect(await texts({ legend: 'none' })).not.toContain('kinase');
  });
});

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The legends and the theme commands headless (_spec/README.md "Legends", "Commands"): legends derived on settle only,
// their rows with accessible names; a legend entry selects its items (clicked again: cleared) and hovers them on the
// `legend` channel; a theme edit is one patch (a channel set and unset, several layers at once, a scale, one value's
// color keeping the others'); the scheme picker's fitting schemes; the reader's substitution redraws without a patch.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { recordDrawLists } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it } from 'vitest';
import type { ThemingBehavior } from '../behaviors/theming.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { legendRows } from '../legends.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Layer = '/views/main/p/domains/d';

function spec(): Node {
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLL' });
  const d = b.data.ranges('domains', {
    target: p,
    rows: {
      start: new Int32Array([0, 6, 14]),
      end: new Int32Array([5, 12, 20]),
      id: ['A', 'B', 'C'],
      category: ['binding', 'linker', 'binding'],
    },
  });
  const v1 = b.data.values('v1', { target: p, values: Float32Array.from({ length: 26 }, (_, i) => i) });
  b.views.scale('s', { label: 'Score' });
  b.views.view('main', (v) => {
    const s = v.section('p', { of: p });
    s.track('letters').layer('l', { type: 'letters', color: theme.scheme('zappo') });
    s.track('domains', { label: 'Domains' }).layer('d', { type: 'ranges', data: d });
    s.track('v1').layer('v', { type: 'values', data: v1 });
  });
  return b.build();
}

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];
afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

async function make() {
  frames = installFrames();
  const viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[viewerItem('viewport.maxPxPerPosition'), 1000]],
  });
  viewers.push(viewer);
  await viewer.load(spec());
  viewer.resize(600);
  await viewer.settle();
  frames.flush();
  return viewer;
}

const layerParams = (v: SeqViewer, address: string) => {
  let n = v.spec.root.value;
  for (const id of address.split('/').filter(Boolean)) n = n.children!.find((c) => c.id === id)!;
  return n.params!;
};

describe('legends', () => {
  it('are derived on settle, merged, with their rows’ accessible names', async () => {
    const v = await make();
    const titles = v.legends.legends.value.map((l) => l.title);
    expect(titles).toContain('Domains');
    const rows = legendRows(v, { of: '/views/main/p/domains' });
    expect(rows).toHaveLength(1);
    expect(rows[0].entries.map((e) => e.name)).toEqual(['Domains: binding', 'Domains: linker']);
    // a frame that does not settle derives none
    const n = v.legends.derived;
    v.view.setRange('/views/main/p', 2, 20);
    frames!.flush();
    expect(v.legends.derived).toBe(n);
  });

  it('an entry selects its items (the only one clicked again clears) and hovers them on the legend channel', async () => {
    const v = await make();
    await v.execute('legend.select', { address: Layer, entry: 'binding' });
    expect(v.selection.value.value.elements.keys.map((k) => k.key).sort()).toEqual(['A', 'C']);
    await v.execute('legend.select', { address: Layer, entry: 'linker', toggle: true });
    expect(v.selection.value.value.elements.keys.map((k) => k.key).sort()).toEqual(['A', 'B', 'C']);
    await v.execute('legend.select', { address: Layer, entry: 'linker', toggle: true });
    await v.execute('legend.select', { address: Layer, entry: 'binding' });
    expect(v.selection.value.value.elements.keys).toEqual([]);
    await v.execute('legend.hover', { address: Layer, entry: 'linker' });
    expect(v.highlight.names.value).toContain('legend');
    await v.execute('legend.hover', {});
    expect(v.highlight.names.value).not.toContain('legend');
  });
});

describe('theme commands', () => {
  it('a channel set and unset, several layers at once: one patch each, undone in one step', async () => {
    const v = await make();
    const before = v.spec.root.value;
    await v.execute('layer.channel', {
      addresses: [Layer, '/views/main/p/v1/v'],
      channel: 'color',
      value: '#ff0000',
    });
    expect(layerParams(v, Layer).color).toBe('#ff0000');
    expect(layerParams(v, '/views/main/p/v1/v').color).toBe('#ff0000');
    await v.execute('layer.channel', { address: Layer, channel: 'color' });
    expect(layerParams(v, Layer).color).toBeUndefined();
    expect(before).not.toBe(v.spec.root.value);
  });

  it('a scale’s params set; one value’s color keeps every other value’s', async () => {
    const v = await make();
    await v.execute('scale.set', { address: '/views/s', scale: { domain: [0, 10] } });
    expect(layerParams(v, '/views/s').domain).toEqual([0, 10]);
    const legend = () => v.legends.legends.value.find((l) => l.title === 'Domains')!;
    await v.settle();
    frames!.flush();
    const old = legend().entries;
    await v.execute('theme.value-color', { address: Layer, entry: 'linker', color: '#123456' });
    await v.settle();
    frames!.flush();
    const now = legend().entries;
    if (old.kind !== 'swatches' || now.kind !== 'swatches') throw new Error('swatches');
    expect(now.items.find((i) => i.label === 'linker')!.color).toBe(0x123456ff);
    expect(now.items.find((i) => i.label === 'binding')!.color).toBe(
      old.items.find((i) => i.label === 'binding')!.color,
    );
  });

  it('one value’s color is refused on a scheme or rules color, never replacing it', async () => {
    const v = await make();
    await v.execute('layer.channel', { address: Layer, value: { kind: 'scheme', params: { name: 'okabe-ito' } } });
    const before = layerParams(v, Layer).color;
    const msg = await v.execute('theme.value-color', { address: Layer, entry: 'linker', color: '#123456' });
    expect(msg).toMatch(/not by a field/);
    expect(layerParams(v, Layer).color).toEqual(before);
  });

  it('hides a value’s items in one patch of the layer’s data, the track counting them; the chip’s click clears it', async () => {
    const v = await make();
    const steps = () => v.spec.revision.value;
    const r0 = steps();
    await v.execute('layer.filter', { address: Layer, entry: 'binding' });
    expect(layerParams(v, Layer).data).toBe('/data/domains/(except :of (item :category binding))');
    expect(steps()).not.toBe(r0);
    await v.execute('layer.filter', { address: Layer, entry: 'linker' });
    expect(layerParams(v, Layer).data).toBe('/data/domains/(except :of (item :category [binding linker]))');
    await v.execute('layer.filter', { address: Layer, entry: 'linker' });
    await v.settle();
    frames!.flush();
    const track = v.frame!.layout.sections[0].tracks.find((t) => t.track.label === 'Domains')!;
    const chip = track.header?.states.find((s) => s.name === 'seq:filtered');
    expect(chip?.count).toBe(2);
    // the chip's click: the input clears the filter
    v.input.click(chip!.x + 1, chip!.box.y + chip!.box.height / 2);
    await v.settle();
    expect(layerParams(v, Layer).data).toBe('/data/domains');
  });

  it('offers the schemes that fit a layer', async () => {
    const v = await make();
    const t = v.behaviors.get('legend') as ThemingBehavior;
    const names = (a: string) => t.fittingSchemes(a).map((s) => s.name);
    expect(names('/views/main/p/letters/l')).toContain('zappo');
    expect(names('/views/main/p/letters/l')).not.toContain('viridis');
    expect(names('/views/main/p/v1/v')).toContain('viridis');
    expect(names(Layer)).toContain('okabe-ito');
  });

  it('the reader’s substitution redraws in safe colors, writing nothing into the spec', async () => {
    const v = await make();
    const root = v.spec.root.value;
    // Zappo's negative residues (D, E) are pure red
    const yellow = () => recordDrawLists(v.frame!.lists, ['data']).includes('#ff0000ff');
    expect(yellow()).toBe(true);
    await v.execute('theme.cvd-safe', { on: true });
    await v.settle();
    frames!.flush();
    expect(yellow()).toBe(false);
    expect(v.spec.root.value).toBe(root);
    expect(v.legends.legends.value.some((l) => l.notes.includes('substituted'))).toBe(true);
  });
});

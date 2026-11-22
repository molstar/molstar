/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section past its cap, its edges (core's _spec/README.md "Render", Section scroll; the viewer's "Input", The one
// rail): a minimal composition (no rows behavior) scrolls it by wheel and rail; the section rail in the right gutter's
// slot, a band's rail left of it, the data's right edge hit as data; a header under the ruler or past the box neither
// hit nor listed (two stacked sections); Esc after a follower's rail drag; the stored top within its content when its
// tracks are hidden while scrolled; the keyboard cursor on a rows band below the cap scrolls the section to its row.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { Rail, SectionLayout } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { HoverProvider, NavigateProvider } from '../behaviors/providers.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Genes = '/views/v/genes';
const Below = '/views/v/below';

/** 30 gene tracks under a ruler (a lanes band of 40 lanes second), then a section below; options as named. */
function spec(o: { follow?: boolean; rows?: boolean } = {}): Node {
  const b = createSeqSpec({ title: 'scroll cases' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(100) });
  b.views.view('v', (v) => {
    const lead = o.follow ? v.section('lead', { of: s }) : undefined;
    lead?.track('ruler', { label: '' }).layer('r', { type: 'ruler' });
    const genes = v.section('genes', { of: s, ...(lead ? { follow: lead } : {}) });
    genes.track('ruler', { label: '' }).layer('r', { type: 'ruler' });
    const tall = b.data.ranges('tall', {
      target: s,
      rows: {
        start: Array.from({ length: 40 }, (_, i) => i),
        end: Array.from({ length: 40 }, () => 150),
        id: Array.from({ length: 40 }, (_, i) => `r${i}`),
      },
    });
    genes.track('tall', { label: 'Tall', height: 120 } as never).layer('l', { type: 'ranges', data: tall });
    for (let i = 0; i < 30; i++) {
      const d = b.data.ranges(`g${i}`, { target: s, rows: { start: [1 + i], end: [400], id: [`gene ${i}`] } });
      genes.track(`t${i}`, { label: `Gene ${i}` }).layer('l', { type: 'ranges', data: d });
    }
    v.section('below', { of: s }).track('letters', { label: 'Below' }).layer('l', { type: 'letters' });
  });
  return b.build();
}

/** An alignment's section: a ruler, 18 tracks of ranges on its columns, then its rows (40), cut by the cap. */
function rowsSpec(): Node {
  const b = createSeqSpec({ title: 'rows below the cap' });
  const labels = Array.from({ length: 40 }, (_, i) => `R${i}`);
  const aln = b.data.alignment('aln', {
    alphabet: 'protein',
    rows: { label: labels, letters: labels.map(() => 'MVLSPA') },
  });
  b.views.view('v', (v) => {
    const s = v.section('s', { of: aln });
    s.track('ruler', { label: '' }).layer('r', { type: 'ruler' });
    for (let i = 0; i < 18; i++) {
      const d = b.data.ranges(`f${i}`, { target: aln, rows: { start: [0], end: [4], id: [`f${i}`] } });
      s.track(`t${i}`, { label: `F${i}` }).layer('l', { type: 'ranges', data: d });
    }
    s.track('rows', { label: '', height: 200 }).layer('r', { type: 'rows', data: aln });
  });
  return b.build();
}

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  frames?.restore();
  frames = undefined;
  vi.useRealTimers();
});

async function make(n: Node, minimal = false) {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  frames = installFrames();
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({
    ...base,
    // the minimal composition: hover and navigation alone (no rows behavior, no cursor)
    ...(minimal ? { behaviors: [{ provider: HoverProvider }, { provider: NavigateProvider }] } : {}),
    ui: undefined,
    config: [[viewerItem('view.jumpMs'), 0]],
  });
  viewer = v;
  await v.load(n);
  v.resize(700);
  await v.settle();
  draw(v);
  v.input.mac = false;
  return v;
}

const draw = (v: SeqViewer) => {
  frames!.flush();
  v.draw();
};
const section = (v: SeqViewer, a = Genes): SectionLayout =>
  v.frame!.layout.sections.find((s) => s.section.address === a)!;
const yOf = (r: Rail, u: number) => r.x + (u / r.total) * r.width;
const max = (s: SectionLayout) => s.placement.contentHeight - s.placement.box.height;

describe('a section past its cap in every composition', () => {
  it('the minimal composition scrolls it to its last track by wheel and by its rail', async () => {
    const v = await make(spec(), true);
    const s = section(v);
    expect(s.sectionRail).toBeDefined();
    const gene = s.tracks[10].box;
    const y = s.placement.box.y + s.placement.box.height - 20;
    for (let i = 0; i < 300; i++) v.input.wheel(gene.x + 50, y, { deltaX: 0, deltaY: 50, deltaMode: 0 });
    draw(v);
    expect(section(v).placement.top).toBeCloseTo(max(section(v)), 6);
    // back to the top by the rail: its window dragged up the whole rail
    const r = section(v).sectionRail!;
    const x = r.top + 6;
    const y0 = yOf(r, (r.start + r.end) / 2);
    expect(v.input.down(x, y0)).toBe(true);
    v.input.move(x, y0 - r.width);
    v.input.up();
    draw(v);
    expect(section(v).placement.top).toBe(0);
  });
});

describe('the rails’ places', () => {
  it('the section rail in the right gutter’s slot, a band’s rail left of it; the data’s right edge hits the data', async () => {
    const v = await make(spec());
    const s = section(v);
    const sr = s.sectionRail!;
    const lanes = s.lanesRails![0];
    const tall = s.tracks[1].box;
    // the section's at the page's right edge past the data; the band's inside its band, left of it
    expect(sr.top).toBeGreaterThanOrEqual(tall.x + tall.width);
    expect(lanes.bottom).toBeLessThanOrEqual(sr.top);
    expect(lanes.top).toBeGreaterThanOrEqual(tall.x);
    // a gene track's last data px is its data, never a rail
    const end = s.placement.box.y + s.placement.box.height;
    const gene = s.tracks.find((t) => /\/t\d+$/.test(t.track.address) && t.box.y + 20 < end)!.box;
    const h = v.hit(gene.x + gene.width - 1, gene.y + 10);
    expect(h).toBeDefined();
    expect(h!.navigator).toBeUndefined();
    expect(h!.layer).toMatch(/\/genes\/t\d+\//);
  });
});

describe('a scrolled section’s headers', () => {
  it('a header past the box is neither hit nor listed; the section below answers there', async () => {
    const v = await make(spec());
    const s = section(v);
    const below = section(v, Below).tracks[0];
    // a hidden gene track lies where the section below is: its line is not hit there, the section below's is
    const hidden = s.tracks.find((t) => t.well.y > s.placement.box.y + s.placement.box.height);
    expect(hidden).toBeDefined();
    expect(hidden!.header).toBeUndefined();
    const x = v.frame!.layout.config.inset + 20;
    const h = v.hit(x, below.labelY - 4);
    expect(h).toBeDefined();
    expect(h!.layer.startsWith(Genes)).toBe(false);
    expect(v.headers!.value.some((g) => g.target.track === hidden!.track.address)).toBe(false);
  });

  it('scrolled, a header under the ruler is neither hit nor listed', async () => {
    const v = await make(spec());
    v.view.setTop(Genes, 300);
    draw(v);
    const s = section(v);
    const ruler = s.tracks[0].well;
    const under = s.tracks.filter((t) => t.well.y + t.well.height <= s.placement.scrollY! && t !== s.tracks[0]);
    expect(under.length).toBeGreaterThan(0);
    for (const t of under) expect(t.header).toBeUndefined();
    // the ruler there is the ruler's (its rail's grab zone: a hit, never a hidden band's); its label column no hidden
    // header's
    const d = v.hit(ruler.x + 40 + v.frame!.layout.config.labelWidth, ruler.y + ruler.height / 2);
    expect(d).toBeDefined();
    expect(d!.layer).not.toMatch(/\/genes\/(t\d+|tall)/);
    for (const t of under) expect(v.headers!.value.some((g) => g.target.track === t.track.address)).toBe(false);
    const h = v.hit(v.frame!.layout.config.inset + 20, ruler.y + ruler.height / 2);
    if (h) expect(h.layer).not.toMatch(/\/genes\/(t\d+|tall)/);
  });
});

describe('Esc after a follower’s rail drag', () => {
  it('puts its scroll back where the press found it', async () => {
    const v = await make(spec({ follow: true }));
    expect(v.view.isFollowing(Genes)).toBe(true);
    const r = section(v).sectionRail!;
    const x = r.top + 6;
    const y0 = yOf(r, (r.start + r.end) / 2);
    v.input.down(x, y0);
    v.input.move(x, y0 + 60);
    expect(v.view.get(Genes)?.top ?? 0).toBeGreaterThan(0);
    v.input.key({ key: 'Escape' });
    expect(v.view.get(Genes)?.top ?? 0).toBe(0);
  });
});

describe('the stored top', () => {
  it('stays within the content when tracks are hidden while scrolled', async () => {
    const v = await make(spec());
    v.view.setTop(Genes, 1e6);
    draw(v);
    const before = section(v).placement.top!;
    expect(v.view.get(Genes)?.top).toBeCloseTo(before, 6);
    for (let i = 0; i < 20; i++) await v.execute('track.hide', { address: `${Genes}/t${i}` });
    await v.settle();
    draw(v);
    const s = section(v);
    expect(s.placement.top!).toBeLessThan(before);
    expect(v.view.get(Genes)?.top).toBeCloseTo(s.placement.top!, 6);
    // a wheel up moves at once (no dead travel back from a stale top)
    const gene = s.tracks[3].box;
    expect(
      v.input.wheel(gene.x + 50, s.placement.box.y + s.placement.box.height - 20, {
        deltaX: 0,
        deltaY: -10,
        deltaMode: 0,
      }),
    ).toBe(true);
  });
});

describe('the keyboard cursor below the cap', () => {
  it('a cursor on a rows band below the cap scrolls the section to its row', async () => {
    const v = await make(rowsSpec());
    const s = section(v, '/views/v/s');
    expect(s.sectionRail).toBeDefined();
    expect(s.placement.top).toBe(0);
    v.input.move(s.tracks[3].box.x + 20, s.tracks[3].box.y + 4);
    // the first key places the cursor on the first row shown: the section scrolls to show it
    v.input.key({ key: 'ArrowDown' });
    frames!.flush();
    for (let i = 0; i < 12; i++) {
      v.input.key({ key: 'ArrowDown' });
      draw(v);
    }
    const after = section(v, '/views/v/s');
    expect(v.selection.cursor.value?.ySystem).toBeDefined();
    expect(after.placement.top!).toBeGreaterThan(0);
    // the cursor's row inside the box
    const l = after.tracks[after.tracks.length - 1].layers[0].lay as { rowY(i: number): number; rowScale: number };
    const i = v.selection.cursor.value!.yPosition!;
    expect(l.rowY(i) + l.rowScale).toBeLessThanOrEqual(after.placement.box.y + after.placement.box.height + 1e-6);
  });
});

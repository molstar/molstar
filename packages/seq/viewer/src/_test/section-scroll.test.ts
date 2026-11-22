/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The one rail over every box that scrolls (core's _spec/README.md "Render", Section scroll; the viewer's "Input", The
// one rail): a rows band, a lanes band, a text past its box and a tracks section past its cap answer the same gesture
// the same way (a drag on the window, a press beside it, Escape, Page and Home/End keys with the pointer in the box);
// a section of 30 gene tracks scrolls by wheel, rail and keys below its fixed ruler at 600, a host's 300 and a spec
// `height` of 200; the wheel's chain (a lanes band, then the section, then the page); an export shows every track.

import type { Node } from '@molstar/core/repr/tree';
import { SeqCoreConfig } from '@molstar/seqx-core/context';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { Rail, SectionLayout } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { scrollerOf, type ScrollerKind } from '../behaviors/scrollers.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Genes = '/views/v/genes';
const Lanes = '/views/v/lanes/tall';
const Text = '/views/v/text';
const Rows = '/views/v/aln/rows';

/** A tracks section (a ruler, `n` gene tracks, one of them a lanes band of 40 lanes), a text, an alignment's rows. */
function spec(o: { n?: number; height?: number; all?: boolean } = {}): Node {
  const b = createSeqSpec({ title: 'scrollers' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(100) });
  b.views.view('v', (v) => {
    const genes = v.section('genes', { of: s, ...(o.height !== undefined ? { height: o.height } : {}) });
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
    for (let i = 0; i < (o.n ?? 30); i++) {
      const d = b.data.ranges(`g${i}`, { target: s, rows: { start: [1 + i], end: [200 + i], id: [`gene ${i}`] } });
      genes.track(`t${i}`, { label: `Gene ${i}` }).layer('l', { type: 'ranges', data: d });
    }
    if (!o.all) return;
    v.section('lanes', { of: s, heights: 'fixed' } as never)
      .track('tall', { height: 120 } as never)
      .layer('l', { type: 'ranges', data: tall });
    // a sequence of its own: sections of one system move together
    const t = b.data.sequence('t', { alphabet: 'dna', letters: 'ACGT'.repeat(100) });
    v.section('text', { of: t, layout: 'text', max_width: 200, height: 90 })
      .track('letters')
      .layer('l', { type: 'letters' });
    const aln = b.data.alignment('aln', {
      alphabet: 'protein',
      rows: {
        label: Array.from({ length: 40 }, (_, i) => `a${i}`),
        letters: Array.from({ length: 40 }, () => 'MEEPMEEP'),
      },
    });
    v.section('aln', { of: aln })
      .track('rows', { height: 200 })
      .layer('r', { type: 'rows', data: aln } as never);
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

async function make(o: { n?: number; height?: number; all?: boolean; maxPx?: number } = {}) {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  frames = installFrames();
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [
      [viewerItem('view.jumpMs'), 0],
      ...(o.maxPx !== undefined ? [[SeqCoreConfig.render, { section: { maxPx: o.maxPx } }] as const] : []),
    ],
  });
  viewer = v;
  await v.load(spec(o));
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
const section = (v: SeqViewer, a: string): SectionLayout =>
  v.frame!.layout.sections.find((s) => s.section.address === a)!;

/** Each scroller's rail in the last frame, its kind and key. */
function railOf(v: SeqViewer, kind: ScrollerKind): { rail: Rail; key: string } {
  const ls = v.frame!.layout.sections;
  if (kind === 'rows') return { rail: ls.find((s) => s.rowsRails)!.rowsRails![0], key: Rows };
  if (kind === 'lanes') {
    const r = section(v, '/views/v/lanes').lanesRails![0];
    return { rail: r, key: Lanes };
  }
  const a = kind === 'section' ? Genes : Text;
  return { rail: section(v, a).sectionRail!, key: a };
}

const yOf = (r: Rail, u: number) => r.x + (u / r.total) * r.width;

describe('the one rail answers alike over every box that scrolls', () => {
  // the text's scroller is a section's too: its rail the same, its top moved through its anchor
  const table: [string, ScrollerKind, 'section' | 'text' | undefined][] = [
    ['a rows band', 'rows', undefined],
    ['a lanes band', 'lanes', undefined],
    ['a text past its box', 'section', 'text'],
    ['a tracks section past its cap', 'section', 'section'],
  ];
  for (const [name, kind, which] of table)
    it(`${name}: a drag, Escape, a press beside the window, Page and Home/End keys`, async () => {
      const v = await make({ all: true });
      const { rail, key } = which === 'text' ? railOf(v, 'text' as never) : railOf(v, kind);
      const sc = () => scrollerOf(v.view, v.frame!.layout, kind, key)!;
      expect(rail.active).toBe(true);
      expect(sc().total).toBeGreaterThan(sc().shown);
      const x = rail.top + 6;
      // a drag on the window by a quarter of its room scrolls by it (units: rows, lanes, px)
      const room = sc().total - sc().shown;
      const y0 = yOf(rail, rail.start + (rail.end - rail.start) / 2);
      expect(v.input.down(x, y0)).toBe(true);
      v.input.move(x, y0 + (room / 4 / rail.total) * rail.width);
      expect(Math.abs(sc().top() - room / 4)).toBeLessThan(room * 0.02);
      // Escape puts it back
      v.input.key({ key: 'Escape' });
      draw(v);
      expect(sc().top()).toBeCloseTo(0, 4);
      // a press beside the window centres it there
      const r1 = railOf(v, kind === 'section' && which === 'text' ? ('text' as never) : kind).rail;
      const at = r1.total * 0.75;
      v.input.down(x, yOf(r1, at));
      v.input.up();
      v.input.click(x, yOf(r1, at), { detail: 1 });
      draw(v);
      expect(sc().top()).toBeCloseTo(Math.min(sc().total - sc().shown, at - sc().shown / 2), 3);
      // the keys with the pointer on the rail: Home, Page Down, End
      v.input.move(x, yOf(r1, 0.1));
      expect(v.input.key({ key: 'Home' })).toBe(true);
      draw(v);
      expect(sc().top()).toBeCloseTo(0, 4);
      expect(v.input.key({ key: 'PageDown' })).toBe(true);
      draw(v);
      expect(sc().top()).toBeCloseTo(Math.min(sc().page, sc().total - sc().shown), 3);
      expect(v.input.key({ key: 'End' })).toBe(true);
      draw(v);
      expect(sc().top()).toBeCloseTo(sc().total - sc().shown, 3);
      // at the end the key is the page's
      expect(v.input.key({ key: 'End' })).toBe(false);
      // session state, every rail's (the rows rail's press beside its window too: a scroll, no step)
      expect(v.view.history.value.back).toBe(0);
    });
});

describe('a section of 30 gene tracks past its cap', () => {
  for (const [name, o, box] of [
    ['the default 600', {}, 600],
    ['a host’s 300', { maxPx: 300 }, 300],
    ['a spec `height` of 200', { height: 200 }, 200],
  ] as const)
    it(`${name}: scrolls by wheel, rail and keys below its fixed ruler`, async () => {
      const v = await make(o);
      const s0 = section(v, Genes);
      expect(s0.placement.box.height).toBe(box);
      const ruler = s0.tracks[0].box.y;
      const gene = s0.tracks[10].box;
      // the wheel over a gene track scrolls the section
      expect(v.input.wheel(gene.x + 50, s0.placement.box.y + box - 20, { deltaX: 0, deltaY: 30, deltaMode: 0 })).toBe(
        true,
      );
      draw(v);
      const s1 = section(v, Genes);
      expect(s1.placement.top).toBe(30);
      expect(s1.tracks[0].box.y).toBe(ruler);
      expect(s1.tracks[10].box.y).toBeCloseTo(gene.y - 30, 6);
      // the rail's window follows
      expect(s1.sectionRail!.start).toBe(30);
      // End with the pointer in it, then the wheel down is the page's
      v.input.move(gene.x + 50, s1.placement.box.y + box - 20);
      expect(v.input.key({ key: 'End' })).toBe(true);
      draw(v);
      const s2 = section(v, Genes);
      expect(s2.placement.top).toBeCloseTo(s2.placement.contentHeight - box, 6);
      expect(v.input.wheel(gene.x + 50, s2.placement.box.y + box - 20, { deltaX: 0, deltaY: 30, deltaMode: 0 })).toBe(
        false,
      );
      expect(s2.tracks[0].box.y).toBe(ruler);
    });

  it('the wheel’s chain: a lanes band first, at its end the section, at the section’s the page', async () => {
    const v = await make();
    const s = section(v, Genes);
    const tall = s.tracks[1].box;
    const scroll = () => v.input.wheel(tall.x + 50, tall.y + 40, { deltaX: 0, deltaY: 400, deltaMode: 0 });
    expect(scroll()).toBe(true);
    draw(v);
    // the lanes moved, the section not yet
    expect(v.view.band('/views/v/genes/tall')?.bandTop).toBeGreaterThan(0);
    expect(section(v, Genes).placement.top).toBe(0);
    // the lanes at their end: the section scrolls (the band moved up under the pointer: aim where it is now)
    for (let i = 0; i < 4; i++) {
      scroll();
      draw(v);
    }
    expect(section(v, Genes).placement.top).toBeGreaterThan(0);
    // at the section's end the wheel is the page's
    for (let i = 0; i < 40; i++) {
      scroll();
      draw(v);
    }
    expect(section(v, Genes).placement.top).toBeCloseTo(
      section(v, Genes).placement.contentHeight - section(v, Genes).placement.box.height,
      6,
    );
    const s1 = section(v, Genes);
    const at = s1.tracks.find((t) => t.box.y > s1.placement.scrollY! + 10)!.box;
    expect(v.input.wheel(at.x + 50, at.y + 5, { deltaX: 0, deltaY: 400, deltaMode: 0 })).toBe(false);
  });

  it('an export shows every track, the live scroll untouched', async () => {
    const v = await make();
    v.view.setTop(Genes, 100);
    draw(v);
    const r = await v.renderSvg();
    expect(r.svg).toContain('>Gene 29<');
    expect(r.svg).toContain('>Gene 0<');
    expect(section(v, Genes).placement.top).toBe(100);
  });
});

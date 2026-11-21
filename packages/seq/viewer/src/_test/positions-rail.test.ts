/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The positions rail beside a declared navigator, and its place (_spec/README.md "Composition", Navigators:
// `navigator.positionsRail`, `navigator.railPlace`): headless on a fake clock. `always` keeps the rail beside a
// declared overview, `unlessNavigator` the old replacement, `none` no rail; its row below the ruler, under the last
// track, or by a declared horizontal navigator; a drag on it pans and the overview's window follows in the same frame.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, overviewNavigator } from '@molstar/seqx-core/lang';
import { SlimRailPx } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { installFrames, type Frames } from './frames.js';

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

const Letters = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAA';

/** A 300-column alignment's rows with a ruler over them; a declared overview (`navigator`: its side) when asked. */
function spec(options: { navigator?: 'top' | 'bottom' | 'right'; ruler?: boolean; fill?: boolean } = {}): Node {
  const b = createSeqSpec({ title: 'positions rail' });
  const rows = Array.from({ length: options.fill ? 80 : 6 }, (_, i) => Letters.repeat(5).slice(i % 10, (i % 10) + 300));
  const aln = b.data.alignment('aln', {
    alphabet: 'protein',
    rows: { label: rows.map((_, i) => `r${i}`), letters: rows },
  });
  b.views.view('v', (v) => {
    const s = v.section('s', { of: aln });
    if (options.ruler !== false) s.track('ruler').layer('ruler', { type: 'ruler' });
    s.track('rows', options.fill ? { label: '' } : { label: '', height: 60 }).layer('r', { type: 'rows', data: aln });
    if (options.navigator) v.navigator('nav', { ...overviewNavigator(s), axes: 'both', place: options.navigator });
  });
  return b.build();
}

async function make(root: Node, config: readonly (readonly [unknown, unknown])[] = []) {
  frames = installFrames();
  const viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[viewerItem('view.jumpMs'), 0], ...config] as never,
  });
  viewers.push(viewer);
  await viewer.load(root);
  viewer.resize(800);
  await viewer.settle();
  frames.flush();
  viewer.draw();
  return viewer;
}

const section = (v: SeqViewer) => v.frame!.layout.sections[0];

describe('the positions rail beside a declared navigator', () => {
  it('is drawn beside a declared overview by default (`always`), replaced by it under `unlessNavigator`, none under `none`', async () => {
    const both = await make(spec({ navigator: 'right' }));
    expect(section(both).rail).toBeDefined();
    const old = await make(spec({ navigator: 'right' }), [[viewerItem('navigator.positionsRail'), 'unlessNavigator']]);
    expect(section(old).rail).toBeUndefined();
    const none = await make(spec(), [[viewerItem('navigator.positionsRail'), 'none']]);
    expect(section(none).rail).toBeUndefined();
    // with no navigator declared, `unlessNavigator` still draws it
    const plain = await make(spec(), [[viewerItem('navigator.positionsRail'), 'unlessNavigator']]);
    expect(section(plain).rail).toBeDefined();
  });

  it('a drag on the rail pans, and the overview’s window follows in the same frame', async () => {
    const v = await make(spec({ navigator: 'right' }));
    v.view.setRange(section(v).rail!.system, 1, 100);
    v.draw();
    const nav = v.navigatorModel('nav')!;
    nav.resize(96, 200);
    const before = { ...nav.draw()!.windows[0] };
    const r = section(v).rail!;
    expect(r.active).toBe(true);
    const start0 = r.start;
    const k = r.width / r.total;
    const x = r.x + 50 * k;
    const y = r.bottom - 2;
    expect(v.input.down(x, y)).toBe(true);
    v.input.move(x + 40 * k, y);
    v.input.up();
    v.draw();
    const after = nav.draw()!.windows[0];
    expect(section(v).rail!.start).toBeGreaterThan(start0);
    expect(after.dx0).toBeGreaterThan(before.dx0);
  });
});

describe('where the rail sits (`navigator.railPlace`)', () => {
  const place = (value: string) => [[viewerItem('navigator.railPlace'), value]] as const;
  const below = (s: ReturnType<typeof section>) => s.tracks[0].box.y + s.tracks[0].box.height;

  it('`ruler` (the default): just below the ruler’s bottom border', async () => {
    const s = section(await make(spec()));
    expect(s.rail!.bottom).toBe(below(s) + SlimRailPx);
  });

  it('`bottom`: a row under the last track', async () => {
    const s = section(await make(spec(), place('bottom')));
    const last = s.tracks[s.tracks.length - 1].box;
    expect(s.rail!.slim).toBe(true);
    expect(s.rail!.top).toBeGreaterThanOrEqual(last.y + last.height);
    expect(s.rail!.bottom - s.rail!.top).toBe(SlimRailPx);
  });

  it('`navigator`: under a declared top navigator (above the first track); beside a bottom one at the bottom; else the ruler', async () => {
    const top = section(await make(spec({ navigator: 'top' }), place('navigator')));
    expect(top.rail!.slim).toBe(true);
    expect(top.rail!.bottom).toBeLessThanOrEqual(top.tracks[0].box.y);
    const bottom = section(await make(spec({ navigator: 'bottom' }), place('navigator')));
    const last = bottom.tracks[bottom.tracks.length - 1].box;
    expect(bottom.rail!.top).toBeGreaterThanOrEqual(last.y + last.height);
    // a vertical overview is no horizontal navigator, and none declared: the ruler's place
    const side = section(await make(spec({ navigator: 'right' }), place('navigator')));
    expect(side.rail!.bottom).toBe(below(side) + SlimRailPx);
    const none = section(await make(spec(), place('navigator')));
    expect(none.rail!.bottom).toBe(below(none) + SlimRailPx);
  });

  it('a section with no ruler takes the bottom row under `ruler`', async () => {
    const s = section(await make(spec({ ruler: false })));
    const last = s.tracks[s.tracks.length - 1].box;
    expect(s.rail!.top).toBeGreaterThanOrEqual(last.y + last.height);
  });
});

describe('a list of places (`navigator.railPlace` as an array)', () => {
  const place = (value: readonly string[]) => [[viewerItem('navigator.railPlace'), value]] as const;

  it('`[ruler, bottom]` draws two rails on the one window; a drag on the bottom one moves the top one in the same frame', async () => {
    const v = await make(spec(), place(['ruler', 'bottom']));
    v.view.setRange(section(v).rail!.system, 1, 100);
    v.draw();
    const rails = section(v).rails!;
    expect(rails.length).toBe(2);
    const [ruler, bottom] = rails;
    expect(ruler.slim).toBe(false);
    expect(bottom.slim).toBe(true);
    expect(bottom.top).toBeGreaterThan(ruler.bottom);
    expect(bottom.active).toBe(true);
    const start0 = ruler.start;
    const k = bottom.width / bottom.total;
    const x = bottom.x + 50 * k;
    const y = bottom.bottom - 2;
    expect(v.input.down(x, y)).toBe(true);
    v.input.move(x + 40 * k, y);
    v.input.up();
    v.draw();
    const after = section(v).rails!;
    expect(after[0].start).toBeGreaterThan(start0);
    expect(after[1].start).toBe(after[0].start);
    expect(after[1].end).toBe(after[0].end);
  });

  it('`pinned` is a 12 px row at the frame’s foot at any scroll of the rows, drags the window, and is no row of an export', async () => {
    const tall = [
      [viewerItem('canvas.height'), 300],
      [viewerItem('navigator.railPlace'), ['ruler', 'pinned']],
    ] as const;
    const v = await make(spec({ fill: true }), tall);
    const track = '/views/v/s/rows';
    const pinned = () => section(v).rails!.find((r) => r.pinned)!;
    v.view.setRange(section(v).rail!.system, 1, 100);
    v.view.setBand(track, { bandTop: 0 });
    v.draw();
    expect(section(v).rails!.length).toBe(2);
    expect(pinned().bottom).toBe(300);
    expect(pinned().bottom - pinned().top).toBe(SlimRailPx);
    // the rows keep clear of it
    const rows = section(v).tracks.at(-1)!.box;
    expect(rows.y + rows.height).toBeLessThanOrEqual(pinned().top);
    v.view.setBand(track, { bandTop: 500 });
    v.draw();
    expect(pinned().bottom).toBe(300);
    // a drag on it pans the window, the rail below the ruler follows
    const [ruler] = section(v).rails!;
    const start0 = ruler.start;
    const p = pinned();
    const k = p.width / p.total;
    const x = p.x + 50 * k;
    const y = p.bottom - 2;
    expect(v.input.down(x, y)).toBe(true);
    v.input.move(x + 40 * k, y);
    v.input.up();
    v.draw();
    expect(section(v).rails![0].start).toBeGreaterThan(start0);
  });

  it('an export has no row for a pinned rail', async () => {
    const config = (place: readonly string[]) =>
      [
        [viewerItem('canvas.height'), 300],
        [viewerItem('navigator.railPlace'), place],
      ] as const;
    const pinned = await make(spec({ fill: true }), config(['ruler', 'pinned']));
    const plain = await make(spec({ fill: true }), config(['ruler']));
    expect((await pinned.renderSvg({ navigators: 'include' })).svg.replace(/s[0-9a-z]+-/g, '')).toBe(
      (await plain.renderSvg({ navigators: 'include' })).svg.replace(/s[0-9a-z]+-/g, ''),
    );
  });

  it('`[ruler, ruler]`, and a fallback onto a listed place, give one rail', async () => {
    expect(section(await make(spec(), place(['ruler', 'ruler']))).rails!.length).toBe(1);
    // no ruler: `ruler` takes the bottom row, which `bottom` then repeats
    expect(section(await make(spec({ ruler: false }), place(['ruler', 'bottom']))).rails!.length).toBe(1);
    // no navigator declared: `navigator` is the ruler's place
    expect(section(await make(spec(), place(['navigator', 'ruler']))).rails!.length).toBe(1);
  });
});

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A host-drawn navigator (_spec/README.md "Models", NavigatorModel, the contract): a host's picture registered with the
// viewer, its navigator driven through the model's methods alone (no canvas, no DOM: what a host's own component
// calls): its frame's windows, a drag, a click, the keys, the wheel, Esc, the readout its picture gives, the cursor,
// the sliders' values; and a declared navigator's picture of an extension's layer type drawn as that type.

import type { Node } from '@molstar/core/repr/tree';
import { defineExtension, type LayerType, type NavigatorPicture } from '@molstar/seqx-core/extension';
import { createSeqSpec, stripNavigator } from '@molstar/seqx-core/lang';
import type { NavigatorRenderer } from '@molstar/seqx-core/render/navigator';
import { afterEach, describe, expect, it } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import type { NavigatorModel } from '../navigator.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import type { SeqViewerSpec } from '../viewer-spec.js';
import { type Frames, installFrames } from './frames.js';

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

const Length = 400;
const P = '/views/v/protein';

/** A coverage strip: the values its `data` param names, filling the navigator; its readout says the position. */
const Coverage: NavigatorPicture = {
  name: 'host:coverage',
  axes: ['positions'],
  size: 24,
  fills: true,
  picture: (c) => {
    const data = c.navigator.params.data;
    if (c.zoom || typeof data !== 'string') return undefined;
    const track: Node = {
      kind: 'track',
      id: 'coverage',
      params: { label: '', height: c.height },
      children: [{ kind: 'layer', id: 'coverage', params: { type: 'values', data } }],
    };
    return { tracks: [track] };
  },
  readout: (p) => `coverage ${p + 1}`,
};

/** An extension's layer type that draws nothing (its fallback `none`): placed only where the extension is composed. */
const Dots: LayerType<object, { key: string }> = {
  name: 'demo:dots',
  shape: 'ranges',
  prepare: () => ({}),
  layout: (_p, ctx) => ({ key: ctx.key }),
  draw: () => {},
  hit: () => undefined,
};
const Demo = defineExtension({
  id: 'demo',
  version: '1.0.0',
  requires: '0.5',
  registry: [{ kind: 'layer-type', type: Dots }],
});

/** A protein with a profile and its domains (as `demo:dots` where asked), and a strip of `content`. */
function spec(content?: string, dots = false): Node {
  const b = createSeqSpec({ title: 'host navigator' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'ACDEFGHIKL'.repeat(Length / 10) });
  b.data.values('profile', { target: p, values: Array.from({ length: Length }, (_, i) => (i % 40) / 40) });
  const d = b.data.ranges('domains', { target: p, rows: { start: [10, 200], end: [120, 330] } });
  b.views.view('v', (v) => {
    const s = v.section('protein', { of: p });
    s.track('letters').layer('letters', { type: 'letters' });
    s.track('domains').layer('d', { type: dots ? 'demo:dots' : 'ranges', data: d, fallback: 'none' } as never);
    if (content) v.navigator('strip', stripNavigator(s, { content: content as 'tracks' }));
  });
  return b.build();
}

async function make(root: Node, extra: Partial<SeqViewerSpec> = {}) {
  frames = installFrames();
  const base = DefaultSeqViewerSpec();
  const viewer = createSeqViewer({
    ...base,
    ui: undefined,
    registry: [...(base.registry ?? []), { kind: 'navigator-picture', picture: Coverage }],
    config: [[viewerItem('viewport.maxPxPerPosition'), 32]] as never,
    ...extra,
  });
  viewers.push(viewer);
  await viewer.load(root);
  viewer.resize(800);
  await viewer.settle();
  frames.flush();
  viewer.draw();
  return viewer;
}

/** Draws the navigator over the viewer's last frame, its wants answered: what a host does after each viewer frame. */
async function drawn(v: SeqViewer, m: NavigatorModel) {
  m.draw();
  await v.settle();
  frames!.flush();
  v.draw();
  return m.draw()!;
}

const range = (v: SeqViewer) => [v.view.get(P)!.start, v.view.get(P)!.end] as const;

describe('a host-drawn navigator through NavigatorModel', () => {
  it('draws the host’s picture, and the gestures, keys and readout run through the model alone', async () => {
    const v = await make(spec());
    const m = v.navigator({ id: 'host', of: P, content: 'host:coverage', params: { data: '/data/profile' } });
    m.resize(400, 0);
    v.view.setRange(P, 100, 140);
    v.draw();
    const f = await drawn(v, m);
    // the picture's size (no `size`: its entry's), its tracks, the window over the view's range
    expect([m.width, m.height]).toEqual([400, 24]);
    const pic = (m as unknown as { renderer: NavigatorRenderer }).renderer.pictureSection!;
    expect(pic.tracks.map((t) => t.layers[0].type.name)).toEqual(['values']);
    const w = f.windows[0];
    expect(w.shown).toBe(true);
    const px = (w.x1 - w.x0) / 40;
    // a hover: the picture's readout outside the window, the window's text over it; the cursor by the part
    m.move(w.dx1 + 100 * px, 12);
    expect(m.readout.value).toMatch(/^coverage \d+$/);
    expect(m.cursor).toBe('pointer');
    m.move((w.dx0 + w.dx1) / 2, 12);
    expect(m.readout.value).toContain('101');
    expect(m.cursor).toBe('grab');
    // a drag moves the view with the hand
    expect(m.down((w.dx0 + w.dx1) / 2, 12)).toBe(true);
    m.move((w.dx0 + w.dx1) / 2 + 50 * px, 12);
    m.up();
    expect(range(v)[0]).toBeCloseTo(150, 0);
    // a press that never moves centres at its click
    v.draw();
    const g = await drawn(v, m);
    const at = g.windows[0].x0 + 200 * px;
    m.down(at, 12);
    m.up();
    m.click(at, 12, 1);
    const [s, e] = range(v);
    expect((s + e) / 2).toBeCloseTo(350, -1);
    // the keys, after the host says the part has the focus
    m.setFocus(true);
    v.draw();
    await drawn(v, m);
    const before = range(v)[0];
    expect(m.key('ArrowLeft')).toBe(true);
    expect(range(v)[0]).toBeLessThan(before);
    expect(m.focusNext()).toBe(false);
    // the wheel zooms about the window's centre; Esc with no gesture is nothing
    const mid = (range(v)[0] + range(v)[1]) / 2;
    expect(m.wheel(10, 12, 0, -100, { mod: true })).toBe(true);
    frames!.flush();
    expect((range(v)[0] + range(v)[1]) / 2).toBeCloseTo(mid, 0);
    expect(m.cancel()).toBe(false);
    m.leave();
    expect(m.readout.value).toBeUndefined();
    // the sliders' values once the view rests
    v.draw();
    await drawn(v, m);
    expect(m.values.value.map((x) => x.section)).toEqual([P]);
    expect(m.issues).toEqual([]);
  });

  it('a declared navigator naming a picture the viewer lacks draws the axis, and the view says so', async () => {
    const base = DefaultSeqViewerSpec();
    const v = await make(spec('host:other'), { registry: base.registry });
    const m = v.navigatorModel('strip')!;
    await drawn(v, m);
    const pic = (m as unknown as { renderer: NavigatorRenderer }).renderer.pictureSection!;
    expect(pic.tracks.map((t) => t.layers[0].type.name)).toEqual(['ruler']);
    expect(v.view.resolved.value!.issues.map((i) => i.message)).toEqual([
      'No picture `host:other` here: drawn as the axis',
    ]);
  });

  it('a strip of a track of an extension’s type draws it as that type where the extension is composed', async () => {
    const types = async (v: SeqViewer) => {
      const m = v.navigatorModel('strip')!;
      await drawn(v, m);
      const pic = (m as unknown as { renderer: NavigatorRenderer }).renderer.pictureSection!;
      return pic.tracks.map((t) =>
        t.layers
          .filter((l) => l.layer.placed)
          .map((l) => l.type.name)
          .join(),
      );
    };
    expect(await types(await make(spec('tracks', true), { extensions: [Demo] }))).toEqual(['letters', 'demo:dots']);
    expect(await types(await make(spec('tracks', true)))).toEqual(['letters', '']);
  });
});

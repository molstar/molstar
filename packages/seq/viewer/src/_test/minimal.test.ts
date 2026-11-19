/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A viewer composed from a smaller registry (_spec/README.md "Composition", Injected, opt in): it draws what its
// registry lists and says what it does not (an unknown layer type, an unknown layout, navigators it cannot export, a
// story it cannot enter), the canvas and its SVG from the same registry (two consumers, one composition).

import type { Node } from '@molstar/core/repr/tree';
import { SeqRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec, stripNavigator } from '@molstar/seqx-core/lang';
import { type DrawLists, LettersType, TracksLayout } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { HoverProvider, NavigateProvider } from '../behaviors/providers.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import type { SeqViewerSpec } from '../viewer-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Letters = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDD';

/** A letters layer and a ranges layer on one sequence, a text section, a declared navigator and a story step. */
function spec(): Node {
  const b = createSeqSpec({ title: 'minimal' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: Letters });
  const a = b.data.ranges('a', {
    target: p,
    rows: { start: new Int32Array([2]), end: new Int32Array([9]), id: ['helix'], category: ['helix'] },
  });
  b.views.view('main', (v) => {
    const s = v.section('s', { of: p });
    s.track('letters').layer('l', { type: 'letters' });
    s.track('ranges').layer('r', { type: 'ranges', data: a });
    v.section('t', { of: p, layout: 'text' }).track('letters').layer('l', { type: 'letters' });
    v.navigator('nav', stripNavigator(s));
  });
  return b.build();
}

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
  vi.useRealTimers();
});

async function make(composition: SeqViewerSpec) {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  frames ??= installFrames();
  const v = createSeqViewer({ ...composition, ui: undefined });
  viewers.push(v);
  await v.load(spec());
  v.resize(600);
  await v.settle();
  frames.flush();
  return v;
}

const marks = (lists: DrawLists) => {
  const d = lists.passes.data;
  return d.rects.n + d.glyphs.n + d.paths.n + d.texts.n;
};
const codes = (issues: readonly { code: string }[]) => [...new Set(issues.map((i) => i.code))].sort();
const minimal = (registry: SeqViewerSpec['registry']): SeqViewerSpec => ({
  registry,
  behaviors: [{ provider: HoverProvider }, { provider: NavigateProvider }],
});

describe('a viewer composed from a smaller registry', () => {
  it('from an empty registry draws no layer and says so: each layer type and the layout it was not given', async () => {
    const v = await make(minimal([{ kind: 'format-runtime', runtime: SeqRuntime }]));
    const layout = v.frame!.layout;
    expect(layout.sections.flatMap((s) => s.tracks.flatMap((t) => t.layers))).toEqual([]);
    expect(codes(v.layoutIssues.value)).toEqual(['seq.render.layer.type', 'seq.render.layout']);
    // every layer said, by its address
    const said = v.layoutIssues.value.filter((i) => i.code === 'seq.render.layer.type').map((i) => i.address);
    expect(said).toEqual(expect.arrayContaining([expect.stringMatching(/\/l$/), expect.stringMatching(/\/r$/)]));
    // its layout table is what the registry gives: the text section is laid out as tracks
    expect(Object.keys(v.layouts)).toEqual([]);
    expect(v.context.types).toEqual({});
  });

  it('the canvas and its SVG take the same registry: letters drawn in both, the rest said in both', async () => {
    const v = await make(
      minimal([
        { kind: 'format-runtime', runtime: SeqRuntime },
        { kind: 'layer-type', type: LettersType },
        { kind: 'layout', layout: TracksLayout },
      ]),
    );
    const drawn = v.frame!.layout.sections.flatMap((s) => s.tracks.flatMap((t) => t.layers.map((l) => l.type.name)));
    expect([...new Set(drawn)]).toEqual(['letters']);
    expect(marks(v.frame!.lists)).toBeGreaterThan(0);
    const svg = await v.renderSvg({ navigators: 'include' });
    const svgDrawn = svg.lists;
    expect(marks(svgDrawn)).toBeGreaterThan(0);
    // the same issues: the ranges layer, the text layout; and the navigators the export cannot draw
    expect(codes(v.layoutIssues.value)).toEqual(['seq.render.layer.type', 'seq.render.layout']);
    expect(codes(svg.issues.filter((i) => i.code.startsWith('seq.render.')))).toEqual([
      'seq.render.layer.type',
      'seq.render.layout',
      'seq.render.navigator.export',
    ]);
    expect(svg.svg).not.toContain('data-navigator');
  });

  it('lists no navigator, story or history: no navigator model, no story step, jumps cut', async () => {
    const v = await make(minimal([{ kind: 'format-runtime', runtime: SeqRuntime }]));
    expect(v.navigatorModel('nav')).toBeUndefined();
    expect(v.navigator({ id: 'host', of: '/views/main/s' })).toBeUndefined();
    expect(v.behaviorIssues.value.map((i) => i.address)).toEqual(['navigator']);
    expect(v.spec.enterStep('any')?.code).toBe('seq.viewer.step.none');
    expect(v.spec.history).toBeUndefined();
    // a jump moves the model at once: no plan, no motion
    expect(v.view.jump(() => v.view.setRange('/views/main/s', 10, 20))).toBe('moved');
    expect(v.moving).toBe(false);
  });

  it('the default spec lists everything: the same view draws every layer, both layouts and the navigator', async () => {
    const v = await make(DefaultSeqViewerSpec());
    expect(codes(v.layoutIssues.value)).toEqual([]);
    const drawn = v.frame!.layout.sections.flatMap((s) => s.tracks.flatMap((t) => t.layers.map((l) => l.type.name)));
    expect([...new Set(drawn)].sort()).toEqual(['letters', 'ranges']);
    expect(v.navigatorModel('nav')).toBeDefined();
    const svg = await v.renderSvg({ navigators: 'include' });
    expect(svg.svg).toContain('data-navigator="nav"');
  });
});

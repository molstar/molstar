/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Stories in the viewer (_spec/README.md "Composition", Motion and Stories; "Ports", Animation), headless on the
// animation port's fake clock: a step played as a transition shown in the frame slot (no revision, no history); a full
// state's edits refused; every leading section's viewport in a jump's motion; reduced motion a cut; one plan sampled
// by the viewer and by a second consumer at the same time, equal.

import { diff, patchToNode } from '@molstar/core/repr/patch';
import { type Node, resolveAddress } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { sample, stripStory } from '@molstar/seqx-core/story';
import { afterEach, describe, expect, it } from 'vitest';
import { headlessAnimation } from '../animation.js';
import { DefaultViewerConfig, viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';

const Protein = '/views/main/protein';
const Other = '/views/main/other';
const Letters = ['views', 'main', 'protein', 'letters'];

function base(): Node {
  const b = createSeqSpec({ title: 'story' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDD' });
  const q = b.data.sequence('q', { alphabet: 'protein', letters: 'MKV'.repeat(20) });
  b.views.view('main', (v) => {
    v.section('protein', { of: p }).track('letters', { height: 20 }).layer('l', { type: 'letters' });
    v.section('other', { of: q }).track('letters').layer('l', { type: 'letters' });
  });
  return b.build();
}

/** The base with a story: `tall` a patch (taller letters, a window on each section), `frozen` a full state. */
function storySpec(): Node {
  const spec = base();
  const tall = (root: Node) => {
    const n = resolveAddress(root, Letters)!;
    return { ...n, params: { ...n.params, height: 40 } };
  };
  const tallRoot = replaceAt(spec, Letters, tall(spec));
  const patch = patchToNode({ ops: diff(spec, tallRoot).ops });
  const view = { sections: { [Protein]: { start: 10, end: 20 }, [Other]: { start: 30, end: 45 } } };
  return {
    ...spec,
    children: [
      ...(spec.children ?? []),
      {
        kind: 'snapshots',
        id: 'snapshots',
        children: [
          { kind: 'snapshot', id: 'tall', params: { label: 'Tall', patch, view, duration: 400, easing: 'linear' } },
          { kind: 'snapshot', id: 'frozen', params: { state: stripStory(tallRoot), duration: 400 } },
          {
            kind: 'snapshot',
            id: 'scrolled',
            params: {
              patch: patchToNode({ ops: [] }),
              view: { sections: { [Protein]: { start: 0, end: 49, top: 30 } } },
              duration: 400,
              easing: 'linear',
            },
          },
        ],
      },
    ],
  };
}

function replaceAt(root: Node, at: readonly string[], node: Node): Node {
  if (!at.length) return node;
  return {
    ...root,
    children: root.children!.map((c) => (c.id === at[0] ? replaceAt(c, at.slice(1), node) : c)),
  };
}

let viewers: SeqViewer[] = [];
afterEach(() => {
  for (const v of viewers) v.retire();
  viewers = [];
});

async function make(config: [unknown, unknown][] = []) {
  const { port, clock } = headlessAnimation();
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    animation: port,
    config: [[viewerItem('viewport.maxPxPerPosition'), 1000], ...(config as never[])],
  });
  viewers.push(v);
  await v.load(storySpec());
  v.resize(600);
  await v.settle();
  v.draw();
  return { v, clock };
}

const heightShown = (v: SeqViewer) => resolveAddress(v.spec.shown.value, Letters)?.params?.height;

describe('a story step played on the application’s clock', () => {
  it('shows the transition in the frame slot, timed by the step; the step’s state and view at its end', async () => {
    const { v, clock } = await make();
    const revision = v.spec.revision.value;
    const { plan, issue } = v.playStep('tall');
    expect(issue).toBeUndefined();
    expect(plan!.duration).toBe(400);
    expect(plan!.channels.map((c) => c.kind).sort()).toEqual(['number', 'view', 'view']);
    expect(v.moving).toBe(true);
    // the model is at the step at once; the frame shows the way there
    expect(v.view.get(Protein)).toMatchObject({ start: 10, end: 20 });
    expect(heightShown(v)).toBe(20);
    const drawn0 = v.animation.drawn.value;
    clock.tick(200);
    expect(heightShown(v)).toBe(30);
    expect(v.animation.drawn.value).toBeGreaterThan(drawn0);
    const mid = v.drawnViewports;
    // every leading section moves, each on its own
    expect(mid[Protein].start).toBeGreaterThan(0);
    expect(mid[Protein].start).toBeLessThan(10);
    expect(mid[Other].start).toBeGreaterThan(0);
    expect(mid[Other].start).toBeLessThan(30);
    clock.runAll(16);
    expect(v.moving).toBe(false);
    expect(v.spec.frame).toBeUndefined();
    expect(heightShown(v)).toBe(40);
    expect(v.drawnViewports[Protein]).toMatchObject({ start: 10, end: 20 });
    expect(v.drawnViewports[Other]).toMatchObject({ start: 30, end: 45 });
    // navigation: no revision, no history step
    expect(v.spec.revision.value).toBe(revision);
    expect(v.view.history.value).toEqual({ back: 0, forward: 0 });
    // leaving it is a transition back, timed by the configuration
    expect(v.playStep(undefined).plan!.duration).toBe(240);
    clock.tick(120);
    expect(heightShown(v)).toBe(30);
    clock.runAll(16);
    expect(heightShown(v)).toBe(20);
    expect(v.playStep('nothing').issue?.code).toBe('seq.viewer.step.unknown');
  });

  it('a step entered by its command is the same transition, its view placed', async () => {
    const { v, clock } = await make();
    expect(v.runCommand('step.enter', { step: 'tall' })).toBe('Shown');
    expect(v.moving).toBe(true);
    expect(v.view.get(Protein)).toMatchObject({ start: 10, end: 20 });
    clock.tick(200);
    expect(heightShown(v)).toBe(30);
    clock.runAll(16);
    expect(heightShown(v)).toBe(40);
    expect(v.drawnViewports[Other]).toMatchObject({ start: 30, end: 45 });
    expect(v.runCommand('step.leave', {})).toBeTruthy();
    expect(v.moving).toBe(true);
    clock.runAll(16);
    expect(heightShown(v)).toBe(20);
  });

  it('is a cut under reduced motion: nothing is requested of the clock', async () => {
    const { v, clock } = await make([[viewerItem('view.motion'), 'cut']]);
    v.playStep('tall');
    expect(v.moving).toBe(false);
    expect(clock.running).toBe(0);
    expect(heightShown(v)).toBe(40);
  });

  it('the application drives a plan it scrubs: `progress`, then `clear`', async () => {
    const { v, clock } = await make();
    const { plan } = v.playStep('tall', { run: false });
    expect(clock.running).toBe(0);
    expect(heightShown(v)).toBe(20);
    v.animation.progress(plan!, 0.25);
    expect(heightShown(v)).toBe(25);
    v.animation.progress(plan!, 1);
    expect(heightShown(v)).toBe(40);
    v.animation.clear();
    expect(v.spec.frame).toBeUndefined();
    expect(heightShown(v)).toBe(40);
  });

  it('one plan sampled by the viewer’s port and by a second consumer at the same time: the same picture', async () => {
    const { v, clock } = await make();
    const { v: w } = await make();
    const { plan } = v.playStep('tall');
    clock.tick(100);
    // the second consumer (an export, a structure's camera) samples the same plan at the clock's time
    w.spec.enterStep('tall');
    w.animation.show(plan!, sample(plan!, clock.now / plan!.duration, undefined));
    w.draw();
    const [a, b] = await Promise.all([v.renderSvg({ dpr: 1 }), w.renderSvg({ dpr: 1 })]);
    expect(a.svg).toBe(b.svg);
    expect(heightShown(w)).toBe(heightShown(v));
  });

  it('a full-state step refuses edits, with a hint to edit the step', async () => {
    const { v, clock } = await make();
    v.playStep('frozen');
    clock.runAll(16);
    expect(v.spec.step.value).toMatchObject({ id: 'frozen', form: 'full' });
    const r = v.spec.applyPatch({ ops: [{ op: 'set', at: Letters, key: 'height', value: 5 }] });
    expect(!r.ok && r.issues.map((i) => i.code)).toEqual(['seq.viewer.step.frozen']);
    v.playStep(undefined);
    clock.runAll(16);
    expect(v.spec.applyPatch({ ops: [{ op: 'set', at: Letters, key: 'height', value: 5 }] }).ok).toBe(true);
  });
});

describe('a jump’s motion through the port', () => {
  it('moves every leading section, each its own way, then stops where the model is', async () => {
    const { v, clock } = await make();
    v.view.jump(() => {
      v.view.setRange(Protein, 20, 30);
      v.view.setRange(Other, 40, 50);
      return true;
    });
    expect(clock.running).toBe(1);
    clock.tick(v.config.view.jumpMs / 2);
    const mid = v.drawnViewports;
    expect(mid[Protein].start).toBeGreaterThan(0);
    expect(mid[Protein].start).toBeLessThan(20);
    expect(mid[Other].start).toBeGreaterThan(0);
    expect(mid[Other].start).toBeLessThan(40);
    // the export and `viewport` read the destination: the motion is the eye's way only
    expect(v.viewport).toEqual({ start: 20, end: 30 });
    clock.runAll(16);
    expect(v.moving).toBe(false);
    expect(v.drawnViewports[Protein]).toMatchObject({ start: 20, end: 30 });
    expect(v.drawnViewports[Other]).toMatchObject({ start: 40, end: 50 });
  });

  it('input stops it where it is drawn, every leading section there', async () => {
    const { v, clock } = await make();
    v.view.jump(() => {
      v.view.setRange(Protein, 20, 30);
      v.view.setRange(Other, 40, 50);
      return true;
    });
    clock.tick(v.config.view.jumpMs / 2);
    const mid = v.drawnViewports;
    v.stopMotion();
    expect(clock.running).toBe(0);
    expect(v.view.get(Protein)).toMatchObject({ start: mid[Protein].start, end: mid[Protein].end });
    expect(v.view.get(Other)).toMatchObject({ start: mid[Other].start, end: mid[Other].end });
  });
});

describe('a ring turned by a step', () => {
  it('tweens its rotation the short way round, its labels hidden while it turns and placed once at its end', async () => {
    const { port, clock } = headlessAnimation();
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, animation: port });
    viewers.push(v);
    const b = createSeqSpec({ title: 'ring' });
    const c = b.data.sequence('c', { alphabet: 'dna', letters: 'ACGT'.repeat(250), cyclic: true });
    const f = b.data.ranges('f', {
      target: c,
      rows: { start: new Int32Array([10, 400]), end: new Int32Array([200, 600]), id: ['a', 'b'] },
    });
    const Ring = ['views', 'main', 'ring'];
    b.views.view('main', (w) => {
      w.section('ring', { of: c, layout: 'rings', rotation: 20 }).track('f').layer('f', { type: 'ranges', data: f });
    });
    b.snapshots.snapshot('turned', {
      patch: { ops: [{ op: 'set', at: Ring, key: 'rotation', value: 300 }] },
      duration: 400,
      easing: 'linear',
    });
    await v.load(b.build());
    v.resize(600);
    await v.settle();
    v.draw();
    const { plan } = v.playStep('turned');
    expect(plan!.channels[0]).toMatchObject({ kind: 'number', delta: -80 });
    const layouts = v.counts.layouts;
    clock.tick(200);
    const mid = resolveAddress(v.spec.shown.value, Ring)!.params!;
    expect(mid).toMatchObject({ rotation: -20, labels: 'none' });
    clock.runAll(16);
    expect(resolveAddress(v.spec.shown.value, Ring)!.params).toMatchObject({ rotation: 300 });
    expect(resolveAddress(v.spec.shown.value, Ring)!.params!.labels).toBeUndefined();
    // laid out each frame of the turn, and once more at its end, where its labels are placed
    expect(v.counts.layouts).toBeGreaterThan(layouts);
  });
});

describe('a jump in a text: its anchor and its letter size', () => {
  it('a step zoom tweens the letter size geometrically, the anchor along', async () => {
    const { port, clock } = headlessAnimation();
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, animation: port });
    viewers.push(v);
    const b = createSeqSpec({ title: 'text' });
    const p = b.data.sequence('p', {
      alphabet: 'protein',
      letters: 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDD'.repeat(20),
    });
    b.views.view('main', (w) => {
      w.section('t', { of: p, layout: 'text' }).track('letters').layer('l', { type: 'letters' });
    });
    await v.load(b.build());
    v.resize(600);
    await v.settle();
    v.draw();
    const T = '/views/main/t';
    expect(v.view.zoom(T, 1.5, 0)).toBe(true);
    v.draw();
    const s0 = v.view.get(T)!.scale!;
    expect(v.view.jump(() => v.view.zoom(T, 2, 0), false)).toBe('moved');
    const s1 = v.view.get(T)!.scale!;
    expect(s1).toBeGreaterThan(s0);
    clock.tick(v.config.view.jumpMs / 2);
    // geometric: at the middle of the eased time, the sizes' geometric mean
    expect(v.drawnViewports[T].scale!).toBeCloseTo(Math.sqrt(s0 * s1), 6);
    clock.runAll(16);
    v.draw();
    expect(v.drawnViewports[T].scale).toBe(v.view.get(T)!.scale);
  });
});

// The review's fixes: each test failed on the fix's parent.

/** What is hit at the middle of the protein's letters, against what the window drawn says is there. */
function hitAgrees(v: SeqViewer) {
  const w = v.drawnViewports[Protein];
  const box = v.frame!.layout.sections[0].tracks[0].box;
  const local = box.width * 0.5;
  const hit = v.hit(box.x + local, box.y + 2);
  return { hit: hit?.position, want: Math.floor(w.start + ((w.end - w.start) * local) / box.width) };
}

describe('input during a story’s motion', () => {
  it('a pan by hand ends a step’s transition: drawn at once where the model is, and hits match the picture', async () => {
    const { v, clock } = await make();
    v.playStep('tall');
    clock.tick(200);
    expect(v.moving).toBe(true);
    v.view.panPx(Protein, -60);
    expect(v.moving).toBe(false);
    expect(v.spec.frame).toBeUndefined();
    v.draw();
    const m = v.view.get(Protein)!;
    expect(v.drawnViewports[Protein]).toMatchObject({ start: m.start, end: m.end });
    const { hit, want } = hitAgrees(v);
    expect(hit).toBe(want);
  });

  it('a pan by hand under a host’s scrub clears its frame', async () => {
    const { v } = await make();
    const { plan } = v.playStep('tall', { run: false });
    v.animation.progress(plan!, 0.5);
    v.draw();
    expect(heightShown(v)).toBe(30);
    v.view.panPx(Protein, -60);
    expect(v.spec.frame).toBeUndefined();
    v.draw();
    expect(heightShown(v)).toBe(40);
    const m = v.view.get(Protein)!;
    expect(v.drawnViewports[Protein]).toMatchObject({ start: m.start, end: m.end });
  });

  it('input stops a step’s motion where it is drawn (stopMotion), every leading section there', async () => {
    const { v, clock } = await make();
    v.playStep('tall');
    clock.tick(200);
    const mid = v.drawnViewports;
    v.stopMotion();
    expect(v.moving).toBe(false);
    expect(clock.running).toBe(0);
    expect(v.view.get(Protein)).toMatchObject({ start: mid[Protein].start, end: mid[Protein].end });
    expect(v.view.get(Other)).toMatchObject({ start: mid[Other].start, end: mid[Other].end });
    const { hit, want } = hitAgrees(v);
    expect(hit).toBe(want);
  });
});

describe('the spec changing under a frame', () => {
  it('an edit during a step’s transition ends it: the edit is shown', async () => {
    const { v, clock } = await make();
    v.playStep('tall');
    clock.tick(100);
    expect(v.spec.applyPatch({ ops: [{ op: 'set', at: [], key: 'title', value: 'edited' }] }).ok).toBe(true);
    expect(v.moving).toBe(false);
    expect(v.spec.frame).toBeUndefined();
    expect(v.spec.shown.value.params?.title).toBe('edited');
    expect(heightShown(v)).toBe(40);
  });

  it('an edit under a host’s scrub drops its frame', async () => {
    const { v } = await make();
    const { plan } = v.playStep('tall', { run: false });
    v.animation.progress(plan!, 0.5);
    expect(v.spec.applyPatch({ ops: [{ op: 'set', at: [], key: 'title', value: 'edited' }] }).ok).toBe(true);
    expect(v.spec.frame).toBeUndefined();
    expect(v.spec.shown.value.params?.title).toBe('edited');
  });
});

describe('a step’s top', () => {
  it('tweens in a step’s transition, and in a jump', async () => {
    const { v, clock } = await make();
    const { plan } = v.playStep('scrolled');
    expect(plan!.channels.some((c) => c.kind === 'view' && c.fields.includes('top'))).toBe(true);
    clock.tick(200);
    expect(v.drawnViewports[Protein].top).toBeCloseTo(15, 6);
    clock.runAll(16);
    v.draw();
    expect(v.drawnViewports[Protein].top).toBe(30);
    v.view.jump(() => {
      v.view.setRange(Protein, 10, 20);
      v.view.setTop(Protein, 0);
      return true;
    });
    clock.tick(v.config.view.jumpMs / 2);
    const top = v.drawnViewports[Protein].top!;
    expect(top).toBeGreaterThan(0);
    expect(top).toBeLessThan(30);
  });
});

describe('reduced motion', () => {
  it('a cut plan’s natural duration is 0', async () => {
    const { v } = await make([[viewerItem('view.motion'), 'cut']]);
    const { plan } = v.playStep('tall');
    expect(v.animation.durationOf(plan!)).toBe(0);
  });

  it('`view.motion` decides for stories as for jumps: `platform` cuts both under reduced motion, `animate` neither', async () => {
    const { v } = await make();
    v.reducedMotion = true;
    v.view.jump(() => v.view.setRange(Protein, 20, 30));
    expect(v.moving).toBe(false);
    v.playStep('tall');
    expect(v.moving).toBe(false);
    const w = await make([[viewerItem('view.motion'), 'animate']]);
    w.v.reducedMotion = true;
    w.v.view.jump(() => w.v.view.setRange(Protein, 20, 30));
    expect(w.v.moving).toBe(true);
    w.clock.runAll(16);
    w.v.playStep('tall');
    expect(w.v.moving).toBe(true);
    w.clock.runAll(16);
  });

  it('a legacy `view.motion: auto` reads as `platform`: a cut under reduced motion', async () => {
    const { v } = await make([[viewerItem('view.motion'), 'auto']]);
    expect(v.config.view.motion).toBe('platform');
    v.reducedMotion = true;
    v.playStep('tall');
    expect(v.moving).toBe(false);
  });

  it('`story.reducedMotion` is no item: one item, `view.motion`', () => {
    expect('reducedMotion' in DefaultViewerConfig.story).toBe(false);
    expect(DefaultViewerConfig.view.motion).toBe('platform');
  });
});

describe('a full-state step refuses every edit', () => {
  it('a whole root and a preview too; a new document leaves the step', async () => {
    const { v, clock } = await make();
    v.playStep('frozen');
    clock.runAll(16);
    const before = v.spec.root.value;
    expect(v.spec.apply({ ...before, params: { ...before.params, title: 'x' } })).toBe('frozen-step');
    expect(v.spec.root.value).toBe(before);
    const r = v.spec.setPreview({ ops: [{ op: 'set', at: Letters, key: 'height', value: 5 }] });
    expect(r.ok).toBe(false);
    expect(v.spec.preview.value).toBeUndefined();
    await v.load(base());
    expect(v.spec.step.value).toBeUndefined();
  });
});

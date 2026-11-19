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
import { viewerItem } from '../config.js';
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

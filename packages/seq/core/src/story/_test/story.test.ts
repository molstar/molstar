/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Stories (_spec/README.md "Stories"): a step's state in both forms and as a delta; a plan's channels by kind; sampling
// as properties (the ends exact, monotone, deterministic, reduced motion a cut at the end); the frame's root; the
// encoder's forms within the chain bound.

import { PCG } from '@molstar/core/data/util/hash-functions';
import { diff, patchToNode } from '@molstar/core/repr/patch';
import { hashNode, type Node, resolveAddress } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { checkSeqSpec, createSeqSpec } from '../../lang/index.js';
import {
  ChannelPhase,
  createStoryCache,
  frameDone,
  frameRoot,
  frameView,
  plan,
  sample,
  stateOf,
  storeSteps,
  storySteps,
  stripStory,
  withSteps,
} from '../index.js';
import { Domains, edit, Letters, Ring, storySpec } from './fixture.js';

const heightOf = (root: Node, at = Letters) => resolveAddress(root, at)?.params?.height;
const patchNode = (a: Node, b: Node) => patchToNode({ ops: diff(a, b).ops });

function withStory(spec: Node, steps: readonly Node[]) {
  return withSteps(spec, steps);
}

describe('stateOf', () => {
  const spec = storySpec();
  const taller = edit(spec, { op: 'set', at: Letters, key: 'height', value: 40 });
  const tallest = edit(taller, { op: 'set', at: Domains, key: 'height', value: 30 });

  it('a patch over the spec, a full state, a delta from a step; cached by the spec’s revision and the step', () => {
    const s = withStory(spec, [
      { kind: 'snapshot', id: 'a', params: { label: 'Taller', patch: patchNode(spec, taller), duration: 500 } },
      { kind: 'snapshot', id: 'b', params: { from: '/snapshots/a', patch: patchNode(taller, tallest) } },
      { kind: 'snapshot', id: 'c', params: { state: stripStory(tallest), easing: 'linear', linger: 10 } },
    ]);
    expect(checkSeqSpec(s).filter((i) => i.severity === 'error')).toEqual([]);
    expect(storySteps(s)).toEqual(['a', 'b', 'c']);
    const cache = createStoryCache();
    const a = stateOf(s, 'a', { cache })!;
    expect(a).toMatchObject({ form: 'patch', label: 'Taller', applies: true, chain: 0 });
    expect(a.timing).toEqual({ duration: 500, easing: 'ease-in-out', linger: 2000 });
    expect(heightOf(a.root)).toBe(40);
    const b = stateOf(s, 'b', { cache })!;
    expect(b).toMatchObject({ form: 'patch', from: 'a', chain: 1, applies: true });
    expect(hashNode(stripStory(b.root))).toBe(hashNode(stripStory(tallest)));
    const c = stateOf(s, 'c', { cache })!;
    expect(c).toMatchObject({ form: 'full', chain: 0, timing: { easing: 'linear', linger: 10 } });
    expect(c.root).toBe(resolveAddress(s, ['snapshots', 'c'])!.params!.state);
    expect(stateOf(s, 'b', { cache })).toBe(b);
    expect(stateOf(s, 'nothing')).toBeUndefined();
  });

  it('a start that is missing or leads back to itself is an issue on the step', () => {
    const s = withStory(spec, [
      { kind: 'snapshot', id: 'a', params: { from: '/snapshots/b', patch: patchNode(spec, taller) } },
      { kind: 'snapshot', id: 'b', params: { from: '/snapshots/a', patch: patchToNode({ ops: [] }) } },
      { kind: 'snapshot', id: 'c', params: { from: '/snapshots/gone', patch: patchNode(spec, taller) } },
    ]);
    for (const id of ['a', 'b', 'c']) {
      const st = stateOf(s, id)!;
      expect(st.applies).toBe(false);
      expect(st.issues.map((i) => i.code)).toContain('seq.story.step.from');
    }
    // the missing start's step still shows its own ops, over the spec
    expect(heightOf(stateOf(s, 'c')!.root)).toBe(40);
  });

  it('a step that no longer applies whole shows the ops that still apply, flagged', () => {
    const note: Node = {
      kind: 'track',
      id: 'note',
      children: [{ kind: 'layer', id: 'x', params: { type: 'letters' } }],
    };
    const step = edit(taller, {
      op: 'insert',
      parent: ['views', 'main', 'protein'],
      place: { last: true },
      node: note,
    });
    // the spec gains a track of the id the step inserts: its insert no longer applies, its height still does
    const s0 = edit(spec, { op: 'insert', parent: ['views', 'main', 'protein'], place: { last: true }, node: note });
    const s = withStory(s0, [{ kind: 'snapshot', id: 'a', params: { patch: patchNode(spec, step) } }]);
    const st = stateOf(s, 'a')!;
    expect(st.applies).toBe(false);
    expect(st.skipped).toBe(1);
    expect(st.issues.map((i) => i.code)).toEqual(['seq.story.step.stale']);
    expect(heightOf(st.root)).toBe(40);
  });
});

describe('plan', () => {
  const spec = storySpec();
  const kinds = (a: Node, b: Node) => plan({ root: a }, { root: b }).channels.map((c) => c.kind);

  it('a channel per change by kind: number, color, cut, enter, leave, order', () => {
    expect(kinds(spec, edit(spec, { op: 'set', at: Letters, key: 'height', value: 40 }))).toEqual(['number']);
    expect(kinds(spec, edit(spec, { op: 'set', at: [...Letters, 'l'], key: 'color', value: '#ff0000' }))).toEqual([
      'color',
    ]);
    // a token is no color a plan can blend: a cut
    expect(kinds(spec, edit(spec, { op: 'set', at: [...Letters, 'l'], key: 'color', value: 'token:cat-2' }))).toEqual([
      'cut',
    ]);
    expect(kinds(spec, edit(spec, { op: 'set', at: [...Letters, 'l'], key: 'type', value: 'values' }))).toEqual([
      'cut',
    ]);
    const track: Node = { kind: 'track', id: 't', params: { height: 10 }, children: [] };
    expect(
      kinds(spec, edit(spec, { op: 'insert', parent: Letters.slice(0, -1), place: { last: true }, node: track })),
    ).toEqual(['enter']);
    expect(kinds(spec, edit(spec, { op: 'remove', at: Domains }))).toEqual(['leave']);
    expect(
      kinds(spec, edit(spec, { op: 'move', at: Domains, parent: Domains.slice(0, -1), place: { after: null } })),
    ).toEqual(['order']);
    expect(kinds(spec, spec)).toEqual([]);
  });

  it('a ring’s rotation turns the short way round, its labels held hidden while it turns', () => {
    const to = edit(spec, { op: 'set', at: Ring, key: 'rotation', value: 350 });
    const p = plan({ root: spec }, { root: to });
    const c = p.channels[0];
    expect(c).toMatchObject({ kind: 'number', from: 0, to: 350, delta: -10, hold: { key: 'labels', value: 'none' } });
    const mid = frameRoot(p, sample(p, 0.5));
    expect(resolveAddress(mid, Ring)?.params).toMatchObject({ rotation: -5, labels: 'none' });
    // placed once, at its end: the end state, labels as it has them
    expect(frameRoot(p, sample(p, 1))).toBe(to);
  });

  it('a view channel pans when the windows overlap, else goes out, across and in; a cyclic one the short arc', () => {
    const view = (start: number, end: number) => ({ sections: { s: { start, end } } });
    const pan = plan({ root: spec, view: view(0, 20) }, { root: spec, view: view(10, 30) }).channels[0];
    expect(pan).toMatchObject({ kind: 'view', path: 'pan' });
    const far = plan({ root: spec, view: view(0, 10) }, { root: spec, view: view(80, 90) });
    expect(far.channels[0]).toMatchObject({ kind: 'view', path: 'out-across-in' });
    const mid = sample(far, 0.5).values;
    expect(mid[1] - mid[0]).toBeCloseTo(90);
    const arc = plan({ root: spec, view: view(0, 10) }, { root: spec, view: view(90, 100) }, { cyclic: { s: 100 } });
    expect(arc.channels[0]).toMatchObject({ path: 'pan', shift: -100 });
    const m = sample(arc, 0.5).values;
    expect((m[0] + m[1]) / 2).toBeCloseTo(0);
    expect(Array.from(sample(arc, 1).values)).toEqual([90, 100]);
    // the view shown switches at the cut
    const switched = plan(
      { root: spec, view: { view: '/views/a', sections: {} } },
      { root: spec, view: { view: '/views/b', sections: {} } },
    );
    expect(frameView(switched, sample(switched, 0.4))?.view).toBe('/views/a');
    expect(frameView(switched, sample(switched, 0.6))?.view).toBe('/views/b');
  });

  it('a text’s anchor and size, a band’s top row and row size, a grid’s y: field by field', () => {
    const a = {
      start: 10,
      end: 60,
      top: 0,
      scale: 12,
      bands: { t: { bandTop: 0, rowScale: 10 } },
      y: { start: 0, end: 5 },
    };
    const b = {
      start: 40,
      end: 90,
      top: 20,
      scale: 24,
      bands: { t: { bandTop: 8, rowScale: 20 } },
      y: { start: 2, end: 9 },
    };
    const p = plan({ root: spec, view: { sections: { s: a } } }, { root: spec, view: { sections: { s: b } } });
    const c = p.channels[0];
    expect(c.kind === 'view' && c.fields).toEqual([
      'start',
      'end',
      'top',
      'scale',
      'bands.t.bandTop',
      'bands.t.rowScale',
      'y.start',
      'y.end',
    ]);
    expect(frameView(p, sample(p, 0))?.sections.s).toEqual(a);
    expect(frameView(p, sample(p, 1))?.sections.s).toEqual(b);
    const mid = frameView(p, sample(p, 0.5))!.sections.s;
    // sizes move geometrically: the middle of 12 and 24 is their geometric mean
    expect(mid.scale).toBeCloseTo(Math.sqrt(12 * 24));
    expect(mid.bands!.t.rowScale).toBeCloseTo(Math.sqrt(200));
    expect(mid.top).toBeCloseTo(10);
  });
});

describe('sample: properties', () => {
  const spec = storySpec();
  /** A random pair of states: heights, a color, a rotation, windows. */
  function pair(rng: PCG) {
    const h = () => 4 + Math.floor(rng.float() * 60);
    const hex = () =>
      `#${Math.floor(rng.float() * 0xffffff)
        .toString(16)
        .padStart(6, '0')}`;
    const state = () => ({
      root: edit(
        spec,
        { op: 'set', at: Letters, key: 'height', value: h() },
        { op: 'set', at: Domains, key: 'height', value: h() },
        { op: 'set', at: [...Letters, 'l'], key: 'color', value: hex() },
        { op: 'set', at: Ring, key: 'rotation', value: Math.floor(rng.float() * 360) },
      ),
      view: (() => {
        const s = rng.float() * 40;
        return { sections: { s: { start: s, end: s + 1 + rng.float() * 20, scale: 1 + rng.float() * 20 } } };
      })(),
    });
    return [state(), state()] as const;
  }

  it('sample(plan, 0) is the start and sample(plan, 1) the end, exactly; deterministic', () => {
    const rng = new PCG(7);
    for (let i = 0; i < 50; i++) {
      const [a, b] = pair(rng);
      const p = plan(a, b);
      expect(frameRoot(p, sample(p, 0))).toBe(a.root);
      expect(frameView(p, sample(p, 0))).toEqual(a.view);
      expect(frameDone(p, sample(p, 1))).toBe(true);
      expect(frameRoot(p, sample(p, 1))).toBe(b.root);
      expect(frameView(p, sample(p, 1))).toEqual(b.view);
      const t = rng.float();
      const x = Array.from(sample(p, t).values);
      expect(Array.from(sample(plan(a, b), t).values)).toEqual(x);
    }
  });

  it('every channel is monotone: numbers, each color component, a pan’s centre and width', () => {
    const rng = new PCG(11);
    for (let i = 0; i < 30; i++) {
      const [a, b] = pair(rng);
      const p = plan(a, b, { config: { ...planConfig, view: { path: 'pan' } } });
      // a window's ends are read as its centre and width (each end alone may overshoot while it zooms)
      const read = (t: number) => {
        const v = Array.from(sample(p, t).values);
        for (const c of p.channels)
          if (c.kind === 'view') {
            const [s, e] = [v[c.slot], v[c.slot + 1]];
            v[c.slot] = (s + e) / 2;
            v[c.slot + 1] = e - s;
          }
        return v;
      };
      let prev = read(0);
      const dir = read(1).map((v, k) => Math.sign(v - prev[k]));
      for (let k = 1; k <= 20; k++) {
        const now = read(k / 20);
        for (let j = 0; j < now.length; j++) expect((now[j] - prev[j]) * dir[j]).toBeGreaterThanOrEqual(-1e-9);
        prev = now;
      }
    }
  });

  it('reduced motion: every channel a cut at the end', () => {
    const rng = new PCG(3);
    const [a, b] = pair(rng);
    const p = plan(a, b, { reducedMotion: true });
    for (const t of [0, 0.3, 0.99]) {
      const f = sample(p, t);
      expect(Array.from(f.phase).every((x) => x === ChannelPhase.before)).toBe(true);
      expect(frameRoot(p, f)).toBe(a.root);
    }
    expect(frameRoot(p, sample(p, 1))).toBe(b.root);
  });

  it('a channel waiting on data holds its start and says so', () => {
    const view = (start: number, end: number) => ({ sections: { s: { start, end } } });
    const p = plan({ root: spec, view: view(0, 10) }, { root: spec, view: view(5, 15) }, { pending: () => true });
    const f = sample(p, 0.5);
    expect(f.waiting).toBe(1);
    expect(Array.from(f.values)).toEqual([0, 10]);
  });

  it('a track entering grows from no height; one going shrinks, then goes', () => {
    const track: Node = { kind: 'track', id: 't', params: { height: 10 }, children: [] };
    const to = edit(spec, { op: 'insert', parent: Letters.slice(0, -1), place: { last: true }, node: track });
    const p = plan({ root: spec }, { root: to, view: undefined }, { easing: 'linear' });
    expect(heightOf(frameRoot(p, sample(p, 0.5)), [...Letters.slice(0, -1), 't'])).toBe(5);
    const back = plan({ root: to }, { root: spec }, { easing: 'linear' });
    expect(heightOf(frameRoot(back, sample(back, 0.25)), [...Letters.slice(0, -1), 't'])).toBe(7.5);
    expect(resolveAddress(frameRoot(back, sample(back, 1)), [...Letters.slice(0, -1), 't'])).toBeUndefined();
  });
});

const planConfig = {
  transition: { duration: 240, easing: 'ease-in-out' as const },
  linger: 2000,
  view: { path: 'auto' as const },
  color: { space: 'oklab' as const },
  cut: { at: 0.5 },
  reducedMotion: 'platform' as const,
  stepStart: 'spec' as const,
  store: 'patch' as const,
  maxChain: 8,
};

describe('storeSteps: the encoder', () => {
  const spec = storySpec();
  /** Twelve states, each a small change of the one before. */
  const states = (() => {
    const out: Node[] = [];
    let r = spec;
    for (let i = 0; i < 12; i++) {
      r = edit(r, { op: 'set', at: i % 2 ? Domains : Letters, key: 'height', value: 10 + i });
      out.push(r);
    }
    return out;
  })();
  const inputs = states.map((root, i) => ({ id: `s${i}`, root }));
  const roundTrip = (steps: Node[]) => {
    const s = withSteps(spec, steps);
    expect(checkSeqSpec(s).filter((i) => i.severity === 'error')).toEqual([]);
    states.forEach((root, i) => {
      const st = stateOf(s, `s${i}`)!;
      expect(st.applies).toBe(true);
      expect(hashNode(stripStory(st.root))).toBe(hashNode(stripStory(root)));
    });
    return steps.map((n) => stateOf(s, n.id!)!);
  };

  it('`patch`: a patch over the spec each; `full`: full states', () => {
    expect(roundTrip(storeSteps(spec, inputs)).every((s) => s.form === 'patch' && s.chain === 0)).toBe(true);
    expect(roundTrip(storeSteps(spec, inputs, { store: 'full' })).every((s) => s.form === 'full')).toBe(true);
  });

  it('`smaller`: deltas within `story.maxChain` of a keyframe; `keep` wins', () => {
    const kept = inputs.map((x, i) => (i === 3 ? { ...x, keep: 'full' as const } : x));
    const steps = roundTrip(storeSteps(spec, kept, { store: 'smaller', maxChain: 2 }));
    expect(steps[3].form).toBe('full');
    expect(Math.max(...steps.map((s) => s.chain))).toBeLessThanOrEqual(2);
    expect(steps.some((s) => s.chain > 0)).toBe(true);
  });

  it('the builder writes each form', () => {
    const b = createSeqSpec();
    b.snapshots.snapshot('p', { patch: { ops: [] }, duration: 100 });
    b.snapshots.snapshot('f', { state: stripStory(spec), from: '/snapshots/p' });
    const built = b.build();
    expect(resolveAddress(built, ['snapshots', 'p'])?.params?.duration).toBe(100);
    expect(stateOf(built, 'f')?.form).toBe('full');
  });
});

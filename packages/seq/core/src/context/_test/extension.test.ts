/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Extensions (_spec/README.md "Context", Extensions; "Resolved view", Layer kinds): a namespaced layer type drawn by a
// composition that has it, its fallback by one without it (one draw list, so the canvas and the SVG alike), the issue
// once a namespace, the shape it reads checked, a newer bare type drawn as its fallback, registration atomic over two
// extensions, the version refused; and a spec with an extension reads, prints, patches and diffs the same with and
// without it.

import { applyPatch, diff } from '@molstar/core/repr/patch';
import { parseSexpr, printSexpr } from '@molstar/core/repr/codecs/sexpr';
import type { Issue, Node } from '@molstar/core/repr/tree';
import { composeSeqContext } from '@molstar/seqx-core/context';
import { createSeqContext } from '@molstar/seqx-core/context/default';
import { createSeqSpec, SeqLanguage, SeqValidateOptions } from '@molstar/seqx-core/lang';
import { RangesType } from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import {
  defineExtension,
  type Hit,
  type LayerType,
  requiresCovers,
  type SeqExtension,
} from '@molstar/seqx-core/extension';

/** A dot at each range's middle (inline ranges only): the smallest layer type an extension can ship. */
interface DotsPrep {
  readonly items: readonly { readonly start: number; readonly end: number }[];
  readonly radius: number;
}
interface DotsLay {
  xs: number[];
  y: number;
  r: number;
  color: number;
  key: string;
  section: string;
  system: string;
}

const DotsType: LayerType<DotsPrep, DotsLay> = {
  name: 'demo:dots',
  shape: 'ranges',
  fallback: 'ranges',
  prepare(layer, model) {
    const a = layer.data ? model.annotations.get(layer.data.address) : undefined;
    const radius = layer.params.radius;
    return { items: a?.items ?? [], radius: typeof radius === 'number' ? radius : 3 };
  },
  layout(prep, ctx, _data, _overlay, prev) {
    const lay = prev ?? { xs: [], y: 0, r: 0, color: 0, key: '', section: '', system: '' };
    lay.xs.length = 0;
    for (const it of prep.items) {
      const x = (ctx.axis.x(it.start) + ctx.axis.x(it.end)) / 2;
      if (Number.isFinite(x)) lay.xs.push(x);
    }
    lay.y = ctx.box.y + ctx.box.height / 2;
    lay.r = prep.radius;
    lay.color = ctx.tokens.colors['seq-range'];
    lay.key = ctx.key;
    lay.section = ctx.section.address;
    lay.system = ctx.section.axis.system;
    return lay;
  },
  draw(lay, pass, out) {
    if (pass !== 'data') return;
    for (const x of lay.xs) out.rect(x - lay.r, lay.y - lay.r, 2 * lay.r, 2 * lay.r, lay.color, 0, 0, lay.r);
  },
  hit(lay, at): Hit | undefined {
    const near = lay.xs.findIndex((x) => Math.abs(x - at[0]) <= lay.r && Math.abs(lay.y - at[1]) <= lay.r);
    return near < 0
      ? undefined
      : { layer: lay.key, section: lay.section, system: lay.system, element: { owner: '', key: near } };
  },
};

const Demo: SeqExtension = defineExtension({
  id: 'demo',
  version: '0.2.0',
  requires: '0.6',
  registry: [{ kind: 'layer-type', type: DotsType }],
  params: { 'demo:dots': { radius: { type: 'number', min: 0, doc: 'The dots’ radius, px.' } } },
});

function spec(layers = 1, extra: Record<string, unknown> = {}): Node {
  // the builder writes the spec's `extensions` for the extensions it names
  const b = createSeqSpec({ title: 'dots', extensions: [Demo] });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTACGTACGTACGTACGT' });
  const r = b.data.ranges('r', {
    target: s,
    rows: { start: new Int32Array([0, 10]), end: new Int32Array([4, 16]) },
  });
  b.views.view('v', (v) => {
    const sec = v.section('sec', { of: s });
    for (let i = 0; i < layers; i++)
      sec
        .track(`t${i}`)
        .layer('dots', { type: 'demo:dots', fallback: 'ranges', data: r, radius: 4.5, ...extra } as never);
  });
  const root = b.build();
  expect(root.params?.extensions).toEqual({ demo: '0.2' });
  return root;
}

const codes = (issues: readonly Issue[]) => issues.map((i) => i.code);
/** The data pass's rects of a radius: the dots. */
const dots = (lists: { passes: { data: { rects: { n: number; radius: Float32Array } } } }, r: number) =>
  Array.from(lists.passes.data.rects.radius.subarray(0, lists.passes.data.rects.n)).filter((x) => x === r).length;

describe('extensions: a namespaced layer type', () => {
  it('is drawn by a composition that has it, and as its fallback by one without it', async () => {
    const withIt = createSeqContext({ extensions: [Demo] });
    const without = createSeqContext();
    const root = spec();
    const view = withIt.view(withIt.resolve(root))!;
    const layer = view.sections[0].tracks[0].layers[0];
    expect([layer.type, layer.fallbacks, layer.shape, layer.placed]).toEqual([
      'demo:dots',
      ['demo:dots', 'ranges'],
      'ranges',
      true,
    ]);
    expect(layer.params).toMatchObject({ radius: 4.5, fallback: 'ranges' });
    expect(codes(view.issues)).toEqual([]);
    const other = without.view(without.resolve(root))!;
    expect(other.sections[0].tracks[0].layers[0].placed).toBe(true);
    expect(other.issues.map((i) => [i.code, i.severity, i.message])).toEqual([
      ['seq.view.layer.extension', 'warning', '1 layer(s) need `demo` 0.2: drawn as `ranges`'],
    ]);
    // one draw list for the canvas and the SVG: the dots (rounded, so paths) with it, the ranges' bars without it
    const a = await withIt.renderSvg(root, { size: { width: 400, height: 'auto' } });
    const b = await without.renderSvg(root, { size: { width: 400, height: 'auto' } });
    expect(a.svg).not.toEqual(b.svg);
    expect([dots(a.lists as never, 4.5), dots(b.lists as never, 4.5)]).toEqual([2, 0]);
    expect(codes(a.issues)).not.toContain('seq.view.layer.extension');
    expect(codes(b.issues)).toContain('seq.view.layer.extension');
  });

  it('says what a reader lacks once a namespace and view, and what draws nothing', () => {
    const without = composeSeqContext({ registry: [{ kind: 'layer-type', type: RangesType }] });
    const v = without.view(without.resolve(spec(3)))!;
    expect(v.issues.filter((i) => i.code === 'seq.view.layer.extension').map((i) => i.message)).toEqual([
      '3 layer(s) need `demo` 0.2: drawn as `ranges`',
    ]);
    const none = without.view(without.resolve(spec(1, { fallback: 'none' })))!;
    expect(none.sections[0].tracks[0].layers[0].placed).toBe(false);
    expect(none.issues.map((i) => i.message)).toEqual(['1 layer(s) need `demo` 0.2: not drawn']);
  });

  it('reads `none` and warns when the layer gives no fallback', () => {
    const root = spec(1, { fallback: undefined });
    expect(codes(createSeqContext().check(root))).toEqual(['seq.lang.layer.fallback']);
    const v = createSeqContext().view(createSeqContext().resolve(root))!;
    expect(v.sections[0].tracks[0].layers[0].fallbacks).toEqual(['demo:dots']);
    expect(v.issues.map((i) => i.message)).toEqual(['1 layer(s) need `demo` 0.2: not drawn']);
  });

  it('checks the shape its data has against the shape the type reads', () => {
    const Values: LayerType = { ...DotsType, name: 'demo:dots', shape: 'values' } as LayerType;
    const values = defineExtension({ ...Demo, defined: undefined, registry: [{ kind: 'layer-type', type: Values }] });
    const ctx = createSeqContext({ extensions: [values] });
    const v = ctx.view(ctx.resolve(spec()))!;
    expect(v.sections[0].tracks[0].layers[0].placed).toBe(false);
    expect(codes(v.issues)).toEqual(['seq.view.layer.shape']);
  });

  it('draws a newer bare type as its fallback, a warning and never the enum’s error', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTACGT' });
    const r = b.data.ranges('r', { target: s, rows: { start: new Int32Array([0]), end: new Int32Array([4]) } });
    b.views.view('v', (v) =>
      v.section('sec', { of: s }).track('t').layer('l', { type: 'arcs', fallback: 'ranges', data: r }),
    );
    // a newer type's own param, kept by this reader (its `fallback` given): no error
    const at = ['views', 'v', 'sec', 't', 'l'];
    const set = applyPatch(b.build(), { ops: [{ op: 'set', at, key: 'curvature', value: 1 }] } as never);
    if (!set.ok) throw new Error('patch');
    const root = set.root;
    expect(createSeqContext().check(root)).toEqual([]);
    const ctx = createSeqContext();
    const v = ctx.view(ctx.resolve(root))!;
    const layer = v.sections[0].tracks[0].layers[0];
    expect([layer.fallbacks, layer.placed]).toEqual([['arcs', 'ranges'], true]);
    expect(v.issues.map((i) => [i.code, i.severity])).toEqual([['seq.view.layer.type', 'warning']]);
    expect(layer.params.curvature).toBe(1);
  });

  it('says a navigator picture this reader lacks, and draws the axis', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTACGT' });
    b.views.view('v', (v) => {
      const sec = v.section('sec', { of: s });
      sec.track('t').layer('l', { type: 'letters' });
      v.navigator('nav', { of: sec, content: 'demo:square' } as never);
    });
    const ctx = createSeqContext();
    const v = ctx.view(ctx.resolve(b.build()))!;
    expect(v.navigators[0].content).toBe('axis');
    expect(v.issues.map((i) => [i.code, i.address, i.message])).toEqual([
      ['seq.view.navigator.picture', '/views/v/navigators/nav', 'No picture `demo:square` here: drawn as the axis'],
    ]);
  });

  it('checks its own params where it is registered, and carries them unchecked where not', () => {
    const bad = spec(1, { radius: -1 });
    expect(codes(createSeqContext({ extensions: [Demo] }).check(bad))).toEqual(['repr.schema.param.range']);
    expect(codes(createSeqContext().check(bad))).toEqual([]);
  });
});

describe('extensions: registration', () => {
  const other = defineExtension({ id: 'other', version: '1.0.0', requires: '0.6', registry: [] });
  const runtime = (name: string) => ({ name, formats: ['shared'], open: () => Promise.reject(new Error('none')) });

  it('checks a descriptor where it is defined: every offender in one error', () => {
    const stray = {
      ...other,
      defined: undefined,
      registry: [{ kind: 'layer-type' as const, type: { ...DotsType, name: 'dots' } }],
    };
    expect(() => defineExtension({ ...stray, config: [{ key: 'other.ok' }, { key: 'loose' }] })).toThrow(
      /seq.context.extension.namespace: other: dots, other: loose/,
    );
    expect(() => defineExtension({ ...other, id: 'seq' })).toThrow(/namespace/);
    // a viewer's parts are its maker's: the core's would drop them unread
    expect(() => defineExtension({ ...other, behaviors: [] } as SeqExtension)).toThrow(/viewer parts/);
    // a composition takes only defined descriptors
    expect(() => createSeqContext({ extensions: [stray] })).toThrow(/not made by defineExtension/);
  });

  it('is atomic over two extensions: their entries’ conflicts in one error, nothing registered', () => {
    const a = defineExtension({ ...other, id: 'a', registry: [{ kind: 'format-runtime', runtime: runtime('a:x') }] });
    const b = defineExtension({ ...other, id: 'b', registry: [{ kind: 'format-runtime', runtime: runtime('b:y') }] });
    expect(() => createSeqContext({ extensions: [a, b] })).toThrow(/shared \(a:x, b:y\)/);
    const ctx = createSeqContext({ extensions: [a] });
    const before = ctx.types;
    expect(() => ctx._extend([a, b])).toThrow(/shared/);
    expect([ctx.types === before, ctx.version, ctx.extensions]).toEqual([true, 0, [a]]);
    // another descriptor under one id, and the same one twice (one)
    expect(() => createSeqContext({ extensions: [Demo, defineExtension({ ...Demo })] })).toThrow(/twice: demo/);
    expect(createSeqContext({ extensions: [Demo, Demo] }).extensions).toEqual([Demo]);
  });

  it('reads versions by their major, 0.x by its minor, never mixing them', () => {
    expect(() => defineExtension({ ...other, requires: '0.4' })).toThrow(
      /seq.context.extension.version: other requires 0.4/,
    );
    expect([
      requiresCovers('0.5', '0.5.3'),
      requiresCovers('5', '0.5.0'),
      requiresCovers('0.1', '1.0.0'),
      requiresCovers('1.2', '1.3.0'),
      requiresCovers('1.2', '1.1.0'),
    ]).toEqual([true, false, false, true, false]);
    // a reader at 0.9 of a spec written for 1.0: the old contract, so the fallback, said
    const old = createSeqContext({ extensions: [defineExtension({ ...Demo, version: '0.9.0' })] });
    const root = spec();
    const v = old.view(old.resolve({ ...root, params: { ...root.params, extensions: { demo: '1.0' } } }))!;
    expect(v.sections[0].tracks[0].layers[0].fallbacks).toEqual(['ranges']);
    expect(codes(v.issues)).toEqual(['seq.view.layer.extension']);
  });

  it('swaps the tables whole at run time, their version moved (the viewer’s `extend`)', () => {
    const ctx = createSeqContext();
    const root = spec();
    expect(ctx.view(ctx.resolve(root))!.issues.length).toBe(1);
    ctx._extend([Demo]);
    expect([ctx.version, 'demo:dots' in ctx.types, Object.keys(ctx.kinds.extensions!)]).toEqual([1, true, ['demo']]);
    expect(ctx.view(ctx.resolve(root))!.issues).toEqual([]);
    ctx._extend([]);
    expect([ctx.version, 'demo:dots' in ctx.types]).toEqual([2, false]);
    expect(ctx.view(ctx.resolve(root))!.issues.length).toBe(1);
  });

  it('draws the fallback where the reader has the extension at an older major, and says an older minor', () => {
    const old = createSeqContext({ extensions: [defineExtension({ ...Demo, version: '0.1.4' })] });
    const v = old.view(old.resolve(spec()))!;
    expect(v.sections[0].tracks[0].layers[0].fallbacks).toEqual(['ranges']);
    expect(codes(v.issues)).toEqual(['seq.view.layer.extension']);
    const minor = createSeqContext({ extensions: [defineExtension({ ...Demo, version: '1.0.0' })] });
    const root = spec();
    const newer = { ...root, params: { ...root.params, extensions: { demo: '1.2' } } };
    const m = minor.view(minor.resolve(newer))!;
    expect(m.sections[0].tracks[0].layers[0].fallbacks[0]).toBe('demo:dots');
    expect(m.issues.map((i) => [i.code, i.severity])).toEqual([['seq.view.extension.version', 'info']]);
  });
});

describe('extensions: a spec with one reads, prints, patches and diffs the same with and without it', () => {
  /** A seeded generator (a 32-bit LCG). */
  const seeded = (seed: number) => {
    let s = seed >>> 0;
    return () => (s = (Math.imul(s, 1664525) + 1013904223) >>> 0) / 2 ** 32;
  };
  const withIt = createSeqContext({ extensions: [Demo] });
  const without = createSeqContext();
  const options = { language: SeqLanguage, validate: SeqValidateOptions };

  it('over seeded random edits of its params and its namespaced keys', () => {
    const rand = seeded(11);
    let root = spec(2, { 'demo:tint': 'warm' });
    for (let run = 0; run < 40; run++) {
      const layer = `/views/v/sec/t${Math.floor(rand() * 2)}/dots`.split('/').filter(Boolean);
      const key = ['radius', 'demo:tint', 'demo:size', 'fallback'][Math.floor(rand() * 4)];
      const value = key === 'fallback' ? (rand() < 0.5 ? 'ranges' : 'none') : Math.round(rand() * 8);
      const r = applyPatch(root, { ops: [{ op: 'set', at: layer, key, value }] } as never, options);
      if (!r.ok) throw new Error(JSON.stringify(r.issues));
      const next = r.root;
      // the language's checks: the same with and without (the extension's own params apart)
      const strip = (is: Issue[]) => codes(is).filter((c) => !c.startsWith('repr.schema.param'));
      expect(strip(withIt.check(next))).toEqual(strip(without.check(next)));
      // read the same: every layer's type and params as written, whichever type draws it
      const read = (ctx: typeof withIt) =>
        ctx.view(ctx.resolve(next))!.sections[0].tracks.map((t) => t.layers.map((l) => [l.type, l.params]));
      expect(read(withIt)).toEqual(read(without));
      // printed and read back unchanged; a diff applied gives the same tree
      expect(parseSexpr(printSexpr(next)).node).toEqual(next);
      const d = diff(root, next);
      const back = applyPatch(root, d, options);
      expect(back.ok && printSexpr(back.root)).toBe(printSexpr(next));
      root = next;
    }
  });
});

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Navigator pictures (_spec/README.md "Render", Navigators, Pictures): a host's picture registered and drawn as spec
// tracks, its window placed and hit as any; its readout; its size; a picture the reader lacks drawn as the axis and
// said; a picture's track of an extension's layer type drawn as that type where the extension is composed (resolved
// with the composition's layer kinds); an extension's picture params checked where it is composed, kept where not.

import type { Node } from '@molstar/core/repr/tree';
import { composeSeqContext } from '@molstar/seqx-core/context';
import { createSeqContext } from '@molstar/seqx-core/context/default';
import { defineExtension, type LayerType, type NavigatorPicture } from '@molstar/seqx-core/extension';
import { createSeqSpec, stripNavigator } from '@molstar/seqx-core/lang';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { navigatorBox } from '../navigator.js';
import { BuiltInNavigatorPictures } from '../built-ins.js';
import type { Settled } from './helpers.js';
import { P53 } from './fixtures.js';
import { settleNavigator } from './navigator-fixtures.js';

const Size = { width: 800, height: 200 };
const NavSize = { width: 800, height: 30 };
const Protein = '/views/v/protein';

/** A coverage strip: the values its `data` param names, filling the navigator; its readout names the position. */
const Coverage: NavigatorPicture = {
  name: 'host:coverage',
  axes: ['positions'],
  size: 30,
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
  readout: (p) => `coverage at ${p + 1}`,
};

/** The smallest layer type an extension ships: it lays nothing out, draws nothing. */
const Dots: LayerType<object, { key: string }> = {
  name: 'demo:dots',
  shape: 'ranges',
  fallback: 'ranges',
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
  params: { 'demo:square': { side: { type: 'number', min: 1, doc: 'The side.' } } },
});

/** p53 with its domains (as `demo:dots` where asked) and a conservation profile, and a navigator of `params`. */
function spec(nav: Record<string, unknown>, dots = false): Node {
  const b = createSeqSpec({ title: 'pictures' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const d = b.data.ranges('domains', { target: p53, rows: { start: [0, 100, 300], end: [60, 290, 393] } });
  const values = Array.from({ length: P53.length }, (_, i) => (i % 50) / 50);
  const c = b.data.values('conservation', { target: p53, values });
  b.views.view('v', (v) => {
    const s = v.section('protein', { of: p53 });
    s.track('letters').layer('letters', { type: 'letters' });
    s.track('domains').layer('ranges', {
      type: dots ? 'demo:dots' : 'ranges',
      data: d,
      ...(dots ? { fallback: 'none' } : {}),
    } as never);
    s.track('conservation').layer('values', { type: 'values', data: c });
    v.navigator('nav', { ...stripNavigator(s), ...nav } as never);
  });
  return b.build();
}

/** The navigator renderer over a composition's tables: its pictures, layer types and kinds. */
const composed = (ctx: ReturnType<typeof composeSeqContext>) => (s: Settled) => ({
  navigator: ctx.view(s.model)!.navigators[0],
  types: ctx.types,
  pictures: ctx.pictures,
  kinds: ctx.kinds,
});

describe('a host’s picture', () => {
  const ctx = createSeqContext({ registry: [{ kind: 'navigator-picture', picture: Coverage }] });

  it('is drawn as its tracks, its window placed and hit as any picture’s', async () => {
    const s = await settleNavigator(
      spec({ content: 'host:coverage', data: '/data/conservation' }),
      0,
      [],
      {},
      composed(ctx),
    );
    expect(s.nav.entry).toBe(Coverage);
    const viewport = { start: 100, end: 200 };
    const { nav } = await s.frames({ size: Size, viewport }, { size: NavSize, moves: { [Protein]: 1 } });
    const pic = s.nav.pictureSection!;
    // one track filling the navigator, its layer the values the param names
    expect(pic.tracks.map((t) => [t.track.address.split('/').pop(), t.layers[0].type.name])).toEqual([
      ['coverage', 'values'],
    ]);
    expect(pic.tracks[0].box.height).toBeGreaterThan(20);
    const w = nav.windows[0];
    const axis = pic.axis!;
    expect(w.shown).toBe(true);
    expect(w.x0).toBeCloseTo(axis.xOfUnit(100), 9);
    expect(w.x1).toBeCloseTo(axis.xOfUnit(200), 9);
    expect(s.nav.hit((w.dx0 + w.dx1) / 2, 15)).toMatchObject({ section: Protein, part: 'body' });
    expect(s.nav.hit(w.dx0 - 2, 15)?.part).toBe('start');
    expect(s.nav.sectionUnitAt(Protein, w.x0)).toBeCloseTo(100, 0);
    // the window follows the view, the picture's base kept
    const moved = await s.frames(
      { size: Size, viewport: { start: 250, end: 300 } },
      { size: NavSize, moves: { [Protein]: 1 } },
    );
    expect(moved.nav.windows[0].x0).toBeCloseTo(axis.xOfUnit(250), 9);
    expect(moved.nav.issues).toEqual([]);
    // its readout, and its box
    expect(s.nav.readout(9)).toBe('coverage at 10');
    expect(navigatorBox(s.nav.navigator, 800, 300, ctx.pictures)).toEqual({ width: 800, height: 30 });
  });

  it('a built-in has no readout of its own: the label is the caller’s', async () => {
    const s = await settleNavigator(spec({ content: 'tracks' }));
    expect(s.nav.entry).toBe(BuiltInNavigatorPictures.tracks);
    expect(s.nav.readout(9)).toBeUndefined();
  });

  it('a reader without it draws the axis, and says so; its params are kept', async () => {
    const written = spec({ content: 'host:coverage', data: '/data/conservation' });
    const bare = createSeqContext();
    expect(bare.check(written)).toEqual([]);
    const v = bare.view(bare.resolve(written))!;
    expect([v.navigators[0].content, v.navigators[0].params.data]).toEqual(['axis', '/data/conservation']);
    expect(v.issues.map((i) => [i.code, i.message])).toEqual([
      ['seq.view.layer.extension', 'No picture `host:coverage` here: drawn as the axis'],
    ]);
    const s = await settleNavigator(written, 0, [], {}, composed(bare));
    await s.frames({ size: Size }, { size: NavSize });
    expect(s.nav.pictureSection!.tracks.map((t) => t.layers[0].type.name)).toEqual(['ruler']);
  });

  it('a composition that lists no picture draws every navigator as the axis, said', () => {
    const none = composeSeqContext();
    const v = none.view(none.resolve(spec({ content: 'tracks' })))!;
    expect(v.navigators[0].content).toBe('axis');
    expect(v.issues.map((i) => i.message)).toEqual(['No picture `tracks` here: drawn as the axis']);
    // the axis needs no entry
    const axis = none.view(none.resolve(spec({ content: 'axis' })))!;
    expect(axis.issues).toEqual([]);
  });

  it('two pictures under one name are refused, every one listed', () => {
    const other = { ...Coverage, size: 12 };
    expect(() =>
      composeSeqContext({
        registry: [
          { kind: 'navigator-picture', picture: Coverage },
          { kind: 'navigator-picture', picture: other },
        ],
      }),
    ).toThrow(/navigator pictures under one name: host:coverage/);
  });
});

describe('a picture’s tracks of an extension’s layer type', () => {
  it('are drawn as that type where the extension is composed (a fallback `none`: nothing where not)', async () => {
    const written = spec({ content: 'tracks' }, true);
    const withDemo = createSeqContext({ extensions: [Demo] });
    const s = await settleNavigator(written, 0, [], {}, composed(withDemo));
    await s.frames({ size: Size }, { size: NavSize });
    const types = (n: typeof s.nav) =>
      n.pictureSection!.tracks.map((t) =>
        t.layers.map((l) => `${l.type.name}${l.layer.placed ? '' : ' (not placed)'}`).join(),
      );
    expect(types(s.nav)).toEqual(['letters', 'demo:dots', 'values']);
    expect(s.nav.allIssues.map((i) => i.code)).not.toContain('seq.view.layer.extension');
    const bare = createSeqContext();
    const b = await settleNavigator(written, 0, [], {}, composed(bare));
    await b.frames({ size: Size }, { size: NavSize });
    // its fallback `none`: a reader without the extension draws nothing there
    expect(types(b.nav)).toEqual(['letters', '', 'values']);
  });
});

describe('an extension’s picture params', () => {
  it('are checked where it is composed, and carried unchecked where not', () => {
    const bad = spec({ content: 'demo:square', side: -1 });
    const codes = (issues: readonly { code: string }[]) => issues.map((i) => i.code);
    expect(codes(createSeqContext({ extensions: [Demo] }).check(bad))).toEqual(['repr.schema.param.range']);
    expect(codes(createSeqContext().check(bad))).toEqual([]);
    expect(resolveSpec(bad).spec.children).toBeDefined();
  });
});

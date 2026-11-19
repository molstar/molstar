/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Navigator fixtures: a view settled with a navigator of its own (the navigator's wants answered by the same store),
// a chromosome with its stained bands and a section on a region of it.

import type { Node } from '@molstar/core/repr/tree';
import type { FormatRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec, ideogramNavigator, stripNavigator } from '@molstar/seqx-core/lang';
import {
  bundledTextMetrics,
  type Frame,
  type FrameInput,
  LightTokens,
  type NavigatorFrame,
  type NavigatorFrameInput,
  NavigatorRenderer,
  type NavigatorRendererOptions,
  type RenderConfigOverrides,
  renderConfig,
} from '../index.js';
import { type Settled, settle } from './helpers.js';

export interface SettledNavigator extends Settled {
  readonly nav: NavigatorRenderer;
  /** The view's frame, then the navigator's over it, once both have settled. */
  frames(input: FrameInput, nav: Omit<NavigatorFrameInput, 'view'>): Promise<{ view: Frame; nav: NavigatorFrame }>;
}

export async function settleNavigator(
  spec: Node,
  navigator = 0,
  runtimes: readonly FormatRuntime[] = [],
  config: RenderConfigOverrides = {},
  options: (s: Settled) => Partial<NavigatorRendererOptions> = () => ({}),
): Promise<SettledNavigator> {
  const s = await settle(spec, runtimes, config);
  const nav = new NavigatorRenderer({
    model: s.model,
    view: s.view,
    navigator: s.view.navigators[navigator],
    data: s.store,
    tokens: LightTokens,
    text: bundledTextMetrics(),
    config: renderConfig(config),
    ...options(s),
  });
  return {
    ...s,
    nav,
    async frames(input, navInput) {
      const view = await s.frame(input);
      let f = nav.frame({ ...navInput, view: view.layout });
      for (let i = 0; i < 4; i++) {
        s.store.want(f.wants);
        await s.store.whenSettled();
        const next = nav.frame({ ...navInput, view: view.layout });
        if (!next.rebuilt) break;
        f = next;
      }
      return { view: s.renderer.frame(input), nav: nav.frame({ ...navInput, view: view.layout }) };
    },
  };
}

/** A chromosome's bands (a few, cytoBandIdeo's form), the centromere's two `acen` bands among them. */
export const Bands = {
  start: [0, 200, 400, 520, 600, 700, 900],
  end: [200, 400, 520, 600, 700, 900, 1000],
  id: ['p13', 'p12', 'p11.2', 'p11.1', 'q11.1', 'q12', 'q13'],
  stain: ['gneg', 'gpos50', 'gpos25', 'acen', 'acen', 'gpos100', 'gvar'],
};

/** A 1000-letter chromosome with its bands, a section on `[300, 360)` of it and an ideogram of the chromosome. */
export function ideogramSpec(options: { readonly reversed?: boolean } = {}): Node {
  const b = createSeqSpec({ title: 'ideogram' });
  const chr = b.data.sequence('chr', { alphabet: 'dna', letters: 'ACGT'.repeat(250) });
  const bands = b.data.ranges('bands', { target: chr, rows: Bands });
  b.views.view('v', (v) => {
    const s = v.section('region', {
      of: chr,
      range: [300, 360],
      ...(options.reversed ? { direction: 'reversed' as const } : {}),
    });
    s.track('letters').layer('letters', { type: 'letters' });
    v.navigator('chr', ideogramNavigator(s, bands, { system: chr }));
  });
  return b.build();
}

/**
 * A chromosome of `length` (cyclic or not) with a section on `range` of it (wrapping the origin when the system is
 * cyclic and the range runs past it), other sections on `others` (a range, or `whole`: the chromosome whole), and a
 * navigator of them all on the chromosome: an ideogram of its bands, or a strip of its axis or its tracks.
 */
export function regionSpec(
  options: {
    readonly length?: number;
    readonly range?: readonly [number, number];
    readonly cyclic?: boolean;
    readonly reversed?: boolean;
    readonly content?: 'bands' | 'axis' | 'tracks';
    readonly place?: 'top' | 'bottom' | 'left' | 'right';
    readonly others?: readonly (readonly [number, number] | 'whole')[];
  } = {},
): Node {
  const length = options.length ?? 1000;
  const b = createSeqSpec({ title: 'region' });
  const chr = b.data.sequence('chr', {
    alphabet: 'dna',
    letters: 'ACGT'.repeat(length / 4),
    ...(options.cyclic ? { cyclic: true } : {}),
  });
  const bands = b.data.ranges('bands', {
    target: chr,
    rows: { start: [0, length / 2], end: [length / 2, length], id: ['p1', 'q1'], stain: ['gneg', 'gpos50'] } as never,
  });
  b.views.view('v', (v) => {
    const s = v.section('region', {
      of: chr,
      range: [...(options.range ?? [300, 330])] as [number, number],
      ...(options.reversed ? { direction: 'reversed' as const } : {}),
    });
    s.track('letters').layer('letters', { type: 'letters' });
    const of: string[] = [s.address];
    (options.others ?? []).forEach((o, i) => {
      const t = v.section(`other${i}`, { of: chr, ...(o === 'whole' ? {} : { range: [...o] as [number, number] }) });
      t.track('letters').layer('letters', { type: 'letters' });
      of.push(t.address);
    });
    const content = options.content ?? 'bands';
    const extra = { system: chr, ...(options.place ? { place: options.place } : {}) };
    v.navigator(
      'chr',
      content === 'bands'
        ? ideogramNavigator(of.length > 1 ? of : s, bands, extra)
        : stripNavigator(of.length > 1 ? of : s, { ...extra, content }),
    );
  });
  return b.build();
}

/** A view spec with a navigator node added to its first view. */
export function withNavigator(spec: Node, params: Record<string, unknown>): Node {
  const views = spec.children!.find((c) => c.kind === 'views')!;
  const view = views.children![0];
  const nav: Node = {
    kind: 'navigators',
    id: 'navigators',
    children: [{ kind: 'navigator', id: 'nav', params: params as never }],
  };
  const v2 = { ...view, children: [...(view.children ?? []), nav] };
  return {
    ...spec,
    children: spec.children!.map((c) => (c === views ? { ...views, children: [v2, ...views.children!.slice(1)] } : c)),
  };
}

/** A chromosome of 60,000 positions and its 60 bands (the stains in turn, the centromere's two at its middle). */
export function bandsSpec(): Node {
  const stains = ['gneg', 'gpos25', 'gneg', 'gpos50', 'gneg', 'gpos75', 'gneg', 'gpos100', 'gvar', 'stalk'];
  const n = 60;
  const b = createSeqSpec({ title: 'sixty bands' });
  const chr = b.data.sequence('chr', { alphabet: 'dna', letters: 'ACGT'.repeat(15_000) });
  const bands = b.data.ranges('bands', {
    target: chr,
    rows: {
      start: Array.from({ length: n }, (_, i) => i * 1000),
      end: Array.from({ length: n }, (_, i) => (i + 1) * 1000),
      id: Array.from({ length: n }, (_, i) => `b${i}`),
      stain: Array.from({ length: n }, (_, i) => (i === 29 || i === 30 ? 'acen' : stains[i % stains.length])),
    } as never,
  });
  b.views.view('v', (v) => {
    const s = v.section('region', { of: chr, range: [20_000, 21_000] });
    s.track('letters').layer('letters', { type: 'letters' });
    v.navigator('chr', ideogramNavigator(s, bands, { system: chr }));
  });
  return b.build();
}

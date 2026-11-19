/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Navigator fixtures: a view settled with a navigator of its own (the navigator's wants answered by the same store),
// a chromosome with its stained bands and a section on a region of it.

import type { Node } from '@molstar/core/repr/tree';
import type { FormatRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec, ideogramNavigator } from '@molstar/seqx-core/lang';
import {
  bundledTextMetrics,
  type Frame,
  type FrameInput,
  LightTokens,
  type NavigatorFrame,
  type NavigatorFrameInput,
  NavigatorRenderer,
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

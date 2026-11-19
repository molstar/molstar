/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BuiltInRender } from '@molstar/seqx-core/render/built-ins';
import { DataStore } from '@molstar/seqx-core/data';
import { AlignmentFileRuntime } from '@molstar/seqx-core/data/alignment-readers';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import {
  attachSources,
  bundledTextMetrics,
  DarkTokens,
  LightTokens,
  type LogoLay,
  Renderer,
  renderSvg,
} from '@molstar/seqx-core/render';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { contrastRgba } from '@molstar/seqx-core/theme/index';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import * as globins from '../globins/family.js';
import { nodeFixtureAssets } from '../node-assets.js';

// The Globins seed's logo in both themes: each stacked letter's color is its class's held at 3:1 against the track,
// in the glyph runs a canvas draws and in the SVG's text.

const Size = { width: 900, height: 600 };
const View = { '/views/main/family': { start: 0, end: 40 } };
const runtimes = [AlignmentFileRuntime, RecordsRuntime];
const assets = nodeFixtureAssets().reader;

describe('the Globins logo’s letter colors', () => {
  for (const [theme, tokens] of [
    ['light', LightTokens],
    ['dark', DarkTokens],
  ] as const)
    it(`${theme}: every letter, in the layout and in the SVG, at 3:1 against the track`, async () => {
      const spec = globins.build();
      const model = resolveSpec(spec, { runtimes });
      const store = new DataStore();
      await attachSources(model, store, { assets });
      const renderer = new Renderer({
        ...BuiltInRender,
        model,
        view: resolveViews(model)[0],
        data: store,
        tokens,
        text: bundledTextMetrics(),
      });
      const input = { size: Size, dpr: 1, settled: true, viewports: View } as const;
      let f = renderer.frame(input);
      for (let i = 0; i < 6; i++) {
        store.want(f.wants);
        await store.whenSettled();
        f = renderer.frame(input);
        if (!f.rebuilt) break;
      }
      f = renderer.frame(input);
      let lay: LogoLay | undefined;
      for (const s of f.layout.sections)
        for (const b of s.lays) for (const l of b ?? []) if (l.type.name === 'logo') lay = l.lay as LogoLay;
      expect(lay!.n).toBeGreaterThan(0);
      const surface = tokens.colors['seq-track'];
      const used = new Set<number>();
      for (let i = 0; i < lay!.n; i++) {
        used.add(lay!.color[i]);
        expect(contrastRgba(lay!.color[i], surface)).toBeGreaterThanOrEqual(3);
      }
      const out = await renderSvg(spec, {
        size: Size,
        theme,
        legend: 'none',
        viewports: View,
        runtimes,
        assets,
      } as never);
      for (const c of used) expect(out.svg).toContain(`fill="#${(c >>> 8).toString(16).padStart(6, '0')}"`);
    });
});

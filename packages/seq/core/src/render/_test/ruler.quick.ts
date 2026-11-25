/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ruler's label fit (_spec/README.md "Render", Layer types, `ruler`): majors, then the insertions by priority (the
// nearest a drawn major first), each where it clears what is drawn; its scratch kept between layouts, so a fit allocates
// nothing once warm. By hand, never in CI: `pnpm -C packages/seq/core test:perf`.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { H, Imgt } from '../../resolve/_test/antibody.js';
import { DefaultRenderConfig } from '../config.js';
import { fit, type RulerLay } from '../ruler.js';
import { bundledTextMetrics } from '../text.js';
import { settle } from './helpers.js';
import { layerLay } from './small-helpers.js';

describe('the ruler’s label fit allocates nothing once warm', () => {
  it('IMGT’s insertions at a mid zoom', async ({ task }) => {
    const b = createSeqSpec({ title: 'ruler' });
    const h = b.data.sequence('h', { alphabet: 'protein', letters: H }, (s) =>
      s.numbering('imgt', { method: 'labels', labels: Imgt, scheme: 'imgt' }),
    );
    b.views.view('main', (v) =>
      v.section('imgt', { of: h }).track('r', { label: 'imgt' }).layer('r', { type: 'ruler', numbering: 'imgt' }),
    );
    const s = await settle(b.build());
    const f = await s.frame({ size: { width: 640, height: 60 }, viewport: { start: 0, end: H.length } });
    const lay = layerLay<RulerLay>(f, 0);
    const minors = Array.from({ length: lay.n }, (_, j) => j).filter((j) => !lay.major[j] && lay.labels[j]);
    expect(minors.length).toBeGreaterThan(0);
    const ctx = { config: DefaultRenderConfig, text: bundledTextMetrics() };
    const probe = quickCase(task, 'ruler fit (insertions)', () => fit(lay, ctx), {});
    console.info(`ruler fit: ${probe.r.bytesPerStep.toFixed(0)} bytes a fit (${probe.r.note})`);
    expect(probe.r.bytesPerStep, probe.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});

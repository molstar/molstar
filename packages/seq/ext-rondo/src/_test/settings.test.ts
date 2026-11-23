/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { configItem } from '@molstar/seqx-core/context';
import { createSeqContext } from '@molstar/seqx-core/context/default';
import { BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import { createSeqViewer } from '@molstar/seqx-viewer';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { describe, expect, it } from 'vitest';
import { composedConfig, DefaultRondoConfig, RondoConfigItems } from '../config.js';
import { rondo, rondoWith } from '../index.js';
import type { RondoLay } from '../lay.js';
import { rondoViewer } from '../viewer.js';
import { genomeChromosomes, lay, ringSpec, Size, settle, toyLinks } from './helpers.js';

// The extension's settings readable by its layer types (_spec/README.md, Configuration): the composition's `config`
// (`rondo.*`) reaches `rondo:ribbons` through the render config it is given, over the descriptor's numbers, so
// `rondoWith` is a convenience, not the only door.

describe('settings as composed', () => {
  it('the composition keeps the settings under a composed extension’s namespace, none else', () => {
    const ctx = createSeqContext({
      extensions: [rondo],
      config: [
        [RondoConfigItems.stemPx, 64],
        [configItem<number>('other.stemPx'), 3],
      ],
    });
    expect(ctx.renderConfig.settings).toEqual({ 'rondo.stemPx': 64 });
    expect(createSeqContext({ extensions: [rondo] }).renderConfig.settings).toBeUndefined();
  });

  it('a layout reads them over the descriptor’s, one object a pair (a warm layout makes none)', () => {
    const settings = { 'rondo.stemPx': 64, 'rondo.fade': 0.2 };
    const c = composedConfig(DefaultRondoConfig, { settings });
    expect(c.stemPx).toBe(64);
    expect(c.fade).toBe(0.2);
    expect(c.maxRibbons).toBe(DefaultRondoConfig.maxRibbons);
    expect(composedConfig(DefaultRondoConfig, { settings })).toBe(c);
    expect(composedConfig(DefaultRondoConfig, {})).toBe(DefaultRondoConfig);
  });

  it('the ribbons laid out with `rondo.stemPx` from the composition: wider stems, as with `rondoWith`', async () => {
    const chroms = genomeChromosomes(6);
    const spec = ringSpec({ chroms, links: toyLinks(chroms, 2000, 5) });
    const plain = lay(await (await settle(spec)).frame({ size: Size, dpr: 1 }))!;
    const composed = lay(
      await (await settle(spec, [rondo], {}, false, [[RondoConfigItems.stemPx, 64]])).frame({ size: Size, dpr: 1 }),
    )!;
    const helper = lay(await (await settle(spec, [rondoWith({ stemPx: 64 })])).frame({ size: Size, dpr: 1 }))!;
    expect(composed.cfg.stemPx).toBe(64);
    expect(composed.stemLevel).toBeGreaterThan(plain.stemLevel);
    expect(composed.stemLevel).toBe(helper.stemLevel);
    expect(composed.ns).toBe(helper.ns);
  });
});

describe('settings of a late extension', () => {
  it('`viewer.extend(rondo)` (a descriptor listing its settings): registered, its `rondo.*` applied, gone with its undo', async () => {
    const base = DefaultSeqViewerSpec();
    const v = createSeqViewer({
      ...base,
      registry: [...(base.registry ?? []), { kind: 'format-runtime' as const, runtime: BedRuntime }],
      config: [...(base.config ?? []), [RondoConfigItems.stemPx, 64]],
      ui: undefined,
    });
    try {
      await v.load(ringSpec());
      v.resize(720, 760);
      await v.settle();
      expect(rondo.config?.length).toBeGreaterThan(0);
      const undo = await v.extend(rondo as never);
      await v.settle();
      v.draw();
      expect(v.context.renderConfig.settings).toEqual({ 'rondo.stemPx': 64 });
      const l = v.frame?.layout.sections.flatMap((s) => s.sectionLays).find((x) => x.type.name === 'rondo:ribbons')
        ?.lay as RondoLay | undefined;
      expect(l?.cfg.stemPx).toBe(64);
      undo();
      await v.settle();
      expect(v.context.renderConfig.settings).toBeUndefined();
    } finally {
      v.retire();
    }
  });

  it('`viewer.extend(rondoViewer)`: the `rondo.*` the composition set apply (the render config composed again)', async () => {
    const base = DefaultSeqViewerSpec();
    const v = createSeqViewer({
      ...base,
      registry: [...(base.registry ?? []), { kind: 'format-runtime' as const, runtime: BedRuntime }],
      config: [...(base.config ?? []), [RondoConfigItems.stemPx, 64]],
      ui: undefined,
    });
    try {
      await v.load(ringSpec());
      v.resize(720, 760);
      await v.settle();
      expect(v.context.renderConfig.settings).toBeUndefined();
      const undo = await v.extend(rondoViewer);
      await v.settle();
      v.draw();
      expect(v.context.renderConfig.settings).toEqual({ 'rondo.stemPx': 64 });
      const l = v.frame?.layout.sections.flatMap((s) => s.sectionLays).find((x) => x.type.name === 'rondo:ribbons')
        ?.lay as RondoLay | undefined;
      expect(l?.cfg.stemPx).toBe(64);
      // its undo: the settings gone with it
      undo();
      await v.settle();
      expect(v.context.renderConfig.settings).toBeUndefined();
    } finally {
      v.retire();
    }
  });
});

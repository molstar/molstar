/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { configItem } from '@molstar/seqx-core/context';
import { createSeqContext } from '@molstar/seqx-core/context/default';
import { describe, expect, it } from 'vitest';
import { composedConfig, DefaultRondoConfig, RondoConfigItems } from '../config.js';
import { rondo, rondoWith } from '../index.js';
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

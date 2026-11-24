/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { theme } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { keepHiddenLegends, viewLegends } from '../index.js';
import { settle } from './helpers.js';
import { moleculesSpec } from './molecules-fixtures.js';

// A letters layer's legend row (_spec/README.md "Themes", Legends): titled by the layer ("Sequence letters", never its
// track's name alone), and listing every category the layer draws in view over all its strips: a rule's, the scheme's
// letters of the chains, and what no rule matched (`otherwise`: the chips drawn in the form's else), never the rule alone.

const atp = theme.rules([theme.rule([theme.eq('component', 'ATP')], '#e15759', 'ATP')]);

async function legends(color = atp) {
  const st = await settle(moleculesSpec({ color }));
  const f = await st.frame({ size: { width: 900, height: 400 }, dpr: 1, settled: true });
  return viewLegends(f.layout, Number.POSITIVE_INFINITY).legends;
}

describe('a letters layer’s legend row', () => {
  it('is titled by the layer, and lists the chains’ scheme and the chips’ otherwise beside the rule', async () => {
    const l = (await legends()).find((x) => x.layers.some((k) => k.endsWith('/complex/letters/letters')))!;
    expect(l.title.endsWith(' letters')).toBe(true);
    const labels = l.entries.kind === 'swatches' ? l.entries.items.map((i) => i.label) : [];
    expect(labels[0]).toBe('ATP');
    // the chains' letters in the scheme's groups, and the chips no rule matched (a glycan, an ion, a water)
    expect(labels).toContain('hydrophobic');
    expect(labels).toContain('otherwise');
    expect(new Set(labels).size).toBe(labels.length);
  });

  it('a hidden layer’s row keeps the title the eye names', async () => {
    const st = await settle(moleculesSpec({ color: atp }));
    const f = await st.frame({ size: { width: 900, height: 400 }, dpr: 1, settled: true });
    const shown = viewLegends(f.layout, Number.POSITIVE_INFINITY).legends;
    const names = new Map<string, string>();
    keepHiddenLegends(f.layout, undefined, shown, new Map(), names);
    const k = [...names.keys()].find((x) => x.endsWith('/complex/letters/letters'))!;
    expect(names.get(k)).toBe(shown.find((x) => x.layers.includes(k))!.title);
  });
});

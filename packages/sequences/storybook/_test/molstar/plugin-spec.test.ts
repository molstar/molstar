/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The plugin's composition (plugin.ts): hover, selection and focus, and the loci labels Mol*'s tooltip is made of.

import { DefaultLociLabelProvider, HighlightLoci, SelectLoci } from '@molstar/plugin/behavior/dynamic/representation';
import { expect, it } from 'vitest';
import { pluginSpec } from '../../src/_stories/support/seqstar-molstar-plugin.js';

it('the composition has the highlight, select and loci label behaviors', () => {
  const transformers = (pluginSpec().behaviors ?? []).map((b) => b.transformer);
  expect(transformers).toContain(HighlightLoci);
  expect(transformers).toContain(SelectLoci);
  expect(transformers).toContain(DefaultLociLabelProvider);
});

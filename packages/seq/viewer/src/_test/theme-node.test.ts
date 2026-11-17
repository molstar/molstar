/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SeqCoreConfig } from '@molstar/seqx-core/context';
import { DarkTokens, LightTokens } from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer } from '../viewer.js';

// The headless viewer in Node (_spec/README.md "Configuration", `theme.source`): the config's theme, no DOM read.

describe('the headless viewer', () => {
  it('uses the built-in table of the config’s theme, light by default', () => {
    const light = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    const dark = createSeqViewer({
      ...DefaultSeqViewerSpec(),
      config: [[SeqCoreConfig.theme, 'dark']],
      ui: undefined,
    });
    expect(light.tokens).toBe(LightTokens);
    expect(dark.tokens).toBe(DarkTokens);
    expect(dark.context.tokens).toBe(DarkTokens);
    light.retire();
    dark.retire();
  });
});

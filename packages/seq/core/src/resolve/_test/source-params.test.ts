/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import type { FormatRuntime } from '../../data/protocol.js';
import { resolveSpec } from '../index.js';

// A source's params ride on the entry its runtime opens (_spec/README.md "Resolve", Systems, sources, annotations): a
// GFF's `features` reaches the runtime (it was dropped, so `features: genes` read as transcripts).

describe('a source’s entry', () => {
  it('carries `features` to its runtime', () => {
    const b = createSeqSpec({ title: 'features' });
    b.data.source('genes', { format: 'gff3', uri: 'f.gff3', features: 'genes' } as never);
    const runtime: FormatRuntime = { name: 'gff', formats: ['gff3'], open: () => Promise.reject(new Error('no')) };
    const model = resolveSpec(b.build(), { runtimes: [runtime] });
    const entry = model.sources.get('/data/genes')?.entry as { features?: string } | undefined;
    expect(entry?.features).toBe('genes');
  });
});

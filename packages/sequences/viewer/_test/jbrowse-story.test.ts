/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The JBrowse 2 demo's document (jbrowse-story-fixtures.ts) opened offline, its four sources served by the fixture
// providers in place of the live files: it validates, opens with no issue but one note, and its tracks ask for their
// regions.

import { validate } from '@molstar/seq-core/format';
import { expect, it } from 'vitest';
import { registries } from '../../core/_test/format/load-fixtures.js';
import { GRCH38, hg38Streams } from '../../core/_test/format/stream-provider.js';
import { KernelRenderers, type LayerRenderer } from '../src/index.js';
import type { Streams } from '../src/streams.js';
import { JBROWSE_SECTION, jbrowseDoc } from './jbrowse-story-fixtures.js';
import { opened } from './layout-fixtures.js';

const flush = () => new Promise((ok) => setTimeout(ok, 0));

it("the JBrowse 2 demo's document validates, opens with no issue but a note and asks for every track's regions", async () => {
  // On chr2 (the fixture serves three chromosomes), opening on 300 positions: letters asked too.
  const doc = jbrowseDoc('chr2', 100_000_000, 100_000_300);
  expect(validate(doc, registries()).filter((i) => i.severity === 'error')).toEqual([]);
  const s = hg38Streams({ chromosomes: GRCH38.slice(0, 3), letters: true });
  const by: Record<string, keyof typeof s> = { genome: 'genome', gc: 'phylop', genes: 'genes', clinvar: 'variants' };
  const t = await opened(doc, {
    renderers: [...KernelRenderers] as LayerRenderer[],
    loaders: {
      load: () => Promise.reject(new Error('none')),
      stream: (src) =>
        src.stream && by[src.stream] ? (s[by[src.stream]] as unknown as { provider: never }).provider : null,
    },
    size: { [JBROWSE_SECTION]: [1000, 400] },
  });
  // The one note: the synthetic variants in ClinVar's place send no `category`, neither a vocabulary nor a column
  // (LB-08), so the ClinVar track's items have the missing color and its legend lists nothing; the note says so.
  const notes = ['svs.legend.unlisted'];
  const issues = () => t.ctx.issues.current.value.filter((i) => !(i.severity === 'info' && notes.includes(i.code)));
  const info = t.ctx.issues.current.value.filter((i) => i.severity === 'info');
  expect(info.map((i) => i.code)).toEqual(notes);
  expect(info[0].message).toMatch(/no values of "category".*an item without one has the missing color/);
  expect(issues().map((i) => `${i.code}: ${i.message}`)).toEqual([]);
  t.layout(JBROWSE_SECTION);
  for (let k = 0; k < 3; k++) {
    (t.ctx.streams as Streams).pass();
    await flush();
  }
  for (const p of [s.genome, s.phylop, s.genes, s.variants]) expect(p.requests.length).toBeGreaterThan(0);
  expect(issues().map((i) => `${i.code}: ${i.message}`)).toEqual([]);
  t.ctx.dispose();
});

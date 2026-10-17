/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { KeyTarget, SeqSet } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import {
  createSeqViewer,
  DefaultRegistry,
  DefaultSeqViewerSpec,
  Hit,
  HitRegion,
  HoverKind,
  manualScheduler,
} from '../src/index.js';
import { fixtureMeasurer } from '../src/testing/index.js';

// The position label of an unobserved residue whose block numbering has no number for it (a model's author numbering
// leaves unobserved residues out): the title still places it, by the first numbering that has one or by its position
// in the sequence, and says so (seqspec:viewer#labels).

type Numbered = { numbering(params: unknown): unknown };

async function labelAt(index: number, numbering: 'segments' | 'with-position') {
  const b = createSeqViewSpec();
  const seq = b.data.sequence({ ref: 'p', value: 'YKVLATGHEK' });
  if (numbering === 'with-position') seq.numbering({ ref: 'p.seq', name: 'seq', method: 'offset', first: 1 });
  const c = b.data.complex({ ref: 'c' });
  // Residues 0–2 unobserved; the author numbering covers the observed ones only (index 3 is 104).
  (c.unit({ ref: 'u', sequence: 'p', observed: [[3, 10]] } as never) as never as Numbered).numbering({
    ref: 'u.auth',
    name: 'auth',
    method: 'segments',
    segments: [[3, 104, 7]],
    default: true,
  } as never);
  b.view({ ref: 'v', label: 'V' })
    .block({ unit: 'u' } as never)
    .track({ label: 'L' })
    .letters();
  const ctx = createSeqViewer(
    { ...DefaultSeqViewerSpec(), registry: [DefaultRegistry] },
    { scheduler: manualScheduler(), text: fixtureMeasurer() },
  );
  await ctx.setDocument(b.getDocument());
  await ctx.graph.update();
  const set = ctx.graph.assembled.value?.set as SeqSet;
  const t: KeyTarget = { kind: -1, id: -1, row: -1 };
  set.byName('u', t);
  const unit = t.id;
  const hit = Hit.create();
  Object.assign(hit, {
    section: 0,
    region: HitRegion.data,
    kind: HoverKind.seq,
    layer: 0,
    block: 0,
    shown: index,
    index,
  });
  Object.assign(hit.seq, { set, unit, index });
  ctx.hover.setFromHit(hit, 'pointer');
  return ctx.labels.compute(ctx.hover.live.current);
}

describe('the label of an unobserved residue', () => {
  it('keeps its number from the block numbering where it has one', async () => {
    const label = await labelAt(4, 'segments');
    expect(label?.title).toMatch(/105/);
    expect(JSON.stringify(label)).not.toMatch(/unobserved/);
  });

  it('falls back to its position in the sequence, said in a line, and stays marked unobserved', async () => {
    const label = await labelAt(1, 'segments');
    expect(label?.title).toMatch(/K ?2$/);
    const text = JSON.stringify(label);
    expect(text).toMatch(/position in the sequence/);
    expect(text).toMatch(/not in the model/);
  });

  it("takes another numbering's number before the bare position", async () => {
    const label = await labelAt(0, 'with-position');
    expect(label?.title).toMatch(/1\b/);
    expect(JSON.stringify(label)).toMatch(/numbered by/);
  });
});

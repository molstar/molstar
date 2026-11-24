/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section of a source's complex costs a warm frame nothing: the view is resolved over the opened model once a node
// version, so with a hover on a chain's block a frame drawn again resolves nothing and makes no object; after the host
// drops a chain (one new opened model, one new view) the frame is warm again at the same 0 bytes, quiet within the
// harness's warm-up. The hover is a residue's of the opened model's unit (its system asserted). By hand, never in
// CI: `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.

import { readFileSync } from 'node:fs';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { FakeStructureViewer } from '../../../support/fake-structure/fake-viewer.js';
import { HostTableFile, hostComplexSpec, hostTableOf, unitAddress, withFakeHost } from '../compose.js';

/**
 * The pointer over a residue of the first block (chain A) on the residues track; the hover's target (none: an element,
 * not a header part), its system and position.
 */
function hoverFirstBlock(v: SeqViewer): [string | undefined, string | undefined, number] {
  const s = v.frame!.layout.sections[0];
  const t = s.tracks[1];
  v.input.move(s.axis!.x(50), t.box.y + t.box.height / 2);
  v.draw();
  return [v.hover.current.target, v.hover.current.system, v.hover.current.position];
}

describe('a section of a source’s complex on the frame’s path', () => {
  it('2OCJ’s complex served by the fake host: a warm frame with a hover, before and after a chain is dropped', async ({
    task,
  }) => {
    const fake = new FakeStructureViewer(hostTableOf(readFileSync(`${PackageFolder}${HostTableFile}`, 'utf8')));
    const base = DefaultSeqViewerSpec();
    const v = createSeqViewer({ ...base, ui: undefined, registry: withFakeHost(base.registry ?? [], fake) });
    v.resize(1440);
    await v.load(hostComplexSpec());
    await v.settle();
    v.draw();
    expect(v.view.resolved.value?.sections[0].blocks).toHaveLength(8);
    expect(hoverFirstBlock(v)).toEqual([undefined, unitAddress('A'), expect.any(Number)]);
    const before = quickCase(task, 'viewer frame over a host’s complex', () => void v.draw(), {}).r;
    fake.dropChain('D');
    // the reopen is asked for the next frame: drawn, then settled
    for (let i = 0; i < 4 && v.view.resolved.value?.sections[0].blocks?.length !== 6; i++) {
      v.draw();
      await new Promise((r) => setTimeout(r, 20));
      await v.settle();
    }
    v.draw();
    expect(v.view.resolved.value?.sections[0].blocks).toHaveLength(6);
    expect(hoverFirstBlock(v)).toEqual([undefined, unitAddress('A'), expect.any(Number)]);
    const after = quickCase(task, 'viewer frame over a host’s complex, a chain dropped', () => void v.draw(), {}).r;
    console.info(
      `frame over a host's complex: ${before.bytesPerStep.toFixed(0)} bytes (${before.note}); after a drop ${after.bytesPerStep.toFixed(0)} (${after.note})`,
    );
    v.retire();
    fake.dispose();
    expect(before.bytesPerStep, before.note).toBeLessThanOrEqual(inject('perfThreshold'));
    expect(after.bytesPerStep, after.note).toBeLessThanOrEqual(inject('perfThreshold'));
    // warm in the harness's own warm-up: no function the frame calls compiled again a step (a record keyed by the
    // rails' addresses deoptimized on each new view: 3.5 B and 35 B of code a step, never quiet)
    expect([after.settled, after.codePerStep], after.note).toEqual(['quiet', 0]);
  });
});

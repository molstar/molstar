/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section of a source's complex costs a warm frame nothing: the view is resolved over the opened model once a node
// version, so with a hover on a chain's block a frame drawn again resolves nothing and makes no object; after the host
// drops a chain (one new opened model, one new view) the frame is warm again at the same 0 bytes. By hand, never in
// CI: `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.

import { readFileSync } from 'node:fs';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { FakeStructureViewer } from '../../../support/fake-structure/fake-viewer.js';
import { HostTableFile, hostComplexSpec, hostTableOf, withFakeHost } from '../compose.js';

/** The pointer over the first block's residues (its track, 30 px in). */
function hoverFirstBlock(v: SeqViewer): void {
  const s = v.frame!.layout.sections[0];
  const t = s.tracks[1];
  v.input.move(s.axis!.sx0 + 30, t.box.y + Math.min(6, t.box.height / 2));
  v.draw();
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
    hoverFirstBlock(v);
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
    hoverFirstBlock(v);
    // the drop's new view and renderer send the frame's path through the compiler again (JIT code, not objects): warmed
    // past the harness's own cap before it is read
    for (let i = 0; i < 100_000; i++) v.draw();
    const after = quickCase(task, 'viewer frame over a host’s complex, a chain dropped', () => void v.draw(), {}).r;
    console.info(
      `frame over a host's complex: ${before.bytesPerStep.toFixed(0)} bytes (${before.note}); after a drop ${after.bytesPerStep.toFixed(0)} (${after.note})`,
    );
    v.retire();
    fake.dispose();
    expect(before.bytesPerStep, before.note).toBeLessThanOrEqual(inject('perfThreshold'));
    expect(after.bytesPerStep, after.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});

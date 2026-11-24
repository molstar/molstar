/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @module-tag medium

import { location } from '@molstar/seqx-core/query';
import { ComplexViews } from '@molstar/seqx-examples';
import { ComplexAddress } from '@molstar/seqx-examples/complex/units';
import { afterAll, beforeAll, describe, expect, it, vi } from 'vitest';
import { linkFake } from '../compose.js';
import { fakeClock, openView, tableOf } from './open.js';

// The Focus and back story's dataset (PDB 1FIN as a complex), headless on a fake clock, composed as the story composes
// it: every unit of the complex that is an entry's chain has its row in the recorded chain table (the ATP of the second
// copy is chain F), a hover and a click on the ligand, a focus. A whole example view opened and settled: medium; the
// 2OCJ twins are its quick representative.

const table = tableOf('fixtures/complex/1FIN-chains.json');
let viewer: Awaited<ReturnType<typeof openView>>;
beforeAll(async () => {
  viewer = await openView(ComplexViews.find((v) => v.id === 'units')!);
}, 20_000);
afterAll(() => viewer.retire());

describe('the Complex story’s link', { tags: ['medium'] }, () => {
  it('every unit of 1FIN’s chains is a row of the table and nothing else is: the transplants are units of other entries', () => {
    const restore = fakeClock();
    const link = linkFake(viewer, table, ComplexAddress);
    const rows = new Set(link.fake.rows.map((r) => r.asym));
    const chains = link.adapter.facts.units.filter((u) => !u.id.startsWith('t-')).map((u) => u.chain);
    expect(chains.sort()).toEqual(['A', 'B', 'C', 'D', 'E', 'F']);
    for (const c of chains) expect(rows.has(c!), c).toBe(true);
    link.dispose();
    restore();
  });

  it('a hover on CDK2 of the first copy lights its chain alone; a click on its ATP’s chip selects the unit; a focus moves the camera', async () => {
    const restore = fakeClock();
    const link = linkFake(viewer, table, ComplexAddress);
    const { fake } = link;
    const owner = {};
    viewer.hover.setLocation(location({ parts: [{ system: `${ComplexAddress}/A`, ranges: [100, 103] }] }), owner);
    expect(fake.highlight.value.map((s) => [fake.rows[s.row].asym, s.start, s.end])).toEqual([['A', 100, 103]]);
    viewer.hover.clear(owner);
    const atp = fake.rows.find((r) => r.asym === 'E')!;
    fake.click({ row: atp.index, index: 0 });
    expect(viewer.selection.value.value.shown.parts.map((p) => [p.system, ...p.ranges])).toEqual([
      [`${ComplexAddress}/atp-a`, 0, 1],
    ]);
    fake.doubleClick({ row: 0, index: 10 });
    await vi.advanceTimersByTimeAsync(20);
    expect(viewer.view.lastJump?.kind).toBe('focus');
    link.dispose();
    restore();
  });

  it('Back before the first focus names no residue: nothing is sent and the fake’s frame stays (F5, the story says so)', async () => {
    const restore = fakeClock();
    const link = linkFake(viewer, table, ComplexAddress);
    link.fake.doubleClick({ row: 0, index: 10 });
    await vi.advanceTimersByTimeAsync(20);
    const frame = link.fake.camera.value;
    const told = link.log.value.length;
    viewer.view.back();
    await vi.advanceTimersByTimeAsync(20);
    expect(viewer.view.lastJump).toMatchObject({ kind: 'back' });
    expect(viewer.view.lastJump?.location).toBeUndefined();
    expect(link.log.value).toHaveLength(told);
    expect(link.fake.camera.value).toBe(frame);
    link.dispose();
    restore();
  });
});

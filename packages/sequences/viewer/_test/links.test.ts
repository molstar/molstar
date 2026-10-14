/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The default keys against the viewer's own `ctx.commands.run` (CM1) on case 10: every invocation a link makes passes
// the built-ins' params tables (link args arrive as strings), lights and clears its channel, and acts; on a cyclic
// sequence a link's range may run across the origin (CG1b).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { Invocation } from '@molstar/seq-core/format';
import { createLinkRunner, formatLink } from '@molstar/seq-core/links';
import { OrderedSet } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { case10 } from '../../core/_test/format/cases.js';
import { driven } from './input-fixtures.js';

describe('default keys on the viewer', () => {
  it("every default key's invocation passes the check and acts", async () => {
    const t = await driven(case10());
    const runs: Promise<boolean>[] = [];
    const runner = createLinkRunner((inv: Invocation) => {
      const p = t.ctx.commands.run(inv, 'host:links');
      runs.push(p);
      return p;
    });
    const T315 = { beg: '315', end: '315', numbering: 'auth' };
    const hl = formatLink('seq-highlight', '1iep/A', T315);
    expect(runner.run(hl, 'enter', 'd1')).toBe(true);
    expect(runner.run(formatLink('seq-select', '1iep/A', { start: 10, stop: 20 }), 'enter', 'd2')).toBe(true);
    await Promise.all(runs);
    expect([...t.ctx.highlight.channels.value.keys()]).toEqual(['links:d1', 'links:d2']);
    expect(t.ctx.highlight.channels.value.get('links:d1')).toMatchObject({ closure: 'containing' });
    runner.run(hl, 'leave', 'd1');
    runner.release('d2');
    runner.run(formatLink('seq-select', '1iep/A', T315), 'click', 'd3');
    runner.run(formatLink('seq-focus', '1iep/A', { ...T315, section: 'chains.main' }), 'click', 'd4');
    runner.run(formatLink('seq-focus', 'site.gatekeeper', { items: '0' }), 'click', 'd5');
    runner.run('!seq-snapshot=next', 'click', 'd6');
    runner.run('!seq-view=alignment', 'click', 'd7');
    expect(await Promise.all(runs)).toEqual(Array(9).fill(true));
    expect(t.ctx.highlight.channels.value.size).toBe(0);
    expect(t.ctx.selection.current.value.seq).not.toBeNull();
    expect(t.ctx.selection.current.value.origin).toBe('host:links');
    expect(t.ctx.snapshots.current.value).toBe('s1');
    expect(t.ctx.views.current.value).toBe('alignment');
  });

  it('a range across a cyclic origin (start > stop): the link highlights and selects both runs', async () => {
    const b = createSeqViewSpec();
    const p = b.data.sequence({ ref: 'p', value: 'ACGTACGTTA'.repeat(24), alphabet: 'dna', cyclic: true });
    b.view({ ref: 'v', label: 'V' }).block({ sequence: p }).track({ label: 'S' }).letters();
    const t = await driven(b.getDocument());
    const runs: Promise<boolean>[] = [];
    const runner = createLinkRunner((inv: Invocation) => {
      const r = t.ctx.commands.run(inv, 'host:links');
      runs.push(r);
      return r;
    });
    const across = formatLink('seq-select', 'p', { start: 200, stop: 30 });
    expect(runner.run(across, 'enter', 'd1')).toBe(true);
    runner.run(across, 'click', 'd1');
    expect(await Promise.all(runs)).toEqual([true, true]);
    const size = (l: unknown) =>
      ((l as { elements?: { indices: unknown }[] } | null)?.elements ?? []).reduce(
        (n, e) => n + OrderedSet.size(e.indices as OrderedSet),
        0,
      );
    expect(size(t.ctx.highlight.channels.value.get('links:d1')?.target)).toBe(70);
    expect(size(t.ctx.selection.current.value.seq)).toBe(70);
  });
});

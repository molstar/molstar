/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Commands (CM1; seqspec:viewer#commands; seqspec:format#descriptions-and-invocations): dispatch in subscription order
// and its errors, `run` on plain data (an unknown id or invalid params resolve false and run nothing; the caller's
// origin), each built-in from an invocation, closure before the verb, `Coordinates` targets through `ctx.systems`, an
// extension's command registered like any entry, and the built-ins unhandled once the `commands` behavior is removed.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import {
  type Invocation,
  list,
  nullable,
  OptionalField,
  RequiredField,
  SimpleParamsSchema,
  str,
} from '@molstar/seq-core/format';
import { type Loci, OrderedSet, type SeqLoci } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { case3, case10 } from '../../core/_test/format/cases.js';
import { gene } from '../../core/_test/format/project-fixtures.js';
import type { Commands } from '../src/commands.js';
import {
  createSeqViewer,
  DefaultSeqViewerSpec,
  defineBehavior,
  type FocusChange,
  isCoordinates,
  type PatchEvent,
  SeqCommand,
  SeqCommands,
  type SeqIssue,
  SeqRegistrationError,
  type SeqViewerContext,
  type SeqViewerRegistryEntry,
} from '../src/index.js';
import { driven, subscriptions } from './input-fixtures.js';
import { shown } from './projection-fixtures.js';

const T315 = { refs: ['1iep/A'], args: { beg: '315', end: '315', numbering: 'auth' } };
const indices = (l: unknown) => {
  const out: number[] = [];
  for (const e of (l as SeqLoci).elements) OrderedSet.forEach(e.indices, (i: number) => out.push(i), null);
  return out;
};
const keyed = (l: Loci | null) => {
  if (l?.kind !== 'seq-loci') return l;
  const values = (s: OrderedSet) => Array.from({ length: OrderedSet.size(s) }, (_, i) => OrderedSet.getAt(s, i));
  return l.elements.map((e) => [l.set.strings.get(l.set.units.key[e.unit]), ...values(e.indices)]);
};
function issuesOf(ctx: SeqViewerContext) {
  const out: SeqIssue[] = [];
  ctx.events.issue.subscribe((i) => out.push(i));
  return out;
}
const commandsOf = (ctx: SeqViewerContext) => ctx.commands as Commands;

describe('dispatch (ported)', () => {
  it('a dispatch with no subscriber resolves', async () => {
    const ctx = createSeqViewer(DefaultSeqViewerSpec());
    await expect(ctx.commands.dispatch(SeqCommands.SetView, { ref: 'v' })).resolves.toBeUndefined();
  });

  it('subscribers run in subscription order, awaited; the command is callable', async () => {
    const ctx = createSeqViewer(DefaultSeqViewerSpec());
    const Pulse = SeqCommand<{ n: number }>('lab:pulse', SimpleParamsSchema({}));
    const order: string[] = [];
    ctx.commands.subscribe(Pulse, async ({ n }) => {
      await Promise.resolve();
      order.push(`a${n}`);
    });
    ctx.commands.subscribe(Pulse, ({ n }) => order.push(`b${n}`));
    await Pulse(ctx, { n: 1 });
    expect(order).toEqual(['a1', 'b1']);
  });

  it('the first rejection rejects the dispatch, later subscribers do not run, seqv.command.error once', async () => {
    const ctx = createSeqViewer(DefaultSeqViewerSpec());
    const issues = issuesOf(ctx);
    const later: string[] = [];
    const sub = ctx.commands.subscribe(SeqCommands.Focus, () => Promise.reject(new Error('no such ref')));
    ctx.commands.subscribe(SeqCommands.Focus, () => later.push('ran'));
    await expect(ctx.commands.dispatch(SeqCommands.Focus, { target: { refs: ['x'] } })).rejects.toThrow('no such ref');
    expect(issues.map((i) => i.code)).toEqual(['seqv.command.error']);
    expect(issues[0].message).toMatch(/seqv:focus.*no such ref/);
    expect(later).toEqual([]);
    sub.unsubscribe();
    await ctx.commands.dispatch(SeqCommands.Focus, { target: { refs: ['x'] } });
    expect(later).toEqual(['ran']);
  });

  it('dispose rejects pending dispatches and every later one, run included', async () => {
    const ctx = createSeqViewer(DefaultSeqViewerSpec());
    ctx.commands.subscribe(SeqCommands.SetView, () => new Promise(() => {}));
    const pending = ctx.commands.dispatch(SeqCommands.SetView, { ref: 'v' });
    ctx.dispose();
    await expect(pending).rejects.toThrow(/disposed/);
    await expect(ctx.commands.dispatch(SeqCommands.SetView, { ref: 'v' })).rejects.toThrow(/disposed/);
    await expect(ctx.commands.run({ command: 'seqv:setView', params: { ref: 'v' } }, 'api')).rejects.toThrow(
      /disposed/,
    );
  });
});

describe('run: plain data', () => {
  /** A context with no behaviors and a recording subscriber on every built-in. */
  function recording() {
    const ctx = createSeqViewer({ ...DefaultSeqViewerSpec(), behaviors: [] });
    const got: [string, unknown][] = [];
    for (const cmd of Object.values(SeqCommands))
      ctx.commands.subscribe(cmd as SeqCommand<unknown>, (p) => got.push([cmd.id, p]));
    return { ctx, got, issues: issuesOf(ctx) };
  }

  it('the built-ins are registry entries; an unknown id resolves false, runs nothing, reports nothing', async () => {
    const { ctx, got, issues } = recording();
    expect(ctx.registry.list('commands').map((c) => c.id)).toEqual([
      'seqv:highlight',
      'seqv:select',
      'seqv:focus',
      'seqv:reveal',
      'seqv:setSnapshot',
      'seqv:setView',
      'seqv:colorScheme',
      'seqv:colorRange',
      'seqv:legend',
      'seqv:retry',
      'seqv:runAction',
    ]);
    for (const inv of [{ command: 'seqv:zoom' }, { command: 'lab:pulse', params: {} }, { command: 3 }, {}])
      expect(await ctx.commands.run(inv as Invocation, 'host:links')).toBe(false);
    expect(await ctx.commands.run({ command: 'seqv:setView', params: { ref: 'v' } }, 'api')).toBe(true);
    expect([got, issues]).toEqual([[['seqv:setView', { ref: 'v' }]], []]);
  });

  it('invalid params resolve false and run nothing: types, missing, unknown, `loci`, bad runs', async () => {
    const { ctx, got, issues } = recording();
    const bad: Invocation[] = [
      { command: 'seqv:focus', params: { target: 3 } },
      { command: 'seqv:focus', params: {} },
      { command: 'seqv:highlight', params: { target: null } }, // no channel
      { command: 'seqv:highlight', params: { channel: 'c', target: null, style: 'glow' } },
      { command: 'seqv:select', params: { target: null, mode: 'replace' } },
      { command: 'seqv:select', params: { target: { refs: ['a'] }, closure: 'unit' } },
      { command: 'seqv:select', params: { target: { refs: ['a'], args: { start: '3' } } } },
      { command: 'seqv:select', params: { target: { refs: 'a' } } },
      { command: 'seqv:select', params: { target: { loci: {} } } }, // loci only through dispatch
      { command: 'seqv:select', params: { target: { coordinates: { systems: ['k'], runs: [0, 5] } } } },
      { command: 'seqv:select', params: { target: { coordinates: { systems: ['k'], runs: [1, 0, 5] } } } },
      { command: 'seqv:select', params: { target: { coordinates: { systems: ['k'], runs: [0, 5, 4] } } } },
      { command: 'seqv:setSnapshot', params: { ref: 's1', animate: 'no' } },
      { command: 'seqv:setView', params: { ref: 'v', extra: 1 } },
      { command: 'seqv:runAction', params: { target: 't' } },
      { command: 'seqv:setView', params: 'v' as never },
    ];
    for (const inv of bad) expect(await ctx.commands.run(inv, 'api'), JSON.stringify(inv)).toBe(false);
    expect([got, issues]).toEqual([[], []]);
  });

  it("the caller's origin goes to the commands whose table takes one, over any origin given", async () => {
    const { ctx, got } = recording();
    const target = { coordinates: { systems: ['seqstar:a'], runs: [0, 1, 2] } };
    await ctx.commands.run({ command: 'seqv:select', params: { target, origin: 'pointer' } }, 'host:agent');
    await ctx.commands.run({ command: 'seqv:runAction', params: { name: 'track.hide', target: 't' } }, 'host:links');
    await ctx.commands.run({ command: 'seqv:highlight', params: { channel: 'c', target: null } }, 'host:links');
    expect(got).toEqual([
      ['seqv:select', { target, origin: 'host:agent' }],
      ['seqv:runAction', { name: 'track.hide', target: 't', origin: 'host:links' }],
      ['seqv:highlight', { channel: 'c', target: null }],
    ]);
  });

  it('a throwing handler: run rejects and reports seqv.command.error', async () => {
    const ctx = createSeqViewer({ ...DefaultSeqViewerSpec(), behaviors: [] });
    const issues = issuesOf(ctx);
    ctx.commands.subscribe(SeqCommands.SetView, () => {
      throw new Error('no view');
    });
    await expect(ctx.commands.run({ command: 'seqv:setView', params: { ref: 'v' } }, 'api')).rejects.toThrow('no view');
    expect(issues.map((i) => [i.code, i.message])).toEqual([
      ['seqv.command.error', 'command "seqv:setView" failed: no view'],
    ]);
  });

  it('an extension command is registered like any entry: checked by its table, undone, a taken id refused', async () => {
    const ctx = createSeqViewer({ ...DefaultSeqViewerSpec(), behaviors: [] });
    const Pulse = SeqCommand<{ link: string; refs: string[] | null }>(
      'lab:pulse',
      SimpleParamsSchema({ link: RequiredField(str, ''), refs: OptionalField(nullable(list(str)), null, '') }),
    );
    const got: unknown[] = [];
    ctx.commands.subscribe(Pulse, (p) => got.push(p));
    const undo = ctx.register({ commands: [Pulse] });
    expect(await ctx.commands.run({ command: 'lab:pulse', params: { link: 'x', refs: ['a'] } }, 'api')).toBe(true);
    expect(await ctx.commands.run({ command: 'lab:pulse', params: { refs: ['a'] } }, 'api')).toBe(false);
    expect(() => ctx.register({ commands: [SeqCommand('lab:pulse', SimpleParamsSchema({}))] })).toThrow(
      SeqRegistrationError,
    );
    expect(() => ctx.register({ commands: [SeqCommand('seqv:focus', SimpleParamsSchema({}))] })).toThrow(
      /already registered/,
    );
    undo();
    expect(await ctx.commands.run({ command: 'lab:pulse', params: { link: 'y' } }, 'api')).toBe(false);
    await Pulse(ctx, { link: 'z', refs: null }); // dispatch from code needs no registration
    expect(got).toEqual([
      { link: 'x', refs: ['a'] },
      { link: 'z', refs: null },
    ]);
  });
});

describe('the built-ins from invocations on case 10', () => {
  const run = (ctx: SeqViewerContext, command: string, params: Record<string, unknown>) =>
    ctx.commands.run({ command: `seqv:${command}`, params }, 'host:links');

  it('highlight: only its channel, with style and closure; null clears it; a target naming nothing does nothing', async () => {
    const t = await driven(case10());
    const issues = issuesOf(t.ctx);
    t.ctx.highlight.set('host', t.ctx.refs.toLoci('abl1a'));
    expect(await run(t.ctx, 'highlight', { channel: 'links:d1', target: T315, closure: 'containing' })).toBe(true);
    const ch = t.ctx.highlight.channels.value;
    expect([...ch.keys()]).toEqual(['host', 'links:d1']);
    // T315 through the unit's auth numbering: the fixture's chain A starts later than 314.
    const [gk] = indices(t.ctx.refs.toLoci('1iep/A', T315.args));
    expect(ch.get('links:d1')).toMatchObject({ style: 'highlight', closure: 'containing' });
    expect(indices(ch.get('links:d1')?.target)).toEqual([gk]);
    await run(t.ctx, 'highlight', {
      channel: 'links:d2',
      target: { refs: ['1iep/A'], args: { start: 3, stop: 6 } },
      style: 'band',
    });
    expect(t.ctx.highlight.channels.value.get('links:d2')).toMatchObject({ style: 'band', closure: 'element' });
    await run(t.ctx, 'highlight', { channel: 'links:d1', target: null });
    await run(t.ctx, 'highlight', { channel: 'links:d2', target: null });
    expect([...t.ctx.highlight.channels.value.keys()]).toEqual(['host']);
    expect(await run(t.ctx, 'highlight', { channel: 'links:u', target: { refs: ['nothing-here'] } })).toBe(true);
    expect([[...t.ctx.highlight.channels.value.keys()], issues]).toEqual([['host'], []]);
  });

  it("range args outside the ref's domain: nothing done, `seqv.target.range` naming it, run resolves false (GS-12)", async () => {
    const t = await driven(case10());
    const issues = issuesOf(t.ctx);
    const unit = t.ctx.refs.toLoci('1iep/A') as SeqLoci;
    const [lo, hi] = [OrderedSet.min(unit.elements[0].indices), OrderedSet.max(unit.elements[0].indices) + 1];
    // As a document's focus is refused (`svs.focus.range`), never clipped to the unit.
    expect(t.ctx.refs.toLoci('1iep/A', { start: lo + 2, stop: hi + 100 })).toBeNull();
    expect(t.ctx.refs.toLoci('1iep/A', { start: lo + 2, stop: lo + 5 })).not.toBeNull();
    const target = { refs: ['1iep/A'], args: { start: lo + 2, stop: hi + 100 } };
    expect(await run(t.ctx, 'select', { target })).toBe(false);
    expect(t.ctx.selection.current.value.seq).toBeNull();
    expect(await run(t.ctx, 'highlight', { channel: 'links:x', target })).toBe(false);
    expect(t.ctx.highlight.channels.value.has('links:x')).toBe(false);
    expect(issues.map((i) => [i.code, i.severity])).toEqual([
      ['seqv.target.range', 'error'],
      ['seqv.target.range', 'error'],
    ]);
    expect(issues[0].message).toContain(`(${lo}..${hi})`);
  });

  it('label args outside a unit, and a range no complex member holds, are refused too, never clipped (LB-26)', async () => {
    const b = createSeqViewSpec();
    const p = b.data.sequence({ ref: 'p', value: 'MKVLAGWTQRSTEDNLKAAA' });
    p.numbering({
      name: 'auth',
      method: 'labels',
      labels: Array.from({ length: 20 }, (_, i) => String(101 + i)),
    } as never);
    b.data.unit({ ref: 'u', sequence: 'p', range: [5, 15] } as never);
    b.data.complex({ ref: 'c' }).unit({ ref: 'cu', sequence: 'p', range: [0, 10] } as never);
    b.view({ ref: 'v', label: 'V' })
      .block({ unit: 'u' } as never)
      .track({ label: 'L' })
      .letters();
    const t = await driven(b.getDocument());
    const span = (ref: string, args: Record<string, unknown>) => {
      const why: string[] = [];
      const l = t.ctx.refs.toLoci(ref, args as never, (m) => why.push(m)) as SeqLoci | null;
      const e = l?.elements.map((x) => [OrderedSet.min(x.indices), OrderedSet.max(x.indices) + 1]) ?? null;
      return [e, why.length];
    };
    // Labels 107..110 are indices [6, 10) inside u [5, 15); 103..110 start before it, 101..103 lie before it.
    expect(span('u', { beg: '107', end: '110', numbering: 'auth' })).toEqual([[[6, 10]], 0]);
    expect(span('u', { beg: '103', end: '110', numbering: 'auth' })).toEqual([null, 1]);
    expect(span('u', { beg: '101', end: '103', numbering: 'auth' })).toEqual([null, 1]);
    expect(span('u', { beg: '999', numbering: 'auth' })).toEqual([null, 1]);
    // A complex: its members that hold the range, none: refused as its member alone is.
    expect(span('c', { start: 5, stop: 30 })).toEqual([null, 1]);
    expect(span('c', { start: 2, stop: 6 })).toEqual([[[2, 6]], 0]);
    const issues = issuesOf(t.ctx);
    expect(
      await run(t.ctx, 'select', { target: { refs: ['u'], args: { beg: '103', end: '110', numbering: 'auth' } } }),
    ).toBe(false);
    expect(t.ctx.selection.current.value.seq).toBeNull();
    expect(issues.map((i) => i.code)).toEqual(['seqv.target.range']);
    expect(issues[0].message).toContain('(labels 103..110)');
  });

  it('refs of several kinds: each kind selected, lit on a channel of its own, focused together (VS-08)', async () => {
    const t = await driven(case3());
    for (const refs of [
      ['sh2', 'sites'],
      ['sites', 'sh2'],
    ]) {
      expect(await run(t.ctx, 'select', { target: { refs } })).toBe(true);
      const s = t.ctx.selection.current.value;
      expect([s.columns !== null, s.items !== null]).toEqual([true, true]);
      expect(await run(t.ctx, 'highlight', { channel: 'links:two', target: { refs } })).toBe(true);
      const kinds = [...t.ctx.highlight.channels.value]
        .filter(([c]) => c.startsWith('links:two'))
        .map(([, e]) => (e.target as { kind: string }).kind)
        .sort();
      expect(kinds).toEqual(['column-loci', 'item-loci']);
    }
    // Cleared: its parts too.
    await run(t.ctx, 'highlight', { channel: 'links:two', target: null });
    expect([...t.ctx.highlight.channels.value.keys()].filter((c) => c.startsWith('links:two'))).toEqual([]);
    // One kind again: the other part's channel goes.
    await run(t.ctx, 'highlight', { channel: 'links:two', target: { refs: ['sh2', 'sites'] } });
    await run(t.ctx, 'highlight', { channel: 'links:two', target: { refs: ['sites'] } });
    expect([...t.ctx.highlight.channels.value.keys()].filter((c) => c.startsWith('links:two'))).toEqual(['links:two']);
    // Focused together: the hull of both kinds' projections, in either order (the sites alone frame less).
    const S = 'family.main';
    const focused = async (refs: string[]) => {
      expect(await run(t.ctx, 'focus', { target: { refs }, section: S, animate: false })).toBe(true);
      const e = t.ctx.focus.current.value.get(S);
      return [e?.block, e?.start, e?.stop] as number[];
    };
    const [a, b] = [await focused(['sh2']), await focused(['sites'])];
    expect(b).not.toEqual(a);
    const hull = [a[0], Math.min(a[1], b[1]), Math.max(a[2], b[2])];
    expect([await focused(['sites', 'sh2']), await focused(['sh2', 'sites'])]).toEqual([hull, hull]);
  });

  it("a command's part channels clear only their own: a host's `<channel>/1` stays (LB-09)", async () => {
    const t = await driven(case3());
    const channels = () => [...t.ctx.highlight.channels.value.keys()].filter((c) => c.startsWith('links:a')).sort();
    t.ctx.highlight.set('links:a/1', t.ctx.refs.toLoci('sites') as never);
    await run(t.ctx, 'highlight', { channel: 'links:a', target: { refs: ['sites'] } });
    expect(channels()).toEqual(['links:a', 'links:a/1']);
    // Two kinds, then one: the command's own second part goes, then the host's is still there once cleared.
    t.ctx.highlight.clear('links:a/1');
    await run(t.ctx, 'highlight', { channel: 'links:a', target: { refs: ['sh2', 'sites'] } });
    await run(t.ctx, 'highlight', { channel: 'links:a', target: { refs: ['sites'] } });
    expect(channels()).toEqual(['links:a']);
    t.ctx.highlight.set('links:a/1', t.ctx.refs.toLoci('sites') as never);
    await run(t.ctx, 'highlight', { channel: 'links:a', target: null });
    expect(channels()).toEqual(['links:a/1']);
  });

  it('focus: T315 centered in its section; animate: false at once; `section` narrows it', async () => {
    const t = await driven(case10());
    const changes: FocusChange[] = [];
    t.ctx.focus.changed.subscribe((c) => changes.push({ ...c }));
    const [gk] = indices(t.ctx.refs.toLoci('1iep/A', T315.args));
    expect(await run(t.ctx, 'focus', { target: T315 })).toBe(true);
    t.advance(1000);
    const vp = t.vp('chains.main').live.current;
    const center = (vp.scrollX + (vp.width - (t.vp('chains.main').axis?.insetX ?? 0)) / 2) / vp.cell;
    expect(Math.floor(center)).toBe(gk);
    expect(t.ctx.focus.current.value.get('chains.main')).toMatchObject({ start: gk, stop: gk + 1 });
    // The same by indices, run directly, at once and in one section.
    changes.length = 0;
    const byIndex = { refs: ['1iep/A'], args: { start: gk, stop: gk + 1 } };
    await run(t.ctx, 'focus', { target: byIndex, animate: false, section: 'chains.main' });
    expect(changes.map((c) => [c.section, c.entry?.start, c.transitionMs, c.origin])).toEqual([
      ['chains.main', gk, 0, 'host:links'],
    ]);
  });

  it('set-snapshot by ref, next with animate: false, prev and the base; set-view', async () => {
    const t = await driven(case10());
    const S = t.ctx.snapshots;
    const calls: unknown[] = [];
    const [apply, next] = [S.apply.bind(S), S.next.bind(S)];
    Object.assign(S, {
      apply: (ref: string | null, o?: object) => {
        calls.push(['apply', ref, o]);
        return apply(ref, o);
      },
      next: (dir: 1 | -1, o?: object) => {
        calls.push(['next', dir, o]);
        return next(dir, o);
      },
    });
    const current: unknown[] = [];
    for (const params of [{ ref: 's2' }, { ref: 'next', animate: false }, { ref: 'prev' }, { ref: null }]) {
      await run(t.ctx, 'setSnapshot', params);
      current.push(S.current.value);
    }
    expect(current).toEqual(['s2', 's3', 's2', null]);
    // `next` and `prev` pass `animate` on as a ref does (the snapshot manager's test reads its effect).
    expect(calls).toEqual([
      ['apply', 's2', { animate: undefined }],
      ['next', 1, { animate: false }],
      ['apply', 's3', { animate: false }],
      ['next', -1, { animate: undefined }],
      ['apply', 's2', { animate: undefined }],
      ['apply', null, { animate: undefined }],
    ]);
    await run(t.ctx, 'setView', { ref: 'alignment' });
    expect(t.ctx.views.current.value).toBe('alignment');
  });

  it("select: the verbs with the caller's origin; null clears", async () => {
    const t = await driven(case10());
    const sel = () => t.ctx.selection.current.value;
    const range = (start: number, stop: number) => ({ refs: ['1iep/A'], args: { start, stop } });
    await run(t.ctx, 'select', { target: range(3, 6) });
    expect([indices(sel().seq), sel().origin]).toEqual([[3, 4, 5], 'host:links']);
    await run(t.ctx, 'select', { target: range(8, 9), mode: 'add' });
    await run(t.ctx, 'select', { target: range(4, 5), mode: 'remove' });
    expect(indices(sel().seq)).toEqual([3, 5, 8]);
    await run(t.ctx, 'select', { target: range(8, 10), mode: 'toggle' });
    expect(indices(sel().seq)).toEqual([3, 5, 8, 9]);
    await t.ctx.commands.run({ command: 'seqv:select', params: { target: null } }, 'host:agent');
    expect([sel().seq, sel().origin]).toEqual([null, 'host:agent']);
  });

  it("run-action with the caller's origin", async () => {
    const t = await driven(case10());
    const patches: PatchEvent[] = [];
    t.ctx.events.patch.subscribe((e) => patches.push(e));
    await run(t.ctx, 'runAction', { name: 'track.hide', target: 'chains.contacts' });
    expect(patches.map((p) => [p.label, p.origin])).toEqual([['Hide track', 'host:links']]);
  });

  it('without the `commands` behavior the built-ins are unhandled: run resolves, nothing changes', async () => {
    const t = await driven(case10());
    t.ctx.behaviors.get('commands')?.remove();
    expect(commandsOf(t.ctx).subscriberCount()).toBe(0);
    expect(await run(t.ctx, 'highlight', { channel: 'links:s', target: T315 })).toBe(true);
    await run(t.ctx, 'select', { target: T315 });
    await run(t.ctx, 'setView', { ref: 'alignment' });
    await run(t.ctx, 'setSnapshot', { ref: 's2' });
    expect([
      t.ctx.highlight.channels.value.size,
      t.ctx.selection.current.value.seq,
      t.ctx.views.current.value,
      t.ctx.snapshots.current.value,
    ]).toEqual([0, null, 'chains', null]);
  });
});

describe("closure before the verb, and Coordinates through ctx.systems (IN3's items)", () => {
  /** Case 14's gene placed on a chromosome: gene 77 closes to residue 11 (VW5). */
  const geneOnChr = () => shown(gene({ chromosome: true }));
  const run = (ctx: SeqViewerContext, command: string, params: Record<string, unknown>) =>
    ctx.commands.run({ command: `seqv:${command}`, params }, 'api');
  const nt77 = { refs: ['gene/unit'], args: { start: 77, stop: 78 } };

  it('a Select invocation closes its target before the verb; at element it stays the nucleotide', async () => {
    const { ctx } = await geneOnChr();
    await run(ctx, 'select', { target: nt77 });
    expect(keyed(ctx.selection.current.value.seq)).toEqual([['gene/unit', 77]]);
    await run(ctx, 'select', { target: nt77, closure: 'containing' });
    expect(keyed(ctx.selection.current.value.seq)).toEqual([
      ['gene/unit', 77],
      ['protein/unit', 11],
    ]);
  });

  it('Coordinates as data: held by key on a highlight, converted for select (closed) and focus', async () => {
    const { ctx } = await geneOnChr();
    const coordinates = { systems: ['seqstar:gene/unit', 'nothing'], runs: [0, 77, 78, 1, 0, 5] };
    await run(ctx, 'highlight', { channel: 'host:3d', target: { coordinates }, closure: 'containing' });
    const held = ctx.highlight.channels.value.get('host:3d')?.target;
    expect(held && isCoordinates(held) && [held.systems, Array.from(held.runs)]).toEqual([
      coordinates.systems,
      coordinates.runs,
    ]);
    expect(held && isCoordinates(held) && held.runs).toBeInstanceOf(Int32Array);
    await run(ctx, 'select', { target: { coordinates }, closure: 'containing' });
    expect(keyed(ctx.selection.current.value.seq)).toEqual([
      ['gene/unit', 77],
      ['protein/unit', 11],
    ]);
    const changes: FocusChange[] = [];
    ctx.focus.changed.subscribe((c) => changes.push({ ...c }));
    await run(ctx, 'focus', { target: { coordinates: { systems: ['seqstar:protein/unit'], runs: [0, 11, 12] } } });
    expect(changes.length).toBeGreaterThan(0);
    await run(ctx, 'select', { target: { coordinates: { systems: ['nothing'], runs: [0, 0, 1] } } });
    expect(ctx.selection.current.value.seq).toBe(null);
  });
});

describe("an application's command, run by invocation", () => {
  it('lab:pulse sets only links:<link> and its leave clears it; removed behaviors leave no subscription; undo: false', async () => {
    const t = await driven(case10());
    const Pulse = SeqCommand<{ link: string; refs: string[] | null }>(
      'lab:pulse',
      SimpleParamsSchema({ link: RequiredField(str, ''), refs: RequiredField(nullable(list(str)), '') }),
    );
    const lab: SeqViewerRegistryEntry = { commands: [Pulse] };
    const PulseHandler = defineBehavior('lab:pulse', 'commands', 0, {}, '', ({ ctx, add }) => {
      add(
        ctx.commands.subscribe(Pulse, (p) => {
          const loci = p.refs ? ctx.refs.toLoci(p.refs[0]) : null;
          ctx.highlight.set(`links:${p.link}`, loci, { style: 'accent' });
        }),
      );
    });
    const clicks: number[] = [];
    const LogClicks = defineBehavior('lab:log-clicks', 'interaction', 0, {}, '', ({ ctx, on }) =>
      on(ctx.input.events.click, (e) => clicks.push(e.hit.track)),
    );
    const base = [subscriptions(t.ctx), commandsOf(t.ctx).subscriberCount()];
    const undo = t.ctx.register(lab);
    const remove = t.ctx.behaviors.add(LogClicks);
    const removeHandler = t.ctx.behaviors.add(PulseHandler);
    const pulse = (link: string, refs: string[] | null) =>
      t.ctx.commands.run({ command: 'lab:pulse', params: { link, refs } }, 'host:links');
    expect(await pulse('x', ['site.gatekeeper'])).toBe(true);
    expect([...t.ctx.highlight.channels.value]).toMatchObject([['links:x', { style: 'accent' }]]);
    await pulse('x', null);
    expect(t.ctx.highlight.channels.value.size).toBe(0);
    t.click('chains.main', 500, 40);
    expect(clicks).toHaveLength(1);
    remove();
    removeHandler();
    expect([subscriptions(t.ctx), commandsOf(t.ctx).subscriberCount()]).toEqual(base);
    undo();
    expect(await pulse('y', ['site.gatekeeper'])).toBe(false); // inert once undone
  });
});

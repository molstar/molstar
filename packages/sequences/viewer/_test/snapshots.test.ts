/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The snapshot layer and the history (seqspec:viewer#the-document, seqspec:viewer#snapshots-and-playback,
// seqspec:viewer#history; VS2): routing under a snapshot and snapshot-tagged local patches, entering a snapshot's view,
// the history across a host's echo, local entries by identity, leaving a snapshot once, `next` with options, the player
// after dispose, committed patches' `base` and actions' origin; the player on the manual clock (ported).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { Patch, SeqViewSpecDoc } from '@molstar/seq-core/format';
import { describe, expect, it } from 'vitest';
import {
  type DocumentPort,
  type FocusChange,
  manualScheduler,
  type PatchEvent,
  type Routed,
  type SeqIssue,
  type SeqViewerContext,
} from '../src/index.js';
import { patch, smallDoc, viewer } from './graph-fixtures.js';

const set = (ref: string, params: Record<string, unknown>) => patch([{ op: 'set', ref, params }]);
const node = (doc: SeqViewSpecDoc | null, ref: string) => {
  const find = (n: SeqViewSpecDoc['root']): SeqViewSpecDoc['root'] | undefined =>
    n.ref === ref ? n : n.children?.map(find).find((x) => x);
  return doc ? find(doc.root) : undefined;
};
const hiddenIn = (ctx: SeqViewerContext, ref: string) => node(ctx.document.shown.value, ref)?.params?.hidden === true;
function issuesOf(ctx: SeqViewerContext) {
  const out: SeqIssue[] = [];
  ctx.events.issue.subscribe((i) => out.push(i));
  return out;
}
const flush = async () => {
  for (let i = 0; i < 5; i++) await new Promise((r) => setTimeout(r, 0));
};

/**
 * `a` and `dom`, view `v` with a focus [2, 12) on `a` and view `w`. `s1` moves the focus to [10, 14) over 400 ms;
 * `s2` labels `v.dom`; `s3` hides `v.dom`; `ins1` and `ins2` insert the track `v.extra`; `back` sets the default view
 * to `v` (as the base does), `first` to null.
 */
function storyDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  b.data.sequence({ ref: 'a', value: 'MKVLATGHEKMKVLATGHEK' });
  b.data.annotation({ ref: 'dom', geometry: 'ranges', target: 'a', value: { start: [0, 10], stop: [4, 16] } });
  const v = b.view({ ref: 'v', label: 'V' });
  v.focus({ ref: 'v.f', target: 'a', start: 2, stop: 12 });
  v.track({ ref: 'v.letters' }).letters({ ref: 'v.letters.l' });
  v.track({ ref: 'v.dom' }).ranges({ ref: 'v.dom.l', data: 'dom' });
  b.view({ ref: 'w', label: 'W' }).track({ ref: 'w.letters' }).letters();
  b.snapshot({ ref: 's1', transition_ms: 400, duration_ms: 1000, ops: [set('v.f', { start: 10, stop: 14 }).ops[0]] });
  b.snapshot({ ref: 's2', duration_ms: 500, ops: [set('v.dom', { label: 'S2' }).ops[0]] });
  b.snapshot({ ref: 's3', ops: [set('v.dom', { hidden: true }).ops[0]] });
  const extra = {
    kind: 'track',
    ref: 'v.extra',
    children: [{ kind: 'layer', ref: 'v.extra.l', params: { type: 'ranges', data: 'dom' } }],
  };
  for (const ref of ['ins1', 'ins2'])
    b.snapshot({ ref, ops: [{ op: 'insert', parent: 'v.main', node: extra as SeqViewSpecDoc['root'] }] });
  b.snapshot({ ref: 'back', ops: [{ op: 'set', ref: 'root', params: { default_view: 'v' } }] });
  b.snapshot({ ref: 'first', ops: [{ op: 'set', ref: 'root', params: { default_view: null } }] });
  const doc = b.getDocument();
  doc.root.params = { ...doc.root.params, default_view: 'v' };
  return doc;
}
async function open(ports = {}) {
  const scheduler = manualScheduler();
  const v = viewer({ scheduler, ...ports });
  const r = await v.ctx.setDocument(storyDoc());
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  return { ...v, scheduler };
}

describe('routing under a snapshot', () => {
  it('Hide of a snapshot-inserted track goes local, tagged: hidden on return, not over another snapshot', async () => {
    const { ctx } = await open();
    await ctx.snapshots.apply('ins1');
    expect(ctx.actions.run('track.hide', 'v.extra')).toBe(true);
    await ctx.graph.update();
    expect(ctx.document.local.value.map((l) => l.snapshot)).toEqual(['ins1']);
    expect(node(ctx.document.base.value, 'v.extra')).toBeUndefined();
    expect(hiddenIn(ctx, 'v.extra')).toBe(true);
    await ctx.snapshots.apply('ins2');
    expect([node(ctx.document.shown.value, 'v.extra') !== undefined, hiddenIn(ctx, 'v.extra')]).toEqual([true, false]);
    await ctx.snapshots.apply(null);
    await ctx.snapshots.apply('ins1');
    expect(hiddenIn(ctx, 'v.extra')).toBe(true);
    expect(ctx.document.local.value.length).toBe(1);
  });

  it("a host's port may commit under a snapshot; masked patches go local with the default port", async () => {
    const sent: Patch[] = [];
    const port: DocumentPort = {
      route: () => 'commit',
      commit: async (p) => {
        sent.push(p);
        return { ok: true, issues: [] };
      },
      history: false,
    };
    const host = await open({ document: port });
    await host.ctx.snapshots.apply('s2');
    await host.ctx.document.propose(set('v.letters', { hidden: true }), 'Hide', 'pointer');
    expect(sent.map((p) => p.ops)).toEqual([set('v.letters', { hidden: true }).ops]);

    const { ctx } = await open();
    const routes: Routed[] = [];
    ctx.events.patch.subscribe((e) => routes.push(e.routed));
    await ctx.applyPatch(patch([{ op: 'insert', parent: 'v.main', node: { kind: 'track', ref: 'v.mine' } }]), {
      layer: 'local',
    });
    await ctx.applyPatch(set('v.dom', { label: 'X' }), { layer: 'local' });
    await ctx.document.propose(set('v.mine', { hidden: true }), 'Hide', 'pointer'); // a node a local patch adds
    await ctx.document.propose(set('v.dom', { label: 'Y' }), 'Label', 'pointer'); // a field a local patch sets
    await ctx.document.propose(set('v.dom', { hidden: true }), 'Hide', 'pointer'); // another field: committed
    expect(routes).toEqual(['local', 'local', 'commit']);
    expect(node(ctx.document.shown.value, 'v.dom')?.params).toMatchObject({ label: 'Y', hidden: true });
  });

  it('an untagged local patch that does not apply over a snapshot is skipped and kept, reported once', async () => {
    const { ctx } = await open();
    const doc = storyDoc();
    doc.snapshots?.push({ ref: 'drop', ops: [{ op: 'remove', ref: 'v.dom' }] });
    await ctx.setDocument(doc);
    const issues = issuesOf(ctx);
    await ctx.applyPatch(set('v.dom', { label: 'mine' }), { layer: 'local', snapshot: null });
    for (const s of ['drop', 's2', 'drop', null]) await ctx.snapshots.apply(s);
    expect(issues.map((i) => i.code)).toEqual(['seqv.patch.skipped', 'seqv.patch.skipped']);
    expect(ctx.document.local.value.length).toBe(1);
    expect(node(ctx.document.shown.value, 'v.dom')?.params?.label).toBe('mine');
  });
});

describe("entering a snapshot's view", () => {
  it("switches back even when its default_view equals the base's; null is the first view", async () => {
    const { ctx } = await open();
    ctx.views.setCurrent('w');
    await ctx.snapshots.apply('back');
    expect(ctx.views.current.value).toBe('v');
    ctx.views.setCurrent('w');
    await ctx.snapshots.apply('first');
    expect(ctx.views.current.value).toBe('v');
    ctx.views.setCurrent('w');
    await ctx.snapshots.apply('s2'); // sets no view
    expect(ctx.views.current.value).toBe('w');
  });

  it("back from one that set a view to one that sets none returns to the base's view", async () => {
    const doc = storyDoc();
    doc.snapshots?.push({ ref: 'to-w', ops: [{ op: 'set', ref: 'root', params: { default_view: 'w' } }] });
    const { ctx } = viewer({ scheduler: manualScheduler() });
    await ctx.setDocument(doc);
    await ctx.snapshots.apply('s2');
    await ctx.snapshots.apply('to-w');
    expect(ctx.views.current.value).toBe('w');
    await ctx.snapshots.apply('s2');
    expect(ctx.views.current.value).toBe('v');
    await ctx.snapshots.apply('to-w');
    await ctx.snapshots.apply(null);
    expect(ctx.views.current.value).toBe('v');
    // A view the reader picked under it stays.
    await ctx.snapshots.apply('to-w');
    ctx.views.setCurrent('v');
    await ctx.snapshots.apply('s2');
    expect(ctx.views.current.value).toBe('v');
    ctx.views.setCurrent('w');
    await ctx.snapshots.apply('s3');
    expect(ctx.views.current.value).toBe('w');
  });
});

describe("the history across a host's echo", () => {
  function echoing() {
    const sent: Patch[] = [];
    const port: DocumentPort = {
      route: () => 'commit',
      commit: async (p) => {
        sent.push(p);
        return { ok: true, issues: [] };
      },
      history: true,
    };
    return { sent, port };
  }
  const echo = (ctx: SeqViewerContext, patches: readonly Patch[], extra?: Patch) => {
    let doc = ctx.document.base.value as SeqViewSpecDoc;
    for (const p of extra ? [...patches, extra] : patches) {
      const n = node(doc, p.ops[0].op === 'set' ? p.ops[0].ref : '') as { params?: object };
      doc = structuredClone(doc);
      Object.assign(node(doc, (p.ops[0] as { ref: string }).ref) as object, {
        params: { ...n.params, ...(p.ops[0] as { params: object }).params },
      });
    }
    return ctx.setDocument(doc);
  };

  it('two commits then one echo keep the history; undo before the echo is refused', async () => {
    const { sent, port } = echoing();
    const { ctx } = await open({ document: port });
    const issues = issuesOf(ctx);
    await ctx.document.propose(set('v.dom', { hidden: true }), 'Hide', 'pointer');
    await ctx.document.propose(set('v.letters', { label: 'L' }), 'Label', 'pointer');
    expect(await ctx.history.undo()).toBe(false);
    expect(issues.map((i) => i.code)).toEqual(['seqv.history.refused']);
    await echo(ctx, sent);
    expect(ctx.history.canUndo.value).toBe(true);
    expect(await ctx.history.undo()).toBe(true);
    expect(sent.at(-1)?.ops).toEqual(set('v.letters', { label: null }).ops);
  });

  it('an echo the host wrote in another form (a table as rows, keys reordered) is still the echo (VS-26)', async () => {
    const { sent, port } = echoing();
    const { ctx } = await open({ document: port });
    await ctx.document.propose(set('v.dom', { hidden: true }), 'Hide', 'pointer');
    let doc = ctx.document.base.value as SeqViewSpecDoc;
    doc = structuredClone(doc);
    Object.assign(node(doc, 'v.dom') as object, {
      params: { hidden: true, ...(node(doc, 'v.dom') as { params: object }).params },
    });
    const dom = node(doc, 'dom') as unknown as { params: { value: { start: number[]; stop: number[] } } };
    const v = dom.params.value;
    dom.params = { ...dom.params, value: v.start.map((start, i) => ({ stop: v.stop[i], start })) as never };
    void sent;
    await ctx.setDocument(doc);
    expect(ctx.history.canUndo.value).toBe(true);
  });

  it('an echo merged with a foreign edit clears it; so does a foreign setDocument with the default port', async () => {
    const { sent, port } = echoing();
    const { ctx } = await open({ document: port });
    await ctx.document.propose(set('v.dom', { hidden: true }), 'Hide', 'pointer');
    await echo(ctx, sent, set('v.letters', { label: 'theirs' }));
    expect([ctx.history.canUndo.value, ctx.history.canRedo.value]).toEqual([false, false]);

    const own = await open();
    await own.ctx.document.propose(set('v.dom', { hidden: true }), 'Hide', 'pointer');
    expect(own.ctx.history.canUndo.value).toBe(true);
    await own.ctx.setDocument(smallDoc());
    expect(own.ctx.history.canUndo.value).toBe(false);
  });
});

describe('local entries', () => {
  it('Hide (local), Reset view, Undo leaves local empty; an undo of a dropped patch is refused', async () => {
    const port: DocumentPort = { route: () => 'local', commit: async () => ({ ok: true, issues: [] }), history: false };
    const { ctx } = await open({ document: port });
    const issues = issuesOf(ctx);
    await ctx.document.propose(set('v.dom', { hidden: true }), 'Hide', 'pointer');
    ctx.document.clearLocal();
    expect(await ctx.history.undo()).toBe(false);
    expect(ctx.document.local.value).toEqual([]);

    await ctx.document.propose(set('v.dom', { hidden: true }), 'Hide', 'pointer');
    await ctx.applyPatch(patch([{ op: 'remove', ref: 'v.dom' }]));
    expect(ctx.document.local.value).toEqual([]);
    expect(await ctx.history.undo()).toBe(false);
    expect(ctx.document.local.value).toEqual([]);
    const own = issues.filter((i) => i.code.startsWith('seqv.'));
    expect(own.map((i) => i.code)).toEqual(['seqv.patch.dropped', 'seqv.history.refused']);
  });
});

describe('snapshots', () => {
  it('one that leaves the base is left once, with one issue, across three base changes', async () => {
    const { ctx } = await open();
    const issues = issuesOf(ctx);
    await ctx.snapshots.apply('s3');
    expect(ctx.snapshots.current.value).toBe('s3');
    await ctx.applyPatch(patch([{ op: 'remove', ref: 'v.dom' }]));
    await ctx.applyPatch(set('v.letters', { label: 'L' }));
    await ctx.applyPatch(set('v.letters', { hidden: true }));
    expect([ctx.snapshots.current.value, ctx.document.snapshot.value]).toEqual([null, null]);
    expect(issues.filter((i) => i.code.startsWith('svs.snapshot') && i.message.includes('"s3"')).length).toBe(1);
    expect(ctx.snapshots.list.value.find((s) => s.ref === 's3')?.valid).toBe(false);
  });

  it('a focus change carries transition_ms; animate: false enters at once, by apply and by next', async () => {
    const { ctx } = await open();
    const changes: FocusChange[] = [];
    ctx.focus.changed.subscribe((c) => changes.push({ ...c }));
    await ctx.snapshots.apply('s1');
    await ctx.snapshots.apply(null, { animate: false });
    await ctx.snapshots.next(1, { animate: false });
    expect(ctx.snapshots.current.value).toBe('s1');
    expect(changes.map((c) => [c.entry?.start, c.transitionMs, c.origin])).toEqual([
      [10, 400, 'snapshot'],
      [2, 0, 'snapshot'],
      [10, 0, 'snapshot'],
    ]);
  });

  it('the player stays transition + duration on each, on the manual clock, and skips an invalid snapshot', async () => {
    const { ctx, scheduler } = await open();
    const issues = issuesOf(ctx);
    // Only s1 to s3 play here; the others are removed. Then a base patch after which s3 no longer applies.
    const doc = storyDoc();
    doc.snapshots = doc.snapshots?.slice(0, 3);
    await ctx.setDocument(doc);
    await ctx.applyPatch(patch([{ op: 'remove', ref: 'v.dom' }]));
    expect(ctx.snapshots.list.value.map((s) => [s.ref, s.valid])).toEqual([
      ['s1', true],
      ['s2', false],
      ['s3', false],
    ]);
    await ctx.snapshots.apply('s3');
    expect(ctx.snapshots.current.value).toBe(null);
    expect(issues.map((i) => i.code)).toContain('svs.snapshot.invalid');

    ctx.snapshots.play();
    await flush();
    expect([ctx.snapshots.current.value, ctx.snapshots.playing.value]).toEqual(['s1', true]);
    scheduler.advance(1399); // s1: 400 transition + 1000
    await flush();
    expect(ctx.snapshots.current.value).toBe('s1');
    scheduler.advance(1); // s2 and s3 are skipped: the end
    await flush();
    expect([ctx.snapshots.current.value, ctx.snapshots.playing.value]).toEqual(['s1', false]);
    ctx.snapshots.play({ loop: true });
    scheduler.advance(10_000);
    await flush();
    expect(ctx.snapshots.current.value).toBe('s1');
    ctx.snapshots.stop();
    scheduler.advance(10_000);
    await flush();
    expect(scheduler.pendingTimers).toBe(0);
  });

  it('disposing while a step loads leaves no timer', async () => {
    const { ctx, scheduler } = await open();
    ctx.snapshots.play({ loop: true });
    ctx.dispose(); // the first step's apply is still awaiting its graph update
    await flush();
    scheduler.advance(10_000);
    await flush();
    expect(scheduler.pendingTimers).toBe(0);
  });
});

describe('committed patches and actions', () => {
  it('a committed patch carries base, a local one none; actions run with their origin', async () => {
    const sent: Patch[] = [];
    let route: Routed = 'commit';
    const port: DocumentPort = {
      route: () => route,
      commit: async (p) => {
        sent.push(p);
        return { ok: true, issues: [] };
      },
      history: false,
    };
    const v = viewer({ document: port });
    await v.ctx.setDocument({ ...storyDoc(), revision: 'r7' });
    const events: PatchEvent[] = [];
    v.ctx.events.patch.subscribe((e) => events.push(e));
    expect(v.ctx.actions.run('track.hide', 'v.dom', 'keyboard')).toBe(true);
    route = 'local';
    expect(v.ctx.actions.run('track.hide', 'v.letters', 'pointer')).toBe(true);
    expect(v.ctx.actions.run('view.switch', 'w', 'keyboard')).toBe(true);
    await flush();
    expect(sent.map((p) => p.base)).toEqual(['r7']);
    expect(v.ctx.document.local.value.map((l) => l.patch.base)).toEqual([undefined]);
    expect(events.map((e) => [e.label, e.origin])).toEqual([
      ['Hide track', 'keyboard'],
      ['Hide track', 'pointer'],
    ]);
    expect([v.ctx.views.current.value, v.ctx.views.current.origin]).toEqual(['w', 'keyboard']);
    expect(v.ctx.actions.can('history.undo', null)).toBe(false);
  });
});

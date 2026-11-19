/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SpecModel: revisions are the truth, patches derived. Patches applied with the language and the checks of a change,
// a commit each, the recent roots by revision; previews as patches applied again when the root changes; story steps
// shown over the root and the preview; the history through Seq*'s port: one step a batch, an agent's turn as one step,
// a change outside it reported, two participants on one host stack undone in order, and no history at all; a stale
// agent patch refused, then rebased by the host's opt-in.

import { State } from '@molstar/core/reactive';
import { applyPatch, type Patch, patchToNode } from '@molstar/core/repr/patch';
import { hashNode, type Node, resolveAddress } from '@molstar/core/repr/tree';
import { HistoryModel, type HistoryPort } from '@molstar/ui/history/history-model';
import { describe, expect, it } from 'vitest';
import { editingBehavior } from '../behaviors/editing.js';
import { viewerConfig } from '../config.js';
import { type EditOutcome, SpecModel } from '../models/spec.js';
import { HistoryParticipants, hostHistory, type UiHistory, uiHistory } from '../ui-history.js';
import { proteinOnGene } from './fixture.js';

const Letters = ['views', 'main', 'protein', 'letters'];
const set = (at: readonly string[], key: string, value: unknown, base?: string): Patch => ({
  ...(base ? { base } : {}),
  ops: [{ op: 'set', at: [...at], key, value: value as never }],
});
const heightOf = (root: Node) => resolveAddress(root, Letters)?.params?.height;

/** The fixture with a story step: the letters track taller, the gene's letters hidden, a note track it adds. */
function withStep(): Node {
  const step: Patch = {
    ops: [
      { op: 'set', at: Letters, key: 'height', value: 40 },
      { op: 'set', at: ['views', 'main', 'gene', 'letters', 'l'], key: 'hidden', value: true },
      {
        op: 'insert',
        parent: ['views', 'main', 'protein'],
        place: { last: true },
        node: { kind: 'track', id: 'note', children: [{ kind: 'layer', id: 'n', params: { type: 'letters' } }] },
      },
    ],
  };
  const r = applyPatch(proteinOnGene(), {
    ops: [
      {
        op: 'insert',
        parent: [],
        place: { last: true },
        node: {
          kind: 'snapshots',
          id: 'snapshots',
          children: [{ kind: 'snapshot', id: 'zoom', params: { label: 'Zoomed', patch: patchToNode(step) } }],
        },
      },
    ],
  });
  if (!r.ok) throw new Error(r.issues.map((i) => i.code).join());
  return r.root;
}

const sync = (o: EditOutcome | Promise<EditOutcome>) => {
  if (o instanceof Promise) throw new Error('ran through a history');
  return o;
};

describe('SpecModel: revisions and patches', () => {
  it('holds the head commit, its revision as the root’s hash, and the static model of what is shown', () => {
    const spec = new SpecModel(proteinOnGene());
    expect(spec.revision.value).toBe(hashNode(spec.root.value));
    expect(spec.head.value.parents).toEqual([]);
    expect(spec.model.value.systems.has('/data/g')).toBe(true);
    // a whole root: an equal tree is no change; another is one commit
    expect(spec.apply(proteinOnGene())).toBe(false);
    const before = spec.revision.value;
    expect(spec.apply(proteinOnGene({ title: 'other' }))).toBe(true);
    expect(spec.head.value.parents).toEqual([before]);
    expect(spec.apply(proteinOnGene({ title: 'b' }), { against: before })).toBe('stale-revision');
  });

  it('applies a patch as one commit; refuses one the checks of a change refuse, and a stale one (EC-81)', () => {
    const spec = new SpecModel(proteinOnGene());
    const r0 = spec.head.value;
    const r = spec.applyPatch(set(Letters, 'height', 30));
    expect(r.ok && r.changed).toBe(true);
    expect(spec.head.value.parents).toEqual([r0.revision]);
    expect(heightOf(spec.root.value)).toBe(30);
    // a range past the protein: a `seq.lang` error the patch introduced
    const past = spec.applyPatch(set(['views', 'main', 'protein'], 'range', [0, 300]));
    expect(past.ok).toBe(false);
    expect(past.issues.map((i) => i.code)).toContain('seq.lang.section.range');
    // made against the first revision: stale
    const stale = spec.applyPatch(set(Letters, 'height', 50, r0.revision));
    expect(!stale.ok && stale.issues.map((i) => i.code)).toEqual(['repr.patch.stale']);
    expect(heightOf(spec.root.value)).toBe(30);
  });

  it('keeps the recent roots by revision (`revisions.keepRecent`), then asks the history', () => {
    const spec = new SpecModel(proteinOnGene(), {}, viewerConfig({ revisions: { keepRecent: 3 } }));
    const revisions = [spec.revision.value];
    for (const h of [11, 12, 13, 14]) {
      spec.applyPatch(set(Letters, 'height', h));
      revisions.push(spec.revision.value);
    }
    expect(revisions.map((r) => spec.rootOf(r) !== undefined)).toEqual([false, false, true, true, true]);
    expect(hashNode(spec.rootOf(revisions[4])!)).toBe(revisions[4]);
  });

  it('tells a commit’s changes before its root is set: the session follows a move', () => {
    const spec = new SpecModel(proteinOnGene());
    const seen: [readonly unknown[], string][] = [];
    spec.beforeCommit((c) => seen.push([c.moved, spec.revision.value]));
    const before = spec.revision.value;
    spec.applyPatch({
      ops: [
        {
          op: 'move',
          at: ['views', 'main', 'protein'],
          parent: ['views', 'main'],
          place: { after: null },
          rename: 'prot',
        },
      ],
    });
    expect(seen).toEqual([[[{ from: ['views', 'main', 'protein'], to: ['views', 'main', 'prot'] }], before]]);
  });

  it('recompiles what the changes reach: a view edit keeps every data entry of the static model', () => {
    const spec = new SpecModel(proteinOnGene());
    const model = spec.model.value;
    spec.applyPatch(set(['views', 'main', 'protein', 'letters', 'l'], 'hidden', true));
    expect(spec.model.value).not.toBe(model);
    for (const [k, v] of model.systems) expect(spec.model.value.systems.get(k)).toBe(v);
    for (const [k, v] of model.maps) expect(spec.model.value.maps.get(k)).toBe(v);
  });
});

const rename = (at: readonly string[], id: string): Patch => ({
  ops: [{ op: 'move', at: [...at], parent: at.slice(0, -1), place: { last: true }, rename: id }],
});
const ProteinSection = ['views', 'main', 'protein'];
const renamed = [{ from: ProteinSection, to: ['views', 'main', 'prot'] }];
const back = [{ from: ['views', 'main', 'prot'], to: ProteinSection }];

describe('SpecModel: the session follows every move of the head', () => {
  it('a refused edit rolled back: told the move, then its undoing; `changed` never told the commit', () => {
    const spec = new SpecModel(proteinOnGene());
    const head = spec.head.value;
    const seen: unknown[] = [];
    spec.beforeCommit((c) => seen.push(c.moved));
    const o = sync(spec.patch('Two', [rename(ProteinSection, 'prot'), set(['views', 'main', 'nothing'], 'height', 1)]));
    expect(o.ok).toBe(false);
    expect(spec.head.value).toBe(head);
    expect(seen).toEqual([renamed, back]);
    expect(spec.changed.value).toBe(head);
  });

  it('undo and redo of a rename tell the session the moves back and forth', async () => {
    const spec = new SpecModel(proteinOnGene());
    const history = spec.attachHistory(uiHistory()) as UiHistory;
    const seen: unknown[] = [];
    spec.beforeCommit((c) => seen.push(c.moved));
    await spec.patch('Rename', [rename(ProteinSection, 'prot'), set([], 'title', 'x')]);
    await history.stack.undo();
    await history.stack.redo();
    expect(seen).toEqual([renamed, [], back, renamed]);
    spec.retire();
  });

  it('a stale rename rebased tells the session its move', () => {
    const spec = new SpecModel(proteinOnGene(), {}, viewerConfig({ edit: { rebaseStale: true } }));
    const seen: unknown[] = [];
    spec.beforeCommit((c) => seen.push(c.moved));
    const base = spec.revision.value;
    spec.applyPatch(set([], 'title', 'meanwhile'));
    const r = spec.applyPatch({ base, ops: rename(ProteinSection, 'prot').ops });
    expect(r.ok && r.rebased).toBe(true);
    expect(r.ok && r.changes.moved).toEqual(renamed);
    expect(seen).toEqual([[], renamed]);
  });
});

describe('SpecModel: previews', () => {
  it('shows a patch over the root, which cancels to nothing or commits as one patch', () => {
    const spec = new SpecModel(proteinOnGene());
    const root = spec.root.value;
    expect(spec.setPreview(set(Letters, 'height', 50)).ok).toBe(true);
    expect(heightOf(spec.shown.value)).toBe(50);
    expect(spec.root.value).toBe(root);
    expect(spec.model.value.spec).toBe(spec.shown.value);
    spec.cancelPreview();
    expect(spec.shown.value).toBe(root);
    // a drag: each move replaces the last, made from the root
    spec.setPreview(set(Letters, 'height', 41));
    spec.setPreview(set(Letters, 'height', 42));
    const r = sync(spec.edit('Resize', (e) => e.commitPreview()));
    expect(r.ok && r.changed).toBe(true);
    expect(spec.preview.value).toBeUndefined();
    expect(heightOf(spec.root.value)).toBe(42);
    expect(spec.head.value.parents).toEqual([hashNode(root)]);
  });

  it('survives an agent’s patch mid-drag: applied again over the new root', () => {
    const spec = new SpecModel(proteinOnGene());
    spec.setPreview(set(Letters, 'height', 60));
    spec.applyPatch(set([], 'title', 'by an agent'));
    expect(spec.root.value.params?.title).toBe('by an agent');
    expect(heightOf(spec.shown.value)).toBe(60);
    expect(spec.shown.value.params?.title).toBe('by an agent');
    sync(spec.edit('Resize', (e) => e.commitPreview()));
    expect(heightOf(spec.root.value)).toBe(60);
    expect(spec.root.value.params?.title).toBe('by an agent');
  });

  it('is dropped with an issue when it no longer applies, or quietly (`preview.dropped`)', () => {
    for (const dropped of ['issue', 'silent'] as const) {
      const spec = new SpecModel(proteinOnGene(), {}, viewerConfig({ preview: { dropped } }));
      spec.setPreview(set(Letters, 'height', 60));
      // an agent removes the track being resized
      expect(spec.applyPatch({ ops: [{ op: 'remove', at: Letters }] }).ok).toBe(true);
      expect(spec.preview.value).toBeUndefined();
      expect(spec.shown.value).toBe(spec.root.value);
      expect(spec.issues.value.map((i) => i.code)).toEqual(dropped === 'issue' ? ['seq.viewer.preview.dropped'] : []);
      // the gesture ends in nothing
      expect(spec.commitPreview()).toBeUndefined();
    }
  });
});

describe('SpecModel: story steps', () => {
  it('are entered and left as navigation, shown over the root and the preview', () => {
    const spec = new SpecModel(withStep());
    const root = spec.root.value;
    expect(spec.enterStep('nothing')?.code).toBe('seq.viewer.step.unknown');
    expect(spec.enterStep('zoom')).toBeUndefined();
    expect(spec.root.value).toBe(root);
    expect(spec.step.value).toMatchObject({ id: 'zoom', label: 'Zoomed', applies: true });
    expect(heightOf(spec.shown.value)).toBe(40);
    expect(resolveAddress(spec.shown.value, ['views', 'main', 'protein', 'note'])).toBeDefined();
    // a preview is an edit of the root: the step is applied over it
    spec.setPreview(set([], 'title', 'previewed'));
    expect(spec.shown.value.params?.title).toBe('previewed');
    expect(heightOf(spec.shown.value)).toBe(40);
    spec.cancelPreview();
    spec.leaveStep();
    expect(spec.step.value).toBeUndefined();
    expect(spec.shown.value).toBe(root);
  });

  it('edits while shown patch the root; the step applies again over it; a masked edit is reported', () => {
    const spec = new SpecModel(withStep());
    spec.enterStep('zoom');
    const r = spec.applyPatch(set(['views', 'main', 'gene', 'letters'], 'height', 12));
    expect(r.ok && r.changed).toBe(true);
    expect(resolveAddress(spec.root.value, ['views', 'main', 'gene', 'letters'])?.params?.height).toBe(12);
    expect(resolveAddress(spec.shown.value, ['views', 'main', 'gene', 'letters'])?.params?.height).toBe(12);
    // a key the step sets too: it lands in the root, and the step still shows its own
    const masked = spec.applyPatch(set(Letters, 'height', 25));
    expect(masked.ok && masked.issues.map((i) => i.code)).toEqual(['seq.viewer.edit.masked']);
    expect(heightOf(spec.root.value)).toBe(25);
    expect(heightOf(spec.shown.value)).toBe(40);
    // a node only the step makes is not in the root: refused, the hint naming the step
    const only = spec.applyPatch(set(['views', 'main', 'protein', 'note'], 'height', 5));
    expect(!only.ok && only.issues[0].code).toBe('repr.patch.target');
    expect(!only.ok && only.issues[0].hint).toContain('/snapshots/zoom');
  });

  it('is shown as not applying once an edit of the root breaks it', () => {
    const spec = new SpecModel(withStep());
    spec.enterStep('zoom');
    // a removal of what the step names is refused (its addresses are refs of the spec)
    const named = spec.applyPatch({ ops: [{ op: 'remove', at: ['views', 'main', 'gene', 'letters'] }] });
    expect(!named.ok && named.issues.map((i) => i.code)).toContain('repr.patch.dangling');
    // a track of the id the step inserts: the step's insert no longer applies
    const note: Node = {
      kind: 'track',
      id: 'note',
      children: [{ kind: 'layer', id: 'x', params: { type: 'letters' } }],
    };
    expect(
      spec.applyPatch({
        ops: [{ op: 'insert', parent: ['views', 'main', 'protein'], place: { last: true }, node: note }],
      }).ok,
    ).toBe(true);
    expect(spec.step.value).toMatchObject({ id: 'zoom', applies: false });
    // the ops that still apply are shown, flagged: its height, not its insert (the root's track stays)
    expect(heightOf(spec.shown.value)).toBe(40);
    expect(resolveAddress(spec.shown.value, ['views', 'main', 'protein', 'note'])).toBe(
      resolveAddress(spec.root.value, ['views', 'main', 'protein', 'note']),
    );
    expect(spec.step.value?.state.skipped).toBe(1);
    expect(spec.issues.value.map((i) => i.code)).toEqual(['seq.viewer.step.stale']);
  });
});

/** A stand-in structure state: a number a host's own command changes, a participant of the same stack. */
function structureState() {
  const state = new State(0);
  const port: HistoryPort<number> = {
    checkpoint: () => state.value,
    restore: async (c) => state.set(c),
  };
  return { state, port };
}

const person = { kind: 'person' };

describe('SpecModel: the history through Seq*’s port', () => {
  it('without a history: edits apply at once, one batch all or nothing, nothing undoes them', () => {
    const spec = new SpecModel(proteinOnGene());
    expect(spec.history).toBeUndefined();
    const before = spec.root.value;
    const o = sync(
      spec.patch('Two', [set(Letters, 'height', 30), set(['views', 'main', 'protein'], 'range', [0, 300])], {
        origin: person,
      }),
    );
    expect(o.ok).toBe(false);
    // the first applied, the second refused: the first taken back
    expect(spec.root.value).toBe(before);
    expect(sync(spec.patch('One', set(Letters, 'height', 30))).changed).toBe(true);
  });

  it('one step a batch; no step for a refusal or no change; undo restores roots, redo puts them back', async () => {
    const spec = new SpecModel(proteinOnGene());
    const history = spec.attachHistory(uiHistory()) as UiHistory;
    const r0 = spec.root.value;
    await spec.patch('Two', [set(Letters, 'height', 30), set([], 'title', 't')]);
    const r1 = spec.root.value;
    expect(history.stack.steps.value.map((s) => s.label)).toEqual(['Two']);
    const refused = await spec.patch('Too far', set(['views', 'main', 'protein'], 'range', [0, 300]));
    expect(refused.ok).toBe(false);
    expect((await spec.patch('Same', set(Letters, 'height', 30))).changed).toBe(false);
    expect(history.stack.steps.value).toHaveLength(1);
    // a preview ends in one step
    spec.setPreview(set(Letters, 'height', 55));
    await spec.edit('Resize', (e) => e.commitPreview());
    expect((await history.stack.undo()).ok).toBe(true);
    expect(spec.root.value).toBe(r1);
    expect((await history.stack.undo()).ok).toBe(true);
    expect(spec.root.value).toBe(r0);
    await history.stack.redo();
    expect(spec.root.value).toBe(r1);
    // the stack resolves a root it checkpointed
    expect(history.root?.(hashNode(r0))).toBe(r0);
    spec.retire();
  });

  it('a change outside a run is reported to the host: its steps can no longer be undone', async () => {
    const spec = new SpecModel(proteinOnGene());
    const history = spec.attachHistory(uiHistory()) as UiHistory;
    await spec.patch('One', set(Letters, 'height', 30));
    spec.applyPatch(set(Letters, 'height', 31));
    const undo = await history.stack.undo();
    expect(!undo.ok && undo.code).toBe('unsafe');
    expect(heightOf(spec.root.value)).toBe(31);
    spec.retire();
  });

  it('two consumers on one history: an agent’s turn as one step between the user’s edits, undone in order', async () => {
    const spec = new SpecModel(proteinOnGene());
    const history = spec.attachHistory(uiHistory()) as UiHistory;
    await spec.patch('User: taller', set(Letters, 'height', 30));
    const agent = { kind: 'agent', turn: 't1' };
    const end = spec.beginTurn('Agent: tidy', agent);
    await spec.patch('retitle', set([], 'title', 'tidy'), { origin: agent });
    await spec.patch('hide', set(['views', 'main', 'gene', 'letters', 'l'], 'hidden', true), { origin: agent });
    end();
    await spec.patch('User: shorter', set(Letters, 'height', 20));
    expect(history.stack.steps.value.map((s) => s.label)).toEqual(['User: shorter', 'Agent: tidy', 'User: taller']);
    await history.stack.undo();
    expect(heightOf(spec.root.value)).toBe(30);
    expect(spec.root.value.params?.title).toBe('tidy');
    await history.stack.undo();
    expect(spec.root.value.params?.title).toBe('p on g');
    expect(resolveAddress(spec.root.value, ['views', 'main', 'gene', 'letters', 'l'])?.params?.hidden).toBeUndefined();
    expect(heightOf(spec.root.value)).toBe(30);
    spec.retire();
  });

  it('a Seq* spec and a stand-in structure state on one host stack, undone in order', async () => {
    const participants = new HistoryParticipants();
    const stack = new HistoryModel(participants);
    const structure = structureState();
    participants.add('structure', structure.port);
    const spec = new SpecModel(proteinOnGene());
    spec.attachHistory(hostHistory(stack, participants));
    // the host's own step, Seq*'s, the host's again
    await stack.run('Color the structure', person, async () => structure.state.set(1));
    await spec.patch('Taller', set(Letters, 'height', 30));
    await stack.run('Recolor the structure', person, async () => structure.state.set(2));
    expect(stack.steps.value.map((s) => s.label)).toEqual(['Recolor the structure', 'Taller', 'Color the structure']);
    await stack.undo();
    expect([structure.state.value, heightOf(spec.root.value)]).toEqual([1, 30]);
    await stack.undo();
    expect([structure.state.value, heightOf(spec.root.value)]).toEqual([1, undefined]);
    await stack.undo();
    expect([structure.state.value, heightOf(spec.root.value)]).toEqual([0, undefined]);
    await stack.redo();
    await stack.redo();
    expect([structure.state.value, heightOf(spec.root.value)]).toEqual([1, 30]);
    // the viewer goes: the host's stack stays, without Seq*'s participant
    spec.retire();
    expect(participants.checkpoint()).toEqual({ structure: 1 });
    stack.dispose();
  });
});

describe('SpecModel: an edit inside a host’s run', () => {
  it('a host’s step whose work runs a Seq* command: one step, no deadlock', async () => {
    const participants = new HistoryParticipants();
    const stack = new HistoryModel(participants);
    const structure = structureState();
    participants.add('structure', structure.port);
    const spec = new SpecModel(proteinOnGene());
    spec.attachHistory(hostHistory(stack, participants));
    const editing = editingBehavior(spec);
    const hide = editing.commands.find((c) => c.id === 'layer.hide')!;
    const outer = stack.run('Hide and recolor', person, async (within) => {
      structure.state.set(5);
      return hide.run({}, { address: '/views/main/protein/letters/l', within }, new AbortController().signal);
    });
    const done = await Promise.race([outer.then(() => 'done'), new Promise((r) => setTimeout(() => r('stuck'), 200))]);
    expect(done).toBe('done');
    expect(stack.steps.value.map((x) => x.label)).toEqual(['Hide and recolor']);
    await stack.undo();
    expect(structure.state.value).toBe(0);
    expect(
      resolveAddress(spec.root.value, ['views', 'main', 'protein', 'letters', 'l'])?.params?.hidden,
    ).toBeUndefined();
    spec.retire();
    stack.dispose();
  });
});

describe('SpecModel: small contracts', () => {
  it('applies a preview once a move: what is shown is the preview’s own result', () => {
    const spec = new SpecModel(proteinOnGene());
    const r = spec.setPreview(set(Letters, 'height', 50));
    expect(r.ok && spec.edited.value).toBe(r.ok && r.root);
  });

  it('commits a preview in one batch: told the preview is gone, what is shown is the new root', () => {
    const spec = new SpecModel(proteinOnGene());
    spec.setPreview(set(Letters, 'height', 50));
    const seen: boolean[] = [];
    spec.preview.listen(() => seen.push(spec.shown.value === spec.root.value && heightOf(spec.root.value) === 50));
    spec.commitPreview();
    expect(seen).toEqual([true]);
  });

  it('a second history attached disposes the first', () => {
    const spec = new SpecModel(proteinOnGene());
    let disposed = 0;
    const port = () => ({
      run: async <T>(_: unknown, work: (c: { id: number }) => T) => work({ id: 1 }),
      isActive: () => true,
      outside: () => {},
      dispose: () => disposed++,
    });
    spec.attachHistory(port);
    spec.attachHistory(port);
    expect(disposed).toBe(1);
    spec.retire();
    expect(disposed).toBe(2);
  });

  it('two viewers’ specs on one host stack: a participant each by default', async () => {
    const participants = new HistoryParticipants();
    const stack = new HistoryModel(participants);
    const a = new SpecModel(proteinOnGene());
    const b = new SpecModel(proteinOnGene());
    a.attachHistory(hostHistory(stack, participants));
    b.attachHistory(hostHistory(stack, participants));
    await a.patch('A', set(Letters, 'height', 30));
    await b.patch('B', set(Letters, 'height', 40));
    await stack.undo();
    expect([heightOf(a.root.value), heightOf(b.root.value)]).toEqual([30, undefined]);
    a.retire();
    b.retire();
    stack.dispose();
  });

  it('names the step shown in a refusal’s hint only when the step makes the target', () => {
    const spec = new SpecModel(withStep());
    spec.enterStep('zoom');
    const none = spec.applyPatch(set(['views', 'main', 'nothing'], 'height', 5));
    expect(none.ok).toBe(false);
    expect(none.issues[0].hint ?? '').not.toContain('/snapshots/zoom');
  });
});

describe('SpecModel: a stale agent patch', () => {
  it('is refused, then rebased by the host’s opt-in through merge3 from the recent roots', async () => {
    const spec = new SpecModel(proteinOnGene());
    const base = spec.revision.value;
    // the agent reads the base revision; the user edits meanwhile
    spec.applyPatch(set(Letters, 'height', 30));
    const agentPatch = set([], 'title', 'by the agent', base);
    const refused = spec.applyPatch(agentPatch, { requireBase: true });
    expect(!refused.ok && refused.issues.map((i) => i.code)).toEqual(['repr.patch.stale']);
    const rebased = spec.applyPatch(agentPatch, { requireBase: true, rebase: true });
    expect(rebased.ok && rebased.rebased).toBe(true);
    expect(spec.root.value.params?.title).toBe('by the agent');
    expect(heightOf(spec.root.value)).toBe(30);
  });

  it('with conflicts: refused, with the conflicts and the partial merge; an unknown base cannot be rebased', () => {
    const spec = new SpecModel(proteinOnGene(), {}, viewerConfig({ edit: { rebaseStale: true } }));
    const base = spec.revision.value;
    spec.applyPatch(set(Letters, 'height', 30));
    const r = spec.applyPatch(set(Letters, 'height', 50, base));
    expect(r.ok).toBe(false);
    expect(!r.ok && r.conflicts?.map((c) => c.kind)).toEqual(['set']);
    expect(!r.ok && r.partial && heightOf(r.partial)).toBe(30);
    const unknown = spec.applyPatch(set(Letters, 'height', 50, 'f'.repeat(32)));
    expect(!unknown.ok && unknown.issues.map((i) => i.code)).toEqual(['seq.viewer.rebase.base']);
  });
});

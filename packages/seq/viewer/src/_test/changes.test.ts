/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Every change has a kind (_spec/README.md "Commands"): each command of the default viewer, run on a viewer with the
// default spec's history, changes exactly what its declared kind says (the spec and one history step for a patch, the
// preview, the session, what is shown without the spec, or nothing); and undo touches only spec patches.

import { hashNode, type Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it } from 'vitest';
import type { EditingBehavior } from '../behaviors/editing.js';
import { type ChangeKind, changeKindOf } from '../commands.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { type UiHistory, uiHistory } from '../ui-history.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Protein = '/views/main/protein';
const LettersLayer = `${Protein}/letters/l`;
const DomainsTrack = `${Protein}/domains`;

function spec(): Node {
  const b = createSeqSpec({ title: 'changes' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLLPENN' });
  const domains = b.data.ranges('domains', {
    target: p,
    rows: { start: new Int32Array([0, 15]), end: new Int32Array([10, 25]), id: ['A', 'B'] },
  });
  b.refs.ref('b', { to: `${domains}/(item :id B)`, label: 'Domain B' });
  b.views.scale('s');
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p });
    s.track('letters').layer('l', { type: 'letters' });
    s.track('domains').layer('d', { type: 'ranges', data: domains });
    // a second section, so the cursor has one to move to
    v.section('below', { of: p }).track('letters').layer('l', { type: 'letters' });
  });
  b.snapshots.snapshot('first', {
    label: 'Domains taller',
    patch: { ops: [{ op: 'set', at: ['views', 'main', 'protein', 'domains'], key: 'height', value: 40 }] },
  });
  return b.build();
}

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

/** The default viewer, headless but with the mounted viewer's history (`uiHistory`). */
async function make() {
  frames = installFrames();
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    history: uiHistory(),
    config: [[viewerItem('viewport.maxPxPerPosition'), 1000]],
  });
  viewers.push(v);
  await v.load(spec());
  v.resize(600);
  await v.settle();
  frames.flush();
  return v;
}

const stackOf = (v: SeqViewer) => (v.spec.history as UiHistory).stack;
const editing = (v: SeqViewer) => v.behaviors.get('editing') as EditingBehavior;

/** What a change can touch, read so that two reads compare. */
function look(v: SeqViewer) {
  const vp = v.view.get(Protein);
  return {
    patch: `${v.spec.revision.value} ${stackOf(v).steps.value.length}`,
    preview: v.spec.preview.value,
    session: [
      v.selection.value.value,
      v.selection.cursor.value,
      ...v.highlight.names.value.map((n) => v.highlight.get(n)),
      // the reader's color-blind substitution
      v.legends.cvdSafe.value,
    ],
    navigate: `${vp?.start} ${vp?.end} ${v.spec.step.value?.id}`,
  };
}

function changed(a: ReturnType<typeof look>, b: ReturnType<typeof look>): ChangeKind[] {
  const out: ChangeKind[] = [];
  if (a.patch !== b.patch) out.push('patch');
  if (a.preview !== b.preview) out.push('preview');
  if (a.session.length !== b.session.length || a.session.some((x, i) => x !== b.session[i])) out.push('session');
  if (a.navigate !== b.navigate) out.push('navigate');
  return out;
}

const run = (v: SeqViewer, id: string, ctx: unknown) => Promise.resolve(v.execute(id, ctx));

const cell = (position: number) => ({ section: Protein, hit: { layer: LettersLayer, system: '/data/p', position } });
const cursorAt = (section: string) => (v: SeqViewer) =>
  v.selection.setCursor({ section, system: '/data/p', position: 5 });

/** The domains colored by their id: a legend whose entries select items. */
const colorById = (v: SeqViewer) =>
  run(v, 'layer.channel', { address: `${DomainsTrack}/d`, value: { kind: 'field', params: { name: 'id' } } });

/** Each command's context, and what makes it change something first (a show needs a hide). */
const Cases: Record<string, { ctx: unknown; before?: (v: SeqViewer) => unknown }> = {
  'view.fit': { ctx: { section: Protein }, before: (v) => v.view.setRange(Protein, 5, 15) },
  'view.focus': { ctx: { section: Protein, focus: { start: 2, end: 8 } } },
  'view.back': {
    ctx: { section: Protein },
    before: (v) => run(v, 'view.focus', { section: Protein, focus: { start: 2, end: 8 } }),
  },
  'view.forward': {
    ctx: { section: Protein },
    before: async (v) => {
      await run(v, 'view.focus', { section: Protein, focus: { start: 2, end: 8 } });
      await run(v, 'view.back', { section: Protein });
    },
  },
  'view.zoom-in': { ctx: { section: Protein } },
  'view.zoom-out': { ctx: { section: Protein }, before: (v) => v.view.setRange(Protein, 5, 15) },
  'view.pan-left': { ctx: { section: Protein }, before: (v) => v.view.setRange(Protein, 10, 20) },
  'view.pan-right': { ctx: { section: Protein }, before: (v) => v.view.setRange(Protein, 10, 20) },
  'address.focus': { ctx: { address: '/refs/b' } },
  'address.reveal': {
    ctx: { address: '/data/p/(range :from 1 :to 5)' },
    before: (v) => v.view.setRange(Protein, 15, 25),
  },
  'address.select': { ctx: { address: '/refs/b' } },
  'address.highlight': { ctx: { address: '/refs/b', channel: 'paper' } },
  'address.hover': { ctx: { address: '/refs/b' } },
  'address.unhover': { ctx: {}, before: (v) => run(v, 'address.hover', { address: '/refs/b' }) },
  'address.copy-link': { ctx: {}, before: (v) => run(v, 'address.select', { address: '/refs/b' }) },
  'refs.add': { ctx: { ref: { id: 'here', to: '/data/p/(range :from 1 :to 5)' } } },
  'refs.remove': { ctx: { address: '/refs/b' } },
  'refs.label': { ctx: { address: '/refs/b', label: 'Second domain' } },
  'selection.bookmark': { ctx: {}, before: (v) => run(v, 'address.select', { address: '/refs/b' }) },
  'highlight.clear': {
    ctx: {},
    before: (v) => run(v, 'address.highlight', { address: '/refs/b', channel: 'paper' }),
  },
  'section.reverse': { ctx: { address: Protein } },
  'section.collapse': { ctx: { address: Protein } },
  'layer.hide': { ctx: { address: LettersLayer } },
  'layer.show': { ctx: { address: LettersLayer }, before: (v) => run(v, 'layer.hide', { address: LettersLayer }) },
  'layer.color': { ctx: { address: `${DomainsTrack}/d`, color: '#d62728' } },
  'layer.channel': { ctx: { address: `${DomainsTrack}/d`, channel: 'opacity', value: 0.5 } },
  'scale.set': { ctx: { address: '/views/s', scale: { domain: [0, 1] } } },
  'theme.value-color': { ctx: { address: `${DomainsTrack}/d`, entry: 'A', color: '#123456' } },
  'legend.select': { ctx: { address: `${DomainsTrack}/d`, entry: 'B' }, before: colorById },
  'legend.hover': { ctx: { address: `${DomainsTrack}/d`, entry: 'B' }, before: colorById },
  'theme.cvd-safe': { ctx: { on: true } },
  'layer.filter': { ctx: { address: `${DomainsTrack}/d`, entry: 'B' }, before: colorById },
  'track.hide': { ctx: { address: DomainsTrack } },
  'track.show': { ctx: { address: DomainsTrack }, before: (v) => run(v, 'track.hide', { address: DomainsTrack }) },
  'track.move': { ctx: { address: DomainsTrack, after: null } },
  'track.add': {
    ctx: {
      address: Protein,
      after: 'domains',
      node: { kind: 'track', id: 'domains-2', children: [{ kind: 'layer', id: 'd', params: { type: 'letters' } }] },
    },
  },
  'step.enter': { ctx: { step: 'first' } },
  'step.leave': { ctx: {}, before: (v) => run(v, 'step.enter', { step: 'first' }) },
  'selection.set': { ctx: cell(3) },
  'selection.extend': { ctx: cell(6), before: (v) => run(v, 'selection.set', cell(3)) },
  'selection.toggle': { ctx: cell(3) },
  'selection.add': { ctx: { range: { system: '/data/p', start: 8, end: 12 } } },
  'selection.subtract': {
    ctx: { address: '/data/p/(range :from 4 :to 6)' },
    before: (v) => run(v, 'selection.add', { range: { system: '/data/p', start: 0, end: 10 } }),
  },
  'selection.clear': { ctx: {}, before: (v) => run(v, 'selection.set', cell(3)) },
  'cursor.prev': { ctx: { section: Protein }, before: cursorAt(Protein) },
  'cursor.next': { ctx: { section: Protein }, before: cursorAt(Protein) },
  'cursor.up': { ctx: { section: Protein }, before: cursorAt('/views/main/below') },
  'cursor.down': { ctx: { section: Protein }, before: cursorAt(Protein) },
  'cursor.extend-prev': { ctx: { section: Protein }, before: cursorAt(Protein) },
  'cursor.extend-next': { ctx: { section: Protein }, before: cursorAt(Protein) },
  'cursor.add-prev': { ctx: { section: Protein }, before: cursorAt(Protein) },
  'cursor.add-next': { ctx: { section: Protein }, before: cursorAt(Protein) },
  'cursor.toggle': { ctx: { section: Protein }, before: cursorAt(Protein) },
  'cursor.focus': { ctx: { section: Protein }, before: cursorAt(Protein) },
  'preview.cancel': {
    ctx: {},
    before: (v) =>
      editing(v).preview({
        ops: [{ op: 'set', at: ['views', 'main', 'protein', 'letters'], key: 'height', value: 30 }],
      }),
  },
};

describe('every change has a kind', () => {
  it('each command of the default viewer changes exactly what its kind says', async () => {
    const ids = (await make()).commands.map((c) => c.id);
    expect(Object.keys(Cases).sort()).toEqual([...ids].sort());
    for (const id of ids) {
      const v = await make();
      const c = Cases[id];
      await c.before?.(v);
      frames!.flush();
      const before = look(v);
      const message = await run(v, id, c.ctx);
      frames!.flush();
      const kind = changeKindOf(id);
      expect(changed(before, look(v)), `${id}: ${message}`).toEqual(kind === 'none' ? [] : [kind]);
      if (kind === 'patch') expect(stackOf(v).steps.value.length, id).toBe(Number(before.patch.split(' ')[1]) + 1);
    }
  });

  it('commands issued together, without waiting, each apply: one step each, in order', async () => {
    const v = await make();
    const runs = [
      v.execute('layer.hide', { address: LettersLayer }),
      v.execute('track.move', { address: DomainsTrack, after: null }),
      v.execute('layer.color', { address: `${DomainsTrack}/d`, color: '#d62728' }),
      v.execute('refs.add', { ref: { id: 'here', to: '/data/p/(range :from 1 :to 5)' } }),
    ];
    expect(await Promise.all(runs)).toEqual(['Hidden', 'Moved', 'Colored', 'Bookmarked']);
    expect(stackOf(v).steps.value.map((s) => s.label)).toEqual([
      'Add a bookmark',
      'Color a layer',
      'Move a track',
      'Hide a layer',
    ]);
  });

  it('a track hidden and shown: its own `hidden`, its row gone from the view and back, its layers untouched', async () => {
    const v = await make();
    await run(v, 'layer.hide', { address: `${DomainsTrack}/d` });
    const root = v.spec.root.value;
    expect(await run(v, 'track.hide', { address: DomainsTrack })).toBe('Hidden');
    expect(v.view.resolved.value!.sections[0].tracks.map((t) => t.address)).toEqual([`${Protein}/letters`]);
    expect(await run(v, 'track.show', { address: DomainsTrack })).toBe('Shown');
    expect(v.spec.revision.value).toBe(hashNode(root));
    expect(v.view.resolved.value!.sections[0].tracks.map((t) => t.address)).toEqual([
      `${Protein}/letters`,
      DomainsTrack,
    ]);
  });

  it('undo touches only spec patches: the session, the viewports and the step shown stay', async () => {
    const v = await make();
    const root = v.spec.root.value;
    await run(v, 'layer.hide', { address: LettersLayer });
    await run(v, 'address.select', { address: '/refs/b' });
    await run(v, 'address.highlight', { address: '/refs/b', channel: 'paper' });
    v.view.setRange(Protein, 4, 14);
    await run(v, 'step.enter', { step: 'first' });
    const kept = look(v);
    expect((await stackOf(v).undo()).ok).toBe(true);
    expect(v.spec.root.value).toBe(root);
    const after = look(v);
    expect(after.session).toEqual(kept.session);
    expect(after.navigate).toBe(kept.navigate);
    expect(stackOf(v).steps.value).toEqual([]);
  });

  it('the selection bookmarked is one patch and one undo step: labeled by its address, undone and redone whole', async () => {
    const v = await make();
    const root = v.spec.root.value;
    expect(await run(v, 'selection.bookmark', {})).toBe('Nothing selected');
    await run(v, 'address.select', { address: '/data/p/(range :from 3 :to 8)' });
    expect(await run(v, 'selection.bookmark', {})).toBe('Bookmarked');
    const row = v.addresses.refs.value.find((r) => r.id === 'bookmark-1')!;
    expect(row).toMatchObject({ label: '/data/p/(range :from 3 :to 8)', to: ['/data/p/(range :from 3 :to 8)'] });
    // a second one takes the next id, a given label replaces the printed address
    expect(await run(v, 'selection.bookmark', { label: 'Loop' })).toBe('Bookmarked');
    expect(v.addresses.refs.value.find((r) => r.id === 'bookmark-2')?.label).toBe('Loop');
    expect(stackOf(v).steps.value.length).toBe(2);
    await stackOf(v).undo();
    await stackOf(v).undo();
    expect(v.spec.root.value).toBe(root);
    await stackOf(v).redo();
    expect(v.addresses.refs.value.map((r) => r.id)).toEqual(['b', 'bookmark-1']);
    // followed as a link in the same viewer: the ref resolves to what was selected
    expect(v.addresses.resolve('/refs/bookmark-1').status).toBe('resolved');
  });

  it('a bookmark is labeled by a patch and removed by one, each undone whole', async () => {
    const v = await make();
    const root = v.spec.root.value;
    expect(await run(v, 'refs.label', { address: '/refs/b', label: 'Second domain' })).toBe('Labeled');
    expect(v.addresses.refs.value[0].label).toBe('Second domain');
    expect(await run(v, 'refs.label', { address: '/refs/b', label: '' })).toBe('Labeled');
    expect(v.addresses.refs.value[0].label).toBe('b');
    expect(await run(v, 'refs.remove', { address: '/refs/b' })).toBe('Removed');
    expect(v.addresses.refs.value).toEqual([]);
    await stackOf(v).undo();
    expect(v.addresses.refs.value.map((r) => r.id)).toEqual(['b']);
    await stackOf(v).undo();
    await stackOf(v).undo();
    expect(v.spec.root.value).toBe(root);
    await stackOf(v).redo();
    await stackOf(v).redo();
    await stackOf(v).redo();
    expect(v.addresses.refs.value).toEqual([]);
    expect(await run(v, 'refs.remove', { address: '/refs/b' })).toBe('No ref there');
  });

  it('highlight.clear clears one channel, or all of them', async () => {
    const v = await make();
    await run(v, 'address.highlight', { address: '/refs/b', channel: 'paper' });
    await run(v, 'address.highlight', { address: '/refs/b', channel: 'notes' });
    expect(await run(v, 'highlight.clear', { channel: 'paper' })).toBe('Cleared');
    expect(v.highlight.names.value).toEqual(['notes']);
    expect(await run(v, 'highlight.clear', { channel: 'paper' })).toBe('Nothing highlighted');
    await run(v, 'address.highlight', { address: '/refs/b', channel: 'paper' });
    expect(await run(v, 'highlight.clear', {})).toBe('Cleared');
    expect(v.highlight.names.value).toEqual([]);
    expect(await run(v, 'highlight.clear', {})).toBe('Nothing highlighted');
  });

  it('track.add inserts a track after another in one patch and one undo step; a taken id is refused', async () => {
    const v = await make();
    const root = v.spec.root.value;
    const node = Cases['track.add'].ctx as { node: Node };
    expect(await run(v, 'track.add', Cases['track.add'].ctx)).toBe('Added');
    expect(v.view.resolved.value!.sections[0].tracks.map((t) => t.address)).toEqual([
      `${Protein}/letters`,
      DomainsTrack,
      `${Protein}/domains-2`,
    ]);
    expect(stackOf(v).steps.value.map((x) => x.label)).toEqual(['Add a track']);
    expect(await run(v, 'track.add', { address: Protein, node: node.node })).not.toBe('Added');
    expect(stackOf(v).steps.value.length).toBe(1);
    await stackOf(v).undo();
    expect(v.spec.root.value).toBe(root);
    expect(await run(v, 'track.add', { address: DomainsTrack, node: node.node })).toBe('No section there');
  });

  it('the bookmark key is the viewer’s with or without a selection, a plain d is nobody’s', async () => {
    const v = await make();
    const mod = { ctrlKey: true, metaKey: true };
    expect(v.config.view.keys['selection.bookmark']).toEqual(['Mod+d']);
    expect(v.input.key({ key: 'd', ...mod })).toBe(true);
    expect(v.addresses.refs.value.map((r) => r.id)).toEqual(['b']);
    await run(v, 'address.select', { address: '/refs/b' });
    expect(v.input.key({ key: 'd', ...mod })).toBe(true);
    await new Promise((r) => setTimeout(r, 0));
    expect(v.addresses.refs.value.map((r) => r.id)).toEqual(['b', 'bookmark-1']);
    expect(v.input.key({ key: 'd' })).toBe(false);
  });
});

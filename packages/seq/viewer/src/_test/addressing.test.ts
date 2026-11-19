/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { applyPatch } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { createSeqContext } from '@molstar/seqx-core/context';
import { SeqRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec, SeqLanguage } from '@molstar/seqx-core/lang';
import { locationOfResult, locationsEqual, rangesLocation, resolveAddressIn } from '@molstar/seqx-core/query';
import { recordDrawLists } from '@molstar/seqx-core/render';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { afterEach, describe, expect, it } from 'vitest';
import type { AddressingBehavior } from '../behaviors/addressing.js';
import { AddressingProvider } from '../behaviors/providers.js';
import { changeKindOf } from '../commands.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

// Addresses in the viewer, headless: the actions a link, a key, a test and an agent run (focus, reveal, select,
// highlight, hover on its own channel), refs as bookmarks, links both ways, the selection and the focus as addresses,
// and the two consumers of one address: the viewer's selection and an export of the same spec resolving it alike.

const Letters = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENN';
const Protein = '/views/main/protein';
const Gene = '/views/main/gene';

/** A protein of 30 residues with two domains, on its gene through a codon map (the gene's section follows). */
function spec(): Node {
  const b = createSeqSpec({ title: 'addressing' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: Letters });
  const g = b.data.sequence('g', { alphabet: 'dna', letters: 'ACGT'.repeat(30) });
  b.data.map('m', { from: p, to: g, method: 'codon', exons: { start: [10, 70], end: [60, 113] }, strand: '+' });
  const domains = b.data.ranges('domains', {
    target: p,
    rows: { start: new Int32Array([0, 15]), end: new Int32Array([10, 25]), id: ['A', 'B'] },
  });
  b.refs.ref('b', { to: `${domains}/(item :id B)`, label: 'Domain B' });
  b.refs.ref('start', { to: `${p}/(range :from 1 :to 5)` });
  b.views.view('main', (v) => {
    const protein = v.section('protein', { of: p });
    protein.track('letters').layer('l', { type: 'letters' });
    protein.track('domains').layer('d', { type: 'ranges', data: '/refs/domains' });
    v.section('gene', { of: g, follow: protein }).track('letters').layer('l', { type: 'letters' });
  });
  b.refs.ref('domains', { to: domains });
  return b.build();
}

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

async function make(rails = true) {
  frames = installFrames();
  const viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [
      [viewerItem('viewport.maxPxPerPosition'), 1000],
      // a slim rail row is the canvas's, never a figure's: compared with an export, none
      ...(rails ? [] : [[viewerItem('navigator.form'), 'none'] as const]),
    ],
  });
  viewers.push(viewer);
  await viewer.load(spec());
  viewer.resize(600);
  await viewer.settle();
  frames.flush();
  return viewer;
}

const range = (v: SeqViewer, s: string) => {
  const vp = v.view.get(s)!;
  return [vp.start, vp.end];
};

describe('actions on addresses', () => {
  it('focus moves the leading section to what the address names, and the follower follows', async () => {
    const v = await make();
    const gene = range(v, Gene);
    const o = v.addresses.focus('/refs/b');
    expect(o).toMatchObject({ done: true, status: 'resolved', message: 'Focused' });
    expect(range(v, Protein)).toEqual([15, 25]);
    expect(range(v, Gene)).not.toEqual(gene);
    // and through a link, the same action
    v.view.fit(Protein);
    const run = v.runLink('!seq-focus=/refs/start');
    expect(run?.outcomes.map((x) => x.done)).toEqual([true]);
    expect(range(v, Protein)).toEqual([0, 5]);
  });

  it('reveal leaves what is shown alone, and pans to what is not', async () => {
    const v = await make();
    v.view.setRange(Protein, 0, 10);
    expect(v.addresses.reveal('/refs/start').done).toBe(true);
    expect(range(v, Protein)).toEqual([0, 10]);
    expect(v.addresses.reveal('/data/p/(at :label 21)').done).toBe(true);
    expect(range(v, Protein)).toEqual([15.5, 25.5]);
  });

  it('select sets the selection, drawn by the selection pass, read back as its address', async () => {
    const v = await make();
    expect(v.addresses.select('/refs/b').done).toBe(true);
    const sel = v.selection.value.value;
    expect(sel.elements.keys).toEqual([{ owner: '/data/domains', key: 'B' }]);
    expect(locationsEqual(sel.shown, rangesLocation('/data/p', [15, 25]))).toBe(true);
    frames!.flush();
    expect(recordDrawLists(v.frame!.lists, ['selection'])).toContain('rect');
    expect(v.addresses.selectionAddress()).toBe('/refs/b');
    // set from an address of ids: read back as the ref that names exactly it
    v.addresses.setSelection('/data/domains/(item :id B)');
    expect(v.selection.value.value.address).toBe('/data/domains/(item :id B)');
    v.addresses.setSelection(undefined);
    expect(v.addresses.selectionAddress()).toBeUndefined();
    expect(v.addresses.copyLink()).toBeUndefined();
  });

  it('highlight sets a named channel; a link’s hover its own, never the pointer’s hover', async () => {
    const v = await make();
    v.runLink('!seq-highlight=/refs/start&seq-channel=paper');
    expect(v.highlight.names.value).toEqual(['paper']);
    v.addresses.hoverAddress('/refs/b');
    expect(v.highlight.names.value).toEqual(['paper', 'link-hover']);
    expect(v.hover.current.owner).toBeUndefined();
    v.addresses.unhover();
    expect(v.highlight.names.value).toEqual(['paper']);
  });

  it('report an address that names nothing, and do nothing', async () => {
    const v = await make();
    const before = range(v, Protein);
    const o = v.addresses.focus('/data/domains/(item :id Z)');
    expect(o).toMatchObject({ done: false, status: 'missing' });
    expect(o.issues.map((i) => i.code)).toEqual(['seq.address.missing']);
    expect(range(v, Protein)).toEqual(before);
    const bad = v.runLink('!seq-focus=data/p&seq-explode=/a');
    expect(bad?.outcomes).toEqual([]);
    expect(bad?.issues.map((i) => i.code)).toEqual(['repr.tree.address.syntax']);
  });

  it('runs as commands of the addressing behavior, each with its change kind', async () => {
    const v = await make();
    const ids = (v.behaviors.get('addressing') as AddressingBehavior).commands.map((c) => c.id);
    expect(ids).toEqual([
      'address.focus',
      'address.reveal',
      'address.select',
      'address.highlight',
      'address.hover',
      'highlight.clear',
      'address.unhover',
      'address.copy-link',
      'refs.add',
      'selection.bookmark',
    ]);
    for (const id of ids) expect(changeKindOf(id), id).toBeDefined();
    expect(ids.map((id) => changeKindOf(id))).toEqual([
      'navigate',
      'navigate',
      'session',
      'session',
      'session',
      'session',
      'session',
      'none',
      'patch',
      'patch',
    ]);
    expect(v.runCommand('address.focus', { address: '/refs/b' })).toBe('Focused');
    expect(v.runCommand('address.focus', {})).toBe('No address');
    expect(v.behavior(AddressingProvider)).toBeDefined();
  });
});

describe('refs as bookmarks', () => {
  it('lists the refs with what each names now; adding one is one patch', async () => {
    const v = await make();
    expect(v.addresses.refs.value.map((r) => [r.address, r.label, r.status])).toEqual([
      ['/refs/b', 'Domain B', 'resolved'],
      ['/refs/start', 'start', 'resolved'],
      ['/refs/domains', 'domains', 'resolved'],
    ]);
    const revision = v.spec.revision.value;
    // the default spec's history runs it as one step, in turn
    expect(await v.execute('refs.add', { ref: { id: 'here', to: '/data/p/(range :from 10 :to 20)' } })).toBe(
      'Bookmarked',
    );
    expect(v.spec.revision.value).not.toBe(revision);
    expect(v.addresses.refs.value.map((r) => r.address)).toContain('/refs/here');
    // a ref that would dangle is refused, and nothing changes
    const now = v.spec.revision.value;
    expect(await v.addresses.addRef({ id: 'gone', to: '/data/nothing' })).not.toBe(true);
    expect(v.spec.revision.value).toBe(now);
  });
});

describe('the session as addresses', () => {
  it('reads the focus as the address focused until the view moves, then as the region shown', async () => {
    const v = await make();
    v.addresses.setFocus('/refs/b');
    expect(v.addresses.focusAddress()).toBe('/refs/b');
    v.view.setRange(Protein, 2, 12);
    expect(v.addresses.focusAddress()).toBe('/data/p/(range :from 3 :to 12)');
    v.view.fit(Protein);
    expect(v.addresses.focusAddress()).toBe('/data/p/(all)');
  });

  it('copies a link of the selection, else of what the pointer is over', async () => {
    const v = await make();
    v.addresses.select('/data/p/(range :from 1 :to 5)');
    expect(v.addresses.copyLink()).toBe(v.addresses.link('focus', '/data/p/(range :from 1 :to 5)'));
    v.selection.clear();
    v.hover.set({ system: '/data/p', position: 3, layer: `${Protein}/letters/l` }, 'test');
    expect(v.addresses.copyLink()).toBe('!seq-focus=/data/p/%28at%20%3Alabel%204%29');
    expect(v.runCommand('address.copy-link', {})).toBe('!seq-focus=/data/p/%28at%20%3Alabel%204%29');
  });
});

describe('the session follows the spec', () => {
  it('rewrites its addresses through a move, and drops one that names nothing with an issue', async () => {
    const v = await make();
    v.addresses.select('/data/p/(range :from 1 :to 5)');
    v.addresses.highlightAddress('/refs/b', 'paper');
    const r = applyPatch(
      v.spec.root.value,
      { ops: [{ op: 'move', at: ['data', 'p'], parent: ['data'], place: { last: true }, rename: 'q' }] },
      { language: SeqLanguage },
    );
    expect(r.ok).toBe(true);
    if (!r.ok) return;
    v.addresses.followMoves(r.changes.moved);
    v.spec.apply(r.root);
    expect(v.addresses.selectionAddress()).toBe('/data/q/(range :from 1 :to 5)');
    expect(v.addresses.issues.value).toEqual([]);
    // a removal the session is not told of: what it named is gone, dropped with an issue naming it
    const gone = applyPatch(
      v.spec.root.value,
      { ops: [{ op: 'remove', at: ['refs', 'b'] }] },
      { language: SeqLanguage, allowDangling: true },
    );
    expect(gone.ok).toBe(true);
    if (!gone.ok) return;
    v.spec.apply(gone.root);
    expect(v.highlight.get('paper')).toBeUndefined();
    expect(v.addresses.issues.value.map((i) => [i.code, i.address])).toEqual([['seq.address.missing', '/refs/b']]);
    expect(v.addresses.selectionAddress()).toBe('/data/q/(range :from 1 :to 5)');
  });

  it('follows a patch the spec model applies: its moves before the new root, one step of the history', async () => {
    const v = await make();
    v.addresses.select('/data/p/(range :from 1 :to 5)');
    v.addresses.highlightAddress('/refs/b', 'paper');
    const o = await v.spec.patch('Rename the protein', {
      ops: [{ op: 'move', at: ['data', 'p'], parent: ['data'], place: { last: true }, rename: 'q' }],
    });
    expect(o.ok && o.changed).toBe(true);
    expect(v.addresses.selectionAddress()).toBe('/data/q/(range :from 1 :to 5)');
    expect(v.highlight.get('paper')?.address).toBe('/refs/b');
    expect(v.addresses.issues.value).toEqual([]);
  });

  it('does not report a focus done while its section is not placed', async () => {
    frames = installFrames();
    const v = createSeqViewer({
      ...DefaultSeqViewerSpec(),
      ui: undefined,
      config: [[viewerItem('viewport.maxPxPerPosition'), 1000]],
    });
    viewers.push(v);
    await v.load(spec());
    // no width yet: the move is kept as the intent
    expect(v.addresses.focus('/refs/b')).toMatchObject({ done: false, message: 'Kept until the view is placed' });
    v.resize(600);
    await v.settle();
    expect(range(v, Protein)).toEqual([15, 25]);
  });
});

describe('two consumers of one address', () => {
  it('the viewer’s selection and an export of the spec resolve it to the same Location, drawn alike', async () => {
    const v = await make(false);
    v.addresses.select('/refs/b');
    frames!.flush();
    const root = spec();
    const model = resolveSpec(root, { runtimes: [SeqRuntime] });
    const r = resolveAddressIn(root, model, '/refs/b');
    const exported = locationOfResult(root, model, r.result!).location;
    expect(locationsEqual(exported, v.selection.value.value.shown)).toBe(true);

    const mine = await v.renderSvg();
    const ctx = createSeqContext({ registry: [{ kind: 'format-runtime', runtime: SeqRuntime }] });
    const theirs = await ctx.renderSvg(root, {
      size: { width: v.width, height: v.height },
      ...(v.viewport ? { viewport: v.viewport } : {}),
      emphasis: { selection: exported },
    });
    ctx.retire();
    expect(mine.svg).toBe(theirs.svg);
    expect(recordDrawLists(theirs.lists, ['selection'])).toBe(recordDrawLists(v.frame!.lists, ['selection']));
    expect(recordDrawLists(theirs.lists, ['selection'])).toContain('rect');
  });
});

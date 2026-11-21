/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { checkSeqSpec } from '@molstar/seqx-core/lang';
import { keySet, location } from '@molstar/seqx-core/query';
import { type Frame, recordDrawLists } from '@molstar/seqx-core/render';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import type { SeqViewerRegistryEntry } from '@molstar/seqx-viewer/header-entries';
import { openedModel } from '@molstar/seqx-core/resolve';
import { HoverCardProvider } from '@molstar/seqx-viewer/behaviors';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { FakeStructureViewer } from '../../../support/fake-structure/fake-viewer.js';
import {
  blockLines,
  ComplexSection,
  hostComplexSpec,
  HostFile,
  HostSource,
  heldLines,
  issueLines,
  BFactorLegend,
  HostChains,
  HostTableFile,
  hostDataSpec,
  hostTableOf,
  unitAddress,
  withDefaultMmcif,
  withFakeHost,
} from '../compose.js';

// The Host data story, headless on a fake clock, composed as the story composes it (compose.ts): one spec naming PDB
// 2OCJ by its file, a viewer with the fake structure viewer's adapter (it serves 2OCJ from its chain table) and a viewer
// with the default adapter (the mmCIF runtime reading the file). The same spec gives the same nodes and the same
// drawing, the fake's run reads no asset, and "Drop chain" in the fake is followed by the sequences through `changes`
// and `reattach`: the other chains' nodes kept by key, the selection kept by key.

const viewers: SeqViewer[] = [];
const fakes: FakeStructureViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  for (const f of fakes.splice(0)) f.dispose();
  vi.restoreAllMocks();
  vi.useRealTimers();
});

const table = () => hostTableOf(readFileSync(`${PackageFolder}${HostTableFile}`, 'utf8'));

/** A settled viewer of the spec, by the fake host or by the default adapter, and what its assets port read. */
async function open(by: 'fake' | 'default', spec = hostDataSpec()) {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  // a frame is a timer 16 ms on: what a host's edit is reopened by
  vi.spyOn(DomHost, 'requestAnimationFrame').mockImplementation((fn) => {
    setTimeout(() => fn(performance.now()), 16);
    return undefined as never;
  });
  const fake = new FakeStructureViewer(table());
  fakes.push(fake);
  const base = DefaultSeqViewerSpec();
  const reads: string[] = [];
  // how often the source is opened (a reopen is an opening): the fake adapter's, counted
  const opens: string[] = [];
  const counted = (registry: ReturnType<typeof withFakeHost>) =>
    registry.map((e) =>
      e.kind === 'format-runtime' && 'runtime' in e && e.runtime.name === 'fake-structure'
        ? {
            ...e,
            runtime: {
              ...e.runtime,
              open: (...a: Parameters<typeof e.runtime.open>) => (opens.push(a[0].address), e.runtime.open(...a)),
            },
          }
        : e,
    );
  const read = (uri: string) => (reads.push(uri), readFileSync(`${PackageFolder}${uri}`));
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    registry:
      by === 'fake'
        ? (counted(withFakeHost(base.registry ?? [], fake)) as SeqViewerRegistryEntry[])
        : withDefaultMmcif(base.registry ?? []),
    ports: {
      assets: () => ({
        reader: {
          readText: async (uri) => read(uri).toString('utf8'),
          readBytes: async (uri) => new Uint8Array(read(uri)),
        },
      }),
    },
  });
  viewers.push(v);
  v.resize(1100);
  await v.load(spec);
  await vi.advanceTimersByTimeAsync(500);
  await v.settle();
  v.draw();
  return { v, fake, reads, opens };
}

const nodesOf = (v: SeqViewer) => v.store.store.coverage(HostSource)?.nodes ?? [];
const unitsOf = (v: SeqViewer) =>
  nodesOf(v)
    .find((n) => n.key === 'complex')!
    .node.children!.filter((c) => c.kind === 'unit')
    .map((u) => u.id);
const drawn = (v: SeqViewer) => recordDrawLists((v.frame as Frame).lists, ['data']);

/** The viewer moved on: the frame the reopen was asked for, the tiles read again, drawn. */
async function follow(v: SeqViewer) {
  await vi.advanceTimersByTimeAsync(500);
  await v.settle();
}

describe('Host data: the fake serving 2OCJ', () => {
  it('the fake reads no asset, the default adapter the file once; equal nodes, and the B-factor tracks drawn the same', async () => {
    expect(checkSeqSpec(hostDataSpec())).toEqual([]);
    const host = await open('fake');
    const file = await open('default');
    expect(host.reads).toEqual([]);
    expect(file.reads).toEqual([HostFile]);
    expect(host.v.modelIssues.value).toEqual([]);
    expect(file.v.modelIssues.value).toEqual([]);
    expect(host.v.layoutIssues.value).toEqual([]);
    const a = nodesOf(file.v);
    const b = nodesOf(host.v);
    expect(a.map((n) => n.key)).toEqual(['1', '2', 'compounds', 'complex']);
    expect(b.map((n) => n.key)).toEqual(a.map((n) => n.key));
    for (let i = 0; i < a.length; i++) expect(b[i].node).toEqual(a[i].node);
    // the tracks read the source's channel: bars (not the gray base alone) on the modeled residues, and the same bars
    const lists = drawn(host.v);
    expect(lists.match(/^rect .* #5b7c99ff$/gm)?.length ?? 0).toBeGreaterThan(100);
    expect(lists).toBe(drawn(file.v));
  });

  it('the legend’s words about the ends are the table’s: chains A and D have no B-factor at residues 1–2 and 197–219 (F12)', () => {
    const fake = new FakeStructureViewer(table());
    for (const asym of HostChains) {
      const row = fake.rows.find((r) => r.asym === asym)!;
      const seen = row.residues.flatMap((r, i) => (r.observed ? [i + 1] : []));
      expect([seen[0], seen[seen.length - 1], row.residues.length]).toEqual([3, 196, 219]);
    }
    expect(BFactorLegend).toContain('residues 1–2 and 197–219');
  });

  it('Drop chain: the sequences follow a frame later through `changes` and `reattach`; the other nodes kept, no asset read', async () => {
    const { v, fake, reads } = await open('fake');
    const sequence = v.store.store.getNode(HostSource, '1');
    const molecule = v.store.store.getNode(HostSource, '2');
    const complex = v.store.store.getNode(HostSource, 'complex')!;
    expect(unitsOf(v)).toEqual(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']);
    expect(heldLines(v)).toEqual([
      'sequence 1: 219 residues',
      'molecule 2: ZINC ION',
      'compounds compounds',
      'complex complex: 8 units A B C D E F G H',
    ]);
    const barsBefore = drawn(v)
      .split('\n')
      .filter((l) => l.startsWith('rect') && l.endsWith('#5b7c99ff'));
    expect(fake.chainToDrop()).toBe('D');
    expect(fake.dropChain('D')).toBe(true);
    // nothing yet: the reopen is asked for the next frame
    expect(unitsOf(v)).toEqual(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']);
    await follow(v);
    expect(unitsOf(v)).toEqual(['A', 'B', 'C', 'E', 'F', 'G']);
    expect(heldLines(v)[3]).toBe('complex complex: 6 units A B C E F G');
    expect(v.store.store.getNode(HostSource, '1')).toBe(sequence);
    expect(v.store.store.getNode(HostSource, '2')).toBe(molecule);
    expect(v.store.store.getNode(HostSource, 'complex')!.version).toBeGreaterThan(complex.version);
    expect(reads).toEqual([]);
    // chain A's track is read again from the reopened source and its bars are as they were, D's gone
    const bars = (t: string) => t.split('\n').filter((l) => l.startsWith('rect') && l.endsWith('#5b7c99ff'));
    expect(bars(drawn(v))).toEqual(barsBefore.slice(0, barsBefore.length / 2));
  });

  it('Drop chain: chain D’s track is blank a few frames later with no further input, chain A’s bars as they were', async () => {
    const { v, fake } = await open('fake');
    const bars = (t: string) => t.split('\n').filter((l) => l.startsWith('rect') && l.endsWith('#5b7c99ff'));
    const before = bars(drawn(v));
    expect(before.length).toBeGreaterThan(200);
    fake.dropChain('D');
    // the clock only: no draw, no settle, no input
    await vi.advanceTimersByTimeAsync(300);
    const after = bars(drawn(v));
    expect(after.length).toBe(before.length / 2);
    expect(after).toEqual(before.slice(0, after.length));
    // and the viewer says why chain D has nothing to draw: its unit is gone from the source
    expect(v.modelIssues.value.map((i) => [i.code, i.address])).toEqual([
      ['seq.resolve.system.unnamed', expect.any(String)],
    ]);
    expect(v.layoutIssues.value.map((i) => [i.code, i.address])).toEqual([
      ['seq.render.axis.pending', '/views/main/chain-D'],
    ]);
    expect(issueLines(v).join('\n')).toContain('seq.render.axis.pending');
  });

  it('changes inside one frame reopen once; a selection on a unit that stays is kept, one on a unit that goes is dropped, said once', async () => {
    const { v, fake, opens } = await open('fake');
    expect(opens).toEqual([HostSource]);
    const select = (chain: string, from: number, to: number) => {
      const l = location({ parts: [{ system: unitAddress(chain), ranges: [from, to] }] });
      v.selection.set({ elements: keySet([]), location: l, shown: l });
    };
    const selected = () => v.selection.value.value.location.parts.map((p) => p.system);
    select('A', 10, 20);
    fake.dropChain('D');
    fake.dropChain('C');
    await follow(v);
    // two edits inside one frame: one opening
    expect(opens).toEqual([HostSource, HostSource]);
    expect(unitsOf(v)).toEqual(['A', 'B', 'E', 'F']);
    expect(selected()).toEqual([unitAddress('A')]);
    select('B', 30, 40);
    fake.dropChain('B');
    await follow(v);
    expect(unitsOf(v)).toEqual(['A', 'E']);
    expect(selected()).toEqual([]);
    expect(v.store.issues.value.filter((i) => i.code === 'seq.viewer.selection.dropped')).toHaveLength(1);
  });
});

describe('Host data: the host’s complex as a section (a section of a source’s complex)', () => {
  /** The residues' track drawn a block (its clip: a chain's band, an ion's chip). */
  const residueBands = (v: SeqViewer) => drawn(v).match(/^clip [\d.]+ 62 /gm)?.length ?? 0;
  it('pending before the source opens (no issue), a block a unit once it holds the complex; the residues drawn', async () => {
    vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
    const fake = new FakeStructureViewer(table());
    fakes.push(fake);
    const base = DefaultSeqViewerSpec();
    const v = createSeqViewer({ ...base, ui: undefined, registry: withFakeHost(base.registry ?? [], fake) });
    viewers.push(v);
    v.resize(1100);
    void v.load(hostComplexSpec());
    expect(blockLines(v)).toEqual(['pending: the source holds no complex yet']);
    expect(v.view.resolved.value?.issues).toEqual([]);
    expect(checkSeqSpec(hostComplexSpec())).toEqual([]);
    const { v: open2 } = await open('fake', hostComplexSpec());
    expect(open2.view.resolved.value?.issues).toEqual([]);
    expect(open2.modelIssues.value).toEqual([]);
    expect(open2.layoutIssues.value).toEqual([]);
    expect(blockLines(open2).map((l) => l.split(': ')[1])).toEqual(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']);
    // the residues drawn: a band a chain's block (zoomed out), a chip an ion's (its name)
    expect(residueBands(open2)).toBe(8);
    expect(drawn(open2).match(/^text "ZINC ION"/gm)).toHaveLength(4);
  });

  it('Drop chain and Restore: the section reflows a frame later, one resolve a node version, none a frame', async () => {
    const { v, fake, reads } = await open('fake', hostComplexSpec());
    const units = () => blockLines(v).map((l) => l.split(': ')[1]);
    expect(units()).toEqual(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']);
    let resolves = 0;
    let versions = 0;
    const offs = [v.view.resolved.listen(() => resolves++), v.store.store.nodeVersion.listen(() => versions++)];
    fake.dropChain('D');
    await follow(v);
    expect(units()).toEqual(['A', 'B', 'C', 'E', 'F', 'G']);
    expect(residueBands(v)).toBe(6);
    expect(v.view.resolved.value?.issues).toEqual([]);
    expect(v.layoutIssues.value).toEqual([]);
    expect(resolves).toBeGreaterThan(0);
    expect(resolves).toBeLessThanOrEqual(versions);
    // frames and a hover read the view resolved for the version: no new one
    const shown = v.view.resolved.value;
    const layouts = v.counts.layouts;
    for (let i = 0; i < 3; i++) v.draw();
    expect(v.view.resolved.value).toBe(shown);
    expect(v.counts.layouts).toBe(layouts);
    expect(fake.restore()).toBe(true);
    await follow(v);
    expect(units()).toEqual(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']);
    expect(residueBands(v)).toBe(8);
    expect(reads).toEqual([]);
    for (const off of offs) off();
    expect(ComplexSection).toBe(v.view.resolved.value?.sections[0].address);
  });
});

describe('Host data: a held unit’s letters, names and copy (the opened model read everywhere the view is)', () => {
  /** The letters drawn on the residues track (its glyph runs). */
  const glyphs = (v: SeqViewer) =>
    (drawn(v).match(/^glyphs .*$/gm) ?? []).join(' ').match(/ [A-Z] [\d.]+ #/g)?.length ?? 0;
  /** The residues track's fills (cells, ticks, bars) by colour. */
  const fills = (v: SeqViewer) =>
    new Set((drawn(v).match(/^rect [\d.]+ 6[2-9] .* (#[0-9a-f]{8})/gm) ?? []).map((l) => l.split(' ')[5]));
  /** A block's bar under a pixel a residue, and a cell whose letter is not held (`seq-missing`). */
  const Bars = '#c4c9cfff';
  const Missing = '#dedad3ff';

  it('the residues track draws the chain’s letters zoomed in, the same before a drop and after; zoomed out its colours, never the missing fill', async () => {
    const { v, fake } = await open('fake', hostComplexSpec());
    const zoomed = async () => {
      v.view.setRange(ComplexSection, 84, 104);
      await follow(v);
      v.draw();
      return glyphs(v);
    };
    const out = async () => {
      v.view.fit(ComplexSection);
      await follow(v);
      v.draw();
      return [...fills(v)];
    };
    const before = await zoomed();
    expect(before).toBeGreaterThan(10);
    // eight blocks at the fit: under a pixel a residue, a bar a block (presence, no data read: the density step)
    expect(await out()).toEqual([Bars]);
    fake.dropChain('D');
    await follow(v);
    // six blocks: over a pixel a residue, the residues' colours, never the missing fill (letters not held)
    const after = await out();
    expect(after).not.toContain(Missing);
    expect(after.filter((c) => c !== Bars).length).toBeGreaterThan(2);
    expect(await zoomed()).toBe(before);
  });

  it('the hover card names a held unit’s residue as its chain (A 109), never its address', async () => {
    const { v } = await open('fake', hostComplexSpec());
    v.view.setRange(ComplexSection, 100, 120);
    await follow(v);
    v.draw();
    const s = v.frame!.layout.sections[0];
    const t = s.tracks[1];
    v.input.move((s.axis!.x(108) + s.axis!.x(109)) / 2, t.box.y + t.box.height / 2);
    await vi.advanceTimersByTimeAsync(600);
    const card = v.behavior(HoverCardProvider)!.shown.value!.card;
    expect(card.title).not.toContain('(node');
    expect(card.title).toMatch(/^A\b/);
  });

  it('copies a held unit’s letters (selection.copyLetters over the opened model)', async () => {
    const copied: string[] = [];
    vi.stubGlobal('navigator', { clipboard: { writeText: async (t: string) => void copied.push(t) } });
    const { v } = await open('fake', hostComplexSpec());
    const letters = String(nodesOf(v).find((n) => n.key === '1')!.node.params!.letters);
    await v.execute('selection.add', { range: { system: unitAddress('A'), start: 10, end: 25 } });
    expect(await v.execute('selection.copy-letters', {})).toBe('Copied 15 letters');
    expect(copied.at(-1)).toBe(letters.slice(10, 25));
  });

  it('the view’s overlay and the frame’s are one: keyed by the node version and the runtimes alike', async () => {
    const { v } = await open('fake', hostComplexSpec());
    const opened = () => openedModel(v.spec.model.value, v.store.currentOverlay(v.spec.model.value));
    expect(v.view.model.value).toBe(opened());
    const before = v.view.model.value;
    // another set of runtimes (a late extension's): another builder, the view re-keyed with it
    v.store.setRuntimes([...v.context.runtimes]);
    expect(v.view.model.value).not.toBe(before);
    expect(v.view.model.value).toBe(opened());
    // the frame's runtimes (the context's) back: one overlay again, the frame's and the view's
    v.draw();
    expect(v.view.model.value).toBe(opened());
  });
});

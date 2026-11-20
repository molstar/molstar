/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// 1FIN as a complex (the Complex's AsAComplex stories), headless on a fake clock, composed as the stories compose them
// (the editing view's viewer with `links.hover` and the host's stack; the atoms view and its SVG export): the blocks and
// chips, every contact (per pair of units) an arc, ATP's chip colored on its own, ATP lighting its CDK2 residues and
// back, a contact's ends on the copy it is in, the author numbers, ATP switched to its atoms by one patch and undone.

import { readFileSync } from 'node:fs';
import { DomHost } from '@molstar/ui/dom';
import {
  type BlockPlacement,
  DarkTokens,
  type LettersLay,
  type PairsLay,
  recordDrawLists,
} from '@molstar/seqx-core/render';
import { ComplexViews } from '@molstar/seqx-examples';
import { Data, entryLigands, ligand } from '@molstar/seqx-examples/complex/data';
import { transplantUnit } from '@molstar/seqx-examples/complex/units';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { LinksHoverChannel, LinksHoverProvider } from '@molstar/seqx-viewer/behaviors';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { composeEditing, type EditingSession } from '../../features/editing/compose.js';
import { atpForm, toggleAtpForm } from '../complex/compose.js';
import { openExample } from './draw-checks.js';

const LinksHoverOn = [{ provider: LinksHoverProvider, params: undefined }] as const;
const view = (id: string) => ComplexViews.find((v) => v.id === id)!;
const Section = '/views/main/complex';

let viewer: SeqViewer | undefined;
let session: EditingSession | undefined;
const raf = DomHost.requestAnimationFrame;

beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
});

afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  viewer?.retire();
  viewer = undefined;
  session?.dispose();
  session = undefined;
  vi.useRealTimers();
});

const fixtures = () => ({
  reader: {
    readText: async (uri: string) => readFileSync(`${PackageFolder}${uri}`, 'utf8'),
    readBytes: async (uri: string) => new Uint8Array(readFileSync(`${PackageFolder}${uri}`)),
  },
});

async function edit(): Promise<EditingSession> {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  const v = view('units');
  const s = composeEditing(v, { ports: { assets: fixtures }, delayMs: 0, behaviors: LinksHoverOn });
  session = s;
  s.viewer.resize(900);
  await s.viewer.load(v.build());
  await vi.advanceTimersByTimeAsync(200);
  s.viewer.draw();
  return s;
}

const placement = (v: SeqViewer) => v.frame!.layout.sections[0].placement as BlockPlacement;
const arcs = (v: SeqViewer) =>
  v.frame!.layout.sections[0].sectionLays.reduce(
    (n, l) => n + (l.type.name === 'pairs' ? (l.lay as PairsLay).n : 0),
    0,
  );
/** The contacts: across each copy, each ATP's residues, each transplant's residues. */
const atpResidues = (chain: string) =>
  [...new Set(Data.ligands.find((l) => l.comp === 'ATP' && l.chain === chain)!.contacts.map(([r]) => r))].sort(
    (a, b) => a - b,
  );
const Contacts =
  Data.contacts.length +
  atpResidues('A').length +
  atpResidues('C').length +
  entryLigands().reduce((n, l) => n + new Set(l.contacts.map(([r]) => r)).size, 0);

describe('1FIN as a complex, headless', () => {
  it('the chains ×2, ATP a chip ×2 colored on its own, the transplants chips; every contact an arc; ATP’s atoms one patch, undone', async () => {
    const s = await edit();
    const v = s.viewer;
    expect(v.layoutIssues.value).toEqual([]);
    const p = placement(v);
    expect(p.strips.axes.map((a, i) => [a.system, p.strips.chip[i]])).toEqual([
      ['/data/1fin/A', 0],
      ['/data/1fin/B', 0],
      ['/data/1fin/atp-a', 1],
      ...entryLigands().map((l) => [`/data/1fin/${transplantUnit(l.entry)}`, 1]),
    ]);
    expect(arcs(v)).toBe(Contacts);
    // ATP's pill its own color, the transplants' the layer's default
    const lays = v.frame!.layout.sections[0].lays;
    const pill = (system: string) => {
      const i = p.strips.axes.findIndex((a) => a.system === system);
      for (const band of lays.slice(i * p.strips.bands, (i + 1) * p.strips.bands))
        for (const l of band ?? []) if (l.type.name === 'letters') return (l.lay as LettersLay).fill[0];
      return undefined;
    };
    expect(pill('/data/1fin/atp-a')).toBe(0xe15759ff);
    expect(pill(`/data/1fin/${transplantUnit(entryLigands()[0].entry)}`)).not.toBe(0xe15759ff);
    expect(await toggleAtpForm(s)).toBe('Atoms');
    await vi.advanceTimersByTimeAsync(50);
    v.draw();
    const q = placement(v);
    const atp = q.strips.axes.findIndex((a) => a.system === '/data/1fin/atp-a/(atoms)');
    expect([atp >= 0, q.strips.chip[atp]]).toEqual([true, 0]);
    expect(v.layoutIssues.value).toEqual([]);
    await s.stack.undo();
    expect(atpForm(s)).toBe('whole');
  });

  it('ATP’s chip lights the CDK2 residues it contacts (`links.hover`); one of them lights its ATP', async () => {
    const s = await edit();
    const v = s.viewer;
    const owner = {};
    const lit = () => v.highlight.get(LinksHoverChannel)?.location.parts.map((q) => [q.system, Array.from(q.ranges)]);
    v.hover.set({ section: Section, system: '/data/1fin/atp-a', position: 0 }, owner);
    const runs = lit()!;
    expect(runs.map((r) => r[0])).toEqual(['/data/1fin/A']);
    const at: number[] = [];
    const r = runs[0][1] as number[];
    for (let k = 0; k < r.length; k += 2) for (let i = r[k]; i < r[k + 1]; i++) at.push(i);
    expect(at).toEqual(atpResidues('A'));
    expect(ligand('ATP').atoms).toHaveLength(31);
    v.hover.clear(owner);
  });

  it('ATP by its atoms: each atom its element; the SVG draws the canvas’s lists in both themes', async () => {
    viewer = await openExample(view('units-atoms'), { fixtures: true, behaviors: LinksHoverOn });
    const v = viewer;
    expect(v.layoutIssues.value).toEqual([]);
    v.setRailsShown(false);
    for (const tokens of [undefined, DarkTokens]) {
      if (tokens) v.setTokens(tokens);
      v.draw();
      const svg = v.renderSvg({ idPrefix: 'cu-' });
      await vi.advanceTimersByTimeAsync(400);
      const r = await svg;
      expect(r.issues).toEqual([]);
      expect(recordDrawLists(r.lists, ['data'])).toBe(recordDrawLists(v.frame!.lists, ['data']));
    }
  });
});

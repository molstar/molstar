/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Antibody's Fc glycan stories (3AVE), headless on a fake clock, composed as the stories compose them (the editing
// view's viewer with `links.hover` listed and the host's stack; the residues view and its SVG export): the blocks and
// chips, every link with both ends in view an arc, the glycan's chip lighting Asn297 and back, the form switched by one
// patch and undone, the glycan's residues each labeled by its code with its seven glycosidic arcs, both themes alike.

import { readFileSync } from 'node:fs';
import { DomHost } from '@molstar/ui/dom';
import { type BlockPlacement, DarkTokens, type PairsLay, recordDrawLists } from '@molstar/seqx-core/render';
import { AntibodyViews } from '@molstar/seqx-examples';
import {
  Asn297,
  Disulfides,
  FcChains,
  fcPosition,
  GlycanLinks,
  FcLength,
  GlycanResidues,
  ZincHistidines,
} from '@molstar/seqx-examples/antibody/fc-glycan';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { LinksHoverChannel, LinksHoverProvider } from '@molstar/seqx-viewer/behaviors';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { composeEditing, type EditingSession } from '../../features/editing/compose.js';
import { glycanForm, toggleGlycanForm } from '../antibody/compose.js';
import { openExample } from './draw-checks.js';

const LinksHoverOn = [{ provider: LinksHoverProvider, params: undefined }] as const;
const view = (id: string) => AntibodyViews.find((v) => v.id === id)!;

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

/** The fixtures read at once (the file system would not settle on the fake clock). */
const fixtures = () => ({
  reader: {
    readText: async (uri: string) => readFileSync(`${PackageFolder}${uri}`, 'utf8'),
    readBytes: async (uri: string) => new Uint8Array(readFileSync(`${PackageFolder}${uri}`)),
  },
});

async function edit(): Promise<EditingSession> {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  const v = view('fc-glycan');
  const s = composeEditing(v, { ports: { assets: fixtures }, delayMs: 0, behaviors: LinksHoverOn });
  session = s;
  s.viewer.resize(900);
  await s.viewer.load(v.build());
  await vi.advanceTimersByTimeAsync(200);
  s.viewer.draw();
  return s;
}

const placement = (v: SeqViewer) => v.frame!.layout.sections[0].placement as BlockPlacement;
const arcs = (v: SeqViewer, layer: string) =>
  v.frame!.layout.sections[0].sectionLays.reduce(
    (n, l) =>
      n + (l.type.name === 'pairs' && (l.lay as { key: string }).key.includes(`/${layer}`) ? (l.lay as PairsLay).n : 0),
    0,
  );
/** The complex's links: per chain the glycan and two disulfides, and the zinc's two histidines on A. */
const Links = FcChains.length * (1 + Disulfides.length) + ZincHistidines.length;

describe('the Fc glycan stories, headless', () => {
  it('the editing view: the chains ×2, the glycans a chip ×2, the zinc a chip; every link an arc; the form one patch, undone', async () => {
    const s = await edit();
    const v = s.viewer;
    expect(v.layoutIssues.value).toEqual([]);
    const p = placement(v);
    expect(p.strips.axes.map((a, i) => [a.system, p.strips.chip[i]])).toEqual([
      ['/data/3ave/A', 0],
      ['/data/3ave/C', 1],
      ['/data/3ave/E', 1],
    ]);
    expect(arcs(v, 'links')).toBe(Links);
    expect(await toggleGlycanForm(s)).toBe('Residues');
    expect(glycanForm(s)).toBe('residues');
    await vi.advanceTimersByTimeAsync(50);
    v.draw();
    const q = placement(v);
    const glycan = q.strips.axes.findIndex((a) => a.system === '/data/3ave/C');
    expect(q.strips.chip[glycan]).toBe(0);
    expect(arcs(v, 'glycan')).toBeGreaterThanOrEqual(GlycanLinks.length);
    expect(v.layoutIssues.value).toEqual([]);
    await s.stack.undo();
    expect(glycanForm(s)).toBe('whole');
  });

  it('the glycan’s chip lights Asn297 (`links.hover`); Asn297 lights the glycan', async () => {
    const s = await edit();
    const v = s.viewer;
    const owner = {};
    const section = '/views/main/fc';
    const lit = () => v.highlight.get(LinksHoverChannel)?.location.parts.map((q) => [q.system, Array.from(q.ranges)]);
    v.hover.set({ section, system: '/data/3ave/C', position: 0 }, owner);
    expect(lit()).toEqual([['/data/3ave/A', [fcPosition(Asn297), fcPosition(Asn297) + 1]]]);
    v.hover.set({ section, system: '/data/3ave/A', position: fcPosition(Asn297) }, owner);
    expect(lit()).toEqual([['/data/3ave/C', [0, GlycanResidues.length]]]);
    v.hover.clear(owner);
  });

  it('by residues: each residue its code, the seven glycosidic links arcs; the SVG draws the canvas’s lists in both themes', async () => {
    viewer = await openExample(view('fc-glycan-residues'), { fixtures: true, behaviors: LinksHoverOn });
    const v = viewer;
    expect(v.layoutIssues.value).toEqual([]);
    v.setRailsShown(false);
    v.draw();
    // at the fit view the glycan's eight residues are a sliver of the section; its arcs drawn whole
    expect(arcs(v, 'glycan')).toBe(GlycanLinks.length);
    // zoomed onto the glycan's block: each residue labeled by its code
    // the section's units: the Fc block's 223, then the glycan's eight
    v.view.setRange('/views/main/fc', FcLength, FcLength + GlycanResidues.length);
    await vi.advanceTimersByTimeAsync(50);
    v.draw();
    const data = recordDrawLists(v.frame!.lists, ['data', 'labels']);
    for (const code of new Set(GlycanResidues)) expect(data).toContain(`text "${code}"`);
    expect(arcs(v, 'glycan')).toBe(GlycanLinks.length);
    v.view.setRange('/views/main/fc', 0, FcLength + GlycanResidues.length + 1);
    await vi.advanceTimersByTimeAsync(50);
    for (const tokens of [undefined, DarkTokens]) {
      if (tokens) v.setTokens(tokens);
      v.draw();
      const svg = v.renderSvg({ idPrefix: 'fc-' });
      await vi.advanceTimersByTimeAsync(400);
      const r = await svg;
      expect(r.issues).toEqual([]);
      expect(recordDrawLists(r.lists, ['data'])).toBe(recordDrawLists(v.frame!.lists, ['data']));
    }
  });
});

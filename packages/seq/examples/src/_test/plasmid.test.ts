/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Plasmid dataset (examples.md, Plasmid): pUC19c, INSDC L09137.2, recorded; its features found in the sequence and
// its sites computed on it, checked against known facts (bla is 861 nt, 286 residues and a stop, MSIQHFRVAL…; EcoRI and
// HindIII cut once, in the MCS); every derived feature keeps its length under every rotation of the record and of its
// reverse complement; the views drawn against the data (each feature where its positions are, across the origin one
// item; the arrows by strand on the piece holding the end they point to, bent through the ring's warp); the scale
// variant and its frame times (performance.md, a ring).

import { createHash } from 'node:crypto';
import { existsSync, readFileSync, statSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { applyPatch } from '@molstar/core/repr/patch';
import { printAddress } from '@molstar/core/repr/tree';
import { createSeqContext } from '@molstar/seqx-core/context/default';
import { encodeSeq, inlineSeq, SeqRuntime } from '@molstar/seqx-core/data';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import {
  DarkTokens,
  type FrameInput,
  hitView,
  LightTokens,
  type RangesLay,
  type RingPlacement,
  recordDrawLists,
} from '@molstar/seqx-core/render';
import { BoundedTries, itWithin, type Within } from '@molstar/test-harness/timing';
import { describe, expect, it } from 'vitest';
import { nodeFixtureAssets, PackageFolder } from '../node-assets.js';
import { ScaleFeatures, ScaleLength, scaleData } from '../plasmid/data.js';
import {
  BlaPrefix,
  deriveFeatures,
  Enzymes,
  featureLetters,
  findSites,
  gff3Text,
  LacZaPrefix,
  PlasmidFiles,
  PlasmidLength,
  PlasmidName,
  recordModule,
  reverseComplement,
  rotate,
  rotatedAnnotations,
  rotateFeature,
  rotateLetters,
  RotatedOrigin,
  seqText,
  translate,
} from '../plasmid/derive.js';
import { PublishedLetters } from '../plasmid/record.js';
import { PlasmidViews } from '../views.js';
import { publishedLetters, seqLetters, settled } from './plasmid-helpers.js';

const update = process.env.SEQ_RECORDS === 'update';
const view = (id: string) => PlasmidViews.find((v) => v.id === id)!;
const letters = PublishedLetters;
const features = deriveFeatures(letters);
const feature = (id: string) => features.find((f) => f.id === id)!;

function compare(file: string, text: string) {
  if (update) writeFileSync(file, text);
  expect(existsSync(file), `no ${file}: run with SEQ_RECORDS=update`).toBe(true);
  expect(text).toBe(readFileSync(file, 'utf8'));
}

const context = () =>
  createSeqContext({
    registry: [
      { kind: 'format-runtime', runtime: SeqRuntime },
      { kind: 'format-runtime', runtime: RecordsRuntime },
    ],
    ports: { assets: () => nodeFixtureAssets() },
  });

describe('the Plasmid fixtures', () => {
  it('the recorded file, the module of its letters, the rotated letters and the GFF3 are what their generators write', () => {
    const fromFile = publishedLetters();
    expect(fromFile).toBe(letters);
    compare(`${PackageFolder}${PlasmidFiles.published}`, seqText(PlasmidName, fromFile));
    compare(`${PackageFolder}${PlasmidFiles.rotated}`, seqText(PlasmidName, rotateLetters(fromFile, RotatedOrigin)));
    compare(fileURLToPath(new URL('../plasmid/record.ts', import.meta.url)), recordModule(fromFile));
    const rotated = rotatedAnnotations(fromFile);
    compare(`${PackageFolder}${PlasmidFiles.rotatedFeatures}`, gff3Text(rotated.features, rotated.sites));
  });

  it('the seq file is what the format’s encoder writes for the record', () => {
    const data = inlineSeq(
      { name: PlasmidName, length: PlasmidLength, alphabet: 'dna', cyclic: true },
      { letters: { form: 'letters', start: 0, letters } },
    );
    expect(JSON.parse(readFileSync(`${PackageFolder}${PlasmidFiles.published}`, 'utf8'))).toEqual(
      JSON.parse(encodeSeq(data)),
    );
  });

  it('sources.json lists every file of the folder: the record’s entry, its licence, the bounds and the checksums', () => {
    const rows = JSON.parse(readFileSync(`${PackageFolder}fixtures/plasmid/sources.json`, 'utf8')) as {
      file: string;
      url: string;
      licence: string[];
      attribution: string;
      derived: string;
      sha256Fixture: string;
      bytesFixture: number;
      bytesDownloaded: number;
    }[];
    expect(rows.map((r) => r.file).sort()).toEqual([
      'plasmid/L09137.2.seq',
      'plasmid/pUC19c-rotated-features.gff3',
      'plasmid/pUC19c-rotated.seq',
    ]);
    const at = (file: string) => rows.find((r) => r.file === `plasmid/${file}`)!;
    const [record, rotated, gff] = [at('L09137.2.seq'), at('pUC19c-rotated.seq'), at('pUC19c-rotated-features.gff3')];
    expect(record.url).toBe('https://www.ebi.ac.uk/ena/browser/api/fasta/L09137.2');
    expect(record.licence).toEqual(['LicenseRef-NCBI-public-data']);
    expect(record.attribution).toContain('L09137.2');
    for (const g of [rotated, gff]) {
      expect(g.licence).toEqual(['LicenseRef-seqx-generated']);
      expect(g.derived).toContain('derive.ts');
    }
    let total = 0;
    for (const r of rows) {
      const bytes = readFileSync(`${PackageFolder}fixtures/${r.file}`);
      expect(bytes.byteLength).toBe(r.bytesFixture);
      expect(createHash('sha256').update(bytes).digest('hex')).toBe(r.sha256Fixture);
      expect(statSync(`${PackageFolder}fixtures/${r.file}`).size).toBeLessThan(1_000_000);
      total += bytes.byteLength;
    }
    expect(total).toBeLessThan(2_000_000);
  });
});

describe('the Plasmid’s anchors', () => {
  it('the record: 2,686 nt of ACGT, the INSDC sequence’s first letters', () => {
    expect(letters.length).toBe(2686);
    expect(letters).toMatch(/^[ACGT]+$/);
    expect(letters.slice(0, 20)).toBe('TCGCGCGTTTCGGTGATGAC');
  });

  it('bla is 861 nt on the minus strand, 286 residues and a stop, MSIQHFRVAL…', () => {
    const bla = feature('bla');
    expect([bla.strand, bla.length, bla.type]).toEqual(['-', 861, 'CDS']);
    const protein = translate(featureLetters(letters, bla));
    expect(protein.length).toBe(287);
    expect(protein.startsWith(BlaPrefix)).toBe(true);
    expect(protein.endsWith('*')).toBe(true);
    expect(protein.indexOf('*')).toBe(286);
  });

  it('lacZα is a 324 nt open reading frame on the minus strand, MTMITPSLHA…', () => {
    const lac = feature('lacZa');
    expect([lac.strand, lac.length]).toEqual(['-', 324]);
    const protein = translate(featureLetters(letters, lac));
    expect(protein.startsWith(LacZaPrefix)).toBe(true);
    expect(protein.indexOf('*')).toBe(protein.length - 1);
  });

  it('EcoRI and HindIII each cut once, in the MCS, which spans from the one to the other', () => {
    const sites = findSites(letters);
    const once = (name: string) => sites.filter((s) => s.enzyme === name);
    expect(once('EcoRI').length).toBe(1);
    expect(once('HindIII').length).toBe(1);
    const mcs = feature('mcs');
    expect(mcs.start).toBe(once('EcoRI')[0].start);
    expect(mcs.start + mcs.length).toBe(once('HindIII')[0].start + 6);
    for (const e of Enzymes.filter((x) => x.mcs)) {
      const s = once(e.name);
      expect(s.length, e.name).toBe(1);
      expect(s[0].start).toBeGreaterThanOrEqual(mcs.start);
      expect(s[0].start + s[0].length).toBeLessThanOrEqual(mcs.start + mcs.length);
    }
    for (const e of Enzymes) expect(once(e.name).length, e.name).toBe(1);
  });

  it('the origin is the 589 bp pMB1 stretch on the minus strand, between lacZα and bla', () => {
    const ori = feature('ori');
    expect([ori.strand, ori.length, ori.type]).toEqual(['-', 589, 'rep_origin']);
    expect(ori.start).toBeGreaterThan(feature('lacZa').start + feature('lacZa').length);
    expect(ori.start + ori.length).toBeLessThanOrEqual(feature('bla').start);
  });

  it('none of the published record’s features crosses its origin; the rotated variant’s bla does, and its primer site', () => {
    for (const f of features) expect(f.start + f.length, f.id).toBeLessThanOrEqual(PlasmidLength);
    const rotated = rotatedAnnotations(letters);
    const crossing = rotated.features
      .filter((f) => f.start + f.length > PlasmidLength)
      .map((f) => `${f.id}${f.strand}`);
    expect(crossing.sort()).toEqual(['bla-', 'primer+']);
  });
});

describe('every derived feature keeps its length under rotation', () => {
  const rc = reverseComplement(letters);
  /** The record at `origin` (or its reverse complement): features found in it are the published ones, moved. */
  function check(origin: number, reverse: boolean) {
    const source = reverse ? rc : letters;
    const turned = rotateLetters(source, origin);
    const found = deriveFeatures(turned);
    expect(found.map((f) => f.id).sort()).toEqual(features.map((f) => f.id).sort());
    for (const f of found) {
      const base = feature(f.id);
      expect(f.length, `${f.id} at ${origin}`).toBe(base.length);
      // the MCS is drawn on the plus strand of whatever the record's orientation is
      if (f.id !== 'mcs') expect(f.strand).toBe(reverse ? (base.strand === '+' ? '-' : '+') : base.strand);
      // the feature's own letters, read on its strand, are the published ones
      // (the MCS has no strand of its own: read on the plus strand, the reverse complement's is the record's complement)
      const own = featureLetters(letters, base);
      expect(featureLetters(turned, f)).toBe(reverse && f.id === 'mcs' ? reverseComplement(own) : own);
      if (!reverse) expect(f.start).toBe(rotateFeature(base, origin).start);
    }
    // the sites: each one's enzyme, start and cut where the published one's are moved to the new origin; on the reverse
    // complement the same palindromic site is on the same bases, so its start maps through the reflection
    const n = PlasmidLength;
    const sites = findSites(turned);
    const base = findSites(letters);
    expect(sites.length).toBe(base.length);
    for (const b of base) {
      const e = Enzymes.find((x) => x.name === b.enzyme)!;
      const found = sites.filter((s) => s.enzyme === b.enzyme);
      expect(found.length, b.enzyme).toBe(1);
      const start = reverse ? rotate(n - b.start - e.site.length, origin) : rotate(b.start, origin);
      expect(found[0].start, `${b.enzyme} at ${origin}`).toBe(start);
      // a palindromic site reads the same on the reverse complement's top strand: the cut is `e.cut` bases in
      expect(found[0].cut).toBe((start + e.cut) % n);
      if (!reverse) expect(found[0].cut).toBe(rotate(b.cut, origin));
    }
  }

  it('a sample of offsets and of the offsets next to each feature’s ends, on the record and its reverse complement', () => {
    const offsets = new Set<number>();
    for (let o = 0; o < PlasmidLength; o += 211) offsets.add(o);
    for (const f of features)
      for (const p of [f.start, f.start + f.length]) for (const d of [-1, 0, 1]) offsets.add(rotate(p + d, 0));
    for (const o of offsets) {
      check(o, false);
      check(o, true);
    }
  });

  it('all 2,686 offsets, on the record and its reverse complement', { tags: ['medium'] }, () => {
    for (let o = 0; o < PlasmidLength; o++) {
      check(o, false);
      check(o, true);
    }
  });

  it('every enzyme’s cut offset is its published one (a hand table, the top strand’s cut)', () => {
    // NEB: the base the top strand is cut after, counted from the site's 5' end
    const hand: Record<string, number> = {
      EcoRI: 1,
      SacI: 5,
      KpnI: 5,
      SmaI: 3,
      BamHI: 1,
      XbaI: 1,
      SalI: 1,
      PstI: 5,
      SphI: 5,
      HindIII: 1,
      NdeI: 2,
      ScaI: 3,
      AatII: 5,
    };
    expect(Object.fromEntries(Enzymes.map((e) => [e.name, e.cut]))).toEqual(hand);
  });

  it('a site across the origin is found, once, with its cut wrapped, for every enzyme and every cut of the site', () => {
    const n = PlasmidLength;
    for (const e of Enzymes) {
      const at = findSites(letters).find((s) => s.enzyme === e.name)!.start;
      for (let inside = 1; inside < e.site.length; inside++) {
        const turned = rotateLetters(letters, at + inside);
        const found = findSites(turned).filter((s) => s.enzyme === e.name);
        expect(found.length, `${e.name} ${inside}`).toBe(1);
        expect(found[0].length).toBe(e.site.length);
        expect(found[0].start).toBe(n - inside);
        expect(found[0].start + found[0].length).toBeGreaterThan(n);
        // the cut wraps (start + cut) % n: before the origin when the cut is inside the first `inside` bases
        expect(found[0].cut).toBe((n - inside + e.cut) % n);
        expect(found[0].cut).toBe(e.cut >= inside ? e.cut - inside : n - inside + e.cut);
      }
    }
  });

  it('the complement handles the IUPAC codes and a translation an ambiguous codon', () => {
    expect(reverseComplement('ACGTRYKMBVDHSWN')).toBe('NWSDHBVKMRYACGT');
    expect(translate('ATGNNNRAAGGG')).toBe('MXXG');
  });
});

/** The marks of a ranges layer by item: each mark's tip (1 on the left, 2 on the right) and its x extent. */
function tips(lay: RangesLay) {
  const out = new Map<string, { tip: number; x0: number; x1: number }[]>();
  for (let j = 0; j < lay.m; j++) {
    const key = String(lay.keys[lay.mitem[j]]);
    out.set(key, [...(out.get(key) ?? []), { tip: lay.mtip[j], x0: lay.mx0[j], x1: lay.mx1[j] }]);
  }
  return out;
}

describe('the Plasmid views drawn against the data', () => {
  const size = { width: 720, height: 900 };

  const rendersInBothThemes = async (views: typeof PlasmidViews) => {
    for (const v of views) {
      const ctx = context();
      const out: string[] = [];
      for (const theme of ['light', 'dark'] as const) {
        const r = await ctx.renderSvg(v.build(), { size, tokens: theme === 'dark' ? DarkTokens : LightTokens });
        expect(r.issues, `${v.id} ${theme}`).toEqual([]);
        expect(r.settled).toBe(true);
        out.push(r.svg);
      }
      expect(out[1], v.id).not.toBe(out[0]);
      ctx.retire();
    }
  };

  // the ring view is the quick representative of the views, every one is the medium lane's
  it('the ring view renders settled in both themes with no issue; dark is not light', () =>
    rendersInBothThemes(PlasmidViews.filter((v) => v.id === 'ring')));

  it('every view renders settled in both themes with no issue; dark is not light', { tags: ['medium'] }, () =>
    rendersInBothThemes(PlasmidViews),
  );

  /** The ring's keys under a position of a track (a hit through the warp). */
  function keysAt(
    layout: Awaited<ReturnType<Awaited<ReturnType<typeof settled>>['frame']>>['layout'],
    track: number,
    pos: number,
  ) {
    const p = layout.sections[0].placement as RingPlacement;
    const keys = new Set<string>();
    for (let y = p.strips.bandY[track] + 1; y < p.strips.bandY[track] + p.strips.bandH[track]; y += 2) {
      p.warp.map(p.frameX(pos + 0.5), y);
      const h = hitView(layout, p.warp.px, p.warp.py);
      if (h?.element) keys.add(String(h.element.key));
    }
    return [...keys].sort();
  }

  it('the ring draws each feature where its positions are, across the origin one item; one label a feature', async () => {
    const s = await settled(view('ring').build());
    const f = await s.frame({ size });
    const p = f.layout.sections[0].placement as RingPlacement;
    const bla = rotateFeature(feature('bla'), RotatedOrigin);
    // bla on the minus ring (1) both sides of the rotated origin: the same item
    expect(keysAt(f.layout, 1, bla.start + 10)).toContain('bla');
    expect(keysAt(f.layout, 1, 0)).toContain('bla');
    expect(keysAt(f.layout, 1, bla.start + bla.length - PlasmidLength - 10)).toContain('bla');
    // the primer site on the plus ring (0) across it
    expect(keysAt(f.layout, 0, PlasmidLength - 5)).toContain('primer');
    expect(keysAt(f.layout, 0, 5)).toContain('primer');
    // the sites in the third ring
    const cut = findSites(letters).find((x) => x.enzyme === 'EcoRI')!.cut;
    expect(keysAt(f.layout, 2, rotate(cut, RotatedOrigin))).toContain('EcoRI');
    const labels = Array.from({ length: p.labels.m }, (_, k) => p.labels.text[p.labels.placed[k]]);
    for (const name of ['bla (AmpR)', 'primer site', 'lacZα', 'MCS', 'ori', 'EcoRI', 'HindIII'])
      expect(labels.filter((t) => t === name).length, name).toBe(1);
    expect([p.title, p.subtitle]).toEqual([PlasmidName, '2,686 bp']);
  });

  // a ring of 2,686 positions shows its features as items (and so arrows) from a device pixel a position: 1,200 CSS px
  // at the device ratio 2 (below it a ring draws bars, which have no ends: reported as a gap of the layer)
  const Big = { size: { width: 1200, height: 1200 }, dpr: 2 } as const;

  it('features are arrows by strand: the plus strand’s at the end, the minus strand’s at the start, on the piece holding it', async () => {
    const s = await settled(view('ring').build());
    // the ring's window over the origin: the detail's two pieces
    const f = await s.frame({ ...Big, viewports: { '/views/main/ring': { start: 2200, end: 3200 } } });
    const bla = rotateFeature(feature('bla'), RotatedOrigin);
    const p = f.layout.sections[0].placement as RingPlacement;
    for (const section of [0, 1]) {
      const sec = f.layout.sections[section];
      const plus = tips(sec.tracks[0].layers[0].lay as RangesLay);
      const minus = tips(sec.tracks[1].layers[0].lay as RangesLay);
      // bla (minus): two pieces, the arrow on the one holding its start, pointing left (to the start)
      expect(
        minus
          .get('bla')!
          .map((m) => m.tip)
          .sort(),
        `section ${section}`,
      ).toEqual([0, 1]);
      // the primer site (plus): two pieces, the arrow on the one holding its end, pointing right
      expect(
        plus
          .get('primer')!
          .map((m) => m.tip)
          .sort(),
        `section ${section}`,
      ).toEqual([0, 2]);
      if (section === 0) {
        // a feature that does not cross the origin has its arrow once: lacZα, the origin (minus) at their starts, the MCS at its end
        expect(minus.get('lacZa')!.map((m) => m.tip)).toEqual([1]);
        expect(minus.get('ori')!.map((m) => m.tip)).toEqual([1]);
        expect(plus.get('mcs')!.map((m) => m.tip)).toEqual([2]);
        // on the ring the strip's x is the position: the minus arrow's piece begins at bla's start, the plus one ends at the primer's end
        const tipped = minus.get('bla')!.find((m) => m.tip)!;
        expect(tipped.x0).toBeCloseTo(p.frameX(bla.start), 0);
        const end = plus.get('primer')!.find((m) => m.tip)!;
        expect(end.x1).toBeCloseTo(p.frameX(25), 0);
      }
    }
  });

  it('the arrowhead is bent through the ring’s warp: its tip at the start of bla, at the ring band’s middle', async () => {
    const s = await settled(view('ring').build());
    const f = await s.frame(Big);
    const p = f.layout.sections[0].placement as RingPlacement;
    const list = recordDrawLists(f.lists, ['data']);
    const vertices: [number, number][] = [];
    for (const line of list.split('\n'))
      if (line.startsWith('path ')) {
        const nums = line.slice(5).split('#')[0].replace(/[MLZ]/g, ' ').trim().split(/\s+/).map(Number);
        for (let i = 0; i + 1 < nums.length; i += 2) vertices.push([nums[i], nums[i + 1]]);
      }
    expect(vertices.length).toBeGreaterThan(200);
    const near = (pos: number, track: number) => {
      // the tip is at the middle of the mark's lane (the features are in the band's first lane)
      const lay = f.layout.sections[0].tracks[track].layers[0].lay as RangesLay;
      p.warp.map(p.frameX(pos), lay.box.y + lay.laneHeight / 2);
      const [x, y] = [p.warp.px, p.warp.py];
      return vertices.some(([vx, vy]) => Math.hypot(vx - x, vy - y) < 1.5);
    };
    const bla = rotateFeature(feature('bla'), RotatedOrigin);
    // minus strand: the tip at bla's start, not at its end (the piece after the origin has no tip)
    expect(near(bla.start, 1)).toBe(true);
    expect(near(bla.start + bla.length - PlasmidLength, 1)).toBe(false);
    // plus strand: the tip at the primer site's end, past the origin, not at its start
    expect(near(25, 0)).toBe(true);
    expect(near(PlasmidLength - 30, 0)).toBe(false);
  });

  it('the ring’s record: its window over the origin, its detail in two pieces', async () => {
    const viewports = { '/views/main/ring': { start: 2600, end: 2780 } };
    const s = await settled(view('ring').build());
    const f = await s.frame({ size: { width: 720, height: 900 }, viewports });
    expect(f.layout.sections[1].axis!.pieces).toEqual([2600, 2686, 0, 94]);
    const ctx = context();
    const r = await ctx.renderSvg(view('ring').build(), { size: { width: 520, height: 720 }, viewports, window: true });
    expect(r.issues).toEqual([]);
    expect(r.settled).toBe(true);
    compare(fileURLToPath(new URL('./records/plasmid-ring.txt', import.meta.url)), recordDrawLists(r.lists));
    ctx.retire();
  });

  it('the crossing view bookmarks bla and the primer site, each an item of its strand’s ranges', async () => {
    const ctx = context();
    const spec = view('crossing').build();
    expect(ctx.check(spec)).toEqual([]);
    const refs = spec.children?.find((c) => c.id === 'refs')?.children?.map((c) => [c.id, c.params?.to]);
    expect(refs).toEqual([
      ['bla', '/data/minus/(item :id bla)'],
      ['primer', '/data/plus/(item :id primer)'],
    ]);
    // before its file is open the rows across the origin wait for the length, and say so; once open, none is left out
    const model = ctx.resolve(spec);
    expect([...new Set(model.issues.map((i) => i.code))]).toEqual(['seq.resolve.length.unknown']);
    const s = await settled(spec);
    expect(s.renderer.options.model.issues).toEqual([]);
    expect(s.renderer.options.model.annotations.get('/data/minus')?.items?.filter((i) => i.key === 'bla')).toHaveLength(
      2,
    );
    ctx.retire();
  });

  it('the rotated linear view: two pieces meeting at the seam, no break, the features whole', async () => {
    const s = await settled(view('rotated-linear').build());
    // 1,386 positions in view: a device pixel each from 1,400 CSS px (below it the features are bars: no arrows)
    const f = await s.frame({ size: { width: 1800, height: 400 } });
    const sec = f.layout.sections[0];
    expect(sec.axis!.pieces).toEqual([2000, 2686, 0, 700]);
    expect(sec.breaks).toEqual([]);
    expect(sec.axis!.total).toBe(686 + 700);
    // the published record: bla on the minus strand is whole (1625 to 2486 is cut at the range's start, 2000: its arrow
    // at its start is not on the cut piece), lacZα whole
    const minus = tips(sec.tracks[1].layers[0].lay as RangesLay);
    expect(minus.get('lacZa')!.map((m) => m.tip)).toEqual([1]);
    expect(minus.get('bla')!.map((m) => m.tip)).toEqual([0]);
  });

  it('the rotated linear view’s start is a section param: a patch moves it, its inverse puts it back', async () => {
    const spec = view('rotated-linear').build();
    const at = ['views', 'main', 'linear'];
    const patch = { ops: [{ op: 'set' as const, at, key: 'range', value: [1000, 1700] }] };
    const r = applyPatch(spec, patch);
    expect(r.ok).toBe(true);
    if (!r.ok) return;
    const s = await settled(r.root);
    const f = await s.frame({ size });
    expect(f.layout.sections[0].axis!.pieces).toEqual([1000, 1700]);
    const back = applyPatch(r.root, r.inverse!);
    expect(back.ok && back.root).toEqual(spec);
  });

  it('the sites layer’s enzymes are a filter on its data: a patch shows two, its inverse all', async () => {
    const spec = view('ring').build();
    const sites = ['views', 'main', 'ring', 'sites', 'sites'];
    const shown = ['EcoRI', 'HindIII'];
    const data = `/data/sites${printAddress([{ kind: 'item', params: { id: shown } }])}`;
    const r = applyPatch(spec, { ops: [{ op: 'set', at: sites, key: 'data', value: data }] });
    expect(r.ok).toBe(true);
    if (!r.ok) return;
    const labelsOf = async (root: typeof spec) => {
      const s = await settled(root);
      const f = await s.frame({ size });
      const p = f.layout.sections[0].placement as RingPlacement;
      return Array.from({ length: p.labels.m }, (_, k) => p.labels.text[p.labels.placed[k]]);
    };
    const some = await labelsOf(r.root);
    for (const e of shown) expect(some).toContain(e);
    expect(some).not.toContain('SacI');
    const back = applyPatch(r.root, r.inverse!);
    expect(back.ok).toBe(true);
    if (back.ok) expect(await labelsOf(back.root)).toContain('SacI');
  });

  it('hiding a strand’s track is a patch, its inverse brings it back', async () => {
    const spec = view('ring').build();
    const minus = ['views', 'main', 'ring', 'minus'];
    const r = applyPatch(spec, { ops: [{ op: 'set', at: minus, key: 'hidden', value: true }] });
    expect(r.ok).toBe(true);
    if (!r.ok) return;
    const labelsOf = async (root: typeof spec) => {
      const s = await settled(root);
      const p = (await s.frame({ size })).layout.sections[0].placement as RingPlacement;
      return Array.from({ length: p.labels.m }, (_, k) => p.labels.text[p.labels.placed[k]]);
    };
    expect(await labelsOf(r.root)).not.toContain('bla (AmpR)');
    const back = applyPatch(r.root, r.inverse!);
    if (back.ok) expect(await labelsOf(back.root)).toContain('bla (AmpR)');
  });

  it('the Files view’s GFF3 holds the same features as the instance: the rotated one, bla ending past the length', () => {
    const text = readFileSync(`${PackageFolder}${PlasmidFiles.rotatedFeatures}`, 'utf8');
    const bla = rotateFeature(feature('bla'), RotatedOrigin);
    expect(text).toContain(`CDS\t${bla.start + 1}\t${bla.start + bla.length}\t.\t-\t.\tID=bla;Name=bla (AmpR)`);
    expect(bla.start + bla.length).toBeGreaterThan(PlasmidLength);
    expect(seqLetters(PlasmidFiles.published)).toBe(letters);
  });
});

describe('the Plasmid at scale', () => {
  it('the generated 200 kb variant: 2,000 features on both strands, a tenth across the origin, from a seed', () => {
    const a = scaleData();
    expect(a.letters.length).toBe(ScaleLength);
    expect(a.features.start.length).toBe(ScaleFeatures);
    expect(a.features.strand.filter((s) => s === '+').length).toBeGreaterThan(800);
    expect(a.features.strand.filter((s) => s === '-').length).toBeGreaterThan(800);
    // across the origin: a range that ends before it starts
    expect(a.features.start.filter((s, i) => s > a.features.end[i]).length).toBe(ScaleFeatures / 10);
    expect(scaleData().letters).toBe(a.letters);
    expect(scaleData(7).letters).not.toBe(a.letters);
  });

  const Size = { width: 1200, height: 1200 };

  // the pan's quick half would take a fifth of a second, so its unbounded run is the medium lane's (itWithin's quick half)
  const panAtScale = async (within: Within) => {
    const s = await settled(view('scale').build());
    const inputs: FrameInput[] = [];
    for (let k = 0; k < 8; k++) {
      const start = 199_000 + k * 400;
      inputs.push({ size: Size, dpr: 2, viewports: { '/views/main/ring': { start, end: start + 3000 } } });
    }
    for (const i of inputs) await s.frame(i);
    const ring = s.renderer.frame(inputs[0]).layout.sections[0];
    const time = () => {
      const t: number[] = [];
      for (let k = 0; k < 40; k++) {
        const t0 = performance.now();
        const f = s.renderer.frame(inputs[(k + 1) % inputs.length]);
        t.push(performance.now() - t0);
        if (!f.rebuilt || f.layout.sections[0] !== ring) throw new Error('a pan lays the detail out, not the ring');
      }
      t.sort((a, b) => a - b);
      return t[t.length >> 1];
    };
    time();
    const median = within(8 * 40, time);
    console.info(`plasmid ring at scale, a pan of its detail: median ${median.toFixed(3)} ms`);
  };

  it('a pan of the detail at scale: the ring drawn from its kept lays', { tags: ['medium'] }, () =>
    panAtScale((_, f) => f()),
  );

  it(
    'a pan of the detail at scale: the ring drawn from its kept lays: within its time bounds',
    { tags: ['heavy'], retry: BoundedTries - 1 },
    () =>
      panAtScale((ms, f) => {
        const t0 = performance.now();
        const r = f();
        expect(performance.now() - t0, `under ${ms} ms`).toBeLessThan(ms);
        return r;
      }),
  );

  itWithin('a hover over the ring at scale: the placement through the warp, the hit', async (within) => {
    const s = await settled(view('scale').build());
    const input: FrameInput = { size: Size, dpr: 2, viewports: { '/views/main/ring': { start: 100, end: 4000 } } };
    const layout = (await s.frame(input)).layout;
    const p = layout.sections[0].placement as RingPlacement;
    const points: { x: number; y: number }[] = [];
    for (let u = 0; u < p.total; u += 997) {
      p.warp.map(p.frameX(u + 0.5), p.strips.bandY[0] + 3);
      points.push({ x: p.warp.px, y: p.warp.py });
    }
    const hover = () => {
      const t: number[] = [];
      for (let k = 0; k < 200; k++) {
        const q = points[k % points.length];
        const t0 = performance.now();
        hitView(layout, q.x, q.y);
        t.push(performance.now() - t0);
      }
      t.sort((a, b) => a - b);
      return t[t.length >> 1];
    };
    hover();
    const median = within(2 * 200, hover);
    console.info(`plasmid ring at scale, hover: median ${median.toFixed(3)} ms`);
    expect(points.length).toBeGreaterThan(100);
  });
});

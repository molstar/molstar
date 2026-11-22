/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BuiltInLayerTypes, BuiltInLayouts } from '@molstar/seqx-core/render/built-ins';
import { existsSync, readFileSync, statSync, writeFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { createSeqContext } from '@molstar/seqx-core/context/default';
import { DataStore, type RangesList, SeqRuntime } from '@molstar/seqx-core/data';
import { RecordsRuntime, readRecords } from '@molstar/seqx-core/data/records';
import {
  bundledTextMetrics,
  hitView,
  LightTokens,
  layoutView,
  type RingPlacement,
  recordDrawLists,
  openModel,
} from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import {
  Cdk2,
  DomainFiles,
  type DomainFile,
  FilesFolder,
  manyContent,
  manyRows,
  ManySequencesFile,
  PlasmidGffFile,
  PlasmidSeqId,
} from '../files/data.js';
import { EntryFile, EntryTracks } from '../files/p53-entry.js';
import * as plasmidGff from '../files/plasmid-gff.js';
import { uniprotTransform } from '../files/uniprot.js';
import { nodeFixtureAssets, PackageFolder } from '../node-assets.js';
import { P53 } from '../p53/data.js';
import { deriveFeatures, PlasmidLength, RotatedOrigin, rotate } from '../plasmid/derive.js';
import { publishedLetters } from './plasmid-helpers.js';
import { FilesViews, TableFiles } from '../views.js';

// The Files dataset (files.md): the committed files are what their generators write, p53's UniProt entry and its
// variants are recorded files listed in `sources.json` (under the bounds, with an allowed licence), the entry read by
// the records runtime says what UniProt says (393 residues; R175, R248 and R273 inside the DNA-binding domain), every
// form of the domains table draws the same, and the plasmid's GFF3 draws a feature across the origin as one item.

const update = process.env.SEQ_RECORDS === 'update';
const dir = (file: string) => `${PackageFolder}${FilesFolder}${file}`;

/** A committed file is what its generator writes (`SEQ_RECORDS=update` writes it again). */
function compare(file: string, content: string | Uint8Array) {
  if (update) writeFileSync(file, content);
  expect(existsSync(file), `no ${file}: run with SEQ_RECORDS=update`).toBe(true);
  const text = readFileSync(file);
  expect(Buffer.from(content).equals(text)).toBe(true);
}

const sha256 = (b: Uint8Array | string) => createHash('sha256').update(b).digest('hex');
const context = () =>
  createSeqContext({
    registry: [
      { kind: 'format-runtime', runtime: RecordsRuntime },
      { kind: 'format-runtime', runtime: SeqRuntime },
    ],
    ports: { assets: () => nodeFixtureAssets() },
  });

describe('the fixtures', () => {
  it.each(Object.keys(DomainFiles) as DomainFile[])('%s is what its generator writes', (file) => {
    const content = DomainFiles[file].content();
    compare(dir(file), typeof content === 'string' ? content : content);
  });

  it('many-sequences.tsv is what its generator writes, listed in sources.json as generated, with its checksum', () => {
    compare(dir(ManySequencesFile), manyContent());
    const rows = JSON.parse(readFileSync(`${PackageFolder}${FilesFolder}sources.json`, 'utf8')) as {
      file: string;
      licence: string[];
      sha256Fixture: string;
      bytesFixture: number;
    }[];
    expect(rows.map((r) => r.file)).toEqual([`files/${ManySequencesFile}`]);
    const bytes = readFileSync(dir(ManySequencesFile));
    expect(rows[0].licence).toEqual(['LicenseRef-seqx-generated']);
    expect(rows[0].bytesFixture).toBe(bytes.byteLength);
    expect(rows[0].sha256Fixture).toBe(sha256(bytes));
  });

  it('every recorded file is listed in sources.json, under the bounds, with an allowed licence and its checksum', () => {
    const rows = JSON.parse(readFileSync(`${PackageFolder}fixtures/p53/sources.json`, 'utf8')) as {
      file: string;
      url: string;
      recorded: string;
      licence: string[];
      attribution: string;
      sha256Fixture: string;
      bytesFixture: number;
    }[];
    expect(rows.map((r) => r.file).sort()).toEqual([
      'p53/2OCJ-chains.json',
      'p53/2OCJ.seq',
      'p53/P04637-variation.json',
      'p53/P04637.json',
    ]);
    let total = 0;
    for (const r of rows) {
      const bytes = readFileSync(`${PackageFolder}fixtures/${r.file}`);
      expect(bytes.byteLength).toBe(r.bytesFixture);
      expect(sha256(bytes)).toBe(r.sha256Fixture);
      expect(r.bytesFixture).toBeLessThan(1024 * 1024);
      expect(r.url).toMatch(/^https:\/\//);
      expect(r.recorded).toMatch(/^\d{4}-\d{2}-\d{2}$/);
      expect(r.attribution).not.toBe('');
      for (const l of r.licence) expect(l).toMatch(/^(CC0-1\.0|CC-BY-4\.0|LicenseRef-NCBI-public-data)/);
      total += r.bytesFixture;
    }
    // the dataset's files together: under 2 MB with the seq files beside them
    const seq =
      statSync(`${PackageFolder}fixtures/p53/gene-region.seq`).size +
      statSync(`${PackageFolder}fixtures/p53/tp53-gene.seq`).size;
    expect(total + seq).toBeLessThan(2 * 1024 * 1024);
  });
});

describe('p53’s UniProt entry', () => {
  const text = readFileSync(`${PackageFolder}${EntryFile}`, 'utf8');
  const read = readRecords({ format: 'json', text }, uniprotTransform({ accession: 'P04637' }), {
    system: 'entry',
    limits: { maxRows: 1_000_000, maxBytes: 64 * 1024 * 1024 },
  });

  it('reads 393 residues, the protein as UniProt writes it', () => {
    expect(read.issues).toEqual([]);
    const own = read.data!.pairs.find((p) => p.request.system === 'P04637')!;
    const tables = (own.answer as { tables: { letters: { letters: string } } }).tables;
    expect(tables.letters.letters.length).toBe(393);
    expect(tables.letters.letters).toBe(P53);
  });

  it('puts R175, R248 and R273 inside the DNA-binding domain', () => {
    const own = read.data!.pairs.find((p) => p.request.system === 'P04637')!;
    const ranges = (own.answer as { tables: { ranges: RangesList } }).tables.ranges;
    const category = ranges.fields!.find((f) => f.name === 'category')!.column as { get(i: number): string | null };
    const rows = Array.from(ranges.start, (_, i) => i).filter((i) => category.get(i) === 'DNA binding');
    expect(rows.length).toBe(1);
    // UniProt counts from 1 and names the last position: [102, 292] is [101, 292) here
    expect([ranges.start[rows[0]], ranges.end[rows[0]]]).toEqual([101, 292]);
    for (const residue of [175, 248, 273]) {
      expect(P53[residue - 1]).toBe('R');
      expect(residue - 1).toBeGreaterThanOrEqual(ranges.start[rows[0]]);
      expect(residue - 1).toBeLessThan(ranges.end[rows[0]]);
    }
  });

  it('keeps the ends’ modifiers as fields; every track’s type has features in the entry', () => {
    const own = read.data!.pairs.find((p) => p.request.system === 'P04637')!;
    const ranges = (own.answer as { tables: { ranges: RangesList } }).tables.ranges;
    expect(ranges.fields!.map((f) => f.name)).toEqual(expect.arrayContaining(['start_modifier', 'end_modifier']));
    const category = ranges.fields!.find((f) => f.name === 'category')!.column as { get(i: number): string | null };
    const types = new Set(Array.from(ranges.start, (_, i) => category.get(i)));
    for (const [type] of EntryTracks) expect(types.has(type), type).toBe(true);
  });

  it('the Proteins API’s variants are single-residue changes inside 1 to 393', () => {
    const v = JSON.parse(readFileSync(`${PackageFolder}fixtures/p53/P04637-variation.json`, 'utf8')) as {
      sequence: string;
      features: { begin: string; end: string; wildType: string; xrefs: { name: string }[] }[];
    };
    expect(v.sequence).toBe(P53);
    expect(v.features.length).toBeGreaterThan(1000);
    let misfits = 0;
    for (const f of v.features) {
      expect(f.begin).toBe(f.end);
      expect(Number(f.begin)).toBeGreaterThanOrEqual(1);
      expect(Number(f.begin)).toBeLessThanOrEqual(393);
      if (f.wildType !== P53[Number(f.begin) - 1]) misfits++;
      for (const x of f.xrefs) expect(['UniProt', 'ClinVar']).toContain(x.name);
    }
    // as recorded: one ClinVar record (position 281) names another wild type than the sequence has; kept as published
    expect(misfits).toBeLessThanOrEqual(1);
  });
});

describe('the domains table in every file form', () => {
  it('draws the same domains from TSV, CSV, JSON and MessagePack, in both layouts', async () => {
    const records: string[] = [];
    for (const [file, id] of TableFiles) {
      const view = FilesViews.find((v) => v.id === id)!;
      const ctx = context();
      const r = await ctx.renderSvg(view.build(), { size: { width: 720, height: 240 } });
      expect(r.issues, file).toEqual([]);
      expect(r.settled, file).toBe(true);
      records.push(recordDrawLists(r.lists));
      ctx.retire();
    }
    expect(records.length).toBe(6);
    for (const r of records) expect(r).toBe(records[0]);
    // the domains are drawn: their labels' items are in the record
    expect(records[0].length).toBeGreaterThan(500);
  });
});

describe('one table, two sequences', () => {
  const text = readFileSync(dir(ManySequencesFile), 'utf8');
  const read = readRecords(
    { format: 'tsv', text },
    {},
    {
      system: 'p53',
      limits: { maxRows: 1_000_000, maxBytes: 64 * 1024 * 1024 },
    },
  );
  const rowsOf = (system: string) => {
    const pair = read.data!.pairs.find((p) => p.request.system === system);
    const ranges = (pair!.answer as { tables: { ranges: RangesList } }).tables.ranges;
    const id = ranges.fields!.find((f) => f.name === 'label')!.column as { get(i: number): string | null };
    return Array.from(ranges.start, (_, i) => [id.get(i), ranges.start[i], ranges.end[i]] as const);
  };
  /** The file's rows of a system as the reader holds them: positions from 0, the end exclusive. */
  const expected = (system: string, own = false) =>
    manyRows()
      .filter((r) => r.target === system || (own && r.target === ''))
      .map((r) => [r.label, r.start - 1, r.end] as const)
      .sort((a, b) => a[1] - b[1]);

  it('names two systems, each holding exactly the rows that target it; the empty target is on the source’s own', () => {
    expect(read.issues).toEqual([]);
    expect(read.data!.pairs.map((p) => p.request.system).sort()).toEqual(['cdk2', 'p53']);
    expect(rowsOf('cdk2')).toEqual(expected('cdk2'));
    expect(rowsOf('p53')).toEqual(expected('p53', true));
    // the row with no target is NLS: on p53, the source's own system, and on no other
    expect(rowsOf('p53').map((r) => r[0])).toContain('Nuclear localization');
    expect(rowsOf('cdk2').map((r) => r[0])).not.toContain('Nuclear localization');
    expect(manyRows().filter((r) => r.target === '').length).toBe(1);
    expect(rowsOf('p53').length + rowsOf('cdk2').length).toBe(manyRows().length);
  });

  it('every row fits its sequence', () => {
    for (const [, , end] of rowsOf('p53')) expect(end).toBeLessThanOrEqual(P53.length);
    for (const [, , end] of rowsOf('cdk2')) expect(end).toBeLessThanOrEqual(Cdk2.length);
  });

  it('draws each section with its own rows’ labels only, once each', async () => {
    const view = FilesViews.find((v) => v.id === 'many-sequences')!;
    const ctx = context();
    const r = await ctx.renderSvg(view.build(), { size: { width: 900, height: 360 } });
    expect(r.issues).toEqual([]);
    expect(r.settled).toBe(true);
    const count = (t: string) => r.svg.split(`>${t}<`).length - 1;
    // a row is drawn in its own section alone (a label at most once, never in both), and whole: a label that does not
    // fit its mark is not drawn (a file's rows are not known whole before they are read: no label lane beside them)
    let drawn = 0;
    for (const row of manyRows()) {
      const n = count(row.label as string);
      expect(n, row.id).toBeLessThanOrEqual(1);
      drawn += n;
    }
    expect(drawn).toBeGreaterThan(0);
    ctx.retire();
  });

  it('draws against the data: a feature’s box spans its positions at the fit', async () => {
    const view = FilesViews.find((v) => v.id === 'many-sequences')!;
    const spec = view.build();
    const store = new DataStore();
    const opened = await openModel(spec, store, {
      runtimes: [RecordsRuntime, SeqRuntime],
      assets: nodeFixtureAssets().reader,
    });
    const input = {
      types: BuiltInLayerTypes,
      layouts: BuiltInLayouts,
      model: opened.model,
      view: opened.view!,
      size: { width: 900, height: 360 },
      viewports: {},
      tokens: LightTokens,
      text: bundledTextMetrics(),
      data: store,
    };
    for (let i = 0; i < 4; i++) {
      store.want(layoutView(input).wants);
      await store.whenSettled();
    }
    const layout = layoutView(input);
    expect(layout.sections.length).toBe(2);
    // what is under a position of a section, over the page's height: its system and the items' keys
    const keysAt = (section: number, pos: number) => {
      const system = ['/data/p53', '/data/cdk2'][section];
      const st = layout.sections[section].placement.strips;
      const keys = new Set<string>();
      for (let y = 0; y < input.size.height; y += 2) {
        const h = hitView(layout, st.axes[0].xOfUnit(pos + 0.5 - st.first[0]), y);
        if (h?.element && h.system === system) keys.add(`${h.system}:${String(h.element.key)}`);
      }
      return [...keys].sort();
    };
    // p53: its fitted window starts at 1, the DBD lies over 94..292 and the empty-target NLS over 319..360
    expect(keysAt(0, 150)).toEqual(['/data/p53:DBD']);
    expect(keysAt(0, 319)).toEqual(['/data/p53:NLS']);
    expect(keysAt(0, 319).join()).not.toMatch(/cdk2/);
    // CDK2: the T-loop at 145..172, nothing of p53's at the same positions
    expect(keysAt(1, 160)).toEqual(['/data/cdk2:TLOOP']);
    expect(keysAt(1, 280)).toEqual([]);
  });
});

describe('the plasmid’s GFF3', () => {
  async function laid() {
    const spec = plasmidGff.build();
    const store = new DataStore();
    // the record's file first (the sequence's length), then the model resolved again with it and the rest opened
    const opened = await openModel(spec, store, {
      runtimes: [RecordsRuntime, SeqRuntime],
      assets: nodeFixtureAssets().reader,
    });
    const { model } = opened;
    expect(model.systems.get('/data/plasmid')?.length).toBe(PlasmidLength);
    const view = opened.view!;
    const input = {
      types: BuiltInLayerTypes,
      layouts: BuiltInLayouts,
      model,
      view,
      size: { width: 720, height: 900 },
      viewports: {},
      tokens: LightTokens,
      text: bundledTextMetrics(),
      data: store,
    };
    for (let i = 0; i < 4; i++) {
      store.want(layoutView(input).wants);
      await store.whenSettled();
    }
    return layoutView(input);
  }

  it('every declared track of the story has at least one item (no track is a blank well)', async () => {
    const rows = readFileSync(`${PackageFolder}${PlasmidGffFile}`, 'utf8')
      .split('\n')
      .filter((l) => l && !l.startsWith('#'))
      .map((l) => l.split('\t')[2]);
    for (const [type] of plasmidGff.GffTracks) expect(rows, type).toContain(type);
    // and drawn: each ring track's layer has marks
    const ring = (await laid()).sections[0];
    const marks = ring.tracks.map((t) => t.layers.map((l) => (l.lay as { m?: number } | undefined)?.m));
    for (let j = 0; j < plasmidGff.GffTracks.length; j++)
      expect(
        marks[j].some((m) => m === undefined || m > 0),
        plasmidGff.GffTracks[j][0],
      ).toBe(true);
  });

  it('reads bla and the primer site, ending past the length, as one item across the origin each, in the ring', async () => {
    const layout = await laid();
    const ring = layout.sections[0];
    const p = ring.placement as RingPlacement;
    const at = (track: number, pos: number) => {
      const keys = new Set<string>();
      for (let y = p.strips.bandY[track] + 1; y < p.strips.bandY[track] + p.strips.bandH[track]; y += 2) {
        p.warp.map(p.frameX(pos + 0.5), y);
        const h = hitView(layout, p.warp.px, p.warp.py);
        if (h?.element) keys.add(String(h.element.key));
      }
      return [...keys].sort();
    };
    const bla = deriveFeatures(publishedLetters()).find((f) => f.id === 'bla')!;
    // the CDS track (0): bla on both sides of the rotated origin, the same item
    expect(at(0, rotate(bla.start + 10, RotatedOrigin))).toContain('bla');
    expect(at(0, rotate(bla.start + bla.length - 10, RotatedOrigin))).toContain('bla');
    expect(at(0, 0)).toContain('bla');
    // the primer site (the track after the four types before it) across the origin
    const primer = GffTrackIndex('primer_bind');
    expect(at(primer, PlasmidLength - 5)).toContain('primer');
    expect(at(primer, 5)).toContain('primer');
    // one label a feature
    const labels = Array.from({ length: p.labels.m }, (_, k) => p.labels.text[p.labels.placed[k]]);
    expect(labels.filter((t) => t === 'bla (AmpR)').length).toBe(1);
    expect(labels.filter((t) => t === 'primer site').length).toBe(1);
    expect(p.title).toBe(PlasmidSeqId);
  });

  it('the SVG draws each across-the-origin feature once with its label, and the file’s region row nowhere', async () => {
    const ctx = context();
    const r = await ctx.renderSvg(plasmidGff.build(), { size: { width: 720, height: 900 } });
    expect(r.issues).toEqual([]);
    expect(r.settled).toBe(true);
    const count = (s: string) => r.svg.split(`>${s}<`).length - 1;
    // once on the ring (outside), once in the linear detail (inside one of its two pieces)
    expect(count('bla (AmpR)')).toBe(2);
    expect(count('primer site')).toBeLessThanOrEqual(2);
    // the region row (`Is_circular`) is no track: its Name is drawn nowhere
    expect(count(PlasmidSeqId)).toBeLessThanOrEqual(2);
    ctx.retire();
  });
});

function GffTrackIndex(type: string): number {
  return plasmidGff.GffTracks.findIndex(([t]) => t === type);
}

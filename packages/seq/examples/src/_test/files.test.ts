/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { existsSync, readFileSync, statSync, writeFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { createSeqContext } from '@molstar/seqx-core/context';
import { DataStore, type RangesList, SeqRuntime } from '@molstar/seqx-core/data';
import { RecordsRuntime, readRecords } from '@molstar/seqx-core/data/records';
import {
  attachSources,
  bundledTextMetrics,
  hitView,
  LightTokens,
  layoutView,
  type RingPlacement,
  recordDrawLists,
} from '@molstar/seqx-core/render';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import { DomainFiles, type DomainFile, FilesFolder, PlasmidSeqId } from '../files/data.js';
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
    expect(rows.map((r) => r.file).sort()).toEqual(['p53/P04637-variation.json', 'p53/P04637.json']);
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

describe('the plasmid’s GFF3', () => {
  async function laid() {
    const spec = plasmidGff.build();
    const model = resolveSpec(spec, { runtimes: [RecordsRuntime, SeqRuntime] });
    const view = resolveViews(model)[0];
    const store = new DataStore();
    await attachSources(model, store, { assets: nodeFixtureAssets().reader });
    const input = {
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
    expect(count('bla (AmpR)')).toBe(1);
    expect(count('primer site')).toBe(1);
    // the region row (`Is_circular`) is no track: its Name is drawn nowhere
    expect(count(PlasmidSeqId)).toBeLessThanOrEqual(2);
    ctx.retire();
  });
});

function GffTrackIndex(type: string): number {
  return plasmidGff.GffTracks.findIndex(([t]) => t === type);
}

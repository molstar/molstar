/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BuiltInLayerTypes, BuiltInLayouts } from '@molstar/seqx-core/render/built-ins';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import type { Node } from '@molstar/core/repr/tree';
import { DataStore, decodeSeq, type FormatRuntime } from '@molstar/seqx-core/data';
import { AlignmentFileRuntime, readNewick, readStockholm } from '@molstar/seqx-core/data/alignment-readers';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import {
  attachSources,
  bundledTextMetrics,
  type FrameInput,
  hitView,
  LightTokens,
  Renderer,
  type RowsLay,
} from '@molstar/seqx-core/render';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import {
  ColumnsFile,
  generatedLabel,
  generatedRow,
  generatedRuntime,
  HemoglobinFile,
  NearestBetaRow,
  ReferenceRow,
  SeedFile,
  SeedTree,
  seedProfile,
} from '../globins/data.js';
import * as globins from '../globins/family.js';
import { HemeContacts, HemeHistidines, HemeUnitContacts } from '../globins/hemes.js';
import { HemoglobinLengths } from '../globins/hemoglobin.js';
import { readCif } from '../../scripts/cif.js';
import { hemeContacts, identity, residuesOf } from '../globins/recording.js';
import { neighborJoining } from '../globins/tree.js';
import { nodeFixtureAssets, PackageFolder } from '../node-assets.js';

// The Globins dataset's facts checked against its fixtures (examples.md, Globins: Anchors), and its views drawn against
// their data (the examples' rule): at the fit and across a tile edge every row in the band has its band of cells and
// its label, and every cell at letters is hit as its row and its column.

const text = (f: string) => readFileSync(`${PackageFolder}${f}`, 'utf8');
const seed = readStockholm(text(SeedFile));
const profile = seedProfile(text(SeedFile));
const hb = decodeSeq(text(HemoglobinFile)).data!;
const lettersOf = (system: string) => {
  const t = hb.pairs.find((p) => p.request.system === system)!.answer;
  return (t.ok ? (t.tables.letters as { letters: string }).letters : '') as string;
};
const FamilyRuntimes = [AlignmentFileRuntime, RecordsRuntime];

describe('the Globins fixtures’ anchors', () => {
  it('the seed: 60 rows of 141 columns read with no issue, no human row; the reference the row nearest 4HHB’s α', () => {
    expect(seed.issues).toEqual([]);
    expect(seed.rows.length).toBe(60);
    expect(seed.columns).toBe(141);
    const labels = seed.rows.map((r) => r.label);
    expect(labels.some((l) => l.startsWith('HBA_HUMAN') || l.startsWith('HBB_HUMAN'))).toBe(false);
    const nearest = (chain: string) =>
      seed.rows
        .map((r) => [r.label, identity(residuesOf(r.letters), chain)] as const)
        .sort((a, b) => b[1] - a[1])[0][0];
    expect(nearest(lettersOf('hba'))).toBe(ReferenceRow);
    expect(nearest(lettersOf('hbb'))).toBe(NearestBetaRow);
  });

  it('the reference row numbers as 4HHB does: its entry is the mature chain (no Met), its histidines at 58 and 87', () => {
    const ref = seed.rows.find((r) => r.label === ReferenceRow)!;
    const residues = residuesOf(ref.letters);
    const at = (n: number) => residues[n - ref.first];
    // UniProt's yak α holds 141 residues, the mature chain: its numbers are 4HHB's (human α's P69905 counts the Met,
    // so its His59 and His88 would be these)
    const record = text(ColumnsFile)
      .split('\n')
      .find((l) => l.startsWith(`${ReferenceRow}\t`))!;
    expect(record.split('\t')[3]).toBe(String(HemoglobinLengths.hba));
    expect(at(HemeHistidines.hba.distal)).toBe('H');
    expect(at(HemeHistidines.hba.proximal)).toBe('H');
  });

  it('the tree is its generator’s (neighbor joining on identity), every row a leaf', () => {
    expect(
      neighborJoining(
        seed.rows.map((r) => r.label),
        seed.rows.map((r) => r.letters),
      ),
    ).toBe(SeedTree);
    const tree = readNewick(SeedTree);
    expect(tree.issues).toEqual([]);
    expect(tree.leaves.map((n) => tree.label[n]).sort()).toEqual(seed.rows.map((r) => r.label).sort());
  });

  it('4HHB: α 141 and β 146 residues, eight helices A to H a chain, a heme a position', () => {
    expect(lettersOf('hba').length).toBe(HemoglobinLengths.hba);
    expect(lettersOf('hbb').length).toBe(HemoglobinLengths.hbb);
    expect(lettersOf('hem').length).toBe(HemoglobinLengths.hem);
    for (const system of ['hba', 'hbb']) {
      const a = hb.pairs.find((p) => p.request.system === system)!.answer;
      const ranges = a.ok ? (a.tables.ranges as { id: { dict: readonly string[]; codes: Int32Array } }) : undefined;
      expect([...ranges!.id.codes].map((c) => ranges!.id.dict[c])).toEqual(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']);
    }
  });

  it('every heme contacts its chain’s proximal and distal histidines', () => {
    for (const c of ['hba', 'hbb'] as const) {
      const letters = lettersOf(c);
      for (const h of [HemeHistidines[c].distal, HemeHistidines[c].proximal]) {
        expect(letters[h - 1], `${c} ${h}`).toBe('H');
        expect(HemeContacts[c]).toContain(h);
      }
    }
  });

  it('4HHB’s heme contacts per pair of units, derived from its atoms: each heme touches its own chain only', async () => {
    const cif = await readCif(gunzipSync(readFileSync(`${PackageFolder}fixtures/structures/4HHB.cif.gz`)));
    expect(hemeContacts(cif)).toEqual(HemeUnitContacts);
  });

  it('the columns file: a record a row, keyed by its label; UniProt’s organism and class', () => {
    const [head, ...lines] = text(ColumnsFile).trim().split('\n');
    expect(head.split('\t')).toEqual(['id', 'OS', 'class', 'uniprot_length']);
    const rows = new Map(lines.map((l) => l.split('\t')).map((c) => [c[0], c] as const));
    expect([...rows.keys()].sort()).toEqual(seed.rows.map((r) => r.label).sort());
    expect(rows.get(ReferenceRow)!.slice(1, 3)).toEqual(['Bos grunniens (Wild yak) (Bos mutus grunniens)', 'Mammalia']);
  });
});

async function settled(spec: Node, runtimes: readonly FormatRuntime[] = []) {
  const model = resolveSpec(spec, { runtimes });
  const view = resolveViews(model)[0];
  const store = new DataStore();
  await attachSources(model, store, { assets: nodeFixtureAssets().reader });
  const renderer = new Renderer({
    types: BuiltInLayerTypes,
    layouts: BuiltInLayouts,
    model,
    view,
    data: store,
    tokens: LightTokens,
    text: bundledTextMetrics(),
  });
  return async (input: FrameInput) => {
    let f = renderer.frame(input);
    for (let i = 0; i < 6; i++) {
      store.want(f.wants);
      await store.whenSettled();
      const next = renderer.frame(input);
      f = next;
      if (!next.rebuilt) break;
    }
    return renderer.frame(input);
  };
}

function rowsOf(f: {
  layout: { sections: readonly { lays: readonly (readonly { lay: unknown; type: { name: string } }[])[] }[] };
}) {
  for (const s of f.layout.sections)
    for (const b of s.lays) for (const l of b ?? []) if (l.type.name === 'rows') return l.lay as RowsLay;
  throw new Error('no rows layer');
}

/** Every display row in the band: its band of cells (one y each) and its label; each cell in view hit as its row and column. */
function checkDrawn(
  f: Parameters<typeof rowsOf>[0] & { layout: Parameters<typeof hitView>[0] },
  letters: (row: number) => string,
) {
  const lay = rowsOf(f);
  const ys = new Set<number>();
  for (let k = 0; k < lay.n; k++) ys.add(lay.y[k]);
  let rows = 0;
  for (let i = lay.first; i < lay.last; i++) {
    const e = lay.display.rowAt(i);
    const y = lay.rowY(i);
    if (y < lay.box.y || y + lay.rowScale > lay.box.y + lay.box.height) continue;
    rows++;
    expect(ys.has(y), `row ${i} has its band`).toBe(true);
    if (e < 0) continue;
    const text = letters(e);
    for (let j = 0; j < Math.min(text.length, 40); j++) {
      if (!lay.axis.cell(j)) continue;
      const x = (lay.axis.sx0 + lay.axis.sx1) / 2;
      const h = hitView(f.layout, x, y + lay.rowScale / 2);
      expect(h).toMatchObject({ target: 'cell', position: j, yPosition: e });
    }
  }
  expect(rows).toBeGreaterThan(0);
  return lay;
}

describe('the Globins views against their data', () => {
  it('the scale variant’s rows are a function of the seed and the index', () => {
    expect(generatedRow(1234, 300, 1, profile)).toEqual(generatedRow(1234, 300, 1, profile));
    expect(generatedRow(1234, 300, 1, profile)).not.toEqual(generatedRow(1235, 300, 1, profile));
  });

  it('the family at letters: every row in the band drawn and hit as its row and column', async () => {
    const frame = await settled(globins.build(), FamilyRuntimes);
    const f = await frame({
      size: { width: 900, height: 600 },
      dpr: 1,
      settled: true,
      viewports: { '/views/main/family': { start: 0, end: 40 } },
    });
    const lay = checkDrawn(f, (r) => seed.rows[r].letters);
    expect(lay.foldsApply).toBe(true);
    expect(lay.tree.length).toBeGreaterThan(0);
  });

  it('the 10,000 rows across a tile edge (rows 1,020 to 1,030): each row drawn, labeled and hit', async () => {
    const frame = await settled(globins.scale(), [generatedRuntime()]);
    const f0 = await frame({ size: { width: 900, height: 600 }, dpr: 1, settled: true });
    const l0 = rowsOf(f0);
    // the band's display row of data row 1,020, at the band's top, at letters
    const top = l0.display.displayRow(1020);
    const f = await frame({
      size: { width: 900, height: 600 },
      dpr: 1,
      settled: true,
      viewports: {
        '/views/main/rows': { start: 0, end: 40, bands: { '/views/main/rows/rows': { bandTop: top, rowScale: 20 } } },
      },
    });
    const row = (r: number) => generatedRow(r, 300, 1, profile);
    const lay = checkDrawn(f, (r) => String.fromCharCode(...row(r)));
    for (let i = lay.first; i < lay.last; i++) {
      const e = lay.display.rowAt(i);
      expect(lay.label(e)).toBe(generatedLabel(e, row(e)));
    }
  });

  it('the end numbers agree with each row’s label range, held and lazy, whatever the band', async () => {
    const check = async (spec: Node, tops: number[]) => {
      const frame = await settled(spec, [generatedRuntime()]);
      for (const bandTop of tops) {
        const f = await frame({
          size: { width: 900, height: 600 },
          dpr: 1,
          settled: true,
          viewports: {
            '/views/main/rows': { start: 0, end: 300, bands: { '/views/main/rows/rows': { bandTop, rowScale: 20 } } },
          },
        });
        const lay = rowsOf(f);
        let k = 0;
        for (let i = lay.first; i < lay.last; i++) {
          const e = lay.display.rowAt(i);
          const range = /\/(\d+)-(\d+)$/.exec(lay.label(e))!;
          expect([lay.numberText[k], lay.lastText[k]], lay.label(e)).toEqual([range[1], range[2]]);
          k++;
        }
        expect(k).toBeGreaterThan(10);
      }
    };
    await check(globins.scale(), [0, 1020, 9970]);
    await check(globins.lazy(), [0, 5_000_000]);
  });

  it('the lazy 10,000,000 rows: the rows in view read alone, drawn and hit', async () => {
    const frame = await settled(globins.lazy(), [generatedRuntime()]);
    const f = await frame({
      size: { width: 900, height: 400 },
      dpr: 1,
      settled: true,
      viewports: {
        '/views/main/rows': {
          start: 0,
          end: 40,
          bands: { '/views/main/rows/rows': { bandTop: 5_000_000, rowScale: 20 } },
        },
      },
    });
    const row = (r: number) => generatedRow(r, 300, 2, profile);
    const lay = checkDrawn(f, (r) => String.fromCharCode(...row(r)));
    expect(lay.held).toBeUndefined();
    expect(lay.first).toBe(5_000_000);
    expect(lay.label(5_000_000)).toBe(generatedLabel(5_000_000, row(5_000_000)));
  });
});

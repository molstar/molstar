/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// V4, ST2: tile addressing, a lazy alignment and streamed annotations by `protocol` in compile and resolve (case 9, the
// tiles provider on the reference server), the newer-minor path, `upgrade`, and the CLIs. The client's and the
// provider's own tests are `@molstar/seq-core/stream/tiles`'s.

import { spawnSync } from 'node:child_process';
import { existsSync, mkdtempSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { createSeqViewSpec } from '@molstar/seq-core/builder';
import {
  type Alignment,
  type Annotation,
  type Issue,
  KeyKind,
  type KeyTarget,
  type SeqSet,
} from '@molstar/seq-core/model';
import { tilesStream } from '@molstar/seq-core/stream/tiles';
import { tileServer } from '@molstar/seq-core/stream/tiles/server';
import { describe, expect, it } from 'vitest';
import {
  compile,
  Lanes,
  type Loaders,
  resolveView,
  type SeqViewSpecDoc,
  StreamProp,
  supportedAll,
  Tiles,
  toJSON,
  type UpgradeStep,
  upgrade,
  validate,
} from '../../src/format/index.js';
import { itWithin } from '../model/timing.js';
import { case1, case9 } from './cases.js';
import { fixtureLoaders, PF00042, registries } from './load-fixtures.js';

const reg = registries();
const codes = (issues: Issue[]) => [...new Set(issues.map((i) => i.code))].sort();
const t: KeyTarget = { kind: -1, id: -1, row: -1 };
const find = (set: SeqSet, key: string) => (set.byKey(set.strings.find(key), t) ? { ...t } : null);
type J = { kind: string; ref?: string; params?: Record<string, unknown>; children?: J[] };
const byRef = (n: J, ref: string): J | undefined =>
  n.ref === ref ? n : n.children?.map((c) => byRef(c, ref)).find(Boolean);
const dataOf = (d: SeqViewSpecDoc) => (d.root as J).children?.find((c) => c.kind === 'data') as J;

describe('tile addressing', () => {
  it('addresses tiles by index and level', () => {
    expect([Tiles.tileOf(255, 0), Tiles.tileOf(256, 0), Tiles.tileOf(1000, 2), Tiles.tileOf(1024, 2)]).toEqual([
      0, 1, 0, 1,
    ]);
    expect([Tiles.start(3, 1), Tiles.sampleOf(1000, 2), Tiles.sampleOf(1027, 2), Tiles.count(18000, 0)]).toEqual([
      1536, 250, 0, 71,
    ]);
    expect([Tiles.levelFor(256), Tiles.levelFor(257), Tiles.levelFor(18000), Tiles.levelFor(320, 100)]).toEqual([
      0, 1, 7, 2,
    ]);
  });
});

describe('case 9 on the test server', () => {
  it("a lazy alignment validates without data; its counts and a protocol annotation's length are required", () => {
    expect(validate(case9(), reg)).toEqual([]);
    const d = case9();
    delete byRef(d.root as J, 'pf00042')?.params?.row_count;
    dataOf(d).children?.push(
      { kind: 'source', ref: 'cov.src', params: { uri: 'https://example.org/cov', protocol: 'seqstar-tiles/1' } },
      { kind: 'annotation', ref: 'cov', params: { geometry: 'values', source_ref: 'cov.src' } },
    );
    expect(validate(d, reg).map((i) => [i.code, i.path.split('/').pop()])).toEqual([
      ['svs.params.missing', 'row_count'],
      ['svs.params.missing', 'length'],
      ['svs.params.missing', 'target'],
    ]);
  });

  it('compiles: the named rows and the reference row fetched whole, the rows annotation and the tree fetching none', async () => {
    const loaders = fixtureLoaders();
    const { set, issues } = await compile(case9(), reg, loaders);
    expect(issues).toEqual([]);
    const named = [`rows/${PF00042.hbaRow}`, `rows/${PF00042.hbbRow}`];
    // The rows, then the statistic's manifest (its coarsest tile) when `cons` is computed.
    expect(loaders.tiles).toEqual(['rows/0', ...named, 'stats/conservation/1/0']);
    const a = set.alignments[find(set, 'pf00042')?.id ?? -1] as Alignment;
    expect([a.lazy, a.rowCount, a.columnCount, a.rows.count]).toEqual([true, 18000, 320, 3]);
    expect(set.strings.get(a.labelOf(PF00042.hbaRow as never))).toMatch(/^SEQ4711\//);
    expect(set.strings.get(a.labelOf(PF00042.hbbRow as never))).toMatch(/^SEQ9133\//);
    expect(set.strings.get(a.labelOf(0 as never))).toMatch(/^SEQ0\//);
    expect(a.labelOf(1 as never)).toBe(-1);
    // A declared reference row is the one fetched.
    const d = case9();
    Object.assign(byRef(d.root as J, 'pf00042')?.params ?? {}, { reference_row: 7 });
    const other = fixtureLoaders();
    expect((await compile(d, reg, other)).issues).toEqual([]);
    expect(other.tiles.filter((x) => x.startsWith('rows/'))).toEqual(['rows/7', ...named]);
    expect(a.tree?.leafOrder.length).toBe(18000);
    expect(set.mappings[find(set, 'pf00042-hba')?.id ?? -1]?.method).toBe('pairwise');
    expect(set.mappings[find(set, 'pf00042-hbb')?.id ?? -1]?.method).toBe('pairwise');
    const meta = set.annotations[find(set, 'meta')?.id ?? -1] as Annotation;
    expect([meta.geometry, meta.count]).toEqual(['rows', PF00042.rows]);
    const cons = set.annotations[find(set, 'cons')?.id ?? -1] as Annotation;
    expect([cons.count, cons.streamed]).toEqual([0, true]);
    expect(cons.custom.props.get(StreamProp)).toMatchObject({
      source: { data_kind: 'values', declared: { measure: 'conservation', length: 320, target: 'pf00042' } },
      manifest: { domain: [0, 0.9] },
    });
    expect(a.custom.props.get(StreamProp)).toMatchObject({ source: { data_kind: 'alignment' } });

    const { view, issues: resolved } = resolveView(case9(), set, 'msa', reg, supportedAll(reg));
    expect(resolved).toEqual([]);
    expect(view.sections.ref.length).toBe(2);
  });

  it('statistics not served are unavailable; so are other params and consensus (svs.provider.lazy)', async () => {
    const d = case9();
    Object.assign(byRef(d.root as J, 'cons')?.params ?? {}, { gaps: 'ignore' });
    const { set, issues } = await compile(d, reg, fixtureLoaders());
    expect(codes(issues)).toEqual(['svs.provider.lazy']);
    expect(find(set, 'cons')).toBeNull();
    const server = tileServer({ rowCount: PF00042.rows, columnCount: PF00042.columns });
    const r = await compile(case9(), reg, fixtureLoaders(undefined, server));
    expect(r.issues.map((i) => i.message).join()).toMatch(/serves no stats\/conservation/);
  });

  /** Case 9 with a streamed annotation over the columns, its own source; layers over it in both sections. */
  async function withCoverage(domain: [number, number] | null) {
    const d = case9();
    const COV = 'https://example.org/pfam/PF00042/coverage';
    dataOf(d).children?.push(
      { kind: 'source', ref: 'cov.src', params: { uri: COV, protocol: 'seqstar-tiles/1' } },
      {
        kind: 'annotation',
        ref: 'cov',
        params: { geometry: 'values', source_ref: 'cov.src', target: 'pf00042', length: 320 },
      },
    );
    const values = (ref: string) => ({ kind: 'layer', ref, params: { type: 'values', data: 'cov', domain } });
    const msa = byRef(d.root as J, 'msa.main') as J;
    msa.children?.push({ kind: 'track', ref: 'cov.track', children: [values('cov.main')] });
    byRef(d.root as J, 'msa.chains')?.children?.push({ kind: 'track', ref: 'cov.chains', children: [values('cov.a')] });
    const cov = tileServer({ rowCount: 1, columnCount: 320, length: 320, values: (i) => i / 320 });
    const base = fixtureLoaders();
    const loaders: Loaders = {
      load: base.load,
      stream: (source) =>
        source.uri === COV ? tilesStream(() => cov.transport)(source) : (base.stream?.(source) ?? null),
    };
    const { set, issues } = await compile(d, reg, loaders);
    expect(issues).toEqual([]);
    // Only its manifest: the domain from the coarsest tile.
    expect(cov.calls).toEqual(['values/1/0']);
    return { set, ...resolveView(d, set, 'msa', reg, supportedAll(reg)) };
  }

  it('a streamed annotation reaches the chains through their rows (ST7); its manifest gives its domain', async () => {
    const { set, view, issues } = await withCoverage(null);
    expect(codes(issues)).toEqual([]);
    const layer = (ref: string) => Array.from(view.layers.ref).indexOf(set.strings.find(ref));
    expect(view.layers.placementCount[layer('cov.main')]).toBe(1);
    // On each chain's block through its row's path: the columns asked back through it (the viewer's).
    const a = layer('cov.a');
    const start = view.layers.placementStart[a];
    expect(view.layers.placementCount[a]).toBeGreaterThan(0);
    expect(view.placements.path[start]).toBeGreaterThanOrEqual(0);
    const given = await withCoverage([0, 1]);
    expect(codes(given.issues)).toEqual([]);
  });

  it('served statistics give their domain', async () => {
    const d = case9();
    const layer = (byRef(d.root as J, 'msa.main') as J).children?.flatMap((c) => c.children ?? []);
    const cons = layer?.find((l) => l.params?.data === 'cons') as J;
    delete cons.params?.domain;
    const { set } = await compile(d, reg, fixtureLoaders());
    const { view, issues } = resolveView(d, set, 'msa', reg, supportedAll(reg));
    expect(issues).toEqual([]);
    const l = Array.from(view.layers.ref).indexOf(set.strings.find(cons.ref ?? ''));
    expect([view.layers.domainMin[l], view.layers.domainMax[l]]).toEqual([0, 0.9]);
  });
});

describe("a lazy alignment's lanes: a run, never a lane per row", () => {
  const n = 1_000_000;
  /** A million-row lazy alignment, a subset of three in one section and the whole in another, compiled. */
  async function million() {
    const uri = 'https://example.org/million';
    const b = createSeqViewSpec();
    const aln = b.data.alignment({
      ref: 'm',
      row_count: n,
      column_count: 40,
      source: { uri, protocol: 'seqstar-tiles/1' },
    });
    // A section has one rows layer: the subset and the whole alignment in two sections of one view.
    const v = b.view({ ref: 'v', label: 'V' });
    v.block({ alignment: aln })
      .track({ ref: 'v.subset' })
      .rows({ row_subset: [999_999, 5, 2] });
    v.section({ ref: 'v.all' }).block({ alignment: aln }).track({ ref: 'v.rows' }).rows();
    const server = tileServer({ rowCount: n, columnCount: 40 });
    const loaders: Loaders = {
      load: (r) => Promise.reject(new Error(`no fixture for ${r.uri}`)),
      stream: tilesStream((u) => (u === uri ? server.transport : null)),
    };
    const { set, issues: compiled } = await compile(b.getDocument(), reg, loaders);
    expect(compiled).toEqual([]);
    return { doc: b.getDocument(), set };
  }

  itWithin("a million rows resolve and compose; a row_subset's lanes stay stored", async (within) => {
    const { doc, set } = await million();
    const { view, issues } = within(2000, () => resolveView(doc, set, 'v', reg, supportedAll(reg)));
    expect(issues).toEqual([]);
    const { lanes, tracks } = view;
    // Stored: the subset's three, in data order; then the run of a million, nothing listed for it.
    expect([lanes.count, Array.from(lanes.row), Array.from(lanes.depth)]).toEqual([n + 3, [2, 5, 999_999], [0, 0, 0]]);
    expect([lanes.lazy.count, Array.from(lanes.lazy.start), Array.from(lanes.lazy.length)]).toEqual([1, [3], [n]]);
    expect([Array.from(tracks.laneStart), Array.from(tracks.laneCount)]).toEqual([
      [0, 3],
      [3, n],
    ]);
    expect([0, 1, 2, 3, 500_003, n + 2, n + 3, -1].map((k) => Lanes.rowOf(lanes, k))).toEqual([
      2,
      5,
      999_999,
      0,
      500_000,
      n - 1,
      -1,
      -1,
    ]);
    expect([Lanes.depthOf(lanes, 1), Lanes.depthOf(lanes, n + 2), Lanes.depthOf(lanes, n + 3)]).toEqual([0, 0, -1]);
  });

  it("case 9: the rows track is a run in tree order; a row's placement lands on its rank", async () => {
    // The case's 1A3N secondary structure on the rows track too: it reaches the alignment through HBA's and HBB's rows.
    const d = case9();
    const { set } = await compile(d, reg, fixtureLoaders());
    const { view, issues } = resolveView(d, set, 'msa', reg, supportedAll(reg));
    expect(issues).toEqual([]);
    const order = (set.alignments[find(set, 'pf00042')?.id ?? -1] as Alignment).tree?.leafOrder as Int32Array;
    const { lanes } = view;
    expect([lanes.row.length, lanes.lazy.count, lanes.lazy.length[0], lanes.lazy.order[0] === order]).toEqual([
      0,
      1,
      PF00042.rows,
      true,
    ]);
    expect([0, 17, PF00042.rows - 1].map((k) => Lanes.rowOf(lanes, k))).toEqual([order[0], order[17], order.at(-1)]);
    const rowsTrack = view.layers.track[Array.from(view.layers.type).indexOf(set.strings.find('rows'))];
    let l = -1;
    for (let k = 0; k < view.layers.count; k++)
      if (view.layers.track[k] === rowsTrack && set.strings.get(view.layers.type[k]) === 'secondary') l = k;
    const at = view.layers.placementStart[l];
    expect(view.layers.placementCount[l]).toBe(2);
    const placed = [0, 1].map((p) => Lanes.rowOf(lanes, view.placements.lane[at + p])).sort((x, y) => x - y);
    expect(placed).toEqual([PF00042.hbaRow, PF00042.hbbRow]);
  });
});

describe('versions', () => {
  it('a newer minor validates with warnings, compiles without its unknowns and draws unknown types by fallback', async () => {
    const d = case1() as SeqViewSpecDoc & { root: J };
    d.version = '1.4.0';
    dataOf(d).children?.push({
      kind: 'gizmo',
      ref: 'g',
      children: [{ kind: 'annotation', ref: 'inside', params: { geometry: 'ranges', value: [] } }],
    });
    const p53 = byRef(d.root, 'p53') as J;
    Object.assign(p53.params ?? {}, { glow: true });
    const layer = JSON.stringify(d.root).match(/"kind":"layer","ref":"([^"]+)"/)?.[1] ?? '';
    Object.assign(byRef(d.root, layer)?.params ?? {}, { type: 'sparkle', fallback: 'ranges' });
    const issues = validate(d, reg);
    expect(issues.every((i) => i.severity === 'warning' && i.code === 'svs.version.newer_minor')).toBe(true);
    expect(issues.length).toBe(4);
    const { set, issues: loaded } = await compile(d, reg, fixtureLoaders());
    expect(loaded.filter((i) => i.severity === 'error')).toEqual([]);
    expect(find(set, 'inside')).toBeNull();
    expect(find(set, 'p53')?.kind).toBe(KeyKind.molecule);
    const view = String(d.root.children?.find((c) => c.kind === 'view')?.ref);
    expect(codes(resolveView(d, set, view, reg, supportedAll(reg)).issues)).toContain('svs.render.fallback');
  });

  it('upgrade runs one step per older version, then the document validates', () => {
    const steps: UpgradeStep[] = [
      { from: '0.8', to: '0.9', run: (d) => ({ ...d, metadata: { title: 'renamed' } }) },
      {
        from: '0.9',
        to: '1.0',
        run: (d) => {
          const { title, ...rest } = d as Record<string, unknown>;
          return { ...rest, metadata: { ...(rest.metadata as object), description: title } };
        },
      },
    ];
    const old = { ...case1(), version: '0.8.1', title: 'p53' };
    const r = upgrade(old, steps);
    expect(r.steps).toEqual(['0.8 -> 0.9', '0.9 -> 1.0']);
    expect(r.doc).toMatchObject({ version: '1.0.0', metadata: { title: 'renamed', description: 'p53' } });
    expect(old.version).toBe('0.8.1');
    expect(validate(r.doc, reg)).toEqual([]);
    const current = case1();
    expect(upgrade(current, steps)).toEqual({ doc: current, steps: [] });
    expect(upgrade({ ...current, version: '0.7.0' }, steps).steps).toEqual([]);
    expect(upgrade(42).doc).toBe(42);
  });
});

describe('CLIs', () => {
  const dir = mkdtempSync(join(tmpdir(), 'svs-cli-'));
  const cli = (name: string, ...args: string[]) =>
    spawnSync(process.execPath, [fileURLToPath(new URL(`../../lib/format/cli/${name}.js`, import.meta.url)), ...args], {
      encoding: 'utf8',
    });
  const file = (name: string, text: string) => {
    writeFileSync(join(dir, name), text);
    return join(dir, name);
  };

  // A process each: slow. The print below is the fast tier's one CLI run; the CLIs run built (`tsc -b`: the slow tier
  // builds first), so it waits for a `lib` (packages/sequences/AGENTS.md, Tests).
  it("svs-validate: exit 0 valid, 0 with a newer minor's warnings, 1 with errors, 2 on usage", {
    tags: ['slow'],
  }, () => {
    const ok = cli('svs-validate', file('ok.json', toJSON(case1())));
    expect([ok.status, ok.stdout, ok.stderr.trim()]).toEqual([
      0,
      '',
      `${join(dir, 'ok.json')}: 0 error(s), 0 warning(s), 0 note(s)`,
    ]);
    const newer = cli('svs-validate', file('newer.json', toJSON({ ...case1(), version: '1.2.0' })));
    expect([newer.status, newer.stdout]).toEqual([0, expect.stringContaining('warning svs.version.newer_minor')]);
    const strict = cli(
      'svs-validate',
      file('strict.json', toJSON({ ...case1(), version: '1.2.0', extra: 1 } as never)),
      '--strict',
    );
    expect([strict.status, strict.stdout]).toEqual([1, expect.stringContaining('error svs.params.unknown')]);
    expect(cli('svs-validate', file('bad.json', '{ nope')).status).toBe(1);
    expect(cli('svs-validate', join(dir, 'missing.json')).status).toBe(2);
    expect(cli('svs-validate').status).toBe(2);
  });

  const built = existsSync(fileURLToPath(new URL('../../lib/format/cli/svs-print-schema.js', import.meta.url)));
  it.skipIf(!built)('svs-print-schema prints the reference', () => {
    const r = cli('svs-print-schema');
    expect(r.status).toBe(0);
    expect(r.stdout).toMatch(/^# SeqViewSpec 1\.0\.0 reference/);
    expect(r.stdout).toContain('### `alignment`');
  });
});

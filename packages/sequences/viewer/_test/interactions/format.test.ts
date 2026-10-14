/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// I1: the vocabulary and queries on a small hand-made set; the fingerprint of case 6 over the recorded 1IEP contacts
// (packages/sequences/core/_test/fixtures/1iep.json: the crystal STI and Vina's four poses) through the format's case and
// fixture loaders; the fingerprint over each kind of coordinate system, polymer partners only, the fallbacks; and that
// the kernel never imports this package.

import { readdirSync, readFileSync, statSync } from 'node:fs';
import { createSeqViewSpec } from '@molstar/seq-core/builder';
import {
  compile,
  createRegistries,
  KernelRegistry,
  MvsRegistry,
  type ResolvedView,
  register,
  resolveView,
  type SeqViewSpecDoc,
  supportedAll,
  validate,
} from '@molstar/seq-core/format';
import {
  type Annotation,
  type CoordinateSystem,
  CoordinateSystemKind,
  type Issue,
  type KeyTarget,
  type SeqSet,
  SeqSets,
  type UnitId,
} from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { case6 } from '../../../core/_test/format/cases.js';
import { fixtureLoaders, MolstarFixture, registries, STI_NAMES } from '../../../core/_test/format/load-fixtures.js';
import {
  FingerprintProp,
  Fingerprints,
  fingerprint,
  fingerprintOf,
  Interactions,
  InteractionsRegistry,
  interactionColors,
  interactionGlyphs,
  Kind,
  KindMask,
  type KindName,
  KindStyle,
  kindByName,
} from '../../src/interactions/index.js';
import { alignmentScene, mappingScene, unitScene } from './fixtures.js';

const t: KeyTarget = { kind: -1, id: -1, row: -1 };
const id = (set: SeqSet, key: string) => (set.byName(key, t) ? t.id : -1);
const annotation = (set: SeqSet, key: string) => set.annotations[id(set, key)] as Annotation;
const codes = (issues: Issue[]) => [...new Set(issues.map((i) => i.code))].sort();
const errors = (issues: Issue[]) => issues.filter((i) => i.severity === 'error');
const molecule = (set: SeqSet, key: string) => ({ kind: CoordinateSystemKind.molecule, id: id(set, key), row: -1 });
const unit = (set: SeqSet, key: string) => ({ kind: CoordinateSystemKind.unit, id: id(set, key), row: -1 });

function withExtension() {
  const reg = createRegistries();
  register(reg, [KernelRegistry, MvsRegistry, MolstarFixture, InteractionsRegistry]);
  return reg;
}

/**
 * A protein P (ACDE) and a two-atom ligand L; four pairs, the third of a kind this vocabulary does not know and the
 * fourth with no category (`other` too). A second protein Q (WYW) no path reaches from P.
 */
function small() {
  const { set, issues } = SeqSets.build({
    molecules: [
      { key: 'p', residues: 'ACDE' },
      { key: 'q', residues: 'WYW' },
      { key: 'l', kind: 'ligand', elements: ['C', 'N'] },
    ],
    units: [
      { key: 'P', molecule: 'p' },
      { key: 'L', molecule: 'l' },
    ],
    annotations: [
      {
        key: 'c',
        geometry: 'pairs',
        columns: {
          target: ['P', 'P', 'P', 'P'],
          start: [0, 1, 1, 3],
          b_target: ['L', 'L', 'L', 'L'],
          b_start: [0, 1, 0, 1],
          category: ['hbond', 'hydrophobic', 'weird', null],
        },
      },
    ],
  });
  expect(errors(issues)).toEqual([]);
  return { set, c: annotation(set, 'c'), P: id(set, 'P') as UnitId, L: id(set, 'L') as UnitId };
}

describe('the vocabulary and queries', () => {
  it('kinds and masks', () => {
    expect(KindMask.count(KindMask.of('hbond', 'hydrophobic'))).toBe(2);
    expect(KindMask.names(KindMask.of('pi_stacking', 'contact'))).toEqual(['contact', 'pi_stacking']);
    expect(KindMask.count(KindMask.all)).toBe(11);
    expect([kindByName('salt_bridge'), kindByName('contact'), kindByName('weird')]).toEqual([2, 0, Kind.other]);
  });

  it('a row with no category is `other`, and `contact` only a stated one', () => {
    expect([kindByName(null), kindByName('')]).toEqual([Kind.other, Kind.other]);
    const { set, c } = small();
    expect([0, 1, 2, 3].map((r) => Interactions.kindOf(set, c, r))).toEqual([
      Kind.hbond,
      Kind.hydrophobic,
      Kind.other,
      Kind.other,
    ]);
  });

  it('contactsOf aggregates kinds per residue, filtered by a mask; re-entrant', () => {
    const { set, c, P, L } = small();
    const out: number[][] = [];
    Interactions.contactsOf(set, c, L, KindMask.all, (u, i, m) => out.push([u, i, m]), null);
    expect(out).toEqual([
      [P, 0, KindMask.of('hbond')],
      [P, 1, KindMask.of('hydrophobic', 'other')],
      [P, 3, KindMask.of('other')],
    ]);
    const hb: number[][] = [];
    Interactions.contactsOf(set, c, L, KindMask.of('hbond'), (u, i, m) => hb.push([u, i, m]), null);
    expect(hb).toEqual([[P, 0, KindMask.of('hbond')]]);
    const nested: number[][] = [];
    Interactions.contactsOf(
      set,
      c,
      L,
      KindMask.all,
      (u, i, m) => {
        nested.push([u, i, m]);
        if (i === 0) Interactions.contactsOf(set, c, L, KindMask.all, () => {}, null);
      },
      null,
    );
    expect(nested).toEqual(out);
    expect(Interactions.contactsOf(set, c, P, KindMask.all, (u) => expect(u).toBe(L), null)).toBeUndefined();
  });

  it("atomsContacting: ligand atoms at one residue with each pair's kind, through `Pairs.endUnits`", () => {
    const { set, c, P, L } = small();
    const out: number[][] = [];
    Interactions.atomsContacting(set, c, P, 1 as never, (u, a, k) => out.push([u, a, k]), null);
    expect(out.sort()).toEqual([
      [L, 0, Kind.other],
      [L, 1, Kind.hydrophobic],
    ]);
    const none: number[] = [];
    Interactions.atomsContacting(set, c, P, 2 as never, (u) => none.push(u), null);
    expect(none).toEqual([]);
    // Nested in its own visitor: each level has its state.
    const inner: number[][] = [];
    Interactions.atomsContacting(
      set,
      c,
      P,
      1 as never,
      (_u, a) => {
        if (a === 0) Interactions.atomsContacting(set, c, P, 0 as never, (u, b, k) => inner.push([u, b, k]), null);
      },
      null,
    );
    expect(inner).toEqual([[L, 0, Kind.hbond]]);
  });

  it('a molecule-level endpoint reaches the ligand units over it, not every unit of the set', () => {
    // Two copies of one ligand molecule; a pair whose b side names the molecule touches both, and a third unit of
    // another molecule is not asked.
    const { set, issues } = SeqSets.build({
      molecules: [
        { key: 'p', residues: 'ACDE' },
        { key: 'l', kind: 'ligand', elements: ['C', 'N'] },
        { key: 'm', kind: 'ligand', elements: ['O'] },
      ],
      units: [
        { key: 'P', molecule: 'p' },
        { key: 'L1', molecule: 'l' },
        { key: 'L2', molecule: 'l' },
        { key: 'M', molecule: 'm' },
      ],
      annotations: [
        { key: 'c', geometry: 'pairs', columns: { target: ['P'], start: [1], b_target: ['l'], b_start: [1] } },
      ],
    });
    expect(errors(issues)).toEqual([]);
    const out: number[][] = [];
    Interactions.atomsContacting(
      set,
      annotation(set, 'c'),
      id(set, 'P') as UnitId,
      1 as never,
      (u, a) => out.push([u, a]),
      null,
    );
    expect(out.sort()).toEqual([
      [id(set, 'L1'), 1],
      [id(set, 'L2'), 1],
    ]);
  });
});

describe('fingerprints over a coordinate system', () => {
  /** Two chains P1, P2 over one sequence p, a ligand L and a cofactor C; L contacts both chains and the cofactor. */
  function chains() {
    const { set, issues } = SeqSets.build({
      molecules: [
        { key: 'p', residues: 'ACDEFG' },
        { key: 'l', kind: 'ligand', elements: ['C', 'N'] },
        { key: 'c', kind: 'ligand', elements: ['FE'] },
      ],
      units: [
        { key: 'P1', molecule: 'p' },
        { key: 'P2', molecule: 'p' },
        { key: 'L', molecule: 'l' },
        { key: 'C', molecule: 'c' },
      ],
      annotations: [
        {
          key: 'pairs',
          geometry: 'pairs',
          columns: {
            target: ['P1', 'P1', 'P2', 'P2', 'C'],
            start: [1, 2, 2, 4, 0],
            b_target: ['L', 'L', 'L', 'L', 'L'],
            b_start: [0, 1, 0, 1, 1],
            category: ['hbond', 'hydrophobic', 'hbond', 'salt_bridge', 'metal'],
          },
        },
      ],
    });
    expect(errors(issues)).toEqual([]);
    return { set, pairs: annotation(set, 'pairs'), L: id(set, 'L') };
  }
  const build = (set: SeqSet, pairs: Annotation, L: number, system: CoordinateSystem) => {
    const issues: Issue[] = [];
    return { f: Fingerprints.build(set, [pairs], Int32Array.of(L), system, { issues }), issues };
  };

  it('a sequence system merges every unit over the molecule; the cofactor, a ligand, is no column (quietly)', () => {
    const { set, pairs, L } = chains();
    const { f, issues } = build(set, pairs, L, molecule(set, 'p'));
    expect([[...f.columns], [...f.cells]]).toEqual([
      [1, 2, 4],
      [KindMask.of('hbond'), KindMask.of('hydrophobic', 'hbond'), KindMask.of('salt_bridge')],
    ]);
    expect(f.system).toEqual(molecule(set, 'p'));
    expect(issues).toEqual([]);
  });

  it("a unit system keeps that one chain's contacts; the other chain's have no path: unmapped, warned once", () => {
    const { set, pairs, L } = chains();
    const { f, issues } = build(set, pairs, L, unit(set, 'P1'));
    expect([[...f.columns], [...f.cells]]).toEqual([
      [1, 2],
      [KindMask.of('hbond'), KindMask.of('hydrophobic')],
    ]);
    expect(f.system).toEqual(unit(set, 'P1'));
    expect(issues.map((i) => [i.code, i.severity, i.message])).toEqual([
      ['interactions.fingerprint.unmapped', 'warning', expect.stringContaining('left out: 2 of "P2"')],
    ]);
  });

  it('a residue with no path into the system: interactions.fingerprint.unmapped, left out', () => {
    const { set, c, L } = small();
    const issues: Issue[] = [];
    const f = Fingerprints.build(set, [c], Int32Array.of(L), molecule(set, 'q'), { issues });
    expect([f.columns.length, f.cells.length]).toEqual([0, 0]);
    expect(issues.map((i) => [i.code, i.severity])).toEqual([['interactions.fingerprint.unmapped', 'warning']]);
    const own = Fingerprints.build(set, [c], Int32Array.of(L), molecule(set, 'p'));
    expect([[...own.columns], [...own.cells]]).toEqual([
      [0, 1, 3],
      [KindMask.of('hbond'), KindMask.of('hydrophobic', 'other'), KindMask.of('other')],
    ]);
  });
});

// --- The provider's systems, in documents -----------------------------------------------------------------

async function load(doc: SeqViewSpecDoc, reg = withExtension()) {
  const { set, issues } = await compile(doc, reg, fixtureLoaders());
  expect(errors(issues)).toEqual([]);
  return { set, issues, reg };
}

describe('the provider over each kind of system', () => {
  it("a unit: `system` validates a unit ref, and the columns are the unit's positions", async () => {
    const { doc } = unitScene();
    expect(errors(validate(doc, withExtension()))).toEqual([]);
    const { set, issues } = await load(doc);
    expect(issues).toEqual([]);
    const fp = annotation(set, 'fp');
    const f = fingerprintOf(fp);
    expect([...(f?.columns ?? [])]).toEqual([3, 7, 12, 20]);
    expect(f?.system).toEqual(unit(set, 'P'));
    expect([fp.geometry, fp.matrix?.states.length]).toEqual(['matrix', 4]);
  });

  it('an alignment: contacts of two chains land in the columns their rows give', async () => {
    const { doc } = alignmentScene();
    const { set } = await load(doc);
    const f = fingerprintOf(annotation(set, 'fp'));
    // P0's 3 and 10 are columns 3 and 10; P1's 10 and 15 are columns 11 and 16 (a gap at column 5 of its row).
    expect([...(f?.columns ?? [])]).toEqual([3, 10, 11, 16]);
    const aln = id(set, 'aln');
    expect(f?.system).toEqual({ kind: CoordinateSystemKind.columns, id: aln, row: -1 });
    const at = (row: number, column: number) => f?.cells[row * 4 + (f?.columns.indexOf(column) ?? -1)];
    expect([at(0, 3), at(1, 3), at(2, 3), at(0, 10), at(2, 11), at(2, 16)]).toEqual([
      KindMask.of('hbond'),
      KindMask.of('hbond'),
      0,
      KindMask.of('hydrophobic'),
      KindMask.of('hbond'),
      KindMask.of('pi_stacking'),
    ]);
  });

  it("a sequence reached through a mapping: the unit's contacts are the sequence's positions", async () => {
    const { doc } = mappingScene();
    const { set } = await load(doc);
    expect([...(fingerprintOf(annotation(set, 'fp'))?.columns ?? [])]).toEqual([8, 12, 20, 30]);
  });

  it('a system that is no unit, sequence or alignment is reported at /params/options/system', async () => {
    const doc = structuredClone(unitScene().doc);
    const data = doc.root.children?.find((c) => c.kind === 'data');
    const node = data?.children?.find((c) => c.ref === 'fp');
    if (node?.params) node.params = { ...node.params, options: { ...(node.params.options as object), system: 'nope' } };
    const reg = withExtension();
    expect(errors(validate(doc, reg)).map((i) => i.path)).toEqual([expect.stringContaining('system')]);
  });
});

// --- Case 6: imatinib in ABL ---------------------------------------------------------------------------------------

/** 1a index of an auth residue number (1a's numbers are 1IEP's auth numbers). */
const at = (auth: number) => auth - 1;
const T315 = at(315);
const M318 = at(318);
const LIGANDS = ['1iep/C', 'poses/0/ligand', 'poses/1/ligand', 'poses/2/ligand', 'poses/3/ligand'];

const declared = (v: ResolvedView, set: SeqSet, type: string) =>
  Array.from(v.layers.declared).indexOf(set.strings.find(type));

describe('case 6: the fingerprint of the crystal pose and four Vina poses', () => {
  it('poses 3 and 4 lose the T315 and M318 hydrogen bonds', async () => {
    const { set, issues } = await load(case6());
    expect(issues.filter((i) => i.code === 'svs.provider.unavailable')).toEqual([]);
    const fp = annotation(set, 'fp');
    expect([fp.geometry, fp.matrix?.states.length, set.strings.get(fp.semanticType)]).toEqual([
      'matrix',
      5,
      'interactions:fingerprint',
    ]);
    const f = fingerprintOf(fp);
    expect(f && fp.custom.props.get(FingerprintProp)).toBe(f);
    if (!f) return;
    expect([...f.rows]).toEqual(LIGANDS.map((k) => id(set, k)));
    expect(f.system).toEqual({ kind: CoordinateSystemKind.molecule, id: id(set, 'abl1a'), row: -1 });
    const lostHbonds = (b: number) => {
      const out: number[] = [];
      Fingerprints.diff(f, 0, b, (col, lost) => KindMask.has(lost, Kind.hbond) && out.push(col), null);
      return out;
    };
    for (const pose of [3, 4]) expect(lostHbonds(pose)).toEqual(expect.arrayContaining([T315, M318]));
    for (const pose of [1, 2]) expect(lostHbonds(pose)).not.toContain(T315);
    for (const pose of [1, 2]) expect(lostHbonds(pose)).not.toContain(M318);
    // The matrix: a row per contacted position, the number of kinds per ligand; no contact is missing (NaN).
    const C = f.columns.length;
    expect([...fp.start]).toEqual([...f.columns]);
    expect(fp.matrix && [...fp.matrix.values]).toEqual(
      Array.from({ length: C * 5 }, (_, k) => KindMask.count(f.cells[(k % 5) * C + Math.floor(k / 5)]) || Number.NaN),
    );
    expect(fp.matrix?.values.some((v) => v === 0)).toBe(false);
  });

  it("the queries on the recorded contacts: T315's partners", async () => {
    const { set } = await load(case6());
    const contacts = annotation(set, 'contacts');
    const A = id(set, '1iep/A') as UnitId;
    const res: number[][] = [];
    Interactions.contactsOf(
      set,
      contacts,
      id(set, '1iep/C') as UnitId,
      KindMask.all,
      (u, i, m) => {
        if (u === A && (i === 90 || i === 93)) res.push([i, m]);
      },
      null,
    );
    expect(res).toEqual([
      [90, KindMask.of('hbond')],
      [93, KindMask.of('hbond')],
    ]); // auth 315 and 318 are chain indices 90 and 93
    const atoms: string[] = [];
    Interactions.atomsContacting(
      set,
      contacts,
      A,
      90 as never,
      (u, a, k) => {
        if (k === Kind.hbond) atoms.push(`${set.strings.get(set.units.key[u])} ${STI_NAMES[a]}`);
      },
      null,
    );
    // The recorded hydrogen bonds at T315: crystal atom 12, pose 1 atom 10, pose 2 atom 16 (1iep.json).
    expect(atoms.sort()).toEqual([
      `1iep/C ${STI_NAMES[12]}`,
      `poses/0/ligand ${STI_NAMES[10]}`,
      `poses/1/ligand ${STI_NAMES[16]}`,
    ]);
  });

  it('a viewer without the type draws the fallback matrix of kind counts', async () => {
    const { set, reg } = await load(case6());
    const all = supportedAll(reg);
    const full = resolveView(case6(), set, 'poses-view', reg, all);
    expect(full.issues).toEqual([]);
    const l = declared(full.view, set, 'interactions:fingerprint');
    expect(set.strings.get(full.view.layers.type[l])).toBe('interactions:fingerprint');
    expect(full.view.layers.placementCount[l]).toBeGreaterThan(0);

    const without = new Set([...all].filter((n) => n !== 'interactions:fingerprint'));
    const r = resolveView(case6(), set, 'poses-view', reg, without);
    expect(codes(r.issues)).toEqual(['svs.render.fallback']);
    const k = declared(r.view, set, 'interactions:fingerprint');
    expect(set.strings.get(r.view.layers.type[k])).toBe('matrix');
    expect(r.view.layers.placementCount[k]).toBe(full.view.layers.placementCount[l]);
    const data = r.view.data[r.view.layers.dataStart[k]];
    expect(data).toBe(annotation(set, 'fp').id);
    // The recorded contacts give one kind per residue and pose: count 1, or missing (another pose's residue).
    expect([r.view.layers.domainMin[k], r.view.layers.domainMax[k]]).toEqual([1, 1]);
  });

  it('without the registry entry: the computed node unavailable, nothing cascades', async () => {
    const reg = registries();
    expect(errors(validate(case6(), reg))).toEqual([]);
    const { set, issues } = await load(case6(), reg);
    expect(issues.filter((i) => i.code === 'svs.provider.unavailable').map((i) => i.message)).toEqual([
      expect.stringContaining('"fp" is unavailable'),
    ]);
    expect(set.byName('fp', t)).toBe(false);
    expect(annotation(set, 'pocket').count).toBeGreaterThan(0);
    const all = supportedAll(reg);
    const poses = resolveView(case6(), set, 'poses-view', reg, all);
    expect(errors(poses.issues)).toEqual([]);
    expect(poses.view.layers.placementCount[declared(poses.view, set, 'interactions:fingerprint')]).toBe(0);
    expect(resolveView(case6(), set, 'chains', reg, all).issues).toEqual([]);
  });
});

describe('builder functions', () => {
  it("fingerprint() writes case 6's computed node", () => {
    const doc = case6();
    const data = doc.root.children?.find((c) => c.kind === 'data');
    const fp = data?.children?.find((c) => c.ref === 'fp');
    const b = createSeqViewSpec();
    fingerprint(b.data, { ref: 'fp', pairs: ['contacts'], ligands: LIGANDS, system: 'abl1a' });
    expect(b.getDocument().root.children?.[0].children?.[0]).toEqual(fp);
  });

  it("case 6's glyph map, written out in the core's cases, draws each kind it names as the extension's style does (VT-12)", () => {
    const doc = case6();
    const glyphs: Record<string, string>[] = [];
    const walk = (n: { params?: Record<string, unknown>; children?: unknown[] }) => {
      const g = n.params?.glyph as { theme?: string; params?: { map?: Record<string, string> } } | undefined;
      if (g?.theme === 'field_map' && g.params?.map) glyphs.push(g.params.map);
      for (const c of n.children ?? []) walk(c as never);
    };
    walk(doc.root);
    expect(glyphs.length).toBeGreaterThan(0);
    for (const map of glyphs)
      for (const [kind, shape] of Object.entries(map))
        expect([kind, shape]).toEqual([kind, KindStyle[kind as KindName][2]]);
  });

  it('theme references whose fallbacks are kernel themes with the colors and glyphs written out', async () => {
    const colors = interactionColors();
    const glyphs = interactionGlyphs();
    expect(colors.fallback?.theme).toBe('field_palette');
    const palette = colors.fallback?.params?.palette as { colors: Record<string, string>; missing_color: string };
    expect(palette.colors.hbond).toBe('#4e79a7');
    // A row with no category is `other`: the palette's missing_color and the map's default are `other`'s.
    expect(palette.missing_color).toBe(palette.colors.other);
    const map = glyphs.fallback?.params?.map as Record<string, unknown>;
    expect(map.hbond).toEqual({ shape: 'diamond', label: 'hydrogen bond' });
    expect(glyphs.fallback?.params?.default).toEqual(map.other);
    const doc = case6();
    const layer = doc.root.children?.[1].children?.[0].children?.find((c) => c.ref === 'chains.contacts')
      ?.children?.[0];
    if (layer?.params) Object.assign(layer.params, { color: colors, glyph: glyphs });
    for (const reg of [registries(), withExtension()]) expect(validate(doc, reg)).toEqual([]);
    const { set, reg } = await load(doc);
    const all = supportedAll(reg);
    expect(resolveView(doc, set, 'chains', reg, all).issues).toEqual([]);
    const without = new Set([...all].filter((n) => n !== 'interactions:kind'));
    const r = resolveView(doc, set, 'chains', reg, without);
    expect(codes(r.issues)).toEqual(['svs.render.theme']);
    const l = declared(r.view, set, 'pairs');
    expect(r.view.layers.themes[l].map((x) => x?.theme)).toEqual([
      'field_palette',
      'field_map',
      undefined,
      undefined,
      undefined,
    ]); // then its size channels
  });
});

describe("the layer type's channels", () => {
  it("`color` the kinds' and `counts` the kind counts' (Blues from 0), and no glyph", () => {
    const type = InteractionsRegistry.layerTypes?.[0];
    expect(type?.channels.map((c) => [c.name, c.kind])).toEqual([
      ['color', 'color'],
      ['counts', 'color'],
    ]);
    expect(type?.channels[0].default?.theme).toBe('interactions:kind');
    expect(type?.channels[1].default).toMatchObject({
      theme: 'field_palette',
      params: { palette: { kind: 'continuous', colors: 'Blues', value_domain: [0, null] } },
    });
    // The provider declares no `defaultColor`: resolve would give it to the counts channel too (routed note).
    expect(InteractionsRegistry.providers?.[0].defaultColor).toBeUndefined();
  });
});

describe('dependencies', () => {
  it('the kernel imports nothing from this package', () => {
    const files = (dir: URL): URL[] =>
      readdirSync(dir).flatMap((name) => {
        const url = new URL(name, dir);
        return statSync(url).isDirectory() ? files(new URL(`${name}/`, dir)) : name.endsWith('.ts') ? [url] : [];
      });
    const kernel = ['model', 'format'].flatMap((p) => files(new URL(`../../../core/src/${p}/`, import.meta.url)));
    expect(kernel.length).toBeGreaterThan(20);
    expect(kernel.filter((f) => readFileSync(f, 'utf8').includes('@molstar/seq-viewer/interactions'))).toEqual([]);
  });
});

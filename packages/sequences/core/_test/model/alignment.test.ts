/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  type Alignment,
  type AlignmentInput,
  Alignments,
  type ColumnIndex,
  Fasta,
  type Issue,
  Newick,
  NO_LABEL,
  type RowIndex,
  Stockholm,
  StringTable,
  type TreeInput,
  Trees,
} from '../../src/model/index.js';
import { pf00042, SH2_A3M, SH2_NO_RF, SH2_STOCKHOLM } from './alignment-fixtures.js';

const codes = (issues: Issue[]) => issues.map((i) => `${i.code} ${i.path}`);
const row = (r: number) => r as RowIndex;
const col = (c: number) => c as ColumnIndex;

// The set arrives with MD2 (the keys of rows, ids across versions, key clashes: set.test); until then an alignment is
// built directly, at the path the set gives it, and `set.strings` is its string table.
function buildOne(input: AlignmentInput): { a: Alignment | null; set: { strings: StringTable }; issues: Issue[] } {
  const strings = new StringTable();
  const issues: Issue[] = [];
  return { a: Alignments.build(input, 0, strings, issues, '/alignments/0'), set: { strings }, issues };
}
const flagged = (a: Alignment | null) => Array.from(a?.insertColumn ?? []).flatMap((x, c) => (x ? [c] : []));
const residuesOf = (a: Alignment, r: number) => {
  const s = a.rows.slotOf(row(r));
  return String.fromCharCode(...a.rows.residues.subarray(a.rows.residueStart[s], a.rows.residueStart[s + 1]));
};

/** The naive model: walk the gapped string (A3M: lower case and `.` are not columns). */
function naive(gapped: string, a3m: boolean) {
  const colToPos: number[] = [];
  const posToCol: number[] = [];
  for (const ch of gapped.replace(/\s/g, '')) {
    const insert = a3m && (ch === '.' || (ch >= 'a' && ch <= 'z'));
    if (a3m && ch === '.') continue;
    const gap = ch === '-' || ch === '.';
    if (!gap) posToCol.push(insert ? -1 : colToPos.length);
    if (!insert) colToPos.push(gap ? -1 : posToCol.length - 1);
  }
  return { colToPos, posToCol };
}

function checkAgainstNaive(a: Alignment, gapped: string[], a3m: boolean) {
  for (let r = 0; r < gapped.length; r++) {
    const n = naive(gapped[r], a3m);
    expect(Array.from({ length: a.columnCount }, (_, c) => a.rows.colToPos(row(r), col(c)))).toEqual(n.colToPos);
    const out = new Int32Array(a.rows.lengthOf(row(r)));
    expect(a.rows.posToCol(row(r), out)).toBe(true);
    expect(Array.from(out)).toEqual(n.posToCol);
    expect(Array.from(a.rows.posMap(row(r)) ?? [])).toEqual(n.posToCol);
    // colOf, the allocation-free read a hover takes (IN3): the map at every position, -1 past the end.
    expect(Array.from({ length: n.posToCol.length + 1 }, (_, p) => a.rows.colOf(row(r), p))).toEqual([
      ...n.posToCol,
      -1,
    ]);
  }
}

describe('alignments', () => {
  it('takes PF00017-shaped insert columns from RF, and the same from case without RF', () => {
    const [rec] = Stockholm.parse(SH2_STOCKHOLM).records;
    const fromRF = buildOne(Stockholm.alignment(rec, 'sh2'));
    expect(flagged(fromRF.a)).toEqual([10, 11, 22]);
    expect(fromRF.issues.map((i) => i.message)).toEqual([
      `Alignment "sh2": 3 insert columns are given (Stockholm's RF); letters are upper case.`,
    ]);
    const [noRF] = Stockholm.parse(SH2_NO_RF).records;
    const fromCase = buildOne(Stockholm.alignment(noRF, 'sh2'));
    expect(flagged(fromCase.a)).toEqual([10, 11, 22]);
    expect(codes(fromCase.issues)).toEqual(['model.alignment.inserts /alignments/0']);

    const a = fromRF.a as Alignment;
    expect([a.rowCount, a.columnCount, a.alphabet, a.lazy, a.referenceRow]).toEqual([4, 33, 'protein', false, 0]);
    // Insert columns stay columns: PTN11's `kd` are its residues 10 and 11, GRB2's column 11 is a gap.
    expect([a.rows.colToPos(row(1), col(10)), a.rows.colToPos(row(1), col(12))]).toEqual([10, 12]);
    expect([a.rows.colToPos(row(2), col(11)), a.rows.colToPos(row(2), col(12))]).toEqual([-1, 11]);
    expect(residuesOf(a, 1)).toBe('WFHPNITGVEKDAENLLLTRGVDGSFLARPS');
    expect(Array.from(a.rows.rowFirst)).toEqual([151, 6, 60, 333]);
    expect(fromRF.set.strings.get(a.labelOf(row(3)))).toBe('P85A_HUMAN/333-360');
    expect(a.labelOf(row(4))).toBe(-1);
    checkAgainstNaive(a, rec.rows, false);
  });

  it('leaves all-lower-case input alone: lower-case letters, not inserts', () => {
    const { a, issues } = buildOne({ key: 'x', flavor: 'aligned', gapped: ['acg-u', 'ac-gu'] });
    expect([a?.insertColumn, a?.alphabet]).toEqual([null, 'rna']);
    expect(issues.map((i) => i.message)).toEqual([
      `Alignment "x": lower case marks no insert columns; letters are upper case.`,
    ]);
  });

  it('reads A3M: lower case is an insert, not a column, and a row equals its ungapped input', () => {
    const { records } = Fasta.parse(SH2_A3M);
    const { input } = Fasta.alignment(records, 'sh2', 'a3m');
    const { a, issues } = buildOne(input);
    expect(codes(issues)).toEqual(['model.alignment.inserts /alignments/0']);
    if (!a) throw new Error('not built');
    expect([a.columnCount, a.insertColumn]).toEqual([30, null]);
    input.gapped.forEach((g, r) => {
      expect(residuesOf(a, r)).toBe(g.replace(/[-.]/g, '').toUpperCase());
    });
    // PTN11: `kd` inserted after column 9; GRB2: `e` after 9 and `g` after 19.
    const { inserts, insertStart } = a.rows;
    expect(Array.from(inserts.subarray(insertStart[1], insertStart[3]))).toEqual([9, 2, 9, 1, 19, 1]);
    expect([a.rows.colToPos(row(1), col(10)), a.rows.colToPos(row(2), col(20))]).toEqual([12, 22]);
    checkAgainstNaive(a, input.gapped, true);
  });

  it('matches the naive model across checkpoints, on PF00042 and on A3M rows with inserts everywhere', () => {
    const rows = pf00042();
    const { a } = buildOne({
      key: 'globin',
      flavor: 'aligned',
      gapped: rows.map((r) => r.gapped),
      labels: rows.map((r) => `${r.label}/${r.start}-${r.start + 1000}`),
    });
    if (!a) throw new Error('not built');
    expect(a.columnCount).toBeGreaterThan(64);
    expect(Array.from(a.rows.rowFirst)).toEqual(rows.map((r) => r.start));
    checkAgainstNaive(
      a,
      rows.map((r) => r.gapped),
      false,
    );
    // Seeded A3M rows: inserts at the start, at checkpoint borders, in runs and at the end.
    let seed = 7;
    const rand = () => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed / 0x7fffffff;
    };
    const gapped = Array.from({ length: 12 }, () => {
      let s = '';
      for (let c = 0; c < 200; c++) {
        while (rand() < 0.15) s += rand() < 0.2 ? '.' : 'acdefgh'[Math.floor(rand() * 7)];
        s += rand() < 0.3 ? '-' : 'ACDEFGH'[Math.floor(rand() * 7)];
      }
      return `${s}kl`;
    });
    const b = buildOne({ key: 'r', flavor: 'a3m', gapped }).a as Alignment;
    checkAgainstNaive(b, gapped, true);
  });

  it('holds only named rows when lazy: -1 and false for rows not held, labels from the source', () => {
    const input: AlignmentInput = {
      key: 'big',
      flavor: 'aligned',
      gapped: [],
      lazy: {
        rowCount: 1_000_000,
        columnCount: 6,
        rows: new Map([
          [70_000, { label: 'HBB_HUMAN/1-5', gapped: 'MV-HLT' }],
          [3, { label: 'HBA_HUMAN/2-6', gapped: '-VLSPA' }],
        ]),
        labelOf: (r) => (r === 5 ? 'SEEN/1-9' : null),
      },
    };
    const { a, set, issues } = buildOne(input);
    expect(issues).toEqual([]);
    if (!a) throw new Error('not built');
    expect([a.lazy, a.rowCount, a.rows.count, a.insertColumn]).toEqual([true, 1_000_000, 2, null]);
    expect([a.rows.slotOf(row(3)), a.rows.slotOf(row(70_000)), a.rows.slotOf(row(4))]).toEqual([0, 1, -1]);
    expect([a.rows.colToPos(row(70_000), col(3)), a.rows.colToPos(row(3), col(0))]).toEqual([2, -1]);
    expect(a.rows.colToPos(row(4), col(1))).toBe(-1);
    expect(a.rows.posToCol(row(4), new Int32Array(10))).toBe(false);
    expect([a.rows.posMap(row(4)), a.rows.lengthOf(row(4))]).toEqual([null, -1]);
    expect(Array.from(a.rows.rowFirst)).toEqual([2, 1]);
    const label = (r: number) => set.strings.get(a.labelOf(row(r)));
    expect([label(3), a.labelOf(row(5)), a.labelOf(row(6)), a.labelOf(row(-1))]).toEqual(['HBA_HUMAN/2-6', -1, -1, -1]);
    // A row not held has its label as text, from the source, never interned.
    const size = set.strings.size;
    const texts = [3, 5, 6, -1].map((r) => a.labelText(row(r)));
    expect([texts, set.strings.size]).toEqual([['HBA_HUMAN/2-6', 'SEEN/1-9', null, null], size]);
    // The row refs (`big/999999`, by key and by text, nothing interned) are the set's: set.test with MD2.
  });

  it('rejects a gapped row with a character outside ASCII, and reads a nucleotide look-alike as protein', () => {
    // `Ł` (U+0141) was truncated to `A` in the row's residues, and `ł` to `a`: the row loaded with no complaint.
    for (const row of ['AŁ-', 'A-ł', 'A😀', 'AC\u0141']) {
      const { a, issues } = buildOne({ key: 'x', flavor: 'aligned', gapped: ['AC-', row.padEnd(3, 'G').slice(0, 3)] });
      expect(a).toBeNull();
      expect(codes(issues)).toEqual(['model.residues.invalid /alignments/0/gapped/1']);
    }
    const lazy = {
      rowCount: 2,
      columnCount: 3,
      rows: new Map([[1, { label: 'a', gapped: 'AŁG' }]]),
      labelOf: () => null,
    };
    expect(codes(buildOne({ key: 'x', flavor: 'aligned', gapped: [], lazy }).issues)).toEqual([
      'model.residues.invalid /alignments/0/lazy/rows/1',
    ]);
    // Rows that are 90% ACGTUN with an E are protein with a warning, not DNA that then rejects the E.
    const peptide = buildOne({ key: 'tag', flavor: 'aligned', gapped: ['AAAAAAAAAGE', 'AAAAAAAA-GE'] });
    expect(peptide.a?.alphabet).toBe('protein');
    expect(codes(peptide.issues)).toEqual(['model.alphabet.detected /alignments/0']);
    expect(buildOne({ key: 'n', flavor: 'aligned', gapped: ['ACGTACGTAC', 'ACGTAC-TAC'] }).a?.alphabet).toBe('dna');
  });

  it('reports ragged rows, bad letters, shapes, descending headers, parents and the reference row', () => {
    const base = { key: 'x', flavor: 'aligned' as const };
    expect(codes(buildOne({ ...base, gapped: ['AC-', 'ACGT'] }).issues)).toEqual([
      'model.alignment.ragged /alignments/0/gapped/1',
    ]);
    expect(codes(buildOne({ ...base, gapped: [] }).issues)).toEqual(['model.alignment.shape /alignments/0/gapped']);
    // A bad letter names its row, as a non-ASCII character does (the row is found on the error path only).
    expect(codes(buildOne({ ...base, gapped: ['AC%'] }).issues)).toEqual([
      'model.residues.invalid /alignments/0/gapped/0',
    ]);
    const third = buildOne({ ...base, gapped: ['MKV', 'MK-', '-K%'] }).issues;
    expect(codes(third)).toEqual(['model.residues.invalid /alignments/0/gapped/2']);
    expect(third[0].message).toContain('row 2: "%" is not a protein letter');
    expect(codes(buildOne({ ...base, gapped: ['AC'], labels: ['a', 'b'] }).issues)).toEqual([
      'model.alignment.shape /alignments/0/labels',
    ]);
    expect(codes(buildOne({ ...base, gapped: ['AC'], insertColumns: new Uint8Array(3) }).issues)).toEqual([
      'model.alignment.shape /alignments/0/insertColumns',
    ]);
    expect(codes(buildOne({ ...base, gapped: ['AC'], referenceRow: 1 }).issues)).toEqual([
      'model.alignment.shape /alignments/0/referenceRow',
    ]);
    const lazy = {
      rowCount: 2,
      columnCount: 3,
      rows: new Map([[1, { label: 'a', gapped: 'AC' }]]),
      labelOf: () => null,
    };
    expect(codes(buildOne({ ...base, gapped: [], lazy }).issues)).toEqual([
      'model.alignment.ragged /alignments/0/lazy/rows/1',
    ]);
    const outside = { ...lazy, rows: new Map([[2, { label: 'a', gapped: 'ACG' }]]) };
    expect(codes(buildOne({ ...base, gapped: [], lazy: outside }).issues)).toEqual([
      'model.alignment.shape /alignments/0/lazy/rows',
    ]);

    const r = buildOne({
      ...base,
      gapped: ['AC', 'AC', 'AC'],
      labels: ['a/10-9', 'b/3-4', 'c/9-5'],
      parentRows: [null, 0, 2],
      rowNumberingFirst: undefined,
    });
    expect(codes(r.issues)).toEqual([
      'model.numbering.descending /alignments/0',
      'model.alignment.parent /alignments/0/parentRows/2',
    ]);
    expect(Array.from(r.a?.rows.rowFirst ?? [])).toEqual([NO_LABEL, 3, NO_LABEL]);
    expect(Array.from(r.a?.parent ?? [])).toEqual([-1, 0, -1]);
    const given = buildOne({ ...base, gapped: ['AC', 'AC'], labels: ['a/5-6', 'b'], rowNumberingFirst: [null, 40] });
    expect(Array.from(given.a?.rows.rowFirst ?? [])).toEqual([NO_LABEL, 40]);
  });
});

describe('trees', () => {
  const labels = ['SRC_HUMAN/151-180', 'PTN11_HUMAN/6-36', 'GRB2_HUMAN/60-90', 'P85A_HUMAN/333-360'];
  const parse = (text: string) => {
    const { tree } = Newick.parse(text);
    if (!tree) throw new Error('not parsed');
    return tree;
  };
  const withTree = (tree: TreeInput, rowLabels = labels) =>
    buildOne({ key: 'sh2', flavor: 'aligned', gapped: rowLabels.map(() => 'AC'), labels: rowLabels, tree });

  it('matches leaves by label or the part before `/start-end`, and visits clades', () => {
    const parsed = parse('((SRC_HUMAN:0.1,GRB2_HUMAN/60-90:0.2)90:0.05,(PTN11_HUMAN,P85A_HUMAN):0.1);');
    const leafRow = Trees.matchLeaves(parsed, labels, 'labels');
    expect(Array.from(leafRow)).toEqual([-1, -1, 0, 2, -1, 1, 3]);
    const { a, set, issues } = withTree({ ...parsed, leafRow });
    expect(issues).toEqual([]);
    const tree = a?.tree;
    if (!tree) throw new Error('no tree');
    expect(Array.from(tree.leafOrder)).toEqual([0, 2, 1, 3]);
    expect(set.strings.get(tree.label[1])).toBe('90');
    const visit = (node: number) => {
      const rows: number[] = [];
      Trees.cladeRows(tree, node, (r, out) => out.push(r), rows);
      return rows;
    };
    expect([visit(0), visit(1), visit(4), visit(3)]).toEqual([[0, 2, 1, 3], [0, 2], [1, 3], [2]]);
    const byOrder = Trees.matchLeaves(parse('((w,x),(y,z));'), labels, 'order');
    expect(Array.from(byOrder)).toEqual([-1, -1, 0, 1, -1, 2, 3]);
  });

  it('puts nodes given in any order into preorder', () => {
    // Root 3 with children 0 (a leaf, row 1) and 2, whose children are 1 (row 0) and 4 (row 2).
    const input: TreeInput = {
      parent: Int32Array.of(3, 2, 3, -1, 2),
      branchLength: Float64Array.of(1, 2, 3, Number.NaN, 4),
      label: ['b', 'a', null, 'root', 'c'],
      leafRow: Int32Array.of(1, 0, -1, -1, 2),
    };
    const { a, set } = withTree(input, ['a', 'b', 'c']);
    const tree = a?.tree;
    if (!tree) throw new Error('no tree');
    expect(Array.from(tree.parent)).toEqual([-1, 0, 0, 2, 2]);
    expect(Array.from(tree.label, (l) => set.strings.get(l))).toEqual(['root', 'b', '', 'a', 'c']);
    expect(Array.from(tree.branchLength).slice(1)).toEqual([1, 3, 2, 4]);
    expect(Array.from(tree.leafOrder)).toEqual([1, 0, 2]);
  });

  it('requires every row to be exactly one leaf', () => {
    const leaves = (text: string, rowLabels = labels) => {
      const parsed = parse(text);
      const r = withTree({ ...parsed, leafRow: Trees.matchLeaves(parsed, rowLabels, 'labels') }, rowLabels);
      expect(r.a?.tree).toBeNull();
      return r.issues.map((i) => `${i.code} ${i.path}: ${i.message}`);
    };
    expect(leaves('(SRC_HUMAN,GRB2_HUMAN,PTN11_HUMAN,P85A_HUMAN,HBA_HUMAN);')).toEqual([
      `model.tree.leaves /alignments/0/tree: Tree: every row must be exactly one leaf: the leaf "HBA_HUMAN" names no row.`,
    ]);
    expect(leaves('(SRC_HUMAN,SRC_HUMAN/151-180,PTN11_HUMAN);')[0]).toMatch(
      /row 0 has several leaves; 2 rows have no leaf\.$/,
    );
    expect(leaves('(A,B);', ['A/1-5', 'A/7-9'])[0]).toMatch(/the leaf "A" names several rows; .*2 rows have no leaf/);
    const shape = (input: TreeInput) => codes(withTree(input, ['a']).issues);
    const one = { branchLength: new Float64Array(2), label: [null, 'a'], leafRow: Int32Array.of(-1, 0) };
    expect(shape({ ...one, parent: Int32Array.of(-1, -1) })).toEqual(['model.tree.shape /alignments/0/tree']);
    expect(shape({ ...one, parent: Int32Array.of(1, 0) })).toEqual(['model.tree.shape /alignments/0/tree']);
    expect(shape({ ...one, parent: Int32Array.of(-1, 0), leafRow: Int32Array.of(0, 0) })).toEqual([
      'model.tree.leaves /alignments/0/tree',
    ]);
  });

  it('keeps the alignment when its tree fails', () => {
    const { a } = withTree({
      parent: Int32Array.of(-1),
      branchLength: new Float64Array(1),
      label: ['x'],
      leafRow: Int32Array.of(-1),
    });
    expect([a?.rowCount, a?.tree]).toEqual([4, null]);
  });
});

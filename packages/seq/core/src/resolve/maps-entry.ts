/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A `map` node and a `numbering` node resolved (_spec/README.md "Resolve", Maps and Numberings): a map's linear part in
// its encoding, its priority and a structure map's facts; a numbering's labels, offset or runs, or a source's plan.
// Apart from model.ts, which places them and keeps them by hash, so neither file grows past its size bound.

import { CategoryColumn, type Issue, issueOf, type Node, type Value } from '@molstar/core/repr/tree';
import { labelsNumbering } from '@molstar/seqx-core/query';
import {
  codingLength,
  codonLinear,
  codonRuns,
  DefaultGeneticCode,
  type GeneticCode,
  GeneticCodes,
  hasUnsplitExon,
} from './codon.js';
import { ResolveIssues } from './issues.js';
import { composeLinear, expand, ImplicitLinear, isEmptyLinear, type Linear, toRuns } from './linear.js';
import { offsetNumbering, segmentsNumbering } from './maps-numbering.js';
import {
  observedDisagrees,
  type StructureMapData,
  segmentRuns,
  segmentsLinear,
  structureData,
} from './maps-segments.js';
import type { MapEntry, NumberingEntry, ResolveConfig, SystemEntry } from './model.js';
import type { Placed } from './model-molecules.js';
import type { UnitEntry } from './molecules.js';
import type { Step } from './paths.js';

const p = (n: Node, key: string): Value | undefined => n.params?.[key];

/**
 * A declared numbering: inline labels, an offset or runs as a numbering of the system (`length` its positions, the
 * open length while unknown); a source's, where the source's rows are.
 */
export function numberingEntry(
  at: Placed,
  system: SystemEntry,
  source: string | undefined,
  length: number,
): NumberingEntry {
  const labels = p(at.node, 'labels');
  const scheme = p(at.node, 'scheme');
  const field = p(at.node, 'field');
  const id = at.node.id ?? '';
  const base = { address: at.address, system: system.address, id, ...(typeof scheme === 'string' ? { scheme } : {}) };
  const method = p(at.node, 'method');
  if (method === 'offset') return { ...base, numbering: offsetNumbering(id, length, Number(p(at.node, 'by') ?? 0)) };
  if (method === 'segments')
    return { ...base, numbering: segmentsNumbering(id, length, p(at.node, 'segments') as Record<string, Value>) };
  if (labels instanceof CategoryColumn) return { ...base, numbering: labelsNumbering(id, labels) };
  if (!source) return base;
  return {
    ...base,
    source,
    field: typeof field === 'string' ? field : 'label',
    // a source read by region, the system by the name the source knows it by
    access: { kind: 'source', source, system: system.name, shape: 'ranges' },
  };
}

type MapCode =
  | 'seq.resolve.map.empty'
  | 'seq.resolve.map.observed'
  | 'seq.resolve.map.codon_length'
  | 'seq.resolve.map.table'
  | 'seq.resolve.length.unknown';

export function mapEntry(
  at: Placed,
  from: SystemEntry,
  to: SystemEntry,
  config: ResolveConfig,
  unit?: UnitEntry,
): { map: MapEntry; issues: Issue[] } {
  const issues: Issue[] = [];
  const report = (code: MapCode, message: string) =>
    issues.push(issueOf(ResolveIssues, code, message, { address: at.address, hint: ResolveIssues.codes[code].hint }));
  const method = p(at.node, 'method');
  let linear: Linear;
  let step: Step = 0;
  let code: GeneticCode | undefined;
  let structure: StructureMapData | undefined;
  if (method === 'segments' || method === 'structure') {
    const table = p(at.node, 'segments') as Record<string, Value> | undefined;
    linear = segmentsLinear(segmentRuns(table), from, to);
    if (method === 'structure') structure = structureData(table, p(at.node, 'status'));
    // unobserved is the unit's own fact: a map onto a unit that says otherwise is reported, the unit winning
    if (structure && unit && observedDisagrees(structure, linear, { length: to.length, observed: unit.observed }))
      report('seq.resolve.map.observed', 'The map’s unobserved residues are not its unit’s');
  } else if (method === 'codon') {
    step = 1;
    const exons = p(at.node, 'exons') as unknown as { start: ArrayLike<number>; end: ArrayLike<number> };
    const strand = p(at.node, 'strand') === '-' ? '-' : '+';
    const size = to.cyclic ? to.length : undefined;
    const unsplit = hasUnsplitExon(exons.start, exons.end, size);
    if (unsplit) report('seq.resolve.length.unknown', 'An exon crosses the origin of a gene of unknown length');
    const coding = codingLength(codonRuns(exons.start, exons.end, strand, size));
    const n = from.length;
    // a coding sequence missing a piece is not judged until it is whole
    if (!unsplit && (coding % 3 !== 0 || (n !== undefined && coding !== 3 * n && coding !== 3 * n + 3)))
      report('seq.resolve.map.codon_length', 'The exons are not three nucleotides a residue');
    const table = p(at.node, 'table');
    code = GeneticCodes[typeof table === 'number' ? table : DefaultGeneticCode];
    if (!code) report('seq.resolve.map.table', 'The genetic code of the map is not built');
    linear = codonLinear(exons.start, exons.end, strand, size);
    // what lies past either sequence maps nowhere (EC-21)
    if (to.length !== undefined) linear = composeLinear(linear, new ImplicitLinear(0, 1, 0, to.length));
    if (n !== undefined) linear = composeLinear(new ImplicitLinear(0, 1, 0, 3 * n), linear);
  } else {
    const by = Number(p(at.node, 'by') ?? 0);
    const lf = from.length ?? Infinity;
    if (to.cyclic && to.length !== undefined && Number.isFinite(lf))
      linear = new ImplicitLinear(by, 1, 0, lf, to.length);
    else linear = new ImplicitLinear(by, 1, Math.max(0, -by), Math.min(lf, (to.length ?? Infinity) - by));
  }
  if (isEmptyLinear(linear)) report('seq.resolve.map.empty', 'The map pairs no position');
  if (config.encoding !== 'cheapest' && linear.runs())
    linear = config.encoding === 'runs' ? toRuns(linear) : expand(linear, linear.lo, linear.hi);
  const priority = p(at.node, 'priority');
  const map: MapEntry = {
    address: at.address,
    kind: method === 'codon' || method === 'segments' || method === 'structure' ? method : 'offset',
    from: from.address,
    to: to.address,
    encoding: linear.encoding,
    step,
    linear,
    ...(code ? { code } : {}),
    priority: typeof priority === 'number' ? priority : 0,
    ...(structure ? { structure } : {}),
  };
  return { map, issues };
}

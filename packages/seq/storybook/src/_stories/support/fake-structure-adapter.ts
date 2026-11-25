/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import {
  abortError,
  type Coverage,
  type FormatRuntime,
  type OpenContext,
  type OpenedSource,
  type RegionAnswer,
  type SourceEntry,
  SourceIssues,
} from '@molstar/seqx-core/data';
import { type StructureDescription, structureNodes } from '@molstar/seqx-core/io/structure';
import type { ChainTable } from './fake-structure/chain-table.js';
import type { FakeStructureViewer } from './fake-structure/fake-viewer.js';

// The fake structure viewer's data adapter (seq-data-ports, The fake host's adapter): the viewer's chain table, and the
// model it holds beside it, read into a `StructureDescription`: the one converter (`structureNodes`) the mmCIF runtime
// uses turns it into the nodes, so a file's reading and the host's give equal nodes. It serves the sources whose `uri`
// (or `ref`) names the entry the viewer holds, answers from memory (a resolved promise, no asset read), serves the
// units' values by `field` from the table's per-residue columns, and says when the viewer's "Drop chain" edits the
// table (`changes`). Any other source, and one naming an `assembly` (the viewer holds the deposited model only), goes
// to the fallback: the default adapter, reading through the assets port.

/** The value channels a table's residues carry, the first one it has the default (`plddt` before `b_factor`). */
const Channels = [
  ['plddt', (r: { plddt?: number }) => r.plddt],
  ['b_factor', (r: { b?: number }) => r.b],
] as const;

/** The deposited model's operator in the table: a chain read without an assembly has none. */
const Deposited = '1_555';

/**
 * The description of a chain table: the entities (the model's, else read from the chains), the components a position,
 * the chain instances (one a copy, under the operator it has), their observed residues with author numbers and values,
 * and from the model its compounds, bonds and a molecule's atoms.
 */
export function descriptionOf(table: ChainTable): StructureDescription {
  const model = table.model;
  const held = new Set(table.chains.map((c) => c.entity));
  const polymerOf = (entity: string) => table.chains.some((c) => c.entity === entity && c.residues[0]?.label != null);
  const entities = model
    ? model.entities.filter((e) => held.has(e.id))
    : [...held].map((id) => ({ id, type: polymerOf(id) ? 'polymer' : 'non-polymer' }));
  const entityRow = new Map(entities.map((e, i) => [e.id, i]));
  const firstChain = (entity: string) => table.chains.find((c) => c.entity === entity)!;

  // components a position: a polymer's from its first chain (every residue, modeled or not), a ligand's one component
  const compEntity: number[] = [];
  const compPosition: number[] = [];
  const compCode: string[] = [];
  const identifiers = { entity: [] as number[], namespace: [] as string[], value: [] as string[] };
  entities.forEach((e, k) => {
    const chain = firstChain(e.id);
    if (e.type === 'polymer') {
      chain.residues.forEach((r, i) => {
        compEntity.push(k);
        compPosition.push(r.label != null ? r.label - 1 : i);
        compCode.push(r.comp);
      });
    } else if (chain.residues[0]) {
      compEntity.push(k);
      compPosition.push(0);
      compCode.push(chain.residues[0].comp);
      if (e.type === 'non-polymer') {
        identifiers.entity.push(k);
        identifiers.namespace.push('ccd');
        identifiers.value.push(chain.residues[0].comp);
      }
    }
  });

  // chain instances: a row a copy; the deposited model's single operator is no operator
  const rows: { chain: (typeof table.chains)[number]; operator: string | undefined }[] = [];
  for (const c of table.chains) {
    if (!entityRow.has(c.entity)) continue;
    const ops = c.operators.length === 1 && c.operators[0] === Deposited ? [undefined] : c.operators;
    for (const operator of ops) rows.push({ chain: c, operator });
  }
  const rowsOf = (asym: string) =>
    rows.flatMap((r, k) => (r.chain.label_asym_id === asym ? [[k, r.operator] as const] : []));

  // observed residues, a chain's in order (a ligand's one residue)
  const resChain: number[] = [];
  const resPosition: number[] = [];
  const resAuth: number[] = [];
  const resIns: string[] = [];
  const values: Record<string, number[]> = {};
  const have = Channels.filter(([, get]) => table.chains.some((c) => c.residues.some((r) => get(r) !== undefined)));
  for (const [name] of have) values[name] = [];
  rows.forEach(({ chain }, k) => {
    const polymer = entities[entityRow.get(chain.entity)!].type === 'polymer';
    chain.residues.forEach((r, i) => {
      if (!r.observed || (!polymer && i > 0)) return;
      const position = polymer ? (r.label != null ? r.label - 1 : i) : 0;
      const auth = Number(r.auth);
      resChain.push(k);
      resPosition.push(position);
      resAuth.push(Number.isNaN(auth) ? position + 1 : auth);
      resIns.push(r.ins);
      for (const [name, get] of have) values[name].push(get(r) ?? Number.NaN);
    });
  });

  // bonds between the chains kept, a row under each of the copies they share
  const link = {
    a_chain: [] as number[],
    a_position: [] as number[],
    b_chain: [] as number[],
    b_position: [] as number[],
    type: [] as string[],
    a_atom: [] as (string | null)[],
    b_atom: [] as (string | null)[],
  };
  for (const b of model?.bonds ?? []) {
    for (const [ra, operator] of rowsOf(b.a.asym)) {
      const rb = rowsOf(b.b.asym).find(([, o]) => o === operator);
      if (!rb) continue;
      link.a_chain.push(ra);
      link.a_position.push(b.a.position);
      link.b_chain.push(rb[0]);
      link.b_position.push(b.b.position);
      link.type.push(b.type);
      link.a_atom.push(b.a.atom ?? null);
      link.b_atom.push(b.b.atom ?? null);
    }
  }

  const used = new Set(compCode);
  const compounds = (model?.compounds ?? []).filter((c) => used.has(c.id));
  const atoms = (model?.atoms ?? []).flatMap((a) =>
    entityRow.has(a.entity) ? [{ k: entityRow.get(a.entity)!, a }] : [],
  );
  return {
    name: table.entry,
    entities: {
      id: entities.map((e) => e.id),
      type: entities.map((e) => e.type as 'polymer'),
      polymerType: entities.map((e) => ('polymerType' in e ? e.polymerType : undefined)),
      name: entities.map((e) => ('name' in e ? e.name : undefined)),
      letters: entities.map((e) => ('letters' in e ? e.letters : undefined)),
    },
    components: { entity: Int32Array.from(compEntity), position: Int32Array.from(compPosition), component: compCode },
    identifiers: {
      entity: Int32Array.from(identifiers.entity),
      namespace: identifiers.namespace,
      value: identifiers.value,
    },
    chains: {
      label_asym_id: rows.map((r) => r.chain.label_asym_id),
      auth_asym_id: rows.map((r) => r.chain.auth_asym_id),
      entity: Int32Array.from(rows, (r) => entityRow.get(r.chain.entity)!),
      operator: rows.map((r) => r.operator),
    },
    residues: {
      chain: Int32Array.from(resChain),
      position: Int32Array.from(resPosition),
      auth_seq_id: Int32Array.from(resAuth),
      ins_code: resIns,
    },
    links: {
      a_chain: Int32Array.from(link.a_chain),
      a_position: Int32Array.from(link.a_position),
      b_chain: Int32Array.from(link.b_chain),
      b_position: Int32Array.from(link.b_position),
      type: link.type,
      a_atom: link.a_atom,
      b_atom: link.b_atom,
    },
    compounds: {
      id: compounds.map((c) => c.id),
      name: compounds.map((c) => c.name),
      type: compounds.map((c) => c.type),
      formula: compounds.map((c) => c.formula),
      parents: compounds.map((c) => (c.parents ? [...c.parents] : undefined)),
    },
    atoms: {
      entity: Int32Array.from(atoms.flatMap(({ k, a }) => a.atoms.map(() => k))),
      residue: Int32Array.from(atoms.flatMap(({ a }) => a.atoms.map((x) => x[0]))),
      name: atoms.flatMap(({ a }) => a.atoms.map((x) => x[1])),
      element: atoms.flatMap(({ a }) => a.atoms.map((x) => x[2])),
    },
    values: Object.fromEntries(Object.entries(values).map(([k, v]) => [k, Float32Array.from(v)])),
  };
}

/** The entry a `uri` names: its last path segment without a query or the `.gz`, `.cif`, `.bcif` endings, lowercase. */
export function entryOfUri(uri: string): string {
  const name = uri.split(/[?#]/)[0].split('/').pop() ?? '';
  return name
    .replace(/\.gz$/i, '')
    .replace(/\.(b?cif)$/i, '')
    .toLowerCase();
}

/** Whether the viewer holds the data a source names: by its `uri`'s entry or its `ref`. */
function holds(viewer: FakeStructureViewer, source: SourceEntry): boolean {
  const entry = viewer.entry.toLowerCase();
  return (
    (source.ref !== undefined && source.ref.toLowerCase() === entry) ||
    (source.uri !== undefined && entryOfUri(source.uri) === entry)
  );
}

/** A source the viewer holds, opened from its table as it is now. */
function serve(viewer: FakeStructureViewer, source: SourceEntry, signal: AbortSignal): OpenedSource {
  const table = viewer.table;
  const at = viewer.edits.value;
  const read = structureNodes(descriptionOf(table), { source: source.address });
  const fields = Channels.map(([name]) => name).filter((f) => read.values(f));
  const lengths = new Map(read.systems.map((s) => [s.name, s.length ?? 0]));
  const coverage: Coverage = {
    systems: read.systems,
    shapes: fields.length
      ? read.systems.map((s) => ({
          system: s.name,
          shape: 'values' as const,
          fields: fields.map((name) => ({ name, type: 'float32' as const })),
        }))
      : [],
    levels: [1],
    nodes: read.nodes,
  };
  const notServed = (why: string): RegionAnswer => ({
    ok: false,
    kind: 'not-served',
    issue: issueOf(SourceIssues, 'seq.source.not_served', why, { address: source.address }),
  });
  return {
    coverage,
    issues: [...read.issues],
    read(request, readSignal): Promise<RegionAnswer> {
      if (readSignal.aborted || signal.aborted) return Promise.reject(abortError());
      const length = lengths.get(request.system);
      if (length === undefined) return Promise.resolve(notServed(`${table.entry} has no unit ${request.system}`));
      if (request.level !== 1) return Promise.resolve(notServed(`${table.entry} serves level 1 only`));
      if (!request.shapes.includes('values')) return Promise.resolve({ ok: true, tables: {} });
      const field = request.field ?? fields[0];
      const values = field === undefined ? undefined : read.values(field)?.get(request.system);
      if (!values)
        return Promise.resolve(notServed(`${table.entry} has no value ${request.field ?? ''} on ${request.system}`));
      const start = Math.max(0, Math.min(length, request.start));
      const end = Math.max(start, Math.min(length, request.end));
      return Promise.resolve({
        ok: true,
        tables: { values: { form: 'positions', start, value: values.slice(start, end) } },
      });
    },
    changes(fn) {
      let on = true;
      const off = viewer.onEdit(() => on && fn());
      // an edit since this reading was opened is told once, soon after
      if (viewer.edits.value !== at) queueMicrotask(() => on && fn());
      return () => {
        on = false;
        off();
      };
    },
    dispose() {},
  };
}

/**
 * The data adapter of a fake structure viewer for `format` (mmCIF's by default): it serves the sources naming the
 * entry the viewer holds from the viewer's table, every other source from `fallback` (the default adapter).
 */
export function fakeStructureAdapter(
  viewer: FakeStructureViewer,
  fallback: FormatRuntime,
  format = 'mmcif',
): FormatRuntime {
  return {
    name: 'fake-structure',
    formats: [format],
    open(source: SourceEntry, ctx: OpenContext): Promise<OpenedSource> {
      if (ctx.signal.aborted) return Promise.reject(abortError());
      // the viewer holds the deposited model: a source asking for an assembly is the file's to read
      const assembly = (source as { readonly assembly?: unknown }).assembly;
      if (assembly !== undefined || !holds(viewer, source)) return fallback.open(source, ctx);
      return Promise.resolve(serve(viewer, source, ctx.signal));
    },
  };
}

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The node kinds and their parameter tables (seqspec:format#nodes-and-refs, seqspec:format#data-nodes,
// seqspec:format#views): the source of the validator, the normalizer and the reference. A kind's own rules
// (`check`) are those a table cannot say.

import type { Geometry } from '../doc.js';
import { seqqlQuery } from '../query.js';
import type { SvsRegistries } from '../registry.js';
import {
  bool,
  dict,
  type FieldType,
  float,
  int,
  isNamespaced,
  isPlainObject,
  list,
  literal,
  nullable,
  ptr,
  ref,
  refined,
  report,
  str,
  tuple,
  union,
} from '../schema/field-schema.js';
import {
  type Fields,
  OptionalField,
  object,
  RequiredField,
  SimpleParamsSchema,
  UnionParamsSchema,
  type ValuesFor,
} from '../schema/params-schema.js';
import { type IndexedNode, nodeList, type TreeContext, TreeSchema } from '../schema/tree-schema.js';
import { KernelLayerTypes } from './kernel.js';
import {
  filter,
  httpsUrl,
  integerLabel,
  keyFilter,
  label,
  orderField,
  provenance,
  themeRef,
  uri,
} from './param-types.js';
import { sectionStyle } from './style.js';
import { labelSides, table, tableHas } from './table.js';

const opt = <V>(type: FieldType<V>, description = '') => OptionalField(nullable(type), null, description);
const req = RequiredField;
const hidden = OptionalField(bool, false, '');
const labelP = opt(str, 'null: from the source, else the ref');
const sourceRef = opt(ref('source'), 'read from a source node');
const prov = opt(provenance);
const alphabet = opt(literal('protein', 'dna', 'rna', 'generic'), 'null: detected');
const range = tuple(int, int);
const direction = literal(1, -1);
/** A share, 0 to 1 (`redundancy`'s threshold: 0.9, not 90). */
const share: FieldType<number> = {
  name: 'float (0..1)',
  check(v, path, ctx) {
    float.check(v, path, ctx);
    if (typeof v === 'number' && (v < 0 || v > 1)) {
      report(ctx, 'svs.params.type', path, `${v} is outside 0..1: write a share such as 0.9, not a percentage`);
    }
  },
};
const dataRefs = union(ref('annotation', 'computed'), list(ref('annotation', 'computed')));
/** One entry of a section's `positions` (seqspec:format#sections): a data ref, or its rows after a layer's `filter`. */
const positionEntry = union(
  ref('annotation', 'computed'),
  object({
    data: req(ref('annotation', 'computed'), 'what it covers'),
    filter: req(filter, "only the rows the filter keeps, as a layer's `filter` does"),
  }),
);
const positionsOf = union(positionEntry, list(positionEntry));
/** A section's location set shown persistently: a highlight channel of the viewer's each. */
const highlightSet = object({
  data: req(ref('annotation', 'computed'), "its items' extents light"),
  items: opt(list(str), 'only these item ids; null: every item'),
  style: OptionalField(literal('band', 'highlight', 'accent'), 'band', 'band: full height behind the tracks'),
});
const simple = <F extends Fields>(fields: F) => SimpleParamsSchema(fields);
const molecule = <F extends Fields>(value: F) =>
  simple({
    ...value,
    source_ref: sourceRef,
    record: opt(union(str, int), 'a record id or 0-based index; required when the source has several'),
    chem_comp_ids: opt(list(nullable(str)), 'a chemical component per position'),
    label: labelP,
    identifiers: OptionalField(
      list(object({ namespace: req(str, ''), value: req(str, ''), version: OptionalField(nullable(str), null, '') })),
      [],
      'uniprot:P04637, chembl:CHEMBL941; never imply a mapping',
    ),
    provenance: prov,
  });

const node = (ctx: TreeContext, r: unknown) => (typeof r === 'string' ? ctx.nodes.get(r) : undefined);
const children = (n: IndexedNode, kind: string) => nodeList(n.node.children).filter((c) => c?.kind === kind);
const at = (n: IndexedNode, key: string) => ptr(ptr(n.path, 'params'), key);
/**
 * A source by `protocol`: its provider has no manifest of its own, so its readers declare the shape
 * (seqspec:sources#the-tiles-protocol).
 */
const isProtocol = (ctx: TreeContext, sourceRef: unknown) => node(ctx, sourceRef)?.params.protocol != null;
function needs(n: IndexedNode, ctx: TreeContext, key: string, when: string) {
  if (n.params[key] == null) report(ctx, 'svs.params.missing', at(n, key), `"${key}" is required ${when}`);
}

/** The geometry of an annotation or computed node; undefined when not known statically. */
function geometryOf(n: IndexedNode | undefined, reg: SvsRegistries): Geometry | undefined {
  if (n?.kind === 'annotation') return n.params.geometry as Geometry;
  if (n?.kind !== 'computed') return undefined;
  const p = n.params.provider;
  if (p === 'column_stats') return n.params.measure === 'consensus' ? 'ranges' : 'values';
  if (p === 'profile') return 'matrix';
  if (p === 'coverage' || p === 'paired' || p === 'query') return 'ranges';
  if (p === 'row_stats' || p === 'redundancy') return 'rows';
  return typeof p === 'string' ? reg.get('providers', p)?.geometry : undefined;
}
/** The data refs of a section's `positions`: each entry's ref or its `data`. */
const positionRefs = (v: unknown): string[] =>
  (Array.isArray(v) ? v : v == null ? [] : [v]).map((e) => (typeof e === 'string' ? e : String(e?.data)));
const refsOf = (v: unknown): string[] => (typeof v === 'string' ? [v] : Array.isArray(v) ? v : []);

const annotationCommon = {
  target: opt(ref('sequence', 'ligand', 'unit', 'alignment', 'row'), "every row's target without a target column"),
  numbering: opt(str, 'required with label columns: their numbering'),
  source_ref: sourceRef,
  uri: opt(uri, 'read from a URI directly; format required with it'),
  format: opt(str, 'a registered format'),
  from_source: opt(ref('sequence', 'ligand', 'complex', 'alignment'), 'read from the file that node came from'),
  category_name: opt(str, 'CIF: the category; null: the first'),
  block_header: opt(str, 'CIF: the data block'),
  block_index: OptionalField(int, 0, "CIF: the data block's index when block_header is null"),
  tag: opt(str, 'Stockholm: GC:<tag>, GR:<tag>, GS:<tag>; FASTA: header'),
  fields: OptionalField(dict(nullable(str)), {}, 'column name -> name in the file, null to ignore one'),
  label: opt(str),
  semantic_type: opt(str, 'open, namespaced: uniprot:Domain, so:0001067'),
  length: opt(int, 'required for a protocol source'),
  links: OptionalField(
    list(object({ label: req(str, ''), url: req(httpsUrl, 'an https:// template; {field} is a column of the row') })),
    [],
    'per-row links',
  ),
  provenance: prov,
};
const annotationCase = <F extends Fields>(g: Geometry, more: F = {} as F) =>
  simple({ ...annotationCommon, value: opt(table(g), 'the table'), ...more });

const side = ref('sequence', 'ligand', 'unit', 'row');
/** A provider's endpoints: any a mapping takes, and alignments (`codon_rows` pairs their rows). */
const providerSide = ref('sequence', 'ligand', 'unit', 'row', 'alignment');
const mappingCase = <F extends Fields>(f: F) =>
  simple({ from: req(side, ''), to: req(side, ''), ...f, provenance: prov });

const stats = {
  alignment: req(ref('alignment'), 'the alignment'),
  row_subset: opt(list(int), 'a subgroup of rows'),
  gaps: OptionalField(literal('ignore', 'count'), 'count', ''),
};

const layerCommon = {
  data: opt(dataRefs, 'annotations or computed nodes, drawn as one; absent with annotation children'),
  blocks: opt(list(ref('block')), 'only these blocks'),
  color: opt(themeRef('color'), "null: the type's default"),
  glyph: opt(themeRef('glyph'), "null: the type's default"),
  domain_group: opt(str, 'layers of the view in one group share their normalized domains (the union, per field)'),
  opacity: OptionalField(float, 1, ''),
  filter: opt(filter, 'only the rows the filter keeps'),
  fallback: opt(str, "a kernel type or omit; null: the type's"),
  options: OptionalField(object(), {}, "an extension type's params"),
};
/** The params every layer type has, typed: the builder's (B1). */
export type LayerCommonParams = ValuesFor<SimpleParamsSchema<typeof layerCommon>>;

export const SvsTreeSchema = TreeSchema({
  rootKind: 'root',
  nodes: {
    root: {
      description: 'the document',
      parent: [],
      ref: 'implied',
      params: simple({ default_view: opt(ref('view'), 'null: the first view') }),
    },
    data: { description: 'everything views show', parent: ['root'], ref: 'implied', count: 'one', params: simple({}) },
    source: {
      description: 'where data comes from, fetched once however many nodes read it',
      parent: ['data'],
      exclusive: [['uri', 'service', 'host_file', 'value', 'stream']],
      params: simple({
        uri: opt(uri, 'an absolute https:// URL'),
        service: opt(str, 'a registered service'),
        id: opt(str, "the service's record"),
        resource: opt(str, 'what to read from the service'),
        host_file: opt(str, 'a name the host resolves'),
        value: opt(str, 'a whole file in a registered text format'),
        format: opt(str, 'a registered format; required with uri, host_file and value unless protocol is set'),
        protocol: opt(literal('seqstar-tiles/1'), 'a streamed source by a registered protocol, at uri'),
        stream: opt(str, "a streamed source: an id the host's stream provider resolves"),
        data_kind: opt(literal('letters', 'items', 'values', 'alignment'), "a streamed source's; required with stream"),
        coordinate_system: opt(
          ref('sequence', 'unit', 'complex', 'alignment', 'row'),
          'items, values: null, the target',
        ),
        max_depth: OptionalField(int, 4, 'how deeply its patches may nest streamed sources'),
        max_bytes: opt(int, 'a size bound; for value, null means 1 MiB'),
      }),
      check(n, ctx) {
        const kind = n.params.data_kind;
        if (n.params.stream != null) needs(n, ctx, 'data_kind', 'with stream');
        if (n.params.coordinate_system != null && (kind === 'letters' || kind === 'alignment'))
          report(ctx, 'svs.params.exclusive', at(n, 'coordinate_system'), `a ${kind} source's reader is its system`);
        if (n.params.service != null) needs(n, ctx, 'id', 'with service');
        if (n.params.service == null && n.params.protocol == null && n.params.stream == null)
          needs(n, ctx, 'format', 'without service or protocol');
        if (n.params.protocol != null) needs(n, ctx, 'uri', 'with protocol');
        const max = n.params.max_bytes;
        if (typeof max === 'number' && max <= 0)
          report(ctx, 'svs.params.type', at(n, 'max_bytes'), 'max_bytes is above 0');
      },
    },
    sequence: {
      description: "a polymer's residues",
      parent: ['data'],
      exclusive: [['value', 'source_ref']],
      params: molecule({
        value: opt(str, 'one-letter codes; no gaps'),
        alphabet,
        cyclic: OptionalField(bool, false, 'a cyclic coordinate system: positions modulo its length; ranges may wrap'),
      }),
    },
    ligand: {
      description: "a ligand's heavy atoms",
      parent: ['data'],
      exclusive: [['value', 'source_ref']],
      params: molecule({ value: opt(list(str), 'the element symbol of each heavy atom') }),
    },
    numbering: {
      description: 'a numbering of its parent',
      parent: ['sequence', 'ligand', 'unit'],
      params: UnionParamsSchema('method', 'how labels are given', {
        offset: simple({
          name: req(str, "unique among its parent's numberings; ordinal and row are reserved"),
          first: OptionalField(int, 1, 'the label of index 0'),
          skip_zero: OptionalField(bool, false, 'the numbering has no 0'),
          default: OptionalField(bool, false, 'the numbering rulers use when a section names none'),
        }),
        segments: simple({
          name: req(str, ''),
          segments: req(list(tuple(int, int, int)), '[index_start, label_start, length]'),
          default: OptionalField(bool, false, ''),
        }),
        labels: simple({
          name: req(str, ''),
          labels: req(list(nullable(integerLabel)), 'one per position; null: unlabeled'),
          default: OptionalField(bool, false, ''),
        }),
        names: simple({
          name: req(str, ''),
          labels: req(list(nullable(label)), 'one name per position; null: unlabeled'),
          default: OptionalField(bool, false, ''),
        }),
      }),
      check(n, ctx) {
        const name = n.params.name;
        if (name === 'ordinal' || name === 'row') {
          report(ctx, 'svs.params.type', at(n, 'name'), `"${name}" is reserved (built in); choose another name`);
        }
        const parent = ctx.nodes.get(n.parent ?? '');
        const same = parent ? children(parent, 'numbering').filter((c) => c.params?.name === name) : [];
        if (same.length > 1 && same[0] !== n.node) {
          report(ctx, 'svs.ref.duplicate', at(n, 'name'), `numbering "${name}" is already defined on "${n.parent}"`);
        }
      },
    },
    unit: {
      description: 'an instance of a sequence or a ligand',
      parent: ['data', 'complex'],
      exclusive: [['sequence', 'ligand']],
      params: simple({
        sequence: opt(ref('sequence'), 'the polymer it instantiates'),
        ligand: opt(ref('ligand'), 'the ligand it instantiates'),
        range: opt(range, 'the shown positions; null: all'),
        observed: opt(list(range), 'observed positions; null: all of range'),
        label: opt(str, "the block's name"),
        origin: opt(str, 'experimental, predicted, docked, designed, <ns>:<name>; sequence is reserved'),
        score: opt(float, 'a score the instance carries'),
      }),
    },
    complex: {
      description: "units shown together: a provider's, or its unit children",
      parent: ['data'],
      params: simple({
        provider: opt(str, '<namespace>:<name>, a registered complex provider; or unit children'),
        source_ref: opt(ref('source'), 'required with provider; without, a streamed letters source'),
        options: OptionalField(object(), {}, "the provider's, checked by its registered table"),
        origin: opt(str, "given to every unit; null: the provider's, or each unit's own"),
        label: opt(str, "null: the provider's, else the ref"),
      }),
      owned: { param: 'options', schema: (p, reg) => reg.get('complexProviders', String(p.provider))?.options },
      check(n, ctx) {
        const p = n.params.provider;
        // A complex is a provider's, its unit children, or a letters stream's (`source_ref` without a provider).
        const streamed = p == null && node(ctx, n.params.source_ref)?.params.stream != null;
        if (+(p != null) + +(children(n, 'unit').length > 0) + +streamed !== 1) {
          report(ctx, 'svs.params.exclusive', n.path, 'give one of provider, unit children and a letters stream');
        }
        if (p != null && !isNamespaced(p))
          report(ctx, 'svs.params.type', at(n, 'provider'), 'a provider is <ns>:<name>');
        if (p != null) needs(n, ctx, 'source_ref', 'with provider');
      },
    },
    alignment: {
      description: 'an alignment, a coordinate system',
      parent: ['data'],
      exclusive: [['value', 'source_ref', 'provider']],
      params: simple({
        value: opt(list(str), 'the gapped rows'),
        flavor: OptionalField(literal('aligned', 'a3m'), 'aligned', 'a3m: lower case an insert, not a column'),
        labels: opt(list(str), "value's row labels"),
        source_ref: sourceRef,
        record: opt(int, 'which alignment of a multi-alignment Stockholm file'),
        alphabet,
        row_count: opt(int, 'required for a protocol source'),
        column_count: opt(int, 'required for a protocol source'),
        reference_row: opt(int, 'the query row; null: row 0'),
        parent_rows: opt(list(nullable(int)), 'per row: the row it is a variant of'),
        provider: opt(str, 'a registered node provider (realign, or <namespace>:<name>): computed from inputs'),
        options: OptionalField(object(), {}, "the provider's, checked by its registered table"),
        inputs: opt(
          union(list(ref('sequence', 'unit', 'row')), ref('alignment')),
          "with provider: what it aligns, in order: sequences, units or rows, or an alignment's rows",
        ),
        label: opt(str),
      }),
      owned: { param: 'options', schema: (p, reg) => reg.get('nodeProviders', String(p.provider))?.options },
      check(n, ctx) {
        const { provider, inputs } = n.params;
        if (provider != null) needs(n, ctx, 'inputs', 'with provider');
        else if (inputs != null) report(ctx, 'svs.params.exclusive', at(n, 'inputs'), 'inputs go with a provider');
        if (!isProtocol(ctx, n.params.source_ref)) return;
        needs(n, ctx, 'row_count', 'for a protocol source');
        needs(n, ctx, 'column_count', 'for a protocol source');
      },
    },
    tree: {
      description: "a tree over its alignment's rows",
      parent: ['alignment'],
      count: 'at_most_one',
      exclusive: [['value', 'source_ref', 'provider']],
      params: simple({
        value: opt(str, 'Newick'),
        source_ref: sourceRef,
        leaves: OptionalField(literal('labels', 'order'), 'labels', 'labels: by row label; order: leaf k is row k'),
        provider: opt(
          str,
          'a registered node provider (nj, upgma, or <namespace>:<name>): computed from the alignment',
        ),
        options: OptionalField(object(), {}, "the provider's, checked by its registered table"),
      }),
      owned: { param: 'options', schema: (p, reg) => reg.get('nodeProviders', String(p.provider))?.options },
      check(n, ctx) {
        const aln = ctx.nodes.get(n.parent ?? '');
        const src = node(ctx, aln?.params.source_ref);
        if (src && (src.params.stream != null || src.params.protocol != null) && n.params.leaves !== 'order')
          report(ctx, 'svs.tree.streamed', at(n, 'leaves'), 'a streamed alignment\'s tree needs leaves: "order"');
      },
    },
    mapping: {
      description: "positions of one coordinate system paired with another's",
      parent: ['data'],
      params: UnionParamsSchema('method', 'how positions are paired', {
        identity: mappingCase({ from_range: opt(range, 'the from indices it covers; null: all that fit') }),
        segments: mappingCase({
          segments: req(list(tuple(int, int, nullable(int))), '[from_start, to_start, length]'),
        }),
        pairwise: mappingCase({
          cigar: req(
            refined('str', (s) => /^([0-9]*[MID=X])+$/.test(s), 'svs.params.type', 'CIGAR (M I D = X)'),
            '',
          ),
        }),
        align: mappingCase({
          aligner: OptionalField(str, 'biomsa', 'a registered aligner'),
          options: OptionalField(object(), {}, "the aligner's params, checked by its table; {}: its defaults"),
        }),
        by_numbering: mappingCase({
          numbering: req(str, 'pair equal labels'),
          to_numbering: opt(str, 'null: the same name'),
        }),
        provider: simple({
          from: req(providerSide, ''),
          to: req(providerSide, ''),
          provider: req(str, 'a registered node provider: codon_rows, or <namespace>:<name>'),
          options: OptionalField(object(), {}, "the provider's, checked by its registered table"),
          provenance: prov,
        }),
        codon: mappingCase({
          from_range: opt(range, 'the residues the CDS codes, the first with the first codon; null: from residue 0'),
          cds: opt(list(range), 'the coding ranges on to (exons), in transcript order; null: one range from start'),
          start: OptionalField(int, 0, 'with cds null: the first index of the coding range'),
          strand: OptionalField(literal('+', '-'), '+', "-: the ranges' reverse complement"),
          genetic_code: OptionalField(int, 1, 'an NCBI genetic code id (gc.prt; 1 the standard code)'),
        }),
      }),
      owned: {
        param: 'options',
        schema: (p, reg) =>
          p.method === 'align'
            ? reg.get('aligners', String(p.aligner ?? 'biomsa'))?.options
            : p.method === 'provider'
              ? reg.get('nodeProviders', String(p.provider))?.options
              : undefined,
      },
    },
    annotation: {
      description: "a table of items; under a layer, that layer's own data",
      parent: ['data', 'layer'],
      exclusive: [['value', 'source_ref', 'uri', 'from_source']],
      params: UnionParamsSchema('geometry', "the table's geometry", {
        ranges: annotationCase('ranges'),
        values: annotationCase('values'),
        pairs: annotationCase('pairs', { b_numbering: opt(str, "side b's numbering; null: as numbering") }),
        matrix: annotationCase('matrix', { states: req(list(str), "the matrix's second axis") }),
        grid: annotationCase('grid', { b_numbering: opt(str, "side b's numbering; null: as numbering") }),
        rows: annotationCase('rows'),
      }),
      check(n, ctx) {
        const p = n.params;
        if (p.uri != null) needs(n, ctx, 'format', 'with uri');
        for (const k of isProtocol(ctx, p.source_ref) ? ['length', 'target'] : [])
          needs(n, ctx, k, 'for a protocol source');
        if (p.value == null) return;
        // `fields` renames a file's columns; an inline table names its own.
        if (isPlainObject(p.fields) && Object.keys(p.fields).length > 0) {
          report(ctx, 'svs.params.exclusive', at(n, 'fields'), "fields renames a file's columns: give value or fields");
        }
        const s = labelSides(p.value);
        const msg = 'label columns need a numbering: set numbering or write indices';
        if (s.a && p.numbering == null && !s.numbering) report(ctx, 'svs.position.form', at(n, 'value'), msg);
        if (s.b && p.numbering == null && p.b_numbering == null && !s.numbering && !s.bNumbering) {
          report(ctx, 'svs.position.form', at(n, 'value'), `side b: ${msg}`);
        }
      },
    },
    computed: {
      description: 'an annotation computed by a provider',
      parent: ['data'],
      params: UnionParamsSchema(
        'provider',
        'a kernel provider or <namespace>:<name>',
        {
          coverage: simple({
            of: opt(ref('unit'), 'null: every unit'),
            reference: opt(ref('sequence'), 'conflicts against it through a mapping'),
          }),
          column_stats: simple({
            ...stats,
            measure: req(literal('conservation', 'consensus', 'identity', 'gap_fraction'), ''),
            reference_row: opt(int, "null: the alignment's"),
          }),
          profile: simple(stats),
          row_stats: simple({
            alignment: req(ref('alignment'), 'the alignment'),
            reference_row: opt(int, "null: the alignment's"),
            items: opt(ref('annotation', 'computed'), 'count the items of this annotation on each row'),
            category: opt(str, 'with items: only this category'),
            columns: opt(range, 'with items: only within these columns, [start, stop)'),
            score: opt(str, 'with items: the maximum of this field of the items on the row'),
          }),
          redundancy: simple({
            alignment: req(ref('alignment'), 'the alignment'),
            threshold: req(share, 'an identity from 0 to 1: a row joins the first representative this identical to it'),
            order: OptionalField(orderField, null, 'the order rows are taken in; null: data order'),
          }),
          paired: simple({
            pairs: req(list(ref('annotation')), 'pair annotations'),
            partner: opt(list(ref('unit')), 'null: every endpoint (an interface)'),
            categories: opt(list(str), 'only these pair categories'),
          }),
          query: simple({
            query: req(seqqlQuery, 'a SeqQL query: { seqql, expr, text }'),
            label: opt(str, 'null: the ref'),
          }),
        },
        { open: simple({ options: OptionalField(object(), {}, "the provider's, checked by its registered table") }) },
      ),
      owned: {
        param: 'options',
        schema: (p, reg) => (isNamespaced(p.provider) ? reg.get('providers', p.provider)?.options : undefined),
      },
    },
    view: {
      description: 'one view',
      parent: ['root'],
      count: 'some',
      params: simple({
        label: req(str, 'what a reader picks it by'),
        layout: OptionalField(literal('tracks', 'text'), 'tracks', 'text: the blocks read as text, wrapping'),
      }),
    },
    section: {
      description: 'one axis of the view with its tracks',
      parent: ['view'],
      count: 'some',
      params: simple({
        label: opt(str, 'shown when a view has several sections'),
        arrangement: OptionalField(literal('side', 'stacked'), 'side', 'blocks side by side, or one lane each'),
        copies: OptionalField(literal('one', 'each'), 'one', 'one block per copy group, or one per unit'),
        numbering: opt(str, 'the numbering rulers and position labels use'),
        ligand_form: OptionalField(literal('position', 'atoms'), 'position', 'a ligand as one position or its atoms'),
        positions: opt(positionsOf, 'show only what these cover'),
        scroll_group: opt(str, 'sections of one group share horizontal scroll and zoom'),
        height: OptionalField(literal('content', 'box'), 'box', 'content: as tall as its rows'),
        height_rows: opt(int, 'a hint for box'),
        row_sizing: OptionalField(literal('uniform', 'natural'), 'uniform', ''),
        style: opt(
          sectionStyle,
          'a preset name (tracks, msa, strip, compact, wrapped, <ns>:<name>) or an inline preset',
        ),
        direction: OptionalField(direction, 1, "its blocks' direction, where a block gives none"),
        follow: opt(ref('section'), "this section and that one follow each other's window"),
        follow_anchor: opt(ref('row', 'unit'), "an alignment axis's anchor; null: its reference_row"),
        highlights: OptionalField(list(highlightSet), [], 'persistent highlights, a channel each'),
        hidden,
      }),
      check(n, ctx) {
        const blocks = children(n, 'block');
        const columns = blocks.filter((b) => b.params?.alignment != null).length;
        if (columns > 1 || (columns === 1 && blocks.length > 1)) {
          report(ctx, 'svs.section.axis', n.path, "an axis is blocks of units or one alignment's columns, not both");
        }
        const pos = n.params.positions;
        positionRefs(pos).forEach((r, i) => {
          if (geometryOf(node(ctx, r), ctx.registries) !== 'rows') return;
          const path = Array.isArray(pos) ? ptr(at(n, 'positions'), i) : at(n, 'positions');
          report(ctx, 'svs.section.positions', path, `"${r}" is a rows annotation; it covers no positions`);
        });
      },
    },
    block: {
      description: "one block of the section's axis",
      parent: ['section'],
      exclusive: [['unit', 'sequence', 'ligand', 'complex', 'alignment']],
      params: simple({
        unit: opt(ref('unit')),
        sequence: opt(ref('sequence'), 'shown through its own unit'),
        ligand: opt(ref('ligand'), 'shown through its own unit'),
        complex: opt(ref('complex'), 'its units in its order, a block each'),
        alignment: opt(ref('alignment'), 'its columns'),
        label: opt(str, "null: the unit's, molecule's or alignment's"),
        numbering: opt(str, "overrides the section's"),
        range: opt(
          range,
          'with unit, sequence or ligand: the piece shown, [start, stop) of the unit (start > stop wraps a cyclic origin); null: all',
        ),
        direction: opt(direction, "-1: the piece from its highest index down; null: the section's"),
        origin: OptionalField(int, 0, 'a cyclic unit without range: the index its piece starts at'),
      }),
      check(n, ctx) {
        const r = n.params.range as [number, number] | null;
        const o = (n.params.origin as number | undefined) ?? 0;
        const what = n.params.complex != null ? "a complex's units" : "an alignment's columns";
        const whole = n.params.complex != null || n.params.alignment != null;
        if (o !== 0 && (whole || r != null)) {
          const why = r != null ? 'a range starts where it starts' : `${what} have none`;
          report(ctx, 'svs.params.exclusive', at(n, 'origin'), `origin is a whole cyclic unit's; ${why}`);
        } else if (o < 0) report(ctx, 'svs.position.range', at(n, 'origin'), `origin ${o} is negative`);
        if (r == null) return;
        if (whole) {
          report(ctx, 'svs.params.exclusive', at(n, 'range'), `range is a unit's piece; ${what} take positions`);
        } else if (r[0] < 0 || r[1] < 0 || r[1] === r[0]) {
          // start > stop wraps a cyclic origin: whether the unit is cyclic is known at resolve.
          report(
            ctx,
            'svs.position.range',
            at(n, 'range'),
            `[${r[0]}, ${r[1]}] is empty or negative: give start < stop`,
          );
        }
      },
    },
    group: {
      description: 'an expandable group of tracks',
      parent: ['section', 'group'],
      params: simple({ label: req(str, ''), collapsed: OptionalField(bool, false, 'the state at open'), hidden }),
    },
    track: {
      description: 'one row of the section',
      parent: ['section', 'group'],
      params: simple({
        label: opt(str, "null: the first layer's data label"),
        height: opt(float, 'a hint, in rows'),
        hidden,
        expanded: opt(bool, "a split layer's rows at open; null: expanded for variants, else collapsed"),
        summary: OptionalField(
          bool,
          false,
          "its group's summary: shown under the group's header only while the group is collapsed (a union of what it holds)",
        ),
        for_each: opt(
          object({
            kind: req(
              literal('sequence', 'ligand', 'unit', 'mapping', 'annotation', 'item_system'),
              "the selected nodes' kind (item_system: a streamed source's items' coordinate systems, each made for it)",
            ),
            source: req(ref('source'), "the streamed source they come from (a nested one's derived ref too)"),
            depth: opt(int, '0: its own nodes, 1: those of the sources it nests, ...; null: any'),
            reaching: OptionalField(bool, true, "only nodes a path places on one of the section's blocks"),
            order: OptionalField(
              orderField,
              null,
              "null, data: key order; else a key table's field ({data, field, direction}) or a list of them, over a table keyed by node",
            ),
            filter: OptionalField(list(keyFilter), [], 'instances kept when every filter keeps them'),
          }),
          'a track template: one track per dynamic node it selects',
        ),
      }),
      check(n, ctx) {
        const ls = children(n, 'layer');
        const split = ls.filter((l) => l.params?.split != null).length;
        if (split > 1 || (split && ls.some((l) => l.params?.type === 'rows')))
          report(ctx, 'svs.layer.split', n.path, 'split one layer of a track, in a track without rows; move the rest');
      },
    },
    layer: {
      description: 'what a track shows',
      parent: ['track'],
      count: 'some',
      params: UnionParamsSchema(
        'type',
        'a kernel layer type or <namespace>:<name>',
        Object.fromEntries(
          KernelLayerTypes.map((t) => [
            t.name,
            simple({ ...layerCommon, ...(t.options as SimpleParamsSchema).fields }),
          ]),
        ),
        { open: simple(layerCommon) },
      ),
      owned: {
        param: 'options',
        schema: (p, reg) => (isNamespaced(p.type) ? reg.get('layerTypes', p.type)?.options : undefined),
      },
      check(n, ctx) {
        const p = n.params;
        const own = children(n, 'annotation');
        if (p.data != null && own.length > 0) {
          report(ctx, 'svs.layer.annotation', at(n, 'data'), 'give data or annotation children, not both');
        }
        const kernel = (t: unknown) => KernelLayerTypes.some((k) => k.name === t);
        if (p.fallback != null && p.fallback !== 'omit' && !kernel(p.fallback)) {
          report(ctx, 'svs.layer.fallback', at(n, 'fallback'), `a fallback is a kernel layer type or "omit"`);
        }
        if (p.of != null && p.row != null) {
          report(ctx, 'svs.params.exclusive', at(n, 'of'), "give of (another polymer's letters) or row, not both");
        }
        if (p.type === 'dot_plot' && typeof p.aspect_ratio === 'number' && !(p.aspect_ratio > 0)) {
          const m = `aspect_ratio ${p.aspect_ratio} is not above 0: write the plot's width ÷ height (1: square, 2: twice as wide), or null for the track's own height`;
          report(ctx, 'svs.params.type', at(n, 'aspect_ratio'), m);
        }
        // A split by a column an inline table lacks (other data is checked when resolved); a field provider's name is a field.
        const by0 = (p.split as { by?: string } | null)?.by;
        const by = by0 && ctx.registries.get('fields', by0) ? undefined : by0;
        for (const d of [...refsOf(p.data).map((r) => node(ctx, r)), ...own]) {
          const v = d?.params?.value;
          if (by && v != null && !tableHas(v, by))
            report(ctx, 'svs.layer.split', at(n, 'split'), `"${by}" is not a column of "${d?.ref}": give one it has`);
        }
        if (isNamespaced(p.type) && p.fallback == null) {
          report(ctx, 'svs.layer.fallback', at(n, 'fallback'), `extension type "${p.type}" needs a fallback`);
        }
        // An extension type's data is checked against its fallback whether or not the type is registered, so a
        // document is valid for every reader, and against the registered type too.
        const reg = ctx.registries;
        const own2 = own.map((c) => c.params?.geometry as Geometry | undefined);
        const geometries = [...refsOf(p.data).map((r) => geometryOf(node(ctx, r), reg)), ...own2];
        for (const name of isNamespaced(p.type) ? [p.fallback, p.type] : [p.type]) {
          const type = typeof name === 'string' ? reg.get('layerTypes', name) : undefined;
          if (!type) continue;
          for (const g of geometries) {
            if (type.accepts !== 'none' && (g === undefined || type.accepts.includes(g))) continue;
            const accepts = type.accepts === 'none' ? 'no data' : type.accepts.join(', ');
            report(ctx, 'svs.layer.data', at(n, 'data'), `${type.name} takes ${accepts}; got ${g ?? 'data'}`);
          }
        }
      },
    },
    focus: {
      description: "the section's initial range",
      parent: ['section'],
      exclusive: [['target', 'annotation']],
      params: simple({
        target: opt(ref('unit', 'sequence', 'ligand', 'alignment', 'row'), 'with start, stop or beg, end'),
        start: opt(int, 'an index; columns for an alignment'),
        stop: opt(int, 'an index'),
        beg: opt(label, 'a label in numbering'),
        end: opt(label, 'a label in numbering'),
        numbering: opt(str, "the labels' numbering"),
        annotation: opt(ref('annotation', 'computed'), "its items' extent"),
        items: opt(list(str), 'item ids; null: every item'),
      }),
      check(n, ctx) {
        const p = n.params;
        const labels = p.beg != null || p.end != null;
        if (labels && (p.start != null || p.stop != null)) {
          report(ctx, 'svs.position.form', n.path, 'use start, stop (indices) or beg, end (labels), not both');
        } else if (labels && p.numbering == null) {
          report(ctx, 'svs.position.form', at(n, 'numbering'), 'beg, end need a numbering');
        }
      },
    },
  },
});

export type SvsKind = keyof typeof SvsTreeSchema.nodes;

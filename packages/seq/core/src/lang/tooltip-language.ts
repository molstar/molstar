/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FormatTemplate } from '@molstar/core/util/string-format';

// The `tooltip` kind (_spec/README.md, Language, Tooltips): which words attach to which data, read the same by every
// reader; the application composes the card. Under `views` (every view) or a `view`; two ways to attach, one result:
// by data (`on` a data node or a layer: each item's text on its own marks) and by range from a table (`target` with
// inline `rows` or a `source`: each row's text on every mark its location covers). MolViewSpec's `tooltip`,
// `tooltip_from_uri` and `tooltip_from_source` map onto it (`schema`, `field`, `text`, `block`, `category`). Spread into
// `SeqLanguage`'s kinds; what the description cannot say is `checkTooltip`.

/** How a tooltip table's rows locate: Seq*'s own, and MolViewSpec's schemas, each read as MVS reads it. */
export const TooltipSchemas = [
  'ranges',
  'address',
  'whole_structure',
  'entity',
  'chain',
  'auth_chain',
  'residue',
  'auth_residue',
  'residue_range',
  'auth_residue_range',
  'atom',
  'auth_atom',
  'all_atomic',
] as const;

export type TooltipSchema = (typeof TooltipSchemas)[number];

/**
 * The columns each MVS schema reads (MolViewSpec's annotation schemas): a row's other columns are ignored, as MVS
 * ignores them (a `beg_label_seq_id` under `chain` is whole chains).
 */
export const MvsSchemaColumns: Readonly<Record<Exclude<TooltipSchema, 'ranges' | 'address'>, readonly string[]>> = {
  whole_structure: [],
  entity: ['label_entity_id'],
  chain: ['label_entity_id', 'label_asym_id'],
  auth_chain: ['auth_asym_id'],
  residue: ['label_entity_id', 'label_asym_id', 'label_seq_id'],
  auth_residue: ['auth_asym_id', 'auth_seq_id', 'pdbx_PDB_ins_code'],
  residue_range: ['label_entity_id', 'label_asym_id', 'beg_label_seq_id', 'end_label_seq_id'],
  auth_residue_range: ['auth_asym_id', 'beg_auth_seq_id', 'end_auth_seq_id'],
  atom: ['label_entity_id', 'label_asym_id', 'label_seq_id', 'label_atom_id', 'type_symbol', 'atom_id', 'atom_index'],
  auth_atom: [
    'auth_asym_id',
    'auth_seq_id',
    'pdbx_PDB_ins_code',
    'auth_atom_id',
    'type_symbol',
    'atom_id',
    'atom_index',
  ],
  all_atomic: [
    'label_entity_id',
    'label_asym_id',
    'label_seq_id',
    'label_atom_id',
    'auth_asym_id',
    'auth_seq_id',
    'pdbx_PDB_ins_code',
    'auth_atom_id',
    'type_symbol',
    'atom_id',
    'atom_index',
    'beg_label_seq_id',
    'end_label_seq_id',
    'beg_auth_seq_id',
    'end_auth_seq_id',
  ],
};

/** A tooltip's defaults: `{}` reads the `tooltip` field; the text is Markdown (inline). */
export const TooltipDefaults = { field: 'tooltip', text: '{}', format: 'markdown' } as const;

export const TooltipKinds = {
  tooltip: {
    doc: 'Words on rendered marks: by data (`on` a data node or a layer, each item’s text on its own marks) or by range from a table (`target` with inline `rows` or a `source`, each row’s text on every mark its location covers, in any layer of the view). Under `views` every view’s, under a `view` that view’s; shown by the application’s card (a provider), never on a position no layer marks.',
    parents: ['views', 'view'],
    ids: 'required',
    params: {
      on: {
        type: 'union',
        of: [
          { type: 'ref', doc: 'A data node or a layer (by data), or an address (by range: a ref, a range, Keys).' },
          { type: 'list', of: { type: 'ref', doc: 'An address.' }, 'min-length': 1, doc: 'Addresses, each a row.' },
        ],
        doc: 'By data: the items of a data node or a layer (its marks only, its filter applied) speak their own fields; by range: an address, or a list, each a one-row table (MVS’s `tooltip { text }` on a component).',
      },
      target: {
        type: 'ref',
        kinds: ['sequence', 'molecule', 'unit', 'alignment', 'assembly', 'complex'],
        doc: 'By range: the system a table’s locations are on (an assembly: a row names its sequence; a complex or a molecule: MVS’s schemas).',
      },
      rows: {
        type: 'table',
        columns: {
          start: { type: 'column', of: 'int32', doc: '`ranges`: the first position of each row.' },
          end: { type: 'column', of: 'int32', doc: '`ranges`: the position after the last.' },
          name: { type: 'column', of: 'category', doc: '`ranges` over an assembly: each row’s sequence.' },
          address: {
            type: 'list',
            of: { type: 'string', doc: 'An address.' },
            doc: '`address`: each row’s location, relative to `target`.',
          },
        },
        doc: 'By range: the tooltip’s own table inline (location columns and text fields, a column each; instead of `source`).',
      },
      source: {
        type: 'ref',
        kinds: ['source'],
        doc: 'By range: a table read from a `source` under `data` (MVS’s `tooltip_from_uri`; instead of `rows`).',
      },
      schema: {
        type: 'string',
        enum: TooltipSchemas,
        doc: 'How a table’s rows locate: `ranges` (absent: `start`, `end`, `name` on an assembly), `address` (an `address` column), or MolViewSpec’s schemas, each read as MVS reads it (a row’s other columns ignored).',
      },
      block: {
        type: 'union',
        of: [
          { type: 'string', doc: 'A block’s header.' },
          { type: 'integer', min: 0, doc: 'A block’s index.' },
        ],
        doc: 'A CIF table’s block: MVS’s `block_header` or `block_index`.',
      },
      category: { type: 'string', doc: 'A CIF table’s category (MVS’s `category_name`).' },
      field: { type: 'string', doc: 'What `{}` reads (MVS’s `field_name`): `tooltip` when absent.' },
      text: {
        type: 'string',
        doc: 'An f-string over the row’s fields (MVS’s `text_format`): `"{}"` when absent; a row with a field it names undefined or empty gives no text. Markdown (inline) in the template only: values are escaped.',
      },
      label: { type: 'string', doc: 'The card row’s label; absent, the text is the row.' },
      format: {
        type: 'string',
        enum: ['markdown', 'text'],
        doc: '`markdown` (absent: strong, emphasis, code, in-view `!seq-…` links, breaks) or `text`, shown as written.',
      },
    },
  },
} as const;

/** A tooltip's template parsed, and the fields it reads (`{}` is `field`); or why it cannot be read. */
export type TooltipTemplate =
  | { readonly template: FormatTemplate; readonly fields: readonly string[] }
  | { readonly error: string };

/** Parses a tooltip's `text` (MVS's `text_format`): its f-string, and the fields it names (`{}` as `field`). */
export function tooltipTemplate(text: string, field: string = TooltipDefaults.field): TooltipTemplate {
  try {
    const template = FormatTemplate(text);
    const fields: string[] = [];
    // the names, read by formatting once with a number for each (a numeric spec reads it); `{}` is `field`
    template.format((name) => {
      const f = name || field;
      if (!fields.includes(f)) fields.push(f);
      return '0';
    });
    return { template, fields };
  } catch (e) {
    return { error: e instanceof Error ? e.message : String(e) };
  }
}

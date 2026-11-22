/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

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

/** A tooltip's defaults: `{}` reads the `tooltip` field (the text is Markdown, inline, unless `format` says `text`). */
export const TooltipDefaults = { field: 'tooltip', text: '{}' } as const;

export const TooltipKinds = {
  tooltip: {
    doc: 'Words on drawn marks, on the card: by data (`on`) or by range (`target` with `rows` or a `source`). Under `views` every view’s.',
    parents: ['views', 'view'],
    ids: 'required',
    params: {
      on: {
        type: 'union',
        of: [
          { type: 'ref', doc: 'A data node, a layer, or an address.' },
          { type: 'list', of: { type: 'ref', doc: 'An address.' }, 'min-length': 1, doc: 'Addresses.' },
        ],
        doc: 'By data: the items speak; by range: addresses, a row each.',
      },
      target: {
        type: 'ref',
        kinds: ['sequence', 'molecule', 'unit', 'alignment', 'assembly', 'complex'],
        doc: 'By range: the rows’ system.',
      },
      rows: {
        type: 'table',
        columns: {
          start: { type: 'column', of: 'int32', doc: 'First position.' },
          end: { type: 'column', of: 'int32', doc: 'Past the last.' },
          target: {
            type: 'column',
            of: 'category',
            doc: 'Its sequence (an assembly).',
          },
          address: {
            type: 'list',
            of: { type: 'string', doc: 'An address.' },
            doc: 'Its location under `target`.',
          },
        },
        doc: 'The table inline.',
      },
      source: {
        type: 'ref',
        kinds: ['source'],
        doc: 'A `source`’s table.',
      },
      schema: {
        type: 'string',
        enum: TooltipSchemas,
        doc: 'How rows locate: `ranges`, `address`, or MVS’s.',
      },
      block: {
        type: 'union',
        of: [
          { type: 'string', doc: 'Header.' },
          { type: 'integer', min: 0, doc: 'Index.' },
        ],
        doc: 'A CIF block.',
      },
      category: { type: 'string', doc: 'A CIF category.' },
      field: { type: 'string', doc: 'What `{}` reads (`tooltip`).' },
      text: {
        type: 'string',
        doc: 'An f-string over the row’s fields (`{}`); values escaped.',
      },
      label: { type: 'string', doc: 'The row’s label.' },
      format: {
        type: 'string',
        enum: ['markdown', 'text'],
        doc: '`markdown` or `text`.',
      },
    },
  },
} as const;

/**
 * The fields a tooltip's `text` names (MVS's `text_format`; `{}` as `field`), or why its braces cannot be read: the
 * f-string's structure as `FormatTemplate` reads it (`{{`, `}}` a brace; a tag `{name:spec}`), not its format specs.
 * Why not the parse: the check reads every spec, the formatter loads with the tooltips (`tooltip-door.ts`); a format
 * spec the formatter refuses gives no text.
 */
export function templateFields(text: string, field: string = TooltipDefaults.field): string[] | string {
  const fields: string[] = [];
  let tag = -1;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (tag < 0) {
      if ((c === '{' || c === '}') && text[i + 1] === c) i++;
      else if (c === '{') tag = i + 1;
      else if (c === '}') return 'unmatched "}"';
    } else if (c === '{') return '"{" within a tag';
    else if (c === '}') {
      const t = text.slice(tag, i);
      const colon = t.indexOf(':');
      const name = (colon < 0 ? t : t.slice(0, colon)) || field;
      if (!fields.includes(name)) fields.push(name);
      tag = -1;
    }
  }
  return tag < 0 ? fields : 'unmatched "{"';
}

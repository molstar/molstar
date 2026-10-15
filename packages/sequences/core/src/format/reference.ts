/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The written reference (seqspec:format#versioning), generated from the parameter tables and the
// registries: a section per node kind, the layer type and theme catalogs, custom tables, codes.

import { SvsVersion } from './doc.js';
import type { ChannelKind, SvsRegistries, ThemeEntry } from './registry.js';
import type { Fields, ParamsSchema } from './schema/params-schema.js';
import type { TreeSchema } from './schema/tree-schema.js';
import { GlyphShapes, LineStyles, PaletteSchema } from './tables/param-types.js';
import { SvsTreeSchema } from './tables/svs-tree.js';
import { GeometryColumns } from './tables/table.js';

const cell = (s: string) => s.replaceAll('|', '\\|');
const code = (s: string) => `\`${cell(s)}\``;
const COUNT = { one: 'exactly one', at_most_one: 'at most one', some: 'one or more' } as const;

function fieldsTable(fields: Fields, skip: ReadonlySet<string> = new Set()): string[] {
  const rows = Object.entries(fields).filter(([k]) => !skip.has(k));
  if (rows.length === 0) return ['No params.', ''];
  return [
    '| Param | Type | Default | |',
    '|---|---|---|---|',
    ...rows.map(([k, f]) => {
      const d = f.required ? 'required' : code(JSON.stringify(f.default));
      return `| ${code(k)} | ${code(f.type.name)} | ${d} | ${cell(f.description)} |`;
    }),
    '',
  ];
}

function paramsTable(s: ParamsSchema, common?: Fields): string[] {
  if (s.type === 'simple') return fieldsTable(s.fields);
  const out = [`A union on ${code(s.discriminator)}: ${s.discriminatorDescription}.`, ''];
  if (common) out.push('Every case:', '', ...fieldsTable(common));
  for (const [name, c] of Object.entries(s.cases))
    out.push(`${code(`${s.discriminator}: "${name}"`)}:`, '', ...fieldsTable(c.fields, same(c, common)));
  if (s.open)
    out.push(
      `${code(`${s.discriminator}: "<namespace>:<name>"`)}:`,
      '',
      ...fieldsTable(s.open.fields, same(s.open, common)),
    );
  return out;
}

/** A params table as the reference prints it (a union's common fields once): for another package's reference. */
export function printParamsTable(s: ParamsSchema): string {
  return paramsTable(s, commonFields(s)).join('\n');
}

/** The common fields a case has as they are (a case that overrides one, a rows layer's `filter`, lists its own). */
const same = (c: { fields: Fields }, common?: Fields) =>
  new Set(Object.keys(common ?? {}).filter((k) => c.fields[k] === common?.[k]));

/** Fields every case of a union shares (written once); a case may override one. */
function commonFields(s: ParamsSchema): Fields | undefined {
  if (s.type === 'simple') return undefined;
  const cases = [...Object.values(s.cases), ...(s.open ? [s.open] : [])];
  const first = cases[0]?.fields ?? {};
  const common = Object.entries(first).filter(
    ([k, f]) => cases.every((c) => k in c.fields) && 2 * cases.filter((c) => c.fields[k] === f).length > cases.length,
  );
  return common.length > 0 ? Object.fromEntries(common) : undefined;
}

function themes(title: string, list: readonly ThemeEntry[], kind: ChannelKind): string[] {
  const of = list.filter((t) => t.params[kind]);
  const table = (t: ThemeEntry) => paramsTable(t.params[kind] as ParamsSchema);
  return [
    `## ${title}`,
    '',
    ...of.flatMap((t) => [`### ${code(t.name)}`, '', `${t.description ?? ''}`, '', ...table(t)]),
  ];
}

/**
 * The types the tables name that are not nodes (seqspec:format#positions, seqspec:format#themes): what `Table<g>`,
 * `Palette` and the rest are.
 */
function typesSection(): string[] {
  return [
    '## Types',
    '',
    '### `Table<geometry>`',
    '',
    'Rows of named columns, written as an array of row objects or as an object of column arrays; a column may be `{dict, codes}` or `{json}`. The columns each geometry reads (a target or position column of another geometry is `svs.columns.unknown`; any other column is a field):',
    '',
    '| Geometry | Columns |',
    '|---|---|',
    ...Object.entries(GeometryColumns).map(([g, cols]) => `| ${code(g)} | ${cols.map(code).join(', ')} |`),
    '',
    '### `Palette`',
    '',
    ...paramsTable(PaletteSchema),
    '### `Glyph`',
    '',
    `A shape name (${GlyphShapes.map(code).join(', ')}) or an object with exactly one of ${code('shape')}, ${code('char')} (one grapheme) and ${code('path')} (SVG path data, at most 4 KB, with ${code('box')}: [w, h]), and optionally ${code('label')} and ${code('fallback')} (a shape).`,
    '',
    '### `LineStyle`',
    '',
    `How a range's item is drawn: ${LineStyles.map(code).join(', ')}; ${code('solid')} is the plain bar or underline.`,
    '',
    '### `Color`',
    '',
    "`#rgb`, `#rrggbb`, a Mol* color name, or a split color `a/b`. No alpha: a layer's `opacity` gives it.",
    '',
    '### `Label`',
    '',
    'A position label: an int (`315`), or a string without whitespace: `"27A"`, `"111.1"` in an integer numbering, a name (`"N13"`) in a `names` numbering.',
    '',
    '### `ThemeRef`',
    '',
    '`{theme, params, fallback}`: a registered theme with the params of its table (at most 4 deep); a `<namespace>:<name>` theme needs a kernel `fallback`.',
    '',
  ];
}

/** The reference as markdown. */
export function printReference(reg: SvsRegistries, schema: TreeSchema = SvsTreeSchema): string {
  const out = [`# SeqViewSpec ${SvsVersion} reference`, '', 'Generated from the parameter tables.', '', '## Nodes', ''];
  for (const [kind, ns] of Object.entries(schema.nodes)) {
    const facts = [`Parent: ${ns.parent.length > 0 ? ns.parent.map(code).join(', ') : 'none (the root)'}`];
    if (ns.count) facts.push(`${COUNT[ns.count]} under its parent`);
    facts.push(ns.ref === 'implied' ? `ref: implied (${code(kind)})` : 'ref: required');
    for (const set of ns.exclusive ?? []) facts.push(`exactly one of ${set.map(code).join(', ')}`);
    out.push(
      `### ${code(kind)}`,
      '',
      `${ns.description}. ${facts.join('; ')}.`,
      '',
      ...paramsTable(ns.params, commonFields(ns.params)),
    );
  }
  out.push(
    ...typesSection(),
    '## Layer types',
    '',
    '| Type | Shows | Data | Where | Axes | Sizing | Fallback |',
    '|---|---|---|---|---|---|---|',
  );
  const types = reg.list('layerTypes');
  for (const t of types) {
    const accepts = t.accepts === 'none' ? 'none' : t.accepts.join(', ');
    out.push(
      `| ${code(t.name)} | ${t.description ?? ''} | ${accepts} | ${t.layouts.join(', ')} | ${t.axes.join(', ')} | ${t.sizing} | ${t.fallback ?? 'none'} |`,
    );
  }
  out.push('', ...types.flatMap((t) => [`### ${code(t.name)}`, '', ...paramsTable(t.options)]));
  const all = reg.list('themes');
  out.push(
    ...themes('Color themes', all, 'color'),
    ...themes('Glyph themes', all, 'glyph'),
    ...themes('Line style themes', all, 'line_style'),
  );
  out.push(
    '## Letter maps',
    '',
    'A `letter_map` (`letters`, `rows`): a name below, or a dict of single letters mapped; unlisted letters are kept.',
    'A map listing alphabets leaves a layer with letters of another as stored.',
    '',
  );
  for (const m of reg.list('letterMaps'))
    out.push(
      `- ${code(m.name)}${m.alphabets ? ` (${m.alphabets.join(', ')})` : ''}: ${Object.entries(m.map)
        .map(([k, v]) => `${k}→${v}`)
        .join(' ')}`,
    );
  out.push(
    '',
    '## Color schemes',
    '',
    'A `scheme` color theme names one (`{"theme": "scheme", "params": {"name": "zappo"}}`); CB: colorblind-safe.',
    '',
    '| Scheme | Label | Kind | CB | What it shows |',
    '|---|---|---|---|---|',
  );
  for (const e of reg.list('colorSchemes')) {
    const fits = 'alphabets' in e.scheme ? `${e.scheme.kind}: ${e.scheme.alphabets.join(', ')}` : e.scheme.kind;
    out.push(`| ${code(e.name)} | ${e.label} | ${fits} | ${e.colorblind ? 'yes' : ''} | ${e.description} |`);
  }
  out.push('', '## Custom tables', '');
  for (const c of reg.list('custom'))
    out.push(`### ${code(`${c.namespace}`)} on ${code(c.kind)}`, '', ...paramsTable(c.schema));
  out.push('## Issue codes', '', '| Code | Severity |', '|---|---|');
  for (const c of reg.list('issueCodes')) out.push(`| ${code(c.code)} | ${c.severity} |`);
  return `${out.join('\n')}\n`;
}

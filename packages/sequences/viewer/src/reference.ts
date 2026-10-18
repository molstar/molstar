/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer's generated reference: its commands, behaviors, config items, the built-in registry entries by kind,
// and its own layer types and issue codes, from a viewer spec (the default one unless given). Edge-rate: docs and
// tests only. A test writes it to packages/sequences/_spec/reference/viewer.md.

import { printParamsTable, type SvsRegistryEntry } from '@molstar/seq-core/format';
import { SeqViewerConfig, type SeqViewerConfigItem } from './config.js';
import type { SeqViewerRegistryEntry } from './registry.js';
import { DefaultSeqViewerSpec, SeqViewerOptionalBehaviors } from './default-spec.js';
import type { SeqViewerSpec } from './spec.js';

const cell = (s: string) => s.replaceAll('|', '\\|').replaceAll('\n', ' ');
const code = (s: string) => `\`${cell(s)}\``;
const json = (v: unknown) => code(JSON.stringify(v));

function configItems(): [string, SeqViewerConfigItem][] {
  const out: [string, SeqViewerConfigItem][] = [];
  for (const [group, items] of Object.entries(SeqViewerConfig))
    for (const [name, item] of Object.entries(items as Record<string, SeqViewerConfigItem>))
      out.push([`${group}.${name}`, item]);
  return out;
}

function table(head: string[], rows: string[][]): string[] {
  return [
    `| ${head.join(' | ')} |`,
    `|${head.map(() => '---').join('|')}|`,
    ...rows.map((r) => `| ${r.join(' | ')} |`),
    '',
  ];
}

/** The registry entries of a spec, flattened per kind, in registration order. */
function entries<T>(spec: SeqViewerSpec, pick: (e: SeqViewerRegistryEntry) => readonly T[] | undefined): T[] {
  return (spec.registry ?? []).flatMap((e) => pick(e) ?? []);
}

/** The format halves' own entries, those of the viewer's namespace (the kernel's are in the format's reference). */
function own(spec: SeqViewerSpec): SvsRegistryEntry[] {
  return (spec.registry ?? []).flatMap((e) => (e.svs === undefined ? [] : Array.isArray(e.svs) ? e.svs : [e.svs]));
}

/** The reference as markdown. */
export function printViewerReference(spec: SeqViewerSpec = DefaultSeqViewerSpec()): string {
  const out = ['# Seq* viewer reference', '', "Generated from the viewer's tables; do not edit.", ''];

  out.push('## Commands', '', 'Run from code (`dispatch`) or as data (`run`, an invocation `{command, params}`).', '');
  for (const c of entries(spec, (e) => e.commands))
    out.push(`### ${code(c.id)}`, '', c.description || '', '', printParamsTable(c.params), '');

  out.push(
    '## Behaviors',
    '',
    'In `DefaultSeqViewerSpec` order; params are defaults, set per behavior or per preset.',
    '',
  );
  out.push(
    ...table(
      ['Name', 'Category', 'Priority', 'What it does', 'Params'],
      spec.behaviors.map(({ behavior: b, defaultParams }) => [
        code(b.name),
        b.category,
        String(b.priority),
        cell(b.display.description ?? ''),
        json({ ...b.params, ...defaultParams }),
      ]),
    ),
  );
  out.push('Built in, not in `DefaultSeqViewerSpec` (`SeqViewerOptionalBehaviors`): a host adds them.', '');
  out.push(
    ...table(
      ['Name', 'Category', 'Priority', 'What it does', 'Params'],
      Object.values(SeqViewerOptionalBehaviors).map((b) => [
        code(b.name),
        b.category,
        String(b.priority),
        cell(b.display.description ?? ''),
        json(b.params),
      ]),
    ),
  );

  out.push('## Config', '', 'Items with defaults, read by their managers at use (`ctx.config.get`, `set`).', '');
  out.push(
    ...table(
      ['Item', 'Key', 'Default', 'What it sets'],
      configItems().map(([name, i]) => [code(name), code(i.key), json(i.defaultValue), cell(i.description)]),
    ),
  );

  out.push('## Built-in registry entries', '', 'What `DefaultRegistry` registers, by kind.', '');
  out.push(
    '### Node kinds',
    '',
    ...table(
      ['Name', 'Stage'],
      entries(spec, (e) => e.nodeKinds).map((n) => [code(n.name), n.stage]),
    ),
  );
  out.push(
    '### Renderers',
    '',
    "One per layer type; a geometry it lacks draws the layer's extents as bars and is reported.",
    '',
    ...table(
      ['Layer type', 'Geometries'],
      entries(spec, (e) => e.renderers).map((r) => [code(r.name), r.geometries.join(', ')]),
    ),
  );
  out.push(
    '### Data columns',
    '',
    'An unknown column type is drawn as `text`.',
    '',
    ...table(
      ['Type', 'Width (CSS px)'],
      entries(spec, (e) => e.columns).map((c) => [code(c.name), String(c.width)]),
    ),
  );
  out.push(
    '### Themes',
    '',
    ...table(
      ['Name', 'Label', 'Category', 'Channels'],
      entries(spec, (e) => e.themes).map((t) => [code(t.name), cell(t.label), t.category, t.kinds.join(', ')]),
    ),
  );
  out.push(
    '### Navigators',
    '',
    ...table(
      ['Name'],
      entries(spec, (e) => e.navigators).map((n) => [code(n.name)]),
    ),
  );
  out.push(
    '### Presets',
    '',
    ...table(
      ['Name', 'Geometry', 'Navigator', 'Behaviors changed'],
      entries(spec, (e) => e.presets).map((p) => [
        code(p.name),
        p.geometry,
        typeof p.navigator === 'string' ? code(p.navigator) : p.navigator ? 'parts' : 'none',
        Object.keys(p.behaviors ?? {})
          .map(code)
          .join(', ') || 'none',
      ]),
    ),
  );
  out.push(
    '### Label providers',
    '',
    ...table(
      ['Name', 'Priority'],
      entries(spec, (e) => e.labels).map((l) => [code(l.name), String(l.priority)]),
    ),
  );
  out.push(
    '### Actions',
    '',
    ...table(
      ['Name', 'Title', 'Key'],
      entries(spec, (e) => e.actions).map((a) => {
        const e = a as { name: string; title?: string; key?: string };
        return [code(e.name), cell(e.title ?? ''), e.key ? code(e.key) : ''];
      }),
    ),
  );

  const svs = own(spec);
  const types = svs.flatMap((e) => e.layerTypes ?? []).filter((t) => t.name.startsWith('seqv:'));
  out.push("## The viewer's layer types", '');
  for (const t of types) out.push(`### ${code(t.name)}`, '', t.description ?? '', '', printParamsTable(t.options), '');
  const codes = svs.flatMap((e) => e.issueCodes ?? []).filter((c) => c.code.startsWith('seqv.'));
  out.push(
    "## The viewer's issue codes",
    '',
    ...table(
      ['Code', 'Severity'],
      codes.map((c) => [code(c.code), c.severity]),
    ),
  );
  return `${out
    .join('\n')
    .replace(/\n{3,}/g, '\n\n')
    .trimEnd()}\n`;
}

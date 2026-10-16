/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ToolGroupView, ToolView } from '@molstar/ai/trace';
import { k } from '@molstar/ai/trace';
import { For, onCleanup, ScrollArea, Show } from '@molstar/prism-ui';
import { distinctUntilChanged, map, skip, take } from 'rxjs';
import { classes as layout } from './_generated/inspector.module.css.js';
import { classes } from './_generated/tools.module.css.js';
import type { InspectorModel, ToolsRow } from './model.js';
import { DetailHead, Facts, firstSentence, Json, listKeys, parts, reveal, TextTabs } from './parts.js';

const HeadHeight = 40;

/**
 * The agent's tools as the model is offered them (packages/ai/_spec/trace.md § The Inspector): the core from the first
 * turn, each skill (its guidance and tools, loaded together), the tools loaded by name; what the conversation has
 * loaded, what cannot be used here and why; per tool its description, its input (as fields and as JSON Schema),
 * whether it only reads, asks first or is turned off, and its size.
 */
export function ToolsView(props: { model: InspectorModel }) {
  const { model } = props;
  let viewport: HTMLDivElement | undefined;
  const root = (
    <div class={layout.view}>
      <div class={layout.split}>
        <ScrollArea
          class={layout.list}
          contentClass={layout.listContent}
          viewportRef={(v) => {
            viewport = v;
          }}
          role="listbox"
          aria-label="Tools"
        >
          <For each={model.toolRows} key={(r) => r.key}>
            {(r) => (r.kind === 'group' ? <GroupRow model={model} row={r} /> : <ToolRow model={model} row={r} />)}
          </For>
        </ScrollArea>
        <ScrollArea class={layout.detail} contentClass={layout.detailContent}>
          <Show when={model.toolRow}>
            {(r) =>
              r.kind === 'group' ? (
                <GroupDetail model={model} group={r.group} />
              ) : (
                <ToolDetail model={model} group={r.group} tool={r.tool} />
              )
            }
          </Show>
        </ScrollArea>
      </div>
    </div>
  ) as HTMLElement;
  listKeys(
    () => root,
    (d) => model.moveTool(d),
  );
  // Opens on the selected tool (one asked for from the trace: in the middle); a moved selection is kept in view.
  const first = model.selectedTool
    .pipe(take(1))
    .subscribe((key) => queueMicrotask(() => reveal(viewport, key, 'center')));
  const moves = model.selectedTool
    .pipe(skip(1))
    .subscribe((key) => queueMicrotask(() => reveal(viewport, key, 'nearest', HeadHeight)));
  onCleanup(() => {
    first.unsubscribe();
    moves.unsubscribe();
  });
  return root;
}

type GroupLoaded = 'all' | 'some' | 'none';

/** How much of a group the conversation has loaded. */
function loadedOf(group: ToolGroupView): GroupLoaded {
  const n = group.tools.filter((t) => t.state === 'loaded').length;
  return n === 0 ? 'none' : n === group.tools.length ? 'all' : 'some';
}

const tokensOf = (tools: readonly ToolView[]) => tools.reduce((n, t) => n + t.tokens, 0);

function GroupRow(props: { model: InspectorModel; row: Extract<ToolsRow, { kind: 'group' }> }) {
  const { model, row } = props;
  const g = row.group;
  const n = g.tools.length;
  return (
    <button
      type="button"
      role="option"
      class={layout.head}
      data-row={row.key}
      aria-selected={model.isSelectedTool(row.key)}
      onClick={() => model.selectTool(row.key)}
    >
      <span class={layout.headTitle}>{g.title}</span>
      {g.unavailable ? (
        <span class={layout.status} data-tip={g.unavailable}>
          unavailable here
        </span>
      ) : g.kind === 'skill' && loadedOf(g) === 'all' ? (
        <span class={layout.status}>loaded</span>
      ) : null}
      <span class={layout.headMeta}>
        {n} {n === 1 ? 'tool' : 'tools'} · ~{k(tokensOf(g.tools))} tokens
      </span>
    </button>
  );
}

function ToolRow(props: { model: InspectorModel; row: Extract<ToolsRow, { kind: 'tool' }> }) {
  const { model, row } = props;
  const t = row.tool;
  return (
    <button
      type="button"
      role="option"
      class={layout.row}
      data-row={row.key}
      data-off={!t.enabled || t.unavailable ? '' : undefined}
      aria-selected={model.isSelectedTool(row.key)}
      onClick={() => model.selectTool(row.key)}
    >
      <span class={layout.role}>{t.read ? 'reads' : 'tool'}</span>
      <span class={layout.body}>
        <span class={layout.name}>{t.name}</span>
        <span class={layout.focus}>
          {t.extra ? 'By name. ' : ''}
          {firstSentence(t.description)}
        </span>
      </span>
      <span class={layout.meta}>
        {t.approval ? <span class={layout.flag}>asks</span> : null}
        {!t.enabled ? <span class={layout.flag}>disabled</span> : null}
        {t.state === 'loaded' ? <span>{t.loadedAtTurn ? `turn ${t.loadedAtTurn}` : 'loaded'}</span> : null}
        <span>~{k(t.tokens)}</span>
      </span>
    </button>
  );
}

const Placement: Record<ToolGroupView['kind'], string> = {
  core: 'Core: offered from the first turn.',
  skill: 'Part of a skill: loaded with its guidance by disclose_skills.',
  single: 'Loaded by name when the model asks for it, several at once.',
};

const LoadedWords: Record<GroupLoaded, string> = { all: 'loaded', some: 'partly loaded', none: 'not loaded' };

function GroupDetail(props: { model: InspectorModel; group: ToolGroupView }) {
  const { model, group } = props;
  const guidance = group.guidance ? (group.guidanceTokens ?? Math.round(group.guidance.length / 4)) : 0;
  return (
    <>
      <DetailHead
        kind={group.kind === 'skill' ? 'skill' : group.kind === 'core' ? 'core' : 'by name'}
        where={group.kind === 'skill' ? group.id : undefined}
      />
      <p class={parts.prose}>{group.summary ?? group.title}</p>
      <Facts
        items={[
          ['Tools', String(group.tools.length)],
          [
            'Size',
            `~${k(tokensOf(group.tools))} tokens of definitions${guidance ? `, ~${k(guidance)} of guidance` : ''}`,
          ],
          ['Now', LoadedWords[loadedOf(group)]],
          ...(group.unavailable ? [['Here', group.unavailable] as const] : []),
        ]}
      />
      <p class={parts.note}>{Placement[group.kind]}</p>
      {group.guidance ? (
        <>
          <h4 class={parts.h}>Guidance</h4>
          <p class={parts.note}>Exactly as the model reads it when the skill loads.</p>
          <pre class={parts.code}>{group.guidance}</pre>
        </>
      ) : null}
      <h4 class={parts.h}>Tools</h4>
      <ul class={classes.toolList}>
        {group.tools.map((t) => (
          <li>
            <button type="button" class={`${parts.link} ${layout.name}`} onClick={() => model.openTool(t.name)}>
              {t.name}
            </button>
            <span class={parts.dim}>{firstSentence(t.description)}</span>
          </li>
        ))}
      </ul>
    </>
  );
}

function ToolDetail(props: { model: InspectorModel; group: ToolGroupView; tool: ToolView }) {
  const { model, group, tool } = props;
  const fields = fieldsOf(tool.input);
  const tab = model.tabs.pipe(
    map((t) => t.definition),
    distinctUntilChanged(),
  );
  const unavailable = tool.unavailable ?? group.unavailable;
  return (
    <>
      <DetailHead kind={`${tool.read ? 'tool · reads' : 'tool'}${tool.extra ? ' · by name' : ''}`} where={group.title}>
        <span class={layout.name}>{tool.name}</span>
      </DetailHead>
      <p class={parts.prose}>{tool.description}</p>
      <Facts
        items={[
          ['Placement', Placement[group.kind]],
          ['Reads', tool.read ? 'only: no effect, no checkpoint' : 'no: it changes something'],
          ['Approval', tool.approval ? 'asks the person first' : 'not needed'],
          ['Enabled', tool.enabled ? 'yes' : 'no: not in this turn’s tool set'],
          ...(tool.modules?.length ? [['Modules', tool.modules.join(', ')] as const] : []),
          ['Size', `~${k(tool.tokens)} tokens`],
          [
            'Now',
            tool.state === 'loaded'
              ? tool.loadedAtTurn
                ? `loaded at turn ${tool.loadedAtTurn}`
                : 'loaded from the start'
              : tool.state === 'deferred'
                ? 'declared by name; loads when called'
                : 'not loaded',
          ],
          ...(unavailable ? [['Here', unavailable] as const] : []),
        ]}
      />
      <TextTabs
        label="Definition"
        tabs={[
          { id: 'input', label: 'Input', count: fields.length },
          { id: 'schema', label: 'JSON Schema' },
        ]}
        value={tab}
        onChange={(t) => model.setTab('definition', t)}
      />
      <Show when={tab}>
        {(t) =>
          t === 'schema' ? (
            <Json value={tool.input} />
          ) : fields.length ? (
            <dl class={classes.fields}>
              {fields.map((f) => (
                <div class={classes.field}>
                  <dt>
                    <span class={layout.name}>{f.path}</span>
                    <span class={parts.dim}>
                      {' '}
                      {f.type}
                      {f.required ? '' : ', optional'}
                    </span>
                  </dt>
                  {f.description ? <dd>{f.description}</dd> : null}
                  {f.values ? <dd class={parts.dim}>One of: {f.values.join(', ')}</dd> : null}
                </div>
              ))}
            </dl>
          ) : (
            <p class={parts.note}>No arguments.</p>
          )
        }
      </Show>
    </>
  );
}

export interface SchemaField {
  path: string;
  type: string;
  required: boolean;
  description?: string;
  values?: string[];
}

type Schema = {
  type?: string | string[];
  properties?: Record<string, Schema>;
  required?: string[];
  items?: Schema;
  enum?: unknown[];
  anyOf?: Schema[];
  oneOf?: Schema[];
  description?: string;
};

/** A JSON Schema's fields, as rows: nested objects as `target.ligand`, in lists as `items[].name`. */
export function fieldsOf(input: object, prefix = ''): SchemaField[] {
  const schema = input as Schema;
  const required = new Set(schema.required ?? []);
  return Object.entries(schema.properties ?? {}).flatMap(([name, s]) => {
    const path = `${prefix}${name}`;
    const item = s.items;
    const values = s.enum ?? item?.enum;
    const field: SchemaField = {
      path,
      type: typeOf(s),
      required: required.has(name),
      ...(s.description ? { description: s.description } : {}),
      ...(values ? { values: values.map(String) } : {}),
    };
    const nested = s.properties ? fieldsOf(s, `${path}.`) : item?.properties ? fieldsOf(item, `${path}[].`) : [];
    return [field, ...nested];
  });
}

function typeOf(s: Schema): string {
  const variants = s.anyOf ?? s.oneOf;
  if (variants) return variants.map(typeOf).join(' or ');
  const t = Array.isArray(s.type) ? s.type.join(' or ') : (s.type ?? (s.enum ? 'choice' : 'any'));
  if (t === 'array') return `list of ${s.items ? typeOf(s.items) : 'any'}`;
  return t;
}

/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type {
  AttemptView,
  ModelStepView,
  PrefixBreak,
  StepView,
  ToolStepView,
  TurnStatus,
  TurnView,
  WindowEntryView,
  WindowPart,
  WindowView,
} from '@molstar/ai-trace';
import { cachedShare, formatUsage, k, ms } from '@molstar/ai-trace/format';
import { type Child, For, observableState, onCleanup, ScrollArea, Show } from '@molstar/prism-ui';
import { distinctUntilChanged, map, skip, take } from 'rxjs';
import { classes as layout } from './_generated/inspector.module.css.js';
import { classes } from './_generated/trace.module.css.js';
import type { DetailTabs, InspectorModel, TraceRow } from './model.js';
import { DetailHead, Facts, firstLine, Json, listKeys, Meter, parts, reveal, TextTabs, Usage } from './parts.js';

/** The height the sticky turn header takes over the rows (px): a revealed row goes below it. */
const HeadHeight = 40;

export const TurnStatusWords: Record<TurnStatus, string> = {
  thinking: 'thinking…',
  working: 'working…',
  waiting: 'waiting for you',
  done: 'done',
  stopped: 'stopped',
  failed: 'failed',
};

/**
 * The agent trace (analysis/viewer-mode.md §5): the whole session at a glance (the overview strip), the conversation as
 * rows (each turn's request, model calls and tool calls; nested runs indented), and the selected row in detail (for a
 * model call: exactly what it saw). j/k or the arrows move through the rows.
 */
export function TraceView(props: { model: InspectorModel }) {
  const { model } = props;
  const has = model.rows.pipe(
    map((r) => r.length > 0),
    distinctUntilChanged(),
  );
  return (
    <Show
      when={has}
      fallback={() => (
        <p class={layout.empty}>
          The agent has not run yet. Ask it something (⌘K, then Ask the assistant): each turn shows here, step by step.
        </p>
      )}
    >
      {() => <Trace model={model} />}
    </Show>
  );
}

function Trace(props: { model: InspectorModel }) {
  const { model } = props;
  let viewport: HTMLDivElement | undefined;
  const root = (
    <div class={layout.view}>
      <Overview model={model} />
      <div class={layout.split}>
        <ScrollArea
          class={layout.list}
          contentClass={layout.listContent}
          viewportRef={(v) => {
            viewport = v;
          }}
          role="listbox"
          aria-label="Trace"
        >
          <For each={model.rows} key={(r) => r.key}>
            {(r) => <Row model={model} row={r} />}
          </For>
        </ScrollArea>
        <ScrollArea class={layout.detail} contentClass={layout.detailContent}>
          <Show when={model.row}>{(r) => <Detail model={model} row={r} />}</Show>
        </ScrollArea>
      </div>
    </div>
  ) as HTMLElement;
  listKeys(
    () => root,
    (d) => model.move(d),
  );
  // Opens on the selected row (in the middle); a moved selection is kept in view.
  const first = model.selected.pipe(take(1)).subscribe((key) => queueMicrotask(() => reveal(viewport, key, 'center')));
  const moves = model.selected
    .pipe(skip(1))
    .subscribe((key) => queueMicrotask(() => reveal(viewport, key, 'nearest', HeadHeight)));
  onCleanup(() => {
    first.unsubscribe();
    moves.unsubscribe();
  });
  return root;
}

// ------------------------------------------------------------------------------------------------- overview strip

function duration(s: StepView): number {
  const first = s.attempts[0];
  const last = s.attempts.at(-1);
  if (!first || !last) return 0;
  return (last.end ?? last.start + 800) - first.start;
}

const failedStep = (s: StepView) =>
  s.kind === 'tool' ? s.status === 'failed' : s.attempts.every((a) => a.status === 'failed');
const runningStep = (s: StepView) => s.attempts.some((a) => a.status === 'running');

function stepLabel(s: StepView): string {
  if (s.kind === 'tool') return s.title ?? s.name;
  return s.reply ? 'model: reply' : `model: ${callNames(s) || '…'}`;
}

const callNames = (s: ModelStepView) => (s.toolCalls ?? []).map((c) => c.split('(')[0]).join(', ');

/** Every step the same width (Steps), or widths that follow durations (Time; idle time between turns left out). */
const unit = (scale: 'steps' | 'time', s: StepView) => (scale === 'steps' ? 1 : Math.max(0.25, duration(s) / 700));

function Overview(props: { model: InspectorModel }) {
  const { model } = props;
  return (
    <div class={classes.overview}>
      <div class={classes.lanes} aria-hidden="true">
        <span>model</span>
        <span>tools</span>
      </div>
      <div class={classes.track}>
        <For each={model.turns} key={(t) => t.id}>
          {(turn, index) => <Column model={model} turn={turn} index={index.value} />}
        </For>
      </div>
      <TextTabs
        label="Overview scale"
        tabs={[
          { id: 'steps', label: 'Steps' },
          { id: 'time', label: 'Time' },
        ]}
        value={model.scale}
        onChange={(s) => model.setScale(s)}
      />
    </div>
  );
}

function Column(props: { model: InspectorModel; turn: TurnView; index: number }) {
  const { model, turn } = props;
  const grow = (f: (scale: 'steps' | 'time') => number) => model.scale.pipe(map(f));
  const request = `req:${turn.id}`;
  return (
    <div
      class={classes.column}
      data-selected={model.isSelected(`turn:${turn.id}`).pipe(map((s) => (s ? '' : undefined)))}
      style={{ flexGrow: grow((scale) => turn.steps.reduce((a, s) => a + unit(scale, s), 0.6)) }}
    >
      <button
        type="button"
        class={classes.block}
        data-kind="request"
        aria-label={`Turn ${props.index + 1}: ${turn.prompt}`}
        data-tip={`Turn ${props.index + 1}: ${turn.prompt}`}
        aria-pressed={model.isSelected(request)}
        onClick={() => model.select(request)}
      />
      {turn.steps.map((s) => {
        const key = `step:${s.id}`;
        const label = `${stepLabel(s)} · ${runningStep(s) ? 'running' : ms(duration(s))}`;
        return (
          <button
            type="button"
            class={classes.block}
            data-kind={s.kind}
            data-nested={s.under !== undefined ? '' : undefined}
            data-failed={failedStep(s) ? '' : undefined}
            data-retried={s.attempts.length > 1 ? '' : undefined}
            data-running={runningStep(s) ? '' : undefined}
            aria-label={label}
            data-tip={label}
            aria-pressed={model.isSelected(key)}
            style={{ flexGrow: grow((scale) => unit(scale, s)) }}
            onClick={() => model.select(key)}
          />
        );
      })}
    </div>
  );
}

// ------------------------------------------------------------------------------------------------- rows

function Row(props: { model: InspectorModel; row: TraceRow }) {
  const { model, row } = props;
  const selected = model.isSelected(row.key);
  const select = () => model.select(row.key);
  if (row.kind === 'turn') {
    const t = row.turn;
    const tools = t.steps.filter((s) => s.kind === 'tool').length;
    return (
      <button
        type="button"
        role="option"
        class={layout.head}
        data-row={row.key}
        aria-selected={selected}
        onClick={select}
      >
        <span class={layout.headTitle}>Turn {row.index + 1}</span>
        <span class={layout.status} data-status={t.status}>
          {TurnStatusWords[t.status]}
        </span>
        <span class={layout.headMeta}>
          {[
            t.end !== undefined ? ms(t.end - t.start) : 'running',
            formatUsage(t.tokens),
            `${tools} tool ${tools === 1 ? 'call' : 'calls'}`,
          ].join(' · ')}
        </span>
      </button>
    );
  }
  if (row.kind === 'request')
    return (
      <Line row={row.key} role="request" selected={selected} onSelect={select}>
        <span>{row.turn.prompt}</span>
        {row.turn.focus ? <span class={layout.focus}>Focus: {row.turn.focus}</span> : null}
      </Line>
    );
  const s = row.step;
  if (s.kind === 'model') {
    const running = runningStep(s);
    return (
      <Line
        row={row.key}
        role={s.reply ? 'reply' : 'model'}
        nested={row.nested}
        selected={selected}
        onSelect={select}
        meta={
          <>
            {s.attempts.length > 1 ? (
              <span class={layout.flag} data-tip={retryNote(s.attempts)}>
                retried
              </span>
            ) : null}
            {s.prefixBreak ? (
              <span class={layout.flag} data-tip={prefixNote(s.prefixBreak)}>
                new prefix
              </span>
            ) : null}
            <span class={layout.time}>{running ? 'thinking…' : ms(duration(s))}</span>
            <span class={layout.tokens}>{s.tokens ? `${k(s.tokens.input)} → ${k(s.tokens.output)}` : ''}</span>
          </>
        }
      >
        {running ? (
          <span class={parts.dim}>waiting for the model</span>
        ) : s.reply ? (
          <span>{s.reply}</span>
        ) : (
          <span>
            <span class={parts.dim}>asks for </span>
            <span class={layout.name}>{callNames(s)}</span>
          </span>
        )}
      </Line>
    );
  }
  const changed = uniqueChanged(s);
  return (
    <Line
      row={row.key}
      role={s.status === 'failed' ? 'failed' : s.read ? 'reads' : 'tool'}
      failed={s.status === 'failed'}
      nested={row.nested}
      selected={selected}
      onSelect={select}
      meta={
        <>
          {changed.map((c) => (
            <span class={layout.chip}>{c.label}</span>
          ))}
          <span class={layout.time}>{s.status === 'running' ? 'running' : ms(duration(s))}</span>
          <span class={layout.tokens} />
        </>
      }
    >
      {s.title ? <span>{s.title}</span> : <span class={layout.name}>{s.name}</span>}
      <span class={layout.args}>{argsLine(s)}</span>
      {s.result !== undefined ? <span class={layout.result}>→ {firstLine(s.result)}</span> : null}
    </Line>
  );
}

function Line(props: {
  row: string;
  role: string;
  selected: ReturnType<InspectorModel['isSelected']>;
  onSelect(): void;
  failed?: boolean;
  nested?: boolean;
  meta?: Child;
  children: Child;
}) {
  return (
    <button
      type="button"
      role="option"
      class={layout.row}
      data-row={props.row}
      data-failed={props.failed ? '' : undefined}
      data-nested={props.nested ? '' : undefined}
      aria-selected={props.selected}
      onClick={() => props.onSelect()}
    >
      <span class={layout.role}>{props.role}</span>
      <span class={layout.body}>
        {props.nested ? <span class={layout.nest}>↳</span> : null}
        {props.children}
      </span>
      <span class={layout.meta}>{props.meta}</span>
    </button>
  );
}

/** A call's arguments in one line: `id=1tqn source=pdb` (the row fades what does not fit). */
function argsLine(s: ToolStepView): string {
  return Object.entries(s.args)
    .map(([key, v]) => `${key}=${typeof v === 'string' ? v : JSON.stringify(v)}`)
    .join(' ');
}

const uniqueChanged = (s: ToolStepView) => [...new Map((s.changed ?? []).map((c) => [c.id, c])).values()];

/** Why a request did not extend the previous one, in a sentence. */
function prefixNote(b: PrefixBreak): string {
  return `It does not extend the previous request: ${b.reason}. It kept ${b.kept} of that request’s ${b.of} parts (tools, system, then each message), so the provider’s cache could serve at most those.`;
}

function retryNote(attempts: readonly AttemptView[]): string {
  const failed = attempts.filter((a) => a.status === 'failed');
  return `${attempts.length} attempts: ${failed.map((a) => a.failure ?? 'failed').join('; ')}`;
}

// ------------------------------------------------------------------------------------------------- detail

function Detail(props: { model: InspectorModel; row: TraceRow }) {
  const { model, row } = props;
  if (row.kind === 'turn') return <TurnDetail turn={row.turn} index={row.index} />;
  if (row.kind === 'request') return <RequestDetail turn={row.turn} />;
  return row.step.kind === 'model' ? (
    <ModelDetail model={model} step={row.step} />
  ) : (
    <ToolDetail model={model} step={row.step} />
  );
}

function TurnDetail(props: { turn: TurnView; index: number }) {
  const { turn } = props;
  const models = turn.steps.filter((s): s is ModelStepView => s.kind === 'model' && s.under === undefined);
  const max = Math.max(1, ...models.map((m) => m.tokens?.input ?? 0));
  return (
    <>
      <DetailHead kind="turn" where={`${props.index + 1} · ${turn.id}`} />
      <p class={parts.prose}>{turn.prompt}</p>
      <Facts
        items={[
          ['Status', TurnStatusWords[turn.status]],
          ['Duration', turn.end !== undefined ? ms(turn.end - turn.start) : 'running'],
          ['Tokens', formatUsage(turn.tokens)],
          ['Model calls', String(turn.steps.filter((s) => s.kind === 'model').length)],
          ['Tool calls', String(turn.steps.filter((s) => s.kind === 'tool').length)],
          ...(turn.failure ? [['Failure', turn.failure] as const] : []),
        ]}
      />
      <h4 class={parts.h}>Input size per model call</h4>
      <p class={parts.note}>Each call carries the whole conversation so far: the growth is what a long turn costs.</p>
      <div class={classes.bars}>
        {models.map((m, i) => {
          const share = m.tokens ? cachedShare(m.tokens) : undefined;
          return (
            <div class={classes.bar}>
              <span>{i + 1}</span>
              <Meter share={(m.tokens?.input ?? 0) / max} />
              <span>{m.tokens ? `${k(m.tokens.input)}${share !== undefined ? ` · ${share}%` : ''}` : '…'}</span>
            </div>
          );
        })}
      </div>
      <p class={parts.note}>The share after the size is what the provider’s cache served.</p>
      <h4 class={parts.h}>Milestones</h4>
      <ul class={classes.milestones}>
        {turn.milestones.map((m) => (
          <li>
            <span class={classes.at}>+{ms(m.at - turn.start)}</span> {m.kind}
            {m.detail ? <span class={parts.dim}> · {m.detail}</span> : null}
          </li>
        ))}
      </ul>
    </>
  );
}

function RequestDetail(props: { turn: TurnView }) {
  const { turn } = props;
  return (
    <>
      <DetailHead kind="request" where={turn.id} />
      <p class={parts.prose}>{turn.prompt}</p>
      <Facts items={[['Focus', turn.focus ?? 'none']]} />
      <p class={parts.note}>
        Sent with the state (what is loaded, selected and shown) and the focus (what the person points at), so the model
        works from what is on the screen.
      </p>
    </>
  );
}

/** A detail's tabs, kept in the model per kind of row (a model call opens on What it saw, wherever one was left). */
function DetailTabsOf<K extends keyof DetailTabs>(props: {
  model: InspectorModel;
  kind: K;
  tabs: readonly { id: DetailTabs[K]; label: string; count?: string | number }[];
  children: (tab: DetailTabs[K]) => Child;
}) {
  const { model, kind } = props;
  const tab = model.tabs.pipe(
    map((t) => t[kind]),
    distinctUntilChanged(),
  );
  return (
    <>
      <TextTabs label="Detail" tabs={props.tabs} value={tab} onChange={(t) => model.setTab(kind, t)} />
      <Show when={tab}>{(t) => props.children(t)}</Show>
    </>
  );
}

function ModelDetail(props: { model: InspectorModel; step: ModelStepView }) {
  const { model, step } = props;
  const failed = step.attempts.filter((a) => a.status === 'failed');
  return (
    <>
      <DetailHead kind={step.under !== undefined ? 'model call · nested' : 'model call'} where={step.callId}>
        {step.model ? <span class={parts.dim}>{step.model}</span> : null}
      </DetailHead>
      <DetailTabsOf
        model={model}
        kind="model"
        tabs={[
          { id: 'saw', label: 'What it saw' },
          { id: 'answer', label: 'Answer' },
          { id: 'timing', label: 'Timing', count: step.attempts.length > 1 ? step.attempts.length : '' },
        ]}
      >
        {(tab) =>
          tab === 'saw' ? (
            step.window ? (
              <WindowDetail window={step.window} prefixBreak={step.prefixBreak} />
            ) : (
              <p class={parts.note}>Not captured.</p>
            )
          ) : tab === 'answer' ? (
            <>
              <Facts
                items={[
                  ['Finish', step.finish ?? (runningStep(step) ? 'running' : '—')],
                  ['Tokens', step.tokens ? <Usage tokens={step.tokens} /> : '…'],
                  ...(failed.length ? [['Failed', failed.map((a) => a.failure ?? 'failed').join('; ')] as const] : []),
                ]}
              />
              {step.reply ? <p class={parts.prose}>{step.reply}</p> : null}
              {step.toolCalls?.length ? <pre class={parts.code}>{step.toolCalls.join('\n')}</pre> : null}
            </>
          ) : (
            <Timing attempts={step.attempts} />
          )
        }
      </DetailTabsOf>
    </>
  );
}

function ToolDetail(props: { model: InspectorModel; step: ToolStepView }) {
  const { model, step } = props;
  const failed = step.status === 'failed';
  const changed = uniqueChanged(step);
  const tool = model.tools.value.groups.flatMap((g) => g.tools).find((t) => t.name === step.name);
  const description =
    tool?.description ?? model.trace.value?.preamble.tools.find((t) => t.name === step.name)?.description;
  return (
    <>
      <DetailHead kind={failed ? 'tool call · failed' : 'tool call'} where={step.callId}>
        <span class={layout.name}>{step.name}</span>
      </DetailHead>
      {step.title ? <p class={parts.prose}>{step.title}</p> : null}
      <DetailTabsOf
        model={model}
        kind="tool"
        tabs={[
          { id: 'result', label: 'Result' },
          { id: 'args', label: 'Arguments' },
          { id: 'definition', label: 'Definition' },
          { id: 'timing', label: 'Timing' },
        ]}
      >
        {(tab) =>
          tab === 'result' ? (
            <>
              {changed.length ? (
                <Facts
                  items={[
                    [
                      'Changed',
                      <span class={layout.chips}>
                        {changed.map((c) => (
                          <span class={layout.chip} data-tip={c.id}>
                            {c.label}
                          </span>
                        ))}
                      </span>,
                    ],
                  ]}
                />
              ) : null}
              <pre class={parts.code} data-failed={failed ? '' : undefined}>
                {step.result ?? '…'}
              </pre>
              {failed ? (
                <p class={parts.note}>
                  The failure went back to the model as the call’s result; the next call fixes it.
                </p>
              ) : null}
            </>
          ) : tab === 'args' ? (
            <Json value={step.args} />
          ) : tab === 'definition' ? (
            <>
              <p class={parts.prose}>{description ?? 'No definition captured.'}</p>
              {tool ? (
                <Facts
                  items={[
                    ['Reads', tool.read ? 'only: no effect, no checkpoint' : 'no: it changes something'],
                    ['Approval', tool.approval ? 'asks the person first' : 'not needed'],
                  ]}
                />
              ) : null}
              <button type="button" class={parts.link} onClick={() => model.openTool(step.name)}>
                Its input, where it loads and its size, in Tools
              </button>
            </>
          ) : (
            <Timing attempts={step.attempts} />
          )
        }
      </DetailTabsOf>
    </>
  );
}

/** The same words as the rows: what people asked, what the model said, what tools returned. */
const EntryRoles: Record<WindowEntryView['role'], string> = {
  system: 'guidance',
  tools: 'tools',
  user: 'request',
  assistant: 'model',
  tool: 'tool',
};

const Parts: readonly (readonly [WindowPart, string])[] = [
  ['guidance', 'Guidance'],
  ['tools', 'Tool definitions'],
  ['state', 'State'],
  ['conversation', 'Conversation'],
];

/** A model call's input: how its window is made up, then the entries, the new ones marked and the older ones folded. */
function WindowDetail(props: { window: WindowView; prefixBreak?: PrefixBreak }) {
  const { window: w } = props;
  const total = w.entries.reduce((a, e) => a + e.tokens, 0);
  const fresh = w.entries.filter((e) => e.fresh).reduce((a, e) => a + e.tokens, 0);
  const size = (part: WindowPart) => w.entries.filter((e) => e.part === part).reduce((a, e) => a + e.tokens, 0);
  // What every call starts with (guidance, tools) stays; the conversation before the first new entry folds.
  const lead = w.entries.findIndex((e) => e.part !== 'guidance' && e.part !== 'tools');
  const firstFresh = w.entries.findIndex((e) => e.fresh);
  const foldFrom = lead < 0 ? w.entries.length : lead;
  const foldTo = firstFresh < 0 ? w.entries.length : Math.max(foldFrom, firstFresh);
  const earlier = observableState(foldTo - foldFrom < 2);
  const folded = w.entries.slice(foldFrom, foldTo);
  return (
    <>
      <div class={classes.window}>
        <div class={classes.makeup} data-tip={`≈${k(total)} tokens`}>
          {Parts.map(([part]) => {
            const t = size(part);
            return t ? <span class={classes.part} data-part={part} style={{ flexGrow: t }} /> : null;
          })}
        </div>
        <div class={classes.legend}>
          {Parts.map(([part, name]) => {
            const t = size(part);
            return t ? (
              <span>
                <i class={classes.swatch} data-part={part} />
                {name} {k(t)}
              </span>
            ) : null;
          })}
        </div>
        <p class={parts.note}>
          ≈{k(total)} tokens
          {w.limit ? `, ${((total / w.limit) * 100).toFixed(1)}% of the model’s window` : ''}; {k(fresh)} new since the
          previous call.
        </p>
        {props.prefixBreak ? <p class={parts.note}>{prefixNote(props.prefixBreak)}</p> : null}
      </div>
      <ol class={classes.entries}>
        {w.entries.slice(0, foldFrom).map((e) => (
          <Entry entry={e} />
        ))}
        {folded.length ? (
          <Show
            when={earlier.asReadonly}
            fallback={() => (
              <li>
                <button type="button" class={classes.entry} onClick={() => earlier.set(true)}>
                  <span class={classes.fresh} />
                  <span class={classes.entryRole}>earlier</span>
                  <span class={classes.entryTitle}>
                    {folded.length} {folded.length === 1 ? 'entry' : 'entries'} from before
                  </span>
                  <span class={classes.entryText}>show</span>
                  <span class={classes.entryTokens}>{k(folded.reduce((a, x) => a + x.tokens, 0))}</span>
                </button>
              </li>
            )}
          >
            {() => folded.map((e) => <Entry entry={e} />)}
          </Show>
        ) : null}
        {w.entries.slice(foldTo).map((e) => (
          <Entry entry={e} />
        ))}
      </ol>
    </>
  );
}

/** One entry of the window: its role, title, first line and size; a click shows it whole. */
function Entry(props: { entry: WindowEntryView }) {
  const e = props.entry;
  const open = observableState(false);
  return (
    <li>
      <button
        type="button"
        class={classes.entry}
        data-fresh={e.fresh ? '' : undefined}
        data-failed={e.failed ? '' : undefined}
        aria-expanded={open.asReadonly}
        onClick={() => open.set(!open.value)}
      >
        <span class={classes.fresh} data-tip={e.fresh ? 'New since the previous call' : undefined} />
        <span class={classes.entryRole}>{EntryRoles[e.role]}</span>
        <span class={classes.entryTitle}>{e.title}</span>
        <span class={classes.entryText}>{firstLine(e.text)}</span>
        <span class={classes.entryTokens}>{k(e.tokens)}</span>
      </button>
      <Show when={open.asReadonly}>{() => <pre class={`${parts.code} ${classes.entryBody}`}>{e.text}</pre>}</Show>
    </li>
  );
}

/** Attempts as a small waterfall: waits hatched, until the first output light, the output darker. */
function Timing(props: { attempts: readonly AttemptView[] }) {
  const { attempts } = props;
  const start = Math.min(...attempts.map((a) => a.start - (a.waitBefore?.ms ?? 0)));
  const end = Math.max(...attempts.map((a) => a.end ?? a.start + 400));
  const span = Math.max(1, end - start);
  const pct = (t: number) => `${((t - start) / span) * 100}%`;
  const width = (from: number, to: number) => `${(Math.max(0, to - from) / span) * 100}%`;
  return (
    <div class={classes.timing}>
      {attempts.map((a) => {
        const until = a.firstOutput ?? a.end ?? end;
        return (
          <div class={classes.attempt}>
            <span class={classes.ordinal}>#{a.ordinal}</span>
            <span class={classes.lane}>
              {a.waitBefore ? (
                <i
                  class={classes.segment}
                  data-kind="wait"
                  data-tip={`${a.waitBefore.reason}, ${ms(a.waitBefore.ms)}`}
                  style={{ left: pct(a.start - a.waitBefore.ms), width: width(a.start - a.waitBefore.ms, a.start) }}
                />
              ) : null}
              <i
                class={classes.segment}
                data-kind={a.status === 'failed' ? 'failed' : a.status === 'running' ? 'running' : 'pre'}
                style={{ left: pct(a.start), width: width(a.start, until) }}
              />
              {a.firstOutput !== undefined && a.end !== undefined ? (
                <i
                  class={classes.segment}
                  data-kind="out"
                  style={{ left: pct(a.firstOutput), width: width(a.firstOutput, a.end) }}
                />
              ) : null}
            </span>
            <span class={classes.attemptText} data-failed={a.status === 'failed' ? '' : undefined}>
              {a.status === 'failed'
                ? (a.failure ?? 'failed')
                : a.status === 'running'
                  ? 'running'
                  : `${ms((a.end ?? a.start) - a.start)}${a.firstOutput !== undefined ? `, the first output after ${ms(a.firstOutput - a.start)}` : ''}`}
              {a.waitBefore ? (
                <span class={parts.dim}>
                  {' '}
                  · after {a.waitBefore.reason}, {ms(a.waitBefore.ms)}
                </span>
              ) : null}
            </span>
          </div>
        );
      })}
      <div class={classes.key}>
        <span>
          <i data-kind="wait" />
          waiting
        </span>
        <span>
          <i data-kind="pre" />
          until the first output
        </span>
        <span>
          <i data-kind="out" />
          output
        </span>
      </div>
    </div>
  );
}

/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { TraceData } from '@molstar/ai/trace';
import { domView, onCleanup, type ShellPage, Show } from '@molstar/prism-ui';
import { distinctUntilChanged, map } from 'rxjs';
import { InspectorPageId } from '../inspector-access.js';
import { classes as layout } from './_generated/inspector.module.css.js';
import { LimitsView } from './limits.js';
import { type InspectorModel, type InspectorView, toolCounts } from './model.js';
import { BudgetMeter, TextTabs } from './parts.js';
import { ToolsView } from './tools.js';
import { TraceView } from './trace.js';

/**
 * The Inspector as a page over the whole app (analysis/viewer-mode.md §5): it needs the room. Its header has the three
 * views as words (Trace with its turns, Tools with its tools, Limits), then whether the agent works and its budget; a
 * page focus (`trace`, `tools`, `limits`, `turn:<id>`, `step:<id>`, `tool:<name>`) opens it there.
 */
export function inspectorPage(model: InspectorModel): ShellPage {
  return {
    id: InspectorPageId,
    title: 'Inspect assistant',
    icon: 'sparkles',
    header: domView(() => <InspectorHeader model={model} />),
    render: domView(({ ctx }) => {
      if (ctx.focus) model.focus(ctx.focus);
      // Opened again at another focus while open (the Assistant tab's rows, a command): it stays, and moves there.
      const opened = ctx.shell.page.subscribe((open) => {
        if (open?.page.id === InspectorPageId && open.focus) model.focus(open.focus);
      });
      onCleanup(() => opened.unsubscribe());
      return <InspectorBody model={model} />;
    }),
  };
}

type RunState = 'working' | 'waiting' | 'idle';

function runState(trace: TraceData | undefined): RunState {
  const status = trace?.turns.at(-1)?.status;
  return status === 'waiting' ? 'waiting' : status === 'thinking' || status === 'working' ? 'working' : 'idle';
}

export function InspectorHeader(props: { model: InspectorModel }) {
  const { model } = props;
  const state = model.trace.pipe(map(runState), distinctUntilChanged());
  const turns = model.trace.pipe(map((t) => (t?.turns.length ? t.turns.length : '')));
  const tools = model.tools.pipe(map((t) => toolCounts(t).tools));
  return (
    <>
      <TextTabs<InspectorView>
        label="Inspector"
        tabs={[
          { id: 'trace', label: 'Trace', count: turns },
          { id: 'tools', label: 'Tools', count: tools },
          { id: 'limits', label: 'Limits' },
        ]}
        value={model.view}
        onChange={(v) => model.setView(v)}
      />
      <span class={layout.spacer} />
      <Show when={model.trace.pipe(map((t) => !!t))}>
        {() => (
          <span class={layout.state} data-state={state}>
            <i class={layout.dot} />
            {state}
          </span>
        )}
      </Show>
      <BudgetMeter budget={model.trace.pipe(map((t) => t?.budget))} />
    </>
  );
}

export function InspectorBody(props: { model: InspectorModel }) {
  const { model } = props;
  return (
    <div class={layout.inspector}>
      <Show when={model.view}>
        {(view) =>
          view === 'trace' ? (
            <TraceView model={model} />
          ) : view === 'tools' ? (
            <ToolsView model={model} />
          ) : (
            <LimitsView model={model} />
          )
        }
      </Show>
    </div>
  );
}

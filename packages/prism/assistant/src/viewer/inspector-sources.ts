/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { compactionFor } from '@molstar/ai-session/compaction';
import {
  createTraceCache,
  type LimitRowView,
  type LimitsData,
  readTrace,
  type SkillInput,
  type ToolEntry,
  type ToolsData,
  toolsView,
  type TraceData,
} from '@molstar/ai-trace';
import { observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { Observable } from 'rxjs';
import type { ViewerAgent } from '@molstar/prism-viewer/agent/session';
import type { InspectorSources } from '../inspector/model.js';

// The Inspector over the viewer's agent (agent.md §1): the trace read from the session's harness through
// `@molstar/ai-trace` as its records change, while something watches it (the Inspector, the Assistant tab's links),
// the Tools view from the viewer's toolkit as it is now, the limits from the session's budget and the model's window.
// DOM-free.

export interface ViewerInspectorSources extends InspectorSources {
  /**
   * Reads the trace again now; resolves once the trace (and the tools and limits that follow it) is up to date. The
   * sources read it themselves only while one of them is subscribed.
   */
  refresh(): Promise<void>;
  dispose(): void;
}

export function createInspectorSources(agent: Pick<ViewerAgent, 'session' | 'toolkit'>): ViewerInspectorSources {
  const { session, toolkit } = agent;
  const trace = observableState<TraceData | undefined>(undefined);
  const tools = observableState<ToolsData>(toolsOf(agent, undefined));
  const limits = observableState<LimitsData>(limitsOf(agent, undefined));
  // Heads, entries and blobs are immutable, the event log append-only: a re-read fetches only what is new.
  const cache = createTraceCache();
  let disposed = false;

  const update = (data: TraceData | undefined) => {
    // No turn yet: nothing to show (the Trace's empty state).
    const next = data?.turns.length ? data : undefined;
    const before = trace.value;
    trace.set(next);
    if (before?.loaded !== next?.loaded || before?.declared !== next?.declared) tools.set(toolsOf(agent, next));
    limits.set(limitsOf(agent, next));
  };

  // One read at a time; changes while one runs make one more read after it (a running turn's records come in bursts).
  let reading: Promise<void> | undefined;
  let again = false;
  const read = (): Promise<void> => {
    if (disposed) return Promise.resolve();
    if (reading) {
      again = true;
      return reading;
    }
    reading = (async () => {
      do {
        again = false;
        try {
          const client = await session.client();
          const data = await readTrace(client, {
            cache,
            window: session.contextWindow(),
            toolCall: (call) => ({ read: toolkit.registration(call.toolId)?.read === true }),
          });
          if (!disposed) update(data);
        } catch (e) {
          // A read that fails (the storage closed under it) leaves the last trace; the next change reads again.
          console.warn('[prism-assistant] the trace could not be read:', e);
        }
      } while (again && !disposed);
    })().finally(() => {
      reading = undefined;
    });
    return reading;
  };

  // The harness's seq (the current one first, then each change), watched while the trace, tools or limits are: a read
  // re-groups the whole session's log, so nothing reads it on a turn's every step while nobody looks.
  let unwatch: (() => void) | undefined;
  let watchers = 0;
  const watch = () => {
    if (watchers++ > 0) return;
    void session.client().then((client) => {
      if (disposed || watchers === 0 || unwatch) return;
      unwatch = client.watch().subscribe(() => void read());
    });
  };
  const release = () => {
    if (--watchers > 0) return;
    unwatch?.();
    unwatch = undefined;
  };
  const watched = <T>(state: ReadonlyState<T>): ReadonlyState<T> => {
    const o = new Observable<T>((subscriber) => {
      watch();
      const sub = state.subscribe(subscriber);
      return () => {
        sub.unsubscribe();
        release();
      };
    }) as ReadonlyState<T>;
    Object.defineProperty(o, 'value', { get: () => state.value });
    return o;
  };
  const offTools = toolkit.onChange(() => tools.set(toolsOf(agent, trace.value)));

  return {
    trace: watched(trace.asReadonly),
    tools: watched(tools.asReadonly),
    limits: watched(limits.asReadonly),
    // A read already running goes round once more: the change asked for is in it.
    refresh: read,
    dispose() {
      disposed = true;
      unwatch?.();
      unwatch = undefined;
      offTools();
    },
  };
}

/** The Tools view: the core, each skill, the tools loaded by name; what the conversation has loaded of them. */
function toolsOf(agent: Pick<ViewerAgent, 'toolkit'>, trace: TraceData | undefined): ToolsData {
  const { toolkit } = agent;
  const enabled = new Set(toolkit.toolSet().enabled);
  const entries: ToolEntry[] = toolkit.tools.map((t) => {
    const id = t.definition.tool.id;
    return {
      definition: t.definition,
      placement: toolkit.placementOf(id),
      modules: t.modules.map(String),
      read: t.read === true,
      approval: t.approval === true,
      enabled: enabled.has(id),
    };
  });
  const skills: SkillInput[] = [
    // A lazy skill (its extension not loaded) lists no tools until it is disclosed.
    ...toolkit.groups.map((g) => ({ id: g.id, title: g.title, summary: g.summary, guidance: toolkit.groupText(g) })),
    ...toolkit.entries.map((e) => ({ id: e.id, title: e.id, summary: e.summary, guidance: e.text })),
  ];
  return toolsView({
    tools: entries,
    skills,
    ...(trace ? { loaded: trace.loaded, declared: trace.declared } : {}),
  });
}

/**
 * The limits a turn runs under (the session's, for its default model): the turn budget (each kind's start and how far
 * the session extends it), the context window and when compaction starts. Read-only: the session takes them at its
 * creation and has no setter.
 */
function limitsOf(agent: Pick<ViewerAgent, 'session'>, trace: TraceData | undefined): LimitsData {
  const { session } = agent;
  const budget = session.budget();
  const latest = trace?.turns.at(-1);
  const tokens = latest?.tokens;
  const window = session.contextWindow();
  const trigger = compactionFor(window).triggerInputTokens;
  const row = (
    key: 'toolCalls' | 'uncachedInput' | 'cacheRead' | 'output',
    label: string,
    note: string,
    unit: string,
    used: number | undefined,
  ): LimitRowView[] => {
    const value = budget.initial[key];
    if (value === undefined) return [];
    const ceiling = budget.ceiling[key];
    return [
      {
        key,
        label,
        note,
        unit,
        value,
        default: value,
        ...(used !== undefined ? { used } : {}),
        ...(ceiling !== undefined ? { ceiling } : {}),
        ...(unit === 'tokens' ? { step: 1000 } : {}),
        editable: false,
      },
    ];
  };
  const toolCalls = latest?.steps.filter((s) => s.kind === 'tool' && s.under === undefined).length;
  return {
    groups: [
      {
        id: 'budget',
        title: 'Turn budget',
        note: 'What one turn may use: it starts with these and is extended toward the ceiling while it needs more. The bars are the latest turn. Set by the viewer.',
        rows: [
          ...row('toolCalls', 'Tool calls', 'commands the agent runs in one turn', 'calls', toolCalls),
          ...row(
            'uncachedInput',
            'Uncached input',
            'tokens the provider reads anew',
            'tokens',
            tokens && Math.max(0, tokens.input - (tokens.cached ?? 0) - (tokens.written ?? 0)),
          ),
          ...row(
            'cacheRead',
            'Cache reads',
            'input served from the provider’s cache',
            'tokens',
            tokens && (tokens.cached ?? 0),
          ),
          ...row('output', 'Output', 'what the model writes, reasoning included', 'tokens', tokens?.output),
        ],
      },
      {
        id: 'model',
        title: 'Model',
        note: 'For the default model, as its provider gives it.',
        rows: [
          {
            key: 'contextWindow',
            label: 'Context window',
            note: 'the most the model reads at once',
            unit: 'tokens',
            value: window,
            default: window,
            step: 1000,
            editable: false,
          },
          {
            key: 'compaction',
            label: 'Compaction trigger',
            note: 'older turns are summarized once a request reaches this',
            unit: 'tokens',
            value: trigger,
            default: trigger,
            step: 1000,
            editable: false,
          },
        ],
      },
    ],
  };
}

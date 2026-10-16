/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ModelStep, scriptedModel } from '@molstar/ai-conformance/testing';
import { story } from '@molstar/prism-ui/testing/index';
import type { AssistantUI, PrismViewerSpec } from '@molstar/prism-viewer/spec';
import { ExampleFiles, RealViewerStory } from '@molstar/prism-viewer/testing/real-viewer';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { InspectorPageId } from '../inspector-access.js';
import { assistant } from '../spec.js';
import type { KeyScope, KeyStore } from '../keys.js';

type Clock = NonNullable<NonNullable<Parameters<typeof scriptedModel>[1]>['clock']>;
import type { SettingsPrefs } from '../viewer/settings-model.js';

// The assistant in the real viewer: its tab, the turn card under the status line, the Inspector (Trace, Tools, Limits),
// all over the viewer's agent. Its model is a scripted one (`scriptedModel`, the replies and tool calls each story
// gives), a provider of the story's own (`createAssistantWith({ provided })`), so nothing reaches a model and no key is
// asked for or read; the keys and the choice are kept in memory, never the person's.

const meta: Meta = { title: 'Viewer/Assistant', parameters: { layout: 'fullscreen' } };
export default meta;

/** Keys kept for the story only. */
function memoryKeys(): KeyStore {
  const keys = new Map<string, { key: string; scope: KeyScope }>();
  return {
    get: async (p) => keys.get(p)?.key,
    set: async (p, key, scope) => void keys.set(p, { key, scope }),
    remove: async (p) => void keys.delete(p),
    list: async () => [...keys].map(([provider, { scope }]) => ({ provider, scope })),
    protection: async () => 'scramble',
  };
}

function memoryPrefs(): SettingsPrefs {
  const values = new Map<string, string>();
  return { get: (k) => values.get(k), set: (k, v) => void values.set(k, v) };
}

/** Real time, for a step that waits (a turn still running). */
const realClock: Clock = {
  now: () => Date.now(),
  sleep: (ms: number, signal?: AbortSignal) =>
    new Promise((resolve, reject) => {
      const timer = setTimeout(resolve, ms);
      signal?.addEventListener('abort', () => {
        clearTimeout(timer);
        reject(signal.reason);
      });
    }),
};

const call = (tool: string, args: Record<string, unknown> = {}) => ({ tool, arguments: args as never });

interface AssistantStoryProps {
  /** The scripted model's steps; none: no model is provided (the setup asks for a key). */
  readonly script?: readonly ModelStep[];
  /** Asked once the assistant is in (a turn of its own). */
  readonly ask?: string;
  /** Then: the Assistant tab, or the Inspector at a focus (`trace`, `tools`, `tool:<id>`, `limits`). */
  readonly open?: 'tab' | { inspector: string };
}

/** The real viewer with 1CBS, its assistant over a scripted model, a request asked and a view opened. */
function AssistantStory(props: AssistantStoryProps) {
  let ui: AssistantUI | undefined;
  const spec = (base: PrismViewerSpec): PrismViewerSpec => ({
    ...base,
    agent: {
      ...(base.agent || {}),
      assistant: async () => {
        const model = props.script ? scriptedModel(props.script, { clock: realClock }) : undefined;
        const m = await assistant({
          keys: memoryKeys(),
          prefs: memoryPrefs(),
          ...(model
            ? {
                provided: [{ id: 'scripted', label: 'Scripted (story)', models: ['scripted'], adapter: model.adapter }],
                provider: 'scripted',
              }
            : {}),
        })();
        // The story keeps the assistant's parts, to drive them.
        return { ...m, createAssistant: async (host) => (ui = await m.createAssistant(host)) };
      },
    },
  });
  const setup = async (v: PrismViewer) => {
    // Opening the tab loads the assistant (its own chunk).
    await v.api.view.open({ tab: 'assistant.tab.assistant' });
    while (!ui) await new Promise((r) => setTimeout(r, 50));
    if (props.ask) ui.commandMenu.ask(props.ask);
    const open = props.open;
    // No view asked for: the panel closes (the turn's card is under the status line).
    if (!open) v.managers.navigation.setDock('right', { collapsed: true });
    else if (open !== 'tab') {
      // The turn's first model call in the trace first.
      await new Promise((r) => setTimeout(r, 300));
      v.views.shell?.openPage(InspectorPageId, open.inspector);
    }
  };
  return <RealViewerStory spec={spec} load={[ExampleFiles.cbs]} setup={setup} />;
}

/** A turn that styles the scene: a surface, the ligand focused; done, its card under the status line. */
const Styled: readonly ModelStep[] = [
  { calls: [call('apply_style', { style: 'surface' }), call('focus', { target: 'ligand' })] },
  { text: 'The protein is shown as a molecular surface, and the camera is on the retinoic acid in its pocket.' },
];

/** The agent module with no model yet: the Assistant tab asks for a key; asking opens it. */
export const NotSetUp: StoryObj = {
  render: () => story(() => <AssistantStory open="tab" />),
};

/** A model set up (the story's scripted one): the Assistant tab with the model, the keys and the Inspect rows. */
export const SetUp: StoryObj = {
  render: () => story(() => <AssistantStory script={Styled} open="tab" />),
};

/** A turn done: what it did, its runs, Undo turn; the card under the status line. */
export const TurnDone: StoryObj = {
  render: () => story(() => <AssistantStory script={Styled} ask="show it as a surface and focus the ligand" />),
};

/** A turn waiting for the person's yes: a load from an address the viewer does not know (Allow, Deny). */
export const TurnWaiting: StoryObj = {
  render: () =>
    story(() => (
      <AssistantStory
        script={[{ calls: [call('load_url', { url: 'https://example.org/structures/model.cif' })] }]}
        ask="load my model from example.org"
      />
    )),
};

/** A turn still running: the model takes its time (the spark pulses; Stop ends it). */
export const TurnRunning: StoryObj = {
  render: () =>
    story(() => (
      <AssistantStory script={[{ text: 'Looking at it…', waitMs: 3_600_000 }]} ask="what is this protein?" />
    )),
};

/** A turn that failed: the model's error, said on the card. */
export const TurnFailed: StoryObj = {
  render: () =>
    story(() => (
      <AssistantStory script={[{ fail: { message: 'The model is overloaded; try again.' } }]} ask="color it red" />
    )),
};

/** The Inspector's trace after a turn: its model calls and tool calls, what each saw and answered. */
export const Trace: StoryObj = {
  render: () =>
    story(() => (
      <AssistantStory script={Styled} ask="show it as a surface and focus the ligand" open={{ inspector: 'trace' }} />
    )),
};

/** The trace of a turn still running: its model call pulses, the turn has no end yet. */
export const TraceRunning: StoryObj = {
  render: () =>
    story(() => (
      <AssistantStory
        script={[{ text: 'Looking at it…', waitMs: 3_600_000 }]}
        ask="what is this protein?"
        open={{ inspector: 'trace' }}
      />
    )),
};

/** The tools the agent has: the core ones loaded from the start, the skills' disclosed on use. */
export const Tools: StoryObj = {
  render: () => story(() => <AssistantStory script={Styled} open={{ inspector: 'tools' }} />),
};

/** A tool in detail: its facts, its input as fields, its JSON Schema. */
export const ToolsAtTool: StoryObj = {
  render: () => story(() => <AssistantStory script={Styled} open={{ inspector: 'tool:apply_style' }} />),
};

/** The limits: the turn's budget, the context window, when it compacts. */
export const Limits: StoryObj = {
  render: () => story(() => <AssistantStory script={Styled} open={{ inspector: 'limits' }} />),
};

/** The Assistant tab after a turn: the Inspect rows say what there is (1 turn, the tools, the limits). */
export const Inspect: StoryObj = {
  render: () =>
    story(() => <AssistantStory script={Styled} ask="show it as a surface and focus the ligand" open="tab" />),
};

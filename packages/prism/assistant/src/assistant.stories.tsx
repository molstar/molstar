/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { story } from '@molstar/prism-ui/testing/index';
import { ViewerShellStory } from '@molstar/prism-viewer/testing/index';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { InspectorPageId } from './inspector-access.js';
import { mockAgentAddOn, mockTrace } from './testing/index.js';

const meta: Meta = { title: 'Shell/Assistant', parameters: { layout: 'fullscreen' } };
export default meta;

/**
 * The agent module included and set up (a mock agent stands in for any-harness): ask from the command menu (⌘K, type
 * a request, Ask AI assistant; “color it red” opens the filled sentence, “show the active site”, “load 1w0e and
 * superpose it”, “try something broken” run turns), watch the turn card, Allow or Deny a download, Stop, Undo turn;
 * History at the rail's foot; the Assistant tab with the model and keys.
 */
export const SetUp: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh">
        <ViewerShellStory loaded={true} addOn={mockAgentAddOn({ configured: true })} />
      </div>
    )),
};

/** The agent module included but no model yet: the Assistant tab asks for a key; the Ask row opens it. */
export const NotSetUp: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh">
        <ViewerShellStory loaded={true} addOn={mockAgentAddOn()} />
      </div>
    )),
};

/** The app with the assistant set up and the Inspector open at a focus (a page over the whole app; Esc or ✕ closes it). */
function Inspector(props: { focus: string; running?: boolean }) {
  return story(() => (
    <div style="height: 100vh">
      <ViewerShellStory
        loaded={true}
        addOn={mockAgentAddOn({ configured: true, running: props.running })}
        openPage={InspectorPageId}
        openPageAt={props.focus}
      />
    </div>
  ));
}

/**
 * The agent trace (⌥T, Assistant trace in the command menu): the overview strip (a column per turn; model calls above,
 * tool calls below; the retry hatched, the failed call in the danger color), the rows (j/k or the arrows move), and the
 * selected row in detail. Turn 2 has a rate-limited retry, a skill loading (a new prefix), a failed tool call the model
 * fixes, and a tool that runs a helper model (indented under it).
 */
export const Trace: StoryObj = { render: () => Inspector({ focus: 'trace' }) };

/** The last turn still running: its model call pulses, the header says working, the turn has no end yet. */
export const TraceRunning: StoryObj = { render: () => Inspector({ focus: 'trace', running: true }) };

/**
 * A model call's What it saw: the window's makeup (guidance, tool definitions, state, conversation), its share of the
 * model's window, what is new since the previous call, why the cache prefix broke (the selections skill loaded), and the
 * entries (earlier ones folded; a click shows one whole).
 */
export const TraceWhatItSaw: StoryObj = {
  render: () => {
    const step = mockTrace()
      .turns.flatMap((t) => t.steps)
      .find((s) => s.kind === 'model' && s.prefixBreak);
    return Inspector({ focus: `step:${step?.id ?? ''}` });
  },
};

/**
 * The agent's tools (Assistant tools): the core, each skill (loaded or not; one unavailable here, and why), the tools
 * loaded by name; reads, asks and disabled as quiet words; a group's guidance as the model reads it.
 */
export const Tools: StoryObj = { render: () => Inspector({ focus: 'tools' }) };

/** A tool in detail (as a tool call's Definition opens it): its facts, its input as fields, its JSON Schema. */
export const ToolsAtTool: StoryObj = { render: () => Inspector({ focus: 'tool:bounds.get' }) };

/**
 * What bounds the agent (Assistant limits): the turn budget per kind with what the latest turn used, the model's context
 * window, compaction trigger and thinking level, and credits as information. Type a value and Enter (Esc puts it
 * back); Reset goes back to the default.
 */
export const Limits: StoryObj = { render: () => Inspector({ focus: 'limits' }) };

/** The Assistant tab with the Inspect block under the model and keys: Trace, Tools, Limits, each opening the page. */
export const Inspect: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh">
        <ViewerShellStory loaded={true} addOn={mockAgentAddOn({ configured: true })} openTab="assistant" />
      </div>
    )),
};

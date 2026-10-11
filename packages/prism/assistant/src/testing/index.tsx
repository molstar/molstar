/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView } from '@molstar/prism-ui';
import type { MockShellAddOn } from '@molstar/prism-ui/testing/index';
import { map } from 'rxjs';
import { InspectorPageId, inspectorCommands, inspectorPage, toolCounts } from '../inspector.js';
import { AgentTurnCard } from '../views/agent-turn.js';
import { AssistantTab } from '../views/assistant-tab.js';
import { HistoryTab } from '../views/history-tab.js';
import { MockAgentModel, mockAssistant } from './mock-agent.js';
import { mockInspector } from './mock-inspector.js';

export { MockAgentModel, MockAISettings, mockAssistant } from './mock-agent.js';
export { mockInspector, mockLimits, mockTrace, MockTools } from './mock-inspector.js';

/**
 * The agent module in the mock shell (`ShellStory`'s `addOn`), contributed as an app would: History in place of Log,
 * the Assistant tab at the rail's foot, the command menu's Ask row, the turn card in the canvas's agent slot, and the
 * person's commands recorded in History. A mock agent with scripted turns stands in for any-harness. The Inspector
 * (`@molstar/prism-assistant/inspector`) is a page over the app, on a fixed mock session (its last turn `running`, or
 * none before the agent is set up), opened by its commands (Agent trace on ⌥T, Agent tools, Agent limits), the
 * Assistant tab's Inspect block and the turn card's Trace.
 */
export function mockAgentAddOn(options: { configured?: boolean; running?: boolean } = {}): MockShellAddOn {
  return ({ shell, status, loaded }) => {
    const agent = new MockAgentModel((text, kind) => status.post(text, kind), {
      configured: options.configured,
      history: loaded,
    });
    const inspector = mockInspector({ running: options.running, empty: !options.configured });
    const open = (focus: string) => shell.openPage(InspectorPageId, focus);
    const counts = inspector.tools.pipe(map(toolCounts));
    return {
      tabs: [
        {
          // History: the commands, the assistant's turns, and the messages under Messages.
          id: 'log',
          title: 'History',
          icon: 'scroll-text',
          order: 2.5,
          place: 'foot',
          render: domView(() => <HistoryTab agent={agent} status={status} />),
        },
        {
          id: 'assistant',
          // The rail's tooltip says what it is for; the panel is headed Assistant.
          title: 'Set up assistant',
          heading: 'Assistant',
          icon: 'sparkles',
          order: 2.7,
          place: 'foot',
          render: domView(() => (
            <AssistantTab
              agent={agent}
              inspect={{
                links: [
                  {
                    id: 'trace',
                    label: 'Trace',
                    detail: inspector.trace.pipe(
                      map((t) => (t?.turns.length ? `${t.turns.length} turns` : 'nothing yet')),
                    ),
                  },
                  {
                    id: 'tools',
                    label: 'Tools',
                    detail: counts.pipe(map((c) => `${c.tools} tools, ${c.skills} skills`)),
                  },
                  { id: 'limits', label: 'Limits', detail: 'turn budget, context window, thinking' },
                ],
                open,
              }}
            />
          )),
        },
      ],
      pages: [inspectorPage(inspector)],
      commands: inspectorCommands(open),
      assistant: mockAssistant(agent, () => shell.show('assistant')),
      agent: <AgentTurnCard agent={agent} onTrace={(id) => open(`turn:${id}`)} />,
      commandDone: (message, commandId) => agent.record(message, commandId),
    };
  };
}

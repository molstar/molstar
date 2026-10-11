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
import { AgentTurnCard } from '../views/agent-turn.js';
import { AssistantTab } from '../views/assistant-tab.js';
import { HistoryTab } from '../views/history-tab.js';
import { MockAgentModel, mockAssistant } from './mock-agent.js';

export { MockAgentModel, MockAISettings, mockAssistant } from './mock-agent.js';

/**
 * The agent module in the mock shell (`ShellStory`'s `addOn`), contributed as an app would: History in place of Log,
 * the Assistant tab at the rail's foot, the command menu's Ask row, the turn card in the canvas's agent slot, and the
 * person's commands recorded in History. A mock agent with scripted turns stands in for any-harness.
 */
export function mockAgentAddOn(options: { configured?: boolean } = {}): MockShellAddOn {
  return ({ shell, status, loaded }) => {
    const agent = new MockAgentModel((text, kind) => status.post(text, kind), {
      configured: options.configured,
      history: loaded,
    });
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
          render: domView(() => <AssistantTab agent={agent} />),
        },
      ],
      assistant: mockAssistant(agent, () => shell.show('assistant')),
      agent: <AgentTurnCard agent={agent} />,
      commandDone: (message, commandId) => agent.record(message, commandId),
    };
  };
}

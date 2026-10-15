/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView } from '@molstar/prism-ui';
import type { AssistantFactory, AssistantUI } from '@molstar/prism-viewer/spec';
import { createKeyStore } from '../keys.js';
import { AgentTurnCard } from '../views/agent-turn.js';
import { AssistantTab } from '../views/assistant-tab.js';
import { createViewerAgentModel } from './agent-model.js';
import { createAISettings, localPrefs } from './settings-model.js';

// The assistant in the viewer's page (viewer-spec.md §2.4, plan 4.4): its tab (the model setup), the command menu's Ask
// row, the turn card under the status line; all over the viewer's agent, made with the models the person set up.

/** Until the person chooses, the agent's default model: asking without a key opens the setup. */
const Fallback = { provider: 'anthropic', model: 'claude-sonnet-5-5' } as const;

export const createAssistant: AssistantFactory = async (host): Promise<AssistantUI> => {
  const { notifications } = host.viewer.managers;
  const settings = createAISettings(createKeyStore(), localPrefs(), {
    report: (text) => notifications.post(text, 'warning'),
  });
  await settings.load();
  const agent = host.createAgent({
    adapterFor: (provider) => settings.adapterFor(provider),
    default: settings.setup.value ?? Fallback,
  });
  const openSetup = () => void host.viewer.api.as({ kind: 'person' }).view.open({ tab: 'assistant.tab.assistant' });
  const model = createViewerAgentModel(host.viewer, agent, settings, { openSetup });
  return {
    ui: {
      tabs: [
        {
          id: 'assistant.tab.assistant',
          title: 'Assistant',
          icon: 'bot',
          order: 2.6,
          render: domView(() => <AssistantTab agent={model} />),
        },
      ],
    },
    commandMenu: {
      ready: () => !!settings.setup.value,
      // Every request is a turn for now; reading words into one command comes with the command menu's assist.
      propose: () => undefined,
      ask: (text) => model.ask(text),
      setup: openSetup,
    },
    agentView: () => <AgentTurnCard agent={model} />,
    dispose() {
      // The model first: it stops a turn running or parked, then the agent goes.
      model.dispose();
      agent.dispose();
    },
  };
};

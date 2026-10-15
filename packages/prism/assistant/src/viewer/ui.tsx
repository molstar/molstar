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
import { type InspectorFocus, InspectorPageId, inspectorCommands } from '../inspector-access.js';
import { inspectLinks } from '../inspector/links.js';
import { InspectorModel } from '../inspector/model.js';
import { inspectorPage } from '../inspector/page.js';
import { createKeyStore } from '../keys.js';
import { AgentTurnCard } from '../views/agent-turn.js';
import { AssistantTab } from '../views/assistant-tab.js';
import { createViewerAgentModel } from './agent-model.js';
import { createInspectorSources } from './inspector-sources.js';
import { createAISettings, localPrefs } from './settings-model.js';
import { keysStorageSource } from './storage.js';

// The assistant in the viewer's page (viewer-spec.md §2.4, plan 4.4): its tab (the model setup, the Inspect rows), the
// command menu's Ask row, the turn card under the status line, the Inspector (a page over the app, with its commands);
// all over the viewer's agent, made with the models the person set up.

/** Until the person chooses, the agent's default model: asking without a key opens the setup. */
const Fallback = { provider: 'anthropic', model: 'claude-sonnet-5-5' } as const;

export const createAssistant: AssistantFactory = async (host): Promise<AssistantUI> => {
  const { notifications } = host.viewer.managers;
  const keys = createKeyStore();
  const settings = createAISettings(keys, localPrefs(), {
    report: (text) => notifications.post(text, 'warning'),
  });
  await settings.load();
  const agent = host.createAgent({
    adapterFor: (provider) => settings.adapterFor(provider),
    default: settings.setup.value ?? Fallback,
  });
  const openSetup = () => void host.viewer.api.as({ kind: 'person' }).view.open({ tab: 'assistant.tab.assistant' });
  const model = createViewerAgentModel(host.viewer, agent, settings, { openSetup });
  // The Inspector (Trace, Tools, Limits) over the agent's session and toolkit, a page over the app (agent.md §1).
  const sources = createInspectorSources(agent);
  const inspector = new InspectorModel(sources);
  const { navigation } = host.viewer.managers;
  const openInspector = (focus: string) => {
    navigation.setCard(undefined);
    navigation.setPage({ id: InspectorPageId, focus });
  };
  const inspect = { links: inspectLinks(inspector), open: openInspector };
  // The stored keys in the viewer's Storage view (journal.md §7), by provider and scope, never the key; last, so a
  // failure above leaves nothing registered.
  const unregisterKeys = host.viewer.storage.register(keysStorageSource(keys, settings));
  return {
    ui: {
      tabs: [
        {
          id: 'assistant.tab.assistant',
          title: 'Assistant',
          icon: 'bot',
          order: 2.6,
          render: domView(() => <AssistantTab agent={model} inspect={inspect} />),
        },
      ],
      pages: [inspectorPage(inspector)],
    },
    commands: inspectorCommands((focus: InspectorFocus) => openInspector(focus)),
    commandMenu: {
      ready: () => !!settings.setup.value,
      // Every request is a turn for now; reading words into one command comes with the command menu's assist.
      propose: () => undefined,
      ask: (text) => model.ask(text),
      setup: openSetup,
    },
    agentView: () => <AgentTurnCard agent={model} onTrace={(id) => openInspector(`turn:${id}`)} />,
    dispose() {
      // The model first: it stops a turn running or parked, then the agent goes.
      model.dispose();
      sources.dispose();
      agent.dispose();
      unregisterKeys();
    },
  };
};

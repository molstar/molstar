/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView } from '@molstar/prism-ui';
import type { AssistantFactory, AssistantHost, AssistantUI } from '@molstar/prism-viewer/spec';
import { type InspectorFocus, InspectorPageId, inspectorCommands } from '../inspector-access.js';
import { inspectLinks } from '../inspector/links.js';
import { InspectorModel } from '../inspector/model.js';
import { inspectorPage } from '../inspector/page.js';
import { createKeyStore, type KeyStore } from '../keys.js';
import { AgentTurnCard } from '../views/agent-turn.js';
import { AssistantTab } from '../views/assistant-tab.js';
import { createViewerAgentModel } from './agent-model.js';
import { createInspectorSources } from './inspector-sources.js';
import { listedModels } from '@molstar/ai/models';
import type { AISetup } from '../model.js';
import { createAISettings, localPrefs, type ProvidedModels, type SettingsPrefs } from './settings-model.js';
import { domainsStorageSource, keptDomains, keysStorageSource } from './storage.js';
import { createDomainRules } from '@molstar/prism-viewer/agent/policies';

// The assistant in the viewer's page (viewer-spec.md §2.4): its tab (the model setup, the Inspect rows), the
// command menu's Ask row, the turn card under the status line, the Inspector (a page over the app, with its commands);
// all over the viewer's agent, made with the models the person set up.

/** Until the person chooses, the agent's default model: asking without a key opens the setup. */
const Fallback = { provider: 'anthropic', model: 'claude-sonnet-5-5' } as const;

/** What a host may give the assistant besides the viewer: models of its own, and where the setup is kept. */
export interface AssistantOptions {
  /** Models the host brings with their providers (a story's scripted model), in effect with no key of the person's. */
  readonly provided?: readonly ProvidedModels[];
  /** Where the person's keys are kept; default the browser's (`createKeyStore`). */
  readonly keys?: KeyStore;
  /** Where the choice of model is kept; default localStorage (`localPrefs`). */
  readonly prefs?: SettingsPrefs;
  /**
   * The model until the person chooses one (default Anthropic's `claude-sonnet-5-5`): a provider id (`anthropic`,
   * `google`, `openai`, a provided one's) and a model id. Either may be left out: the provider is then the catalog's for
   * the model (else Anthropic), the model the provider's first listed. In effect while its provider has a key.
   */
  readonly provider?: string;
  readonly model?: string;
}

/** The host's model as `AssistantOptions` names it, with what it leaves out filled in (`Fallback` when nothing is). */
export function preferredModel(options: Pick<AssistantOptions, 'provider' | 'model' | 'provided'>): AISetup {
  const { model } = options;
  if (!options.provider && !model) return Fallback;
  const provider = options.provider ?? listedModels().find((m) => m.id === model)?.provider ?? Fallback.provider;
  const first =
    provider === Fallback.provider
      ? Fallback.model
      : (options.provided?.find((p) => p.id === provider)?.models[0] ?? listedModels(provider)[0]?.id);
  return { provider, model: model ?? first ?? Fallback.model };
}

/** The assistant as the viewer's `AssistantFactory`, with a host's options (a story's scripted model). */
export function createAssistantWith(options: AssistantOptions = {}): AssistantFactory {
  return (host) => makeAssistant(host, options);
}

export const createAssistant: AssistantFactory = createAssistantWith();

async function makeAssistant(host: AssistantHost, options: AssistantOptions): Promise<AssistantUI> {
  const { notifications } = host.viewer.managers;
  const keys = options.keys ?? createKeyStore();
  const prefs = options.prefs ?? localPrefs();
  // Domains the person allowed (agent.md §4): for the session, or kept in this browser (the prefs, as a list).
  const rules = createDomainRules(keptDomains(prefs.get('domains')), (browser) =>
    prefs.set('domains', JSON.stringify(browser)),
  );
  const preferred = preferredModel(options);
  const settings = createAISettings(keys, prefs, {
    report: (text) => notifications.post(text, 'warning'),
    ...(options.provided ? { provided: options.provided } : {}),
    ...(options.provider || options.model ? { preferred } : {}),
  });
  await settings.load();
  const agent = host.createAgent(
    {
      providerFor: (provider) => settings.providerFor(provider),
      default: settings.setup.value ?? preferred,
    },
    { rules },
  );
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
  const unregisterDomains = host.viewer.storage.register(domainsStorageSource(rules));
  return {
    ui: {
      tabs: [
        {
          id: 'assistant.tab.assistant',
          title: 'Assistant',
          icon: 'bot',
          // At the rail's foot, beside Settings (where the assistant's setup is too).
          order: 2.65,
          place: 'foot',
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
      unregisterDomains();
    },
  };
}

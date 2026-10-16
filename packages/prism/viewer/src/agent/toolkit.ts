/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  createToolkit,
  failed,
  type ToolGroup,
  type Toolkit,
  type ToolkitPart,
  type ToolRegistration,
} from '@molstar/ai';
import { type AgentPolicy, mayAsk } from './policies.js';
import { CoreTools, type OperationToolSpec, type SkillDeclaration, Skills } from './skills.js';
import { type AgentViewer, operationTool, tool, toolResult, type ViewerAgentApi } from './tools.js';

// The declarations are plain data in `skills.ts` (a spec extends them without this module and zod; agent.md §3).
export type { SkillDeclaration } from './skills.js';

// The core and the skills (agent.md §3): a small core from the first turn (Stories' is about 24 tools); the rest in
// skills, disclosed on request. A skill over a lazy feature is a lazy skill (plan track A6): in the catalog by its
// summary; disclosing it loads the feature, then its tools. A tool whose operation the viewer does not have (a feature
// left out of the spec) is left out of its skill, and a skill with none is left out.

/** The agent's own tools that are not operations: the scene's outline, asking the person, undoing its last turn. */
function ownTools(): ToolRegistration<ViewerAgentApi>[] {
  return [
    tool({
      id: 'undo_turn',
      description:
        'Undoes your previous turn at once (what it changed comes back as it was). Fails when the person changed something after it: then ask them.',
      input: { kind: 'object', fields: {} },
      modules: ['viewer'],
      run: async ({ api }) => {
        const turn = api.viewer.lastTurn;
        if (!turn) return failed('not-found', 'There is no earlier turn of yours to undo.');
        // Never the person's later steps (`discardLater`): that is theirs to choose.
        return toolResult(await api.viewer.api.run('history.undo-turn', { turn, discardLater: false }));
      },
    }),
    tool({
      id: 'get_scene',
      description:
        "The scene's outline: structures with their components and representations, volumes, measurements. detail: summary (the structures and their components' names), components (as each turn says it), full (the default: representations and hidden rows too).",
      input: {
        kind: 'object',
        fields: { detail: { kind: 'string', enum: ['summary', 'components', 'full'], optional: true } },
      },
      modules: ['account'],
      read: true,
      run: ({ args, api }) => {
        const detail = (args as { detail?: 'summary' | 'components' | 'full' } | undefined)?.detail ?? 'full';
        return { ok: true, did: 'The scene:', body: api.account.outline(400, detail) };
      },
    }),
    tool({
      id: 'ask',
      description:
        'Asks the person a question when the request is unclear (which chain? load a second copy?); end the turn after it: their answer comes as the next request.',
      input: {
        kind: 'object',
        fields: {
          question: { kind: 'string' },
          options: { kind: 'array', item: { kind: 'string' }, optional: true, description: 'Choices to offer.' },
        },
      },
      modules: [],
      run: ({ args }) => ({
        ok: true,
        did: `Asked: ${args.question}${args.options?.length ? ` (${args.options.join(' / ')})` : ''}. End your turn now.`,
      }),
    }),
  ];
}

/** The tools of a skill the viewer has now (a tool whose operation is missing is left out). */
function toolsOf(
  viewer: AgentViewer,
  specs: readonly OperationToolSpec[],
  policy?: AgentPolicy,
): ToolRegistration<ViewerAgentApi>[] {
  return specs.filter((s) => viewer.describe(s.operation)).map((s) => operationTool(viewer, s, policy));
}

export interface ViewerToolkitOptions {
  /** Loads a lazy extension (a lazy skill's load); true when it is loaded. */
  readonly loadExtension: (id: string) => Promise<boolean>;
  /** Whether an extension is lazy and not loaded yet. */
  readonly isStub: (id: string) => boolean;
  /** What goes without asking (agent.md §4). */
  readonly policy?: AgentPolicy;
  /** The core's operation tools (default `CoreTools`; the spec's `agent.core` applied to them). */
  readonly core?: readonly OperationToolSpec[];
  /** The skills: the viewer's (default `Skills`, the spec's `agent.skills` applied) and the extensions' (viewer-spec.md §3). */
  readonly skills?: readonly SkillDeclaration[];
}

/** The viewer's toolkit, and its refresh between turns. */
export interface ViewerToolkit extends Toolkit<ViewerAgentApi> {
  /**
   * The tools over operations described again as the viewer is now (a registry's choices, a feature added or gone:
   * agent.md §2); a part whose definitions changed is replaced, a tool gone is the eviction hint. Between turns only (a
   * turn's tool set is fixed at its start). True when anything changed.
   */
  refresh(): boolean;
}

interface RefreshablePart {
  readonly build: () => ToolkitPart<ViewerAgentApi>;
  signature: string;
  remove: () => void;
}

/** What a part tells the model: its tools' descriptions and schemas, its groups' tools. */
function signatureOf(part: ToolkitPart<ViewerAgentApi>): string {
  return JSON.stringify([
    (part.tools ?? []).map((t) => [t.definition.tool.id, t.definition.description, t.definition.inputSchema]),
    (part.groups ?? []).map((g) => [g.id, g.tools]),
  ]);
}

/** The viewer's toolkit: the core, and the skills (lazy ones over lazy features). */
export function createViewerToolkit(viewer: AgentViewer, options: ViewerToolkitOptions): ViewerToolkit {
  const skills = options.skills ?? Skills;
  const own = ownTools();
  const eager = skills.filter((s) => !s.feature || !options.isStub(s.feature));
  const lazy = skills.filter((s) => s.feature && options.isStub(s.feature));
  const toolkit = createToolkit<ViewerAgentApi>({ tools: own, core: own.map((t) => t.definition.tool.id), groups: [] });
  // Every part over operations is added on its own, so a refresh can replace it.
  const parts: RefreshablePart[] = [];
  const addPart = (build: () => ToolkitPart<ViewerAgentApi>) => {
    const part = build();
    parts.push({ build, signature: signatureOf(part), remove: toolkit.add(part) });
  };
  addPart(() => {
    const core = toolsOf(viewer, options.core ?? CoreTools, options.policy);
    return { tools: core, core: core.map((t) => t.definition.tool.id) };
  });
  const skillPart = (s: SkillDeclaration) => (): ToolkitPart<ViewerAgentApi> => {
    const tools = toolsOf(viewer, s.tools, options.policy);
    // A skill with no tool (its feature left out) is left out.
    return tools.length ? { tools, groups: [{ ...groupOf(s), tools: tools.map((t) => t.definition.tool.id) }] } : {};
  };
  for (const s of eager) addPart(skillPart(s));
  for (const s of lazy) {
    let removeStub = () => {};
    let loading: Promise<void> | undefined;
    removeStub = toolkit.add({
      groups: [
        {
          ...groupOf(s),
          tools: [],
          // The stubs declare their operations' flags (viewer-spec.md §3): a turn knows from its start which may ask,
          // by the rule the real tools use (`mayAsk`).
          pending: s.tools
            .filter((t) => viewer.describe(t.operation))
            .map((t) => ({ id: t.id, approval: mayAsk(viewer.describe(t.operation)!.definition) })),
          // Once (two disclosures share it); a failed load leaves the stub, and a later disclosure tries again.
          load: () =>
            (loading ??= (async () => {
              if (!(await options.loadExtension(s.feature!))) throw new Error(`${s.title} could not be loaded.`);
              removeStub();
              addPart(skillPart(s));
            })().catch((e) => {
              loading = undefined;
              throw e;
            })),
        },
      ],
    });
  }
  return Object.assign(toolkit, {
    refresh(): boolean {
      let changed = false;
      for (const part of parts) {
        const next = part.build();
        const signature = signatureOf(next);
        if (signature === part.signature) continue;
        // Out, then in: the same ids are taken by the new definitions.
        part.remove();
        part.remove = toolkit.add(next);
        part.signature = signature;
        changed = true;
      }
      return changed;
    },
  });
}

function groupOf(s: SkillDeclaration): Omit<ToolGroup<ViewerAgentApi>, 'tools'> {
  return { id: s.id, title: s.title, summary: s.summary, guidance: s.guidance };
}

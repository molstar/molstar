/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { estimatedTokens, type ToolDefinition } from '../messages/index.js';
import type { ToolGroupView, ToolsData, ToolView } from './types.js';

/**
 * The Tools view (analysis/viewer-mode.md §5) from plain inputs, so a host passes the contents of its toolkit without
 * this package knowing it: the tools with where they load, the skills that group them, and what a trace says is loaded.
 */

/**
 * Where a tool loads: with the first turn, with a skill (its tools, or its extras: loaded by name or when called), or by
 * name. Also as `tools/`' toolkit gives it (`placementOf`), so a host passes that as it is: a `group` is a
 * skill, and an `entry`'s tools show under the entry (pass the entry in `skills` for its title).
 */
export type ToolPlacement =
  | { readonly kind: 'core' }
  | { readonly kind: 'skill'; readonly skill: string; readonly extra?: boolean }
  | { readonly kind: 'single' }
  | { readonly kind: 'group'; readonly group: { readonly id: string }; readonly extra: boolean }
  | { readonly kind: 'entry'; readonly entry: { readonly id: string } };

/** The skill a placement puts a tool under, and whether as an extra; undefined for the core and single tools. */
function skillOf(placement: ToolPlacement): { readonly id: string; readonly extra: boolean } | undefined {
  switch (placement.kind) {
    case 'skill':
      return { id: placement.skill, extra: placement.extra === true };
    case 'group':
      return { id: placement.group.id, extra: placement.extra };
    case 'entry':
      return { id: placement.entry.id, extra: false };
    default:
      return undefined;
  }
}

export interface ToolEntry {
  readonly definition: ToolDefinition;
  readonly placement: ToolPlacement;
  /** The API modules it declares. */
  readonly modules?: readonly string[];
  /** It only reads (default false). */
  readonly read?: boolean;
  /** It needs the person's approval (default false). */
  readonly approval?: boolean;
  /** The turn's tool set includes it (default true). */
  readonly enabled?: boolean;
  /** Why this tool cannot be used here, when its skill's reason does not say it. */
  readonly unavailable?: string;
}

export interface SkillInput {
  readonly id: string;
  readonly title: string;
  /** One line: what it is for. */
  readonly summary?: string;
  /** Exactly as the model receives it when it loads the skill. */
  readonly guidance?: string;
  /** Why it cannot be used here. */
  readonly unavailable?: string;
}

export interface ToolsInput {
  readonly tools: readonly ToolEntry[];
  /** The skills, in the order shown; one a tool names but this lacks is shown by its id. */
  readonly skills?: readonly SkillInput[];
  /** The tools the conversation has loaded and the turn that first showed each: a trace's `loaded`. */
  readonly loaded?: readonly { readonly id: string; readonly turn?: number }[];
  /** The tools the conversation declares without loading: a trace's `declared`. */
  readonly declared?: readonly string[];
  readonly titles?: { readonly core?: string; readonly single?: string };
}

/** A definition's size in tokens, as the trace estimates it (characters / 4 of its JSON). */
export const definitionTokens = (definition: ToolDefinition): number => estimatedTokens(JSON.stringify(definition));

export function toolsView(input: ToolsInput): ToolsData {
  const loaded = new Map((input.loaded ?? []).map((t) => [t.id, t.turn]));
  const declared = new Set(input.declared ?? []);
  const view = (entry: ToolEntry, group: { unavailable?: string | undefined }): ToolView => {
    const { definition } = entry;
    const id = definition.tool.id;
    const unavailable = entry.unavailable ?? group.unavailable;
    const turn = loaded.get(id);
    return {
      id,
      name: definition.name,
      description: definition.description,
      input: definition.inputSchema,
      ...(entry.modules ? { modules: entry.modules } : {}),
      read: entry.read === true,
      approval: entry.approval === true,
      enabled: entry.enabled !== false,
      state: loaded.has(id) ? 'loaded' : declared.has(id) ? 'deferred' : 'offered',
      ...(turn !== undefined ? { loadedAtTurn: turn } : {}),
      tokens: definitionTokens(definition),
      ...(skillOf(entry.placement)?.extra ? { extra: true } : {}),
      ...(unavailable ? { unavailable } : {}),
    };
  };

  const core = input.tools.filter((t) => t.placement.kind === 'core');
  const single = input.tools.filter((t) => t.placement.kind === 'single');
  const skillIds = [
    ...(input.skills ?? []).map((s) => s.id),
    ...input.tools.flatMap((t) => skillOf(t.placement)?.id ?? []),
  ];
  const groups: ToolGroupView[] = [];
  if (core.length > 0)
    groups.push({
      id: 'core',
      kind: 'core',
      title: input.titles?.core ?? 'Core tools',
      tools: core.map((t) => view(t, {})),
    });
  for (const id of [...new Set(skillIds)]) {
    const skill = input.skills?.find((s) => s.id === id);
    const own = input.tools.filter((t) => skillOf(t.placement)?.id === id);
    const ordered = [
      ...own.filter((t) => !skillOf(t.placement)?.extra),
      ...own.filter((t) => skillOf(t.placement)?.extra),
    ];
    groups.push({
      id,
      kind: 'skill',
      title: skill?.title ?? id,
      ...(skill?.summary ? { summary: skill.summary } : {}),
      ...(skill?.guidance !== undefined
        ? { guidance: skill.guidance, guidanceTokens: estimatedTokens(skill.guidance) }
        : {}),
      ...(skill?.unavailable ? { unavailable: skill.unavailable } : {}),
      tools: ordered.map((t) => view(t, { unavailable: skill?.unavailable })),
    });
  }
  if (single.length > 0)
    groups.push({
      id: 'single',
      kind: 'single',
      title: input.titles?.single ?? 'Tools loaded by name',
      tools: single.map((t) => view(t, {})),
    });
  return { groups };
}

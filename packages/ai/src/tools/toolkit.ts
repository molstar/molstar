/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Message, ToolDefinition } from '../messages/index.js';
import { z } from 'zod';
import { closestHint, suggest } from './closest.js';
import { DISCLOSE_SKILLS, DISCLOSE_TOOLS, LABEL_DEFINITIONS, LABEL_ENTRY, LABEL_SKILL } from './names.js';
import type { ToolRegistration } from './registry.js';
import { type ToolInput, zodBacked } from './input.js';
import { failed } from './result.js';
import { defineTool } from './tool.js';

/**
 * What the agent sees of its tools, and when (tools-and-session.md §3): a small core from the first turn; everything
 * else loaded on demand, added where the model asks for it, so the prompt's prefix only grows at its end.
 * - A group (a skill) is a family of tools with the know-how to use them: the disclosure tool adds its guidance as a
 *   context entry keyed by the group (`groupKey`: what says it is loaded, for disclosure and compaction alike) and its
 *   tools' definitions; its extras, which its guidance names, load by name or when called.
 * - A tool outside the core and the groups is loaded by name, several at once.
 * - An entry is another text the model may load the same way (a data source's notes, a server's tools), with the
 *   tools that come with it.
 * The toolkit holds the registrations and these declarations, checks that each tool has exactly one place, and
 * generates the disclosure tools whose descriptions are the catalog.
 */

/** What the disclosure tools read of the conversation: the host provides it as the API's `disclosure` module. */
export interface DisclosureModule {
  /** The ids of the tools the conversation has loaded. */
  loaded(): Promise<string[]>;
  /** The keys (`groupKey`) of the groups and entries whose text is in the conversation (`groupsIn`). */
  skills(): Promise<string[]>;
}

/** The module the disclosure tools need of the host's API. */
export interface DisclosureApi {
  disclosure?: DisclosureModule;
}

export interface ToolGroup<Api> {
  readonly id: string;
  /** Its name for people (the Tools view: "Animation skill & tools"). */
  readonly title: string;
  /** One line for the catalog: what it is for. */
  readonly summary: string;
  /** What the person wants it for: a host's own listing may group by these. */
  readonly tags?: readonly string[];
  readonly guidance: string;
  /** Its tools, by id: loaded with its guidance (what nearly every use of it calls). */
  readonly tools: readonly string[];
  /** Tools its guidance names that load by name, or when called. */
  readonly extras?: readonly string[];
  /** The modules `unavailable` reads: the disclosure tool is given them. */
  readonly modules?: readonly (keyof Api)[];
  /** Why it cannot be used here, if it cannot (a module the host does not provide). */
  unavailable?(api: Partial<Api>): { message: string; hint?: string } | undefined;
  /**
   * A lazy skill (tools-and-session.md §3, a dynamic toolkit): in the catalog by its summary before its code is
   * loaded. Disclosing it calls `load`, which brings its tools by replacing this group with the real one
   * (`Toolkit.add`, the stub's part removed); the disclosure then goes on with the group as it is after the load.
   */
  load?(): Promise<void>;
  /**
   * The tools a lazy skill brings, by id, before they are registered: a turn enables them from its start (a turn's
   * tool set is fixed when it begins), preauthorized unless `approval`.
   */
  readonly pending?: readonly { readonly id: string; readonly approval?: boolean }[];
}

/** Another disclosable text: loaded with the groups' disclosure tool, by its id, with the tools that come with it. */
export interface DisclosableEntry {
  readonly id: string;
  readonly summary: string;
  readonly text: string;
  /** Tools loaded with it, never alone. */
  readonly tools?: readonly string[];
}

/** Where a tool is: the core, a group (its tools, or its extras), an entry (never alone), or loaded singly. */
export type Placement<Api> =
  | { kind: 'core' }
  | { kind: 'group'; group: ToolGroup<Api>; extra: boolean }
  | { kind: 'entry'; entry: DisclosableEntry }
  | { kind: 'single' };

export interface ToolkitSpec<Api> {
  /** Every tool but the disclosure tools, which the toolkit makes. */
  tools: readonly ToolRegistration<Api>[];
  /** Tool ids offered from the first turn (the disclosure tools are always among them). */
  core: readonly string[];
  groups: readonly ToolGroup<Api>[];
  entries?: readonly DisclosableEntry[];
  /** The disclosure tools' ids (default `disclose_skills`, `disclose_tools`). */
  names?: { skills?: string; tools?: string };
}

/** What a host adds to a toolkit after it is made (an extension loaded): removed again by `add`'s return. */
export interface ToolkitPart<Api> {
  readonly tools?: readonly ToolRegistration<Api>[];
  /** Of these tools, the ones offered from the first turn. */
  readonly core?: readonly string[];
  readonly groups?: readonly ToolGroup<Api>[];
  readonly entries?: readonly DisclosableEntry[];
}

export interface Toolkit<Api> {
  /** The registrations: the given ones and the disclosure tools (as they are now). */
  readonly tools: readonly ToolRegistration<Api>[];
  readonly groups: readonly ToolGroup<Api>[];
  readonly entries: readonly DisclosableEntry[];
  readonly names: { readonly skills: string; readonly tools: string };
  placementOf(toolId: string): Placement<Api>;
  /** The key of a group's or entry's text in the context: it and a fingerprint of its text. */
  groupKey(group: ToolGroup<Api> | DisclosableEntry): string;
  /** A group's guidance as the model reads it: the text, and which of its tools load by name. */
  groupText(group: ToolGroup<Api>): string;
  /** A group's guidance, or an entry's text, as a context entry: identical text is not appended again. */
  groupEntry(group: ToolGroup<Api> | DisclosableEntry): Message;
  /** The keys of the groups and entries whose text is among these messages (as it is now: a changed text is not). */
  groupsIn(messages: readonly Message[]): string[];
  coreDefinitions(): ToolDefinition[];
  allDefinitions(): ToolDefinition[];
  definitionEntries(definitions: readonly ToolDefinition[], deferred?: boolean): Message[];
  /** The two disclosure registrations (already among `tools`). */
  disclosureTools(): ToolRegistration<Api>[];
  /** A registration by tool id, as the toolkit is now (bindings resolve through it at call time). */
  registration(toolId: string): ToolRegistration<Api> | undefined;
  /**
   * A dynamic toolkit (decision 25 of the prism-ui spec, `.v6/prism-ui/spec/README.md`): adds tools, groups and entries between and during turns; the same checks
   * as at creation (one place for each tool), nothing changed when they fail. Returns the removal.
   */
  add(part: ToolkitPart<Api>): () => void;
  /** Bumped by every change. */
  readonly version: number;
  /** Called after every change. */
  onChange(listener: () => void): () => void;
  /** Tools the toolkit had and no longer has: compaction drops their definitions (the eviction hint). */
  removed(): readonly string[];
  /** What a turn enables and preauthorizes: the registered tools and those lazy skills will bring. */
  toolSet(): { enabled: string[]; preauthorized: string[] };
}

/** A short fingerprint of a text, for a context entry's key (FNV-1a). */
export function fingerprint(text: string): string {
  let h = 0x811c9dc5;
  for (let i = 0; i < text.length; i++) {
    h ^= text.charCodeAt(i);
    h = Math.imul(h, 0x01000193);
  }
  return (h >>> 0).toString(16);
}

/** Loads tools: one system message holding their definitions (never text; harness.md §2). */
export const definitionEntries = (definitions: readonly ToolDefinition[], deferred = false): Message[] =>
  definitions.length === 0
    ? []
    : [
        {
          role: 'system',
          parts: definitions.map((definition) => ({ kind: 'tool-definition', definition, deferred })),
          label: { kind: LABEL_DEFINITIONS },
        },
      ];

/** The first sentence of a description, for a catalog line. */
function firstSentence(text: string): string {
  const m = /^(.+?[.!?])(\s|$)/.exec(text);
  return (m?.[1] ?? text).replace(/\s+/g, ' ');
}

const isGroup = <Api>(x: ToolGroup<Api> | DisclosableEntry): x is ToolGroup<Api> => 'guidance' in x;

/**
 * A toolkit over registrations and declarations. Throws, naming the tool, when a tool has two places or a declaration
 * names a tool that is not registered.
 */
export function createToolkit<Api extends DisclosureApi>(spec: ToolkitSpec<Api>): Toolkit<Api> {
  const names = { skills: spec.names?.skills ?? DISCLOSE_SKILLS, tools: spec.names?.tools ?? DISCLOSE_TOOLS };
  const parts: ToolkitPart<Api>[] = [];
  const of = (list: readonly ToolkitPart<Api>[]): ToolkitSpec<Api> => ({
    tools: [...spec.tools, ...list.flatMap((p) => p.tools ?? [])],
    core: [...spec.core, ...list.flatMap((p) => p.core ?? [])],
    groups: [...spec.groups, ...list.flatMap((p) => p.groups ?? [])],
    entries: [...(spec.entries ?? []), ...list.flatMap((p) => p.entries ?? [])],
    names,
  });
  // The disclosure tools read the toolkit as it is when they run (a lazy skill's load changes it under them).
  let built: Built<Api> = buildToolkit(of(parts), (): Built<Api> => built);
  let version = 0;
  const listeners = new Set<() => void>();
  const removed = new Set<string>();
  const change = (next: Built<Api>) => {
    const before = new Set(built.tools.map((t) => t.definition.tool.id));
    built = next;
    const now = new Set(built.tools.map((t) => t.definition.tool.id));
    for (const id of before) if (!now.has(id)) removed.add(id);
    for (const id of now) removed.delete(id);
    version++;
    for (const l of [...listeners]) l();
  };
  return {
    get tools() {
      return built.tools;
    },
    get groups() {
      return built.groups;
    },
    get entries() {
      return built.entries;
    },
    names,
    placementOf: (id) => built.placementOf(id),
    groupKey: (x) => built.groupKey(x),
    groupText: (g) => built.groupText(g),
    groupEntry: (x) => built.groupEntry(x),
    groupsIn: (m) => built.groupsIn(m),
    coreDefinitions: () => built.coreDefinitions(),
    allDefinitions: () => built.allDefinitions(),
    definitionEntries,
    disclosureTools: () => built.disclosure,
    registration: (id) => built.byId.get(id),
    add(part) {
      // Checked by building it: a conflict throws and leaves the toolkit as it was.
      const next = buildToolkit(of([...parts, part]), () => built);
      parts.push(part);
      change(next);
      return () => {
        const at = parts.indexOf(part);
        if (at < 0) return;
        parts.splice(at, 1);
        change(buildToolkit(of(parts), () => built));
      };
    },
    get version() {
      return version;
    },
    onChange(listener) {
      listeners.add(listener);
      return () => listeners.delete(listener);
    },
    removed: () => [...removed],
    toolSet: () => {
      const pending = built.groups.flatMap((g) => g.pending ?? []);
      return {
        enabled: [...built.byId.keys(), ...pending.map((p) => p.id)],
        preauthorized: [
          ...built.tools.filter((t) => !t.approval).map((t) => t.definition.tool.id),
          ...pending.filter((p) => !p.approval).map((p) => p.id),
        ],
      };
    },
  };
}

/** A toolkit's state for one set of declarations: what `createToolkit` rebuilds on each change. */
interface Built<Api> {
  readonly tools: readonly ToolRegistration<Api>[];
  readonly groups: readonly ToolGroup<Api>[];
  readonly entries: readonly DisclosableEntry[];
  readonly byId: ReadonlyMap<string, ToolRegistration<Api>>;
  readonly disclosure: ToolRegistration<Api>[];
  placementOf(toolId: string): Placement<Api>;
  groupKey(group: ToolGroup<Api> | DisclosableEntry): string;
  groupText(group: ToolGroup<Api>): string;
  groupEntry(group: ToolGroup<Api> | DisclosableEntry): Message;
  groupsIn(messages: readonly Message[]): string[];
  coreDefinitions(): ToolDefinition[];
  allDefinitions(): ToolDefinition[];
}

function buildToolkit<Api extends DisclosureApi>(spec: ToolkitSpec<Api>, current: () => Built<Api>): Built<Api> {
  const names = { skills: spec.names?.skills ?? DISCLOSE_SKILLS, tools: spec.names?.tools ?? DISCLOSE_TOOLS };
  const entries = spec.entries ?? [];
  const disclosureIds = [names.skills, names.tools];

  // ---------------------------------------------------------------- one place for each tool
  const registered = new Set<string>();
  for (const t of spec.tools) {
    const id = t.definition.tool.id;
    if (registered.has(id)) throw new Error(`Tool "${id}" is registered twice.`);
    if (disclosureIds.includes(id)) throw new Error(`Tool "${id}" is a disclosure tool: the toolkit makes it.`);
    registered.add(id);
  }
  const placed = new Map<string, { placement: Placement<Api>; where: string }>();
  const put = (id: string, placement: Placement<Api>, where: string) => {
    if (!registered.has(id) && !disclosureIds.includes(id))
      throw new Error(`Tool "${id}" is listed in ${where} but is not registered.`);
    const had = placed.get(id);
    if (had) throw new Error(`Tool "${id}" has more than one place: ${had.where} and ${where}.`);
    placed.set(id, { placement, where });
  };
  for (const id of new Set([...spec.core, ...disclosureIds])) put(id, { kind: 'core' }, 'the core');
  const owners = new Set<string>();
  for (const group of spec.groups) {
    if (owners.has(group.id)) throw new Error(`Group "${group.id}" is declared twice.`);
    owners.add(group.id);
    for (const p of group.pending ?? [])
      if (registered.has(p.id)) throw new Error(`Tool "${p.id}" is pending in group "${group.id}" but is registered.`);
    for (const id of group.tools) put(id, { kind: 'group', group, extra: false }, `group "${group.id}"`);
    for (const id of group.extras ?? []) put(id, { kind: 'group', group, extra: true }, `group "${group.id}" (extras)`);
  }
  for (const entry of entries) {
    if (owners.has(entry.id)) throw new Error(`"${entry.id}" is declared twice (a group or an entry).`);
    owners.add(entry.id);
    for (const id of entry.tools ?? []) put(id, { kind: 'entry', entry }, `entry "${entry.id}"`);
  }
  const placementOf = (toolId: string): Placement<Api> => placed.get(toolId)?.placement ?? { kind: 'single' };

  // ---------------------------------------------------------------- texts and keys
  const groupText = (group: ToolGroup<Api>): string => {
    const extras = group.extras?.length
      ? `\n\nLoaded with this skill: ${group.tools.join(', ')}. By name when needed (${names.tools}, or call it): ${group.extras.join(', ')}.`
      : '';
    return `${group.guidance}${extras}`;
  };
  const textOf = (x: ToolGroup<Api> | DisclosableEntry) => (isGroup(x) ? groupText(x) : x.text);
  const groupKey = (x: ToolGroup<Api> | DisclosableEntry) =>
    `${isGroup(x) ? 'skill' : 'entry'}:${x.id}@${fingerprint(textOf(x))}`;
  const groupEntry = (x: ToolGroup<Api> | DisclosableEntry): Message => ({
    role: 'user',
    parts: [{ kind: 'text', text: textOf(x) }],
    label: { kind: isGroup(x) ? LABEL_SKILL : LABEL_ENTRY, key: groupKey(x) },
  });
  // By the entry's label, whose key holds its text's fingerprint (an earlier text is not the group's now).
  const groupsIn = (messages: readonly Message[]): string[] => {
    const current = new Set([...spec.groups, ...entries].map(groupKey));
    return messages.flatMap((m) =>
      (m.label?.kind === LABEL_SKILL || m.label?.kind === LABEL_ENTRY) && current.has(m.label.key ?? '')
        ? [m.label.key as string]
        : [],
    );
  };

  // ---------------------------------------------------------------- disclosure
  const disclosure = disclosureTools(spec, names, entries, { placementOf, groupText, groupKey, groupEntry }, current);
  const tools = [...spec.tools, ...disclosure];
  const byId = new Map(tools.map((t) => [t.definition.tool.id, t.definition]));
  const definitionsOf = (ids: readonly string[]) => ids.flatMap((id) => byId.get(id) ?? []);
  return {
    tools,
    groups: spec.groups,
    entries,
    byId: new Map(tools.map((t) => [t.definition.tool.id, t])),
    disclosure,
    placementOf,
    groupKey,
    groupText,
    groupEntry,
    groupsIn,
    coreDefinitions: () => definitionsOf([...placed].filter(([, p]) => p.placement.kind === 'core').map(([id]) => id)),
    allDefinitions: () => tools.map((t) => t.definition),
  };
}

interface Parts<Api> {
  placementOf(toolId: string): Placement<Api>;
  groupText(group: ToolGroup<Api>): string;
  groupKey(x: ToolGroup<Api> | DisclosableEntry): string;
  groupEntry(x: ToolGroup<Api> | DisclosableEntry): Message;
}

/**
 * The disclosure tools, built over the registrations (their catalogs name the tools there). Both read which tools the
 * conversation already has from the API's `disclosure` module.
 */
function disclosureTools<Api extends DisclosureApi>(
  spec: ToolkitSpec<Api>,
  names: { skills: string; tools: string },
  entries: readonly DisclosableEntry[],
  { placementOf, groupKey, groupEntry }: Parts<Api>,
  current: () => Built<Api>,
): ToolRegistration<Api>[] {
  const byId = new Map(spec.tools.map((t) => [t.definition.tool.id, t.definition]));
  const singles = spec.tools.filter((t) => placementOf(t.definition.tool.id).kind === 'single');
  const loaded = async (api: DisclosureApi) => new Set(api.disclosure ? await api.disclosure.loaded() : []);
  /** The keys of the groups' guidance in the context: a group is loaded exactly when its key is there. */
  const keys = async (api: DisclosureApi) => new Set(api.disclosure ? await api.disclosure.skills() : []);
  const enumOf = (list: readonly string[]) => (list.length ? { enum: [...list] } : {});

  const modules = [...new Set<keyof Api>(['disclosure' as keyof Api, ...spec.groups.flatMap((g) => g.modules ?? [])])];
  const groupIds = spec.groups.map((g) => g.id);
  const also = entries.length
    ? ` Other texts, loaded the same way: ${entries.map((e) => `"${e.id}" (${e.summary})`).join('; ')}.`
    : '';

  const discloseSkills = defineTool<Api, keyof Api, ToolInput<{ skills: string[] }>>({
    id: names.skills,
    description: `Load skills: their guidance and their tools, all a request needs in one call, once per conversation. Skills: ${spec.groups.map((s) => `"${s.id}" (${s.summary})`).join('; ')}.${also}`,
    input: zodBacked<{ skills: string[] }>(
      { kind: 'object', fields: { skills: { kind: 'array', item: { kind: 'string' }, description: 'Skill names.' } } },
      z.object({
        // Ids as an enum in the schema (the groups only, so it never changes with entries), but any id parses: an
        // entry's id, or one it does not know, is answered.
        skills: z.array(z.string().meta(enumOf(groupIds))).describe('Skill names.'),
      }),
    ),
    modules,
    run: async ({ args, api }) => {
      const have = await loaded(api);
      const given = await keys(api);
      const load: (ToolGroup<Api> | DisclosableEntry)[] = [];
      const already: string[] = [];
      const refused: { message: string; hint?: string }[] = [];
      const wanted = new Set(args.skills.map((id) => id.trim().toLowerCase()));
      const known = [...groupIds, ...entries.map((e) => e.id)];
      for (const id of wanted)
        if (!known.some((k) => k.toLowerCase() === id))
          refused.push({ message: `No skill "${id}".`, hint: closestHint('skills and other texts', id, known) });
      // A lazy skill loads first: its real group (and tools) replace the stub's, and the rest reads the toolkit after.
      const lazy = spec.groups.filter((g) => g.load && wanted.has(g.id.toLowerCase()));
      for (const g of lazy) {
        try {
          await g.load!();
        } catch (e) {
          refused.push({
            message: `The ${g.id} skill could not be loaded: ${e instanceof Error ? e.message : String(e)}`,
          });
        }
      }
      const now = lazy.length ? current() : undefined;
      const all: (ToolGroup<Api> | DisclosableEntry)[] = now
        ? [...now.groups, ...now.entries]
        : [...spec.groups, ...entries];
      const keyOf = now ? now.groupKey : groupKey;
      const byIdNow = now ? new Map(now.tools.map((t) => [t.definition.tool.id, t.definition])) : byId;
      for (const x of all.filter((x) => wanted.has(x.id.toLowerCase()))) {
        const unavailable = isGroup(x) ? x.unavailable?.(api) : undefined;
        if (unavailable) refused.push(unavailable);
        else if (isGroup(x) && x.load) refused.push({ message: `The ${x.id} skill is not loaded.` });
        else if (given.has(keyOf(x))) already.push(x.id);
        else load.push(x);
      }
      if (!load.length && !already.length) {
        const [first] = refused;
        return failed(
          'skill.unavailable',
          refused.map((r) => r.message).join(' '),
          refused.length === 1 ? first?.hint : undefined,
        );
      }
      const toolIds = [...new Set(load.flatMap((x) => x.tools ?? []))].filter((id) => !have.has(id));
      const said = [
        ...load.map((x) =>
          isGroup(x)
            ? `Skill ${x.id} loaded; its tools: ${x.tools.join(', ')}${x.extras?.length ? `; by name when needed: ${x.extras.join(', ')}` : ''}.`
            : `${x.id} loaded${x.tools?.length ? `; its tools: ${x.tools.join(', ')}` : ''}.`,
        ),
        ...(already.length
          ? [
              `Skill${already.length > 1 ? 's' : ''} ${already.join(', ')} ${already.length > 1 ? 'are' : 'is'} already loaded: the guidance is earlier in this conversation, and the tools are available.`,
            ]
          : []),
      ];
      // The guidance follows this result as its own entry, keyed by the group (after the step's results: the prompt's
      // prefix does not change). The tools' definitions come with it.
      const contributions = [
        ...definitionEntries(toolIds.flatMap((id) => byIdNow.get(id) ?? [])),
        ...load.map(now ? now.groupEntry : groupEntry),
      ];
      return {
        ok: true,
        did: `${said.join(' ')}${load.length ? ' Its guidance follows.' : ''}`,
        // Loading changes the conversation, which the harness keeps, and nothing of the host's: no write hook.
        effect: 'none',
        ...(refused.length ? { warnings: refused.map((r) => `${r.message}${r.hint ? ` ${r.hint}` : ''}`) } : {}),
        ...(contributions.length ? { contributions } : {}),
      };
    },
  });

  const singleIds = singles.map((t) => t.definition.tool.id);
  const catalog = singles.length
    ? ` Tools: ${singles.map((t) => `${t.definition.name} (${firstSentence(t.definition.description)})`).join('; ')}.`
    : '';
  const discloseTools = defineTool<Api, keyof Api, ToolInput<{ tools: string[] }>>({
    id: names.tools,
    description: `Load tools by name, all you need in one call; they stay available for the rest of the conversation.${catalog} Tools that belong to a skill come with ${names.skills}; a skill's guidance names the ones it loads by name.`,
    input: zodBacked<{ tools: string[] }>(
      { kind: 'object', fields: { tools: { kind: 'array', item: { kind: 'string' }, description: 'Tool names.' } } },
      z.object({
        // The names as an enum in the schema, but any name parses: one in a skill, or none, is answered.
        tools: z.array(z.string().meta(enumOf(singleIds))).describe('Tool names.'),
      }),
    ),
    modules: ['disclosure' as keyof Api],
    run: async ({ args, api }) => {
      if (args.tools.length === 0)
        return failed(
          'tools.none',
          `Name the tools to load${singleIds[0] ? `, e.g. tools: ["${singleIds[0]}"]` : ''}.`,
        );
      const have = await loaded(api);
      const load: ToolDefinition[] = [];
      const already: string[] = [];
      const problems: string[] = [];
      const all = [...byId.keys()];
      for (const name of new Set(args.tools.map((t) => t.trim()))) {
        const place = placementOf(name);
        const definition = byId.get(name);
        if (!definition) problems.push(`no tool "${name}" (closest: ${suggest(name, all, 3).join(', ')})`);
        else if (place.kind === 'group' && !place.extra)
          problems.push(`${name} belongs to the ${place.group.id} skill: ${names.skills}(["${place.group.id}"])`);
        else if (place.kind === 'entry' && !have.has(name))
          problems.push(`${name} comes with ${place.entry.id}: ${names.skills}(["${place.entry.id}"])`);
        else if (place.kind === 'core' || have.has(name)) already.push(name);
        else load.push(definition);
      }
      if (load.length === 0 && already.length === 0) return failed('tools.unknown', `${problems.join('; ')}.`);
      const said = [
        load.length ? `Loaded: ${load.map((d) => d.name).join(', ')}.` : '',
        already.length ? `Already available: ${already.join(', ')}.` : '',
      ].filter(Boolean);
      return {
        ok: true,
        did: said.join(' '),
        effect: 'none',
        ...(problems.length ? { warnings: [`Not loaded: ${problems.join('; ')}.`] } : {}),
        ...(load.length ? { contributions: definitionEntries(load) } : {}),
      };
    },
  });

  return [discloseSkills, discloseTools];
}

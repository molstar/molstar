/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Message } from '../../messages/index.js';
import type { DisclosableEntry, ToolGroup, Toolkit } from '../../tools/index.js';
import { LABEL_REQUEST } from '../compose.js';

/**
 * What compaction keeps is the host's (tools-and-session.md §6): what a request is, the core tools, the groups
 * (skills) with their guidance, the memory, the pins. A session builds the domain from its toolkit and the request
 * label; the host adds the rest.
 */

/** A skill or another disclosable text: what compaction needs to know of it. */
export interface CompactionGroup {
  readonly id: string;
  /** The tools that come with it. */
  readonly tools?: readonly string[];
}

export interface CompactionDomain {
  /** Whether a message is a turn's request: a turn starts at each request and runs to the next. */
  isRequest(message: Message): boolean;
  /** The person's own words in a request message (default: the prompt in its `LABEL_REQUEST` label). */
  requestText?(message: Message): string;
  /** The tools every conversation has loaded: always kept. */
  coreToolIds: readonly string[];
  groups: readonly CompactionGroup[];
  /** The ids of the groups whose text is among these messages (as given, whatever its tools). */
  groupsIn(messages: readonly Message[]): string[];
  /** A group's guidance (or an entry's text) as the context entry that gives it. */
  groupEntry(group: CompactionGroup): Message;
  /** The memory as it now stands: the new head starts from it (its events go with the old context). */
  memory?(): string;
  /** Skills and tools to keep whatever the model says (the host's own: a story authored in code). */
  pins?(): { skills: readonly string[]; tools: readonly string[] };
  /** More texts a wish may quote as the person's own (their notes in the memory). */
  wishSources?(): readonly string[];
  /** Tools the host no longer has (a dynamic toolkit's removed ones): their definitions are not kept. */
  evicted?(): readonly string[];
}

const textOf = (m: Message): string => m.parts.map((p) => (p.kind === 'text' ? p.text : '')).join('');

/** The person's request in a turn's user message: its label's prompt (a message without one: all its text). */
export function requestTextOf(m: Message): string {
  const prompt = m.label?.kind === LABEL_REQUEST ? m.label.data?.prompt : undefined;
  return typeof prompt === 'string' ? prompt : textOf(m);
}

/** The domain a toolkit gives: its core, its groups and entries, a request labelled `LABEL_REQUEST`. */
export function toolkitDomain<Api>(toolkit: Toolkit<Api>): CompactionDomain {
  const all: (ToolGroup<Api> | DisclosableEntry)[] = [...toolkit.groups, ...toolkit.entries];
  const byId = new Map(all.map((x) => [x.id, x] as const));
  return {
    isRequest: (m) => m.label?.kind === LABEL_REQUEST,
    requestText: requestTextOf,
    coreToolIds: toolkit.coreDefinitions().map((d) => d.tool.id),
    groups: all.map((x) => ({ id: x.id, ...(x.tools ? { tools: x.tools } : {}) })),
    // The toolkit finds the keys of the texts in the context (`skill:id@…`, `entry:id@…`); a group is loaded when its
    // key is there, whatever its tools.
    groupsIn: (messages) => {
      const keys = toolkit.groupsIn(messages);
      return all
        .filter((x) => keys.some((k) => k.startsWith(`skill:${x.id}@`) || k.startsWith(`entry:${x.id}@`)))
        .map((x) => x.id);
    },
    groupEntry: (group) => {
      const found = byId.get(group.id);
      if (!found) throw new Error(`No group "${group.id}".`);
      return toolkit.groupEntry(found);
    },
  };
}

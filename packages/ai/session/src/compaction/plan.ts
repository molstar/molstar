/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Message, ToolDefinition } from '@molstar/ai-agent';
import type { HeadEdit } from '@molstar/ai-harness';
import { type CompactionDomain, requestTextOf } from './domain.js';
import { CARRY_OVER_HEADER, type CompactionSettings, SIZE_NOTE } from './settings.js';

/**
 * The plan of a compaction (tools-and-session.md §6): what the model's answer says to carry over, checked, and the
 * new head built from it.
 */

export interface CarryOver {
  wishes: string[];
  decisions: string;
  open: string;
  notes: string;
  thisTurn: string;
  retention: { tools: string[]; skills: string[] };
}

const SECTIONS = ['wishes', 'decisions', 'open', 'notes', 'this turn'] as const;

/** The sections of the compaction answer; undefined when none is there. */
export function parseCarryOver(text: string): CarryOver | undefined {
  const section = (name: string) => {
    const m = new RegExp(`\\[${name}\\]([\\s\\S]*?)\\[\\/${name}\\]`, 'i').exec(text);
    return m ? (m[1] ?? '').trim() : undefined;
  };
  const found = SECTIONS.map(section);
  if (found.every((s) => s === undefined)) return undefined;
  const list = (s: string | undefined) =>
    (s ?? '')
      .split(/[,\n]/)
      .map((x) => x.trim().replace(/^["'`]|["'`]$/g, ''))
      .filter((x) => x && !/^(none|nothing|-)$/i.test(x));
  const retention = section('retention') ?? '';
  const inner = (name: string) => new RegExp(`\\[${name}\\]([\\s\\S]*?)\\[\\/${name}\\]`, 'i').exec(retention)?.[1];
  return {
    wishes: (found[0] ?? '')
      .split('\n')
      .map((l) => l.trim())
      .filter((l) => l.startsWith('-')),
    decisions: found[1] ?? '',
    open: found[2] ?? '',
    notes: found[3] ?? '',
    thisTurn: found[4] ?? '',
    retention: { tools: list(inner('tools')), skills: list(inner('skills')) },
  };
}

/**
 * A quote reduced to what any text layer keeps: letters and digits only, accents and case gone. A PDF's text breaks
 * words at line ends ("under-\ngoes") and inside them ("com plex"), splits accents and loses symbols; none of that
 * decides whether the words are there.
 */
export function normalizedForQuotes(s: string): string {
  return s
    .normalize('NFKD')
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '');
}

/** Quoted phrases in a wish line ("…" or “…”). */
function quotesIn(line: string): string[] {
  return [...line.matchAll(/"([^"]+)"|“([^”]+)”/g)].map((m) => (m[1] ?? m[2] ?? '').trim()).filter(Boolean);
}

/**
 * Wishes whose quotes are the person's own words (found in their requests or in wishes carried over before), and the
 * ones dropped. A wish without a quote is dropped; one that also quotes something else (a tool's error, say) keeps
 * only the person's quotes.
 */
export function checkWishes(
  wishes: readonly string[],
  sources: readonly string[],
): { kept: string[]; dropped: string[] } {
  const haystack = sources.map(normalizedForQuotes).join('\n');
  const kept: string[] = [];
  const dropped: string[] = [];
  for (const w of wishes) {
    const quotes = quotesIn(w);
    const theirs = quotes.filter((q) => {
      const n = normalizedForQuotes(q);
      return n.length >= 4 && haystack.includes(n);
    });
    if (theirs.length === 0) dropped.push(w);
    else kept.push(theirs.length === quotes.length ? w : `- ${theirs.map((q) => `"${q}"`).join(' ')}`);
  }
  return { kept, dropped };
}

// ------------------------------------------------------------------ reading the head

/** A head entry as compaction reads it: its id (the message's hash) and the message. */
export interface HeadEntry {
  readonly id: string;
  readonly message: Message;
}

const textOf = (m: Message): string => m.parts.map((p) => (p.kind === 'text' ? p.text : '')).join('');

const isResults = (m: Message) => m.role === 'user' && m.parts.some((p) => p.kind === 'tool-result');
const isCarryOver = (m: Message) => m.role === 'user' && textOf(m).startsWith(CARRY_OVER_HEADER);

/** Wish lines of a carry-over this code wrote earlier. */
function carriedWishes(m: Message): string[] {
  if (!isCarryOver(m)) return [];
  const match = /Wishes:\n([\s\S]*?)(\n\n|$)/.exec(textOf(m));
  return (match?.[1] ?? '').split('\n').filter((l) => l.startsWith('-'));
}

/** An earlier carry-over's text after its header, without its wishes (carried apart). */
function carriedBody(m: Message): string[] {
  if (!isCarryOver(m)) return [];
  const body = textOf(m)
    .slice(CARRY_OVER_HEADER.length)
    .replace(/Wishes:\n[\s\S]*?(\n\n|$)/, '')
    .replace(new RegExp(`${SIZE_NOTE}[^\n]*`, 'g'), '')
    .trim();
  return body ? [body] : [];
}

/** A message's size in tokens, estimated (about 3 characters a token errs high for both providers). */
export const sizeOf = (m: Message) => Math.ceil(JSON.stringify(m).length / 3);

/** A turn of the conversation: its request, the model's messages and the results, in order. */
export interface Turn {
  readonly start: number;
  readonly entries: HeadEntry[];
}

/**
 * The head's conversation, turn by turn: a turn starts at each request and runs to the next; system messages and
 * notes (skills, carry-overs, the memory) stand apart.
 */
function turnsOf(entries: readonly HeadEntry[], domain: CompactionDomain): Turn[] {
  const turns: Turn[] = [];
  for (const [i, e] of entries.entries()) {
    const m = e.message;
    if (m.role === 'system') continue;
    if (domain.isRequest(m)) turns.push({ start: i, entries: [] });
    else if (m.role === 'user' && !isResults(m)) continue;
    turns.at(-1)?.entries.push(e);
  }
  return turns;
}

const toolCallsOf = (m: Message): string[] =>
  m.role === 'assistant' ? m.parts.flatMap((p) => (p.kind === 'tool-call' ? [p.tool.id] : [])) : [];

/** A mechanical digest of turns: the request and reply verbatim (clipped), the tools as a line. */
export function digestOf(turns: readonly Turn[], domain: CompactionDomain): string {
  const clip = (s: string, n: number) => (s.length > n ? `${s.slice(0, n - 1).trimEnd()}…` : s);
  const requestText = domain.requestText ?? requestTextOf;
  return turns
    .map((t) => {
      const messages = t.entries.map((e) => e.message);
      const request = messages.find((m) => domain.isRequest(m));
      const replies = messages
        .filter((m) => m.role === 'assistant')
        .map(textOf)
        .filter(Boolean);
      const failed = new Set(
        messages.flatMap((m) =>
          m.role === 'user'
            ? m.parts.flatMap((p) => (p.kind === 'tool-result' && p.status === 'failed' ? [p.requestId] : []))
            : [],
        ),
      );
      const calls = messages.flatMap((m) =>
        m.role === 'assistant'
          ? m.parts.flatMap((p) =>
              p.kind === 'tool-call' ? [`${p.tool.id}${failed.has(p.requestId) ? ' ✗' : ''}`] : [],
            )
          : [],
      );
      return [
        `- Asked: ${clip(request ? requestText(request) : '', 400)}`,
        calls.length ? `  Did: ${calls.join(', ')}` : '',
        replies.length ? `  Replied: ${clip(replies.at(-1) ?? '', 600)}` : '',
      ]
        .filter(Boolean)
        .join('\n');
    })
    .join('\n');
}

// ------------------------------------------------------------------ the new head

export interface CompactionPlan {
  /** The new head: kept entries by id, new messages in place, in order. */
  edit: HeadEdit;
  /** For the trace and tests. */
  report: CompactionReport;
}

export interface CompactionReport {
  fallback: boolean;
  wishesKept: number;
  wishesDropped: string[];
  toolsKept: string[];
  skillsKept: string[];
  unknownRetention: string[];
  turnsVerbatim: number;
  turnsCondensed: number;
  /** The carry-over call's tokens as the provider reported them (reasoning is part of output). */
  call?: { input: number; cached: number; output: number; reasoning: number };
  /** Why the carry-over call gave no answer (the plan then fell back to a digest). */
  callFailure?: string;
  /** The carried-over context's length, in characters. */
  carriedChars: number;
  /** The new context's size, estimated from the old one's size (tokens). */
  estimatedTokens?: number;
  /** It stays above the target: accepted anyway. */
  overTarget?: boolean;
  /**
   * What no compaction can condense (the guidance, the kept tools, the current request and its last step) is above
   * the target by itself: a settings matter, so no size note asks the model to condense harder.
   */
  targetOutOfReach?: boolean;
  /** That part's estimated size, when the plan is over its target. */
  fixedTokens?: number;
}

export interface PlanInput {
  entries: readonly HeadEntry[];
  domain: CompactionDomain;
  /** The last turn is the one being run (mid-turn, or at its end); at a turn's start, none is. */
  current?: boolean;
  /** The compaction call's answer, or undefined when it failed. */
  answer: string | undefined;
  /** Skills and tools to keep regardless (the domain's pins). */
  pins: { skills: readonly string[]; tools: readonly string[] };
  settings: Pick<CompactionSettings, 'recentTurnsTokens' | 'maxSummaryTokens'>;
  /** The memory as it now stands (the domain's): the new head starts from it. */
  memory?: string;
  /**
   * The plan without it came out above the target while what can be condensed is what is over: ask the next
   * compaction to condense harder, and say how far (`fixed` cannot be condensed).
   */
  sizeNote?: { estimated: number; target: number; fixed: number };
}

const textMessage = (role: 'user' | 'system', text: string): Message => ({ role, parts: [{ kind: 'text', text }] });

export function planCompaction(input: PlanInput): CompactionPlan {
  const { entries, domain } = input;
  const turns = turnsOf(entries, domain);
  const requestText = domain.requestText ?? requestTextOf;
  // The turn being run: the last (at a turn's start its request is not in the head yet: none).
  const current = input.current === false ? undefined : turns.at(-1);
  const past = turns.filter((t) => t !== current);
  const midTurn = (current?.entries ?? []).some((e) => e.message.role === 'assistant');

  // Recent completed turns stay verbatim, newest first, within their budget.
  const recent: Turn[] = [];
  let used = 0;
  for (const t of [...past].reverse()) {
    const n = t.entries.reduce((a, e) => a + sizeOf(e.message), 0);
    if (used + n > input.settings.recentTurnsTokens) break;
    recent.unshift(t);
    used += n;
  }
  const condensed = past.filter((t) => !recent.includes(t));

  // The current turn: whole at a turn's start; mid-turn, its request and its last step (the call and its results
  // together, so the head stays projectable).
  let currentKept: HeadEntry[] = current?.entries ?? [];
  if (current && midTurn) {
    const request = current.entries.find((e) => domain.isRequest(e.message));
    const at = current.entries.findLastIndex((e) => e.message.role === 'assistant');
    const step = at >= 0 ? current.entries.slice(at, at + 2).filter((e, i) => i === 0 || isResults(e.message)) : [];
    currentKept = [request, ...step].filter((e): e is HeadEntry => e !== undefined);
  }
  const verbatim = [...recent.flatMap((t) => t.entries), ...currentKept];

  // What the person said, and wishes carried over before: what a wish may quote.
  const messages = entries.map((e) => e.message);
  const earlierWishes = messages.flatMap(carriedWishes);
  const sources = [
    ...messages.filter((m) => domain.isRequest(m)).map(requestText),
    ...earlierWishes,
    ...(domain.wishSources?.() ?? []),
  ];
  const parsed = input.answer ? parseCarryOver(input.answer) : undefined;
  const fallback = !parsed;
  const wishes = parsed ? checkWishes(parsed.wishes, sources) : { kept: earlierWishes, dropped: [] };

  // Tools: the latest definition of each loaded tool in the head, and of each tool only declared.
  const definitions = new Map<string, ToolDefinition>();
  const declarations = new Map<string, ToolDefinition>();
  for (const m of messages)
    for (const p of m.parts)
      if (p.kind === 'tool-definition')
        (p.deferred ? declarations : definitions).set(p.definition.tool.id, p.definition);
  // A skill is loaded when its guidance is in the head, whatever its tools.
  const loadedSkills = domain.groupsIn(messages);
  const wantedSkills = new Set([
    ...input.pins.skills,
    ...(parsed?.retention.skills ?? []).filter((s) => loadedSkills.includes(s)),
  ]);
  const unknownRetention = [
    ...(parsed?.retention.skills ?? []).filter((s) => !loadedSkills.includes(s)),
    ...(parsed?.retention.tools ?? []).filter((t) => !definitions.has(t)),
  ];
  const wanted = domain.groups.filter((g) => wantedSkills.has(g.id));
  const keepTools = new Set<string>([
    ...domain.coreToolIds,
    ...input.pins.tools,
    ...(parsed?.retention.tools ?? []),
    ...wanted.flatMap((g) => g.tools ?? []),
    ...verbatim.flatMap((e) => toolCallsOf(e.message)),
  ]);
  // What the current turn disclosed stays with it.
  if (current)
    for (const e of entries.slice(current.start))
      for (const p of e.message.parts)
        if (p.kind === 'tool-definition' && !p.deferred) keepTools.add(p.definition.tool.id);
  const keptLoaded = [...definitions.values()].filter((d) => keepTools.has(d.tool.id));
  // A kept tool stays loaded, now up front. Any other declared tool stays declared, to be loaded again where it is
  // disclosed; a loaded one never declared is dropped.
  const keptDeclared = [...declarations.values()].filter((d) => !keptLoaded.some((k) => k.tool.id === d.tool.id));

  // The carried-over text.
  const skillsDropped = loadedSkills.filter((s) => !wantedSkills.has(s));
  const block = (title: string, body: string) => (body.trim() ? `${title}:\n${body.trim()}` : '');
  const carried = [
    CARRY_OVER_HEADER,
    block('Wishes', wishes.kept.join('\n')),
    ...(parsed
      ? [
          block('Decisions', parsed.decisions),
          block('Open', parsed.open),
          block('Notes', parsed.notes),
          block('This turn so far', parsed.thisTurn),
        ]
      : [
          block('Carried over before', messages.flatMap(carriedBody).join('\n\n')),
          block('Earlier turns', digestOf(condensed, domain)),
        ]),
    skillsDropped.length
      ? `Skills loaded earlier and dropped here (disclose again when needed): ${skillsDropped.join(', ')}.`
      : '',
    input.sizeNote ? sizeNoteText(input.sizeNote) : '',
  ]
    .filter(Boolean)
    .join('\n\n');

  const notes: Message[] = [textMessage('user', carried)];
  // The memory as it now stands: the events of the old context fold into it.
  if (input.memory) notes.push(textMessage('user', input.memory));
  // Retained skills come back at their current guidance, as long as they fit the summary budget.
  let budget = input.settings.maxSummaryTokens - Math.ceil(carried.length / 3);
  const skillsKept: string[] = [];
  for (const group of wanted) {
    const entry = domain.groupEntry(group);
    const tokens = sizeOf(entry) + 100;
    if (tokens > budget) continue;
    budget -= tokens;
    skillsKept.push(group.id);
    notes.push(entry);
  }

  // One system message: the latest guidance (an update after a release replaces the earlier one), its text only;
  // the kept tools after it, in one message.
  const guidance = messages.filter((m) => m.role === 'system' && textOf(m).length > 0).at(-1);
  const tools = [
    ...keptLoaded.map((definition) => ({ kind: 'tool-definition' as const, definition, deferred: false })),
    ...keptDeclared.map((definition) => ({ kind: 'tool-definition' as const, definition, deferred: true })),
  ];
  const edit: HeadEdit = {
    entries: [
      ...(guidance ? [textMessage('system', textOf(guidance))] : []),
      ...(tools.length ? [{ role: 'system' as const, parts: tools }] : []),
      ...notes,
      ...verbatim.map((e) => e.id),
    ],
  };
  return {
    edit,
    report: {
      fallback,
      wishesKept: wishes.kept.length,
      wishesDropped: wishes.dropped,
      toolsKept: keptLoaded.map((d) => d.tool.id),
      skillsKept,
      unknownRetention,
      turnsVerbatim: recent.length,
      turnsCondensed: condensed.length,
      carriedChars: carried.length,
    },
  };
}

/** The note a carry-over gets when it, not the part no compaction can condense, is over the target. */
function sizeNoteText(n: { estimated: number; target: number; fixed: number }): string {
  const k = (x: number) => x.toLocaleString('en-US');
  return `${SIZE_NOTE} this condensed conversation is about ${k(n.estimated)} tokens, above its target of ${k(n.target)}; about ${k(n.fixed)} of it cannot be condensed (guidance, tools, the current request), so what is carried over should be at most about ${k(Math.max(0, n.target - n.fixed))}. At the next compaction, condense harder.`;
}

/**
 * A new head's size in tokens, estimated from the old head's size: the same tokens per character. Deferred tool
 * declarations do not count (a provider shows at most their names).
 */
export function estimatedTokens(edit: HeadEdit, old: readonly HeadEntry[], oldTokens: number): number {
  const chars = (m: Message) =>
    JSON.stringify(
      m.role === 'system' ? { ...m, parts: m.parts.filter((p) => p.kind !== 'tool-definition' || !p.deferred) } : m,
    ).length;
  const byId = new Map(old.map((e) => [e.id, e.message]));
  const oldChars = old.reduce((a, e) => a + chars(e.message), 0);
  const newChars = edit.entries.reduce((a, e) => {
    const m = typeof e === 'string' ? byId.get(e) : e;
    return a + (m ? chars(m) : 0);
  }, 0);
  return oldChars === 0 ? 0 : Math.ceil((oldTokens * newChars) / oldChars);
}

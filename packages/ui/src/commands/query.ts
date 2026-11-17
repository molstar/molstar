/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What typed text matches (commands §4.3): recent commands, fuzzy matches with presets, a direct alias with its tail.

import { fuzzyScore } from './match.js';
import type { CommandDefinition, CommandPreset } from './types.js';

export interface CommandMatch<Ctx> {
  command: CommandDefinition<Ctx>;
  /** Why it cannot run now, if it cannot. */
  blocked?: string;
  /** The way in it was found by, if one: opened with `open(id, tail, preset.id)`. */
  preset?: CommandPreset;
}

export interface CommandQuery<Ctx> {
  /** Recently used, most recent first, when nothing is typed. */
  recent: CommandMatch<Ctx>[];
  /** Best first when something is typed; every other command, in order, when nothing is. */
  matches: CommandMatch<Ctx>[];
  /** The text starts with a command's alias and goes on: that command, and the rest as its tail. */
  direct?: CommandMatch<Ctx> & { tail: string };
}

/** What a query reads besides the commands. */
export interface QuerySource<Ctx> {
  recents: readonly string[];
  byId: ReadonlyMap<string, CommandDefinition<Ctx>>;
  match(command: CommandDefinition<Ctx>, preset?: CommandPreset): CommandMatch<Ctx>;
  presets(command: CommandDefinition<Ctx>): readonly CommandPreset[];
}

export function aliases(c: { alias?: string | readonly string[] }): readonly string[] {
  return c.alias === undefined ? [] : typeof c.alias === 'string' ? [c.alias] : c.alias;
}

export function queryCommands<Ctx>(
  all: readonly CommandDefinition<Ctx>[],
  text: string,
  source: QuerySource<Ctx>,
): CommandQuery<Ctx> {
  const q = text.trim();
  const lower = text.toLowerCase().trimStart().replace(/\s+/g, ' ');
  const direct = directAlias(all, lower, source);
  if (direct) {
    const tail = collapse(text)
      .slice(direct.alias.length + 1)
      .trim();
    const match = source.match(direct.command);
    // A blocked one keeps the other ways in sight: what the whole text matches.
    const others = match.blocked ? ranked(all, q, source).filter((m) => m.command !== direct.command) : [];
    return { recent: [], matches: [match, ...others], direct: { ...match, tail } };
  }
  if (!q) {
    const recent: CommandMatch<Ctx>[] = [];
    const seen = new Set<string>();
    for (const id of source.recents) {
      const c = source.byId.get(id);
      if (c && !seen.has(id)) {
        seen.add(id);
        recent.push(source.match(c));
      }
    }
    return { recent, matches: all.filter((c) => !seen.has(c.id)).map((c) => source.match(c)) };
  }
  return { recent: [], matches: ranked(all, q, source) };
}

/** The text with its runs of space as one, its start trimmed (as the alias is compared). */
function collapse(text: string): string {
  return text.trimStart().replace(/\s+/g, ' ');
}

/**
 * The command whose alias the text starts with and goes on past: the longest alias that fits, and among commands
 * sharing it one that can run. None while the text is still the start of a title or an alias ("focus on" of "Focus
 * on…", "add scene" of "add scene"): that is a query, not a command with arguments.
 */
function directAlias<Ctx>(
  all: readonly CommandDefinition<Ctx>[],
  lower: string,
  source: QuerySource<Ctx>,
): { command: CommandDefinition<Ctx>; alias: string } | undefined {
  const typed = lower.trimEnd();
  for (const c of all) {
    if (c.title.toLowerCase().replace(/…$/, '').startsWith(typed)) return undefined;
    for (const a of aliases(c)) if (a.toLowerCase().startsWith(typed)) return undefined;
  }
  let best: { command: CommandDefinition<Ctx>; alias: string } | undefined;
  let bestRuns = false;
  for (const c of all) {
    if (!c.params) continue;
    for (const alias of aliases(c)) {
      const a = alias.toLowerCase();
      if (!lower.startsWith(`${a} `) || typed.length <= a.length) continue;
      const longer = !best || a.length > best.alias.length;
      if (!longer && a.length < (best?.alias.length ?? 0)) continue;
      const runs = !source.match(c).blocked;
      if (longer || (!bestRuns && runs)) {
        best = { command: c, alias: a };
        bestRuns = runs;
      }
    }
  }
  return best;
}

/**
 * Commands whose title or an alias fuzzy-matches the text, best first, and their presets that the text says more about
 * than the command alone (by the preset's title, its words, or both titles).
 */
function ranked<Ctx>(all: readonly CommandDefinition<Ctx>[], q: string, source: QuerySource<Ctx>): CommandMatch<Ctx>[] {
  const score = (text: string, less = 0) => (fuzzyScore(text, q) ?? -1) - less;
  const hits: { c: CommandDefinition<Ctx>; p?: CommandPreset; s: number }[] = [];
  for (const c of all) {
    let s = score(c.title);
    for (const a of aliases(c)) s = Math.max(s, score(a, 0.5));
    if (s >= 0) hits.push({ c, s });
    for (const p of source.presets(c)) {
      let own = score(p.title);
      for (const k of p.keywords ?? []) own = Math.max(own, score(k, 0.5));
      const both = score(`${c.title} ${p.title}`);
      if (own >= 0 || (s < 0 && both >= 0)) hits.push({ c, p, s: Math.max(own, both) });
    }
  }
  // Better matches first by whole points (a title that starts with the text, a word that does, then looser ones); equal
  // ones keep their category together, in the order categories were registered, then go by title.
  const categories = new Map<string, number>();
  for (const c of all) if (!categories.has(c.category)) categories.set(c.category, categories.size);
  const title = (x: (typeof hits)[number]) => (x.p ? `${x.c.title} ${x.p.title}` : x.c.title);
  hits.sort(
    (a, b) =>
      Math.floor(b.s) - Math.floor(a.s) ||
      (categories.get(a.c.category) ?? 0) - (categories.get(b.c.category) ?? 0) ||
      title(a).localeCompare(title(b)),
  );
  return hits.map((x) => source.match(x.c, x.p));
}

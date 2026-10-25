/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CommandDefinition, CommandMenuItem, CommandMenuPrefix } from '@molstar/prism-ui/core';
import type { ViewerApi } from '@molstar/prism-viewer/api/area';
import { operationCommand } from '@molstar/prism-viewer/commands/operations';
import type { CommandPrefixContext, PrismStore } from '@molstar/prism-viewer/spec';
import { applyCompletion, complete } from '../language/completion.js';
import { describeCommand } from '../language/describe.js';
import { formatDiagnostic } from '../language/diagnostics.js';
import { type Dialect, DialectNames, Dialects } from '../language/ir.js';
import { hasErrors } from '../language/script.js';

// The language in Prism's command menu (_spec §9): `pymol ` and `vmd ` prefixes, and `!` for the person's default
// dialect. A prefix's rows are the line's parse: what it will do (Enter runs it as one step) or what does not read,
// then completions that edit the line. `undo` and `redo` alone run the viewer's own history.

export const DialectPreference = 'command-language.dialect';

/** Without preferences (an embedded viewer), the choice lasts for the session: per viewer (its area), not the page's. */
const sessionDialects = new WeakMap<object, Dialect>();

/** The person's default dialect (`!`), from the preferences; none when unset. */
export function defaultDialect(preferences: PrismStore | undefined, viewer: object): Dialect | undefined {
  try {
    const v = preferences?.get(DialectPreference);
    if (Dialects.includes(v as Dialect)) return v as Dialect;
  } catch {
    // An unreadable store: the session's choice, if one was made.
  }
  return sessionDialects.get(viewer);
}

function setDialect(preferences: PrismStore | undefined, viewer: object, dialect: Dialect): void {
  if (!preferences) {
    sessionDialects.set(viewer, dialect);
    return;
  }
  try {
    preferences.set(DialectPreference, dialect);
  } catch {
    sessionDialects.set(viewer, dialect);
  }
}

/** Rows show this many completions at most. */
const MaxCompletions = 8;

type Post = CommandPrefixContext['post'];

/**
 * Runs a line from the menu (which has closed). The viewer says a person's run on the status line, what it did and
 * printed (`did`) or why it failed; its warnings it does not: they are said here.
 */
async function runLine(api: ViewerApi, dialect: Dialect, text: string, post: Post): Promise<void> {
  const r = await api.script?.run({ dialect, text });
  if (r?.ok) for (const w of r.warnings ?? []) post(w, 'warning');
}

/** The rows for a line typed after a prefix. */
export function prefixRows(
  api: ViewerApi,
  prefix: string,
  dialect: Dialect,
  text: string,
  post: Post = () => {},
): CommandMenuItem[] {
  const area = api.script;
  if (!area) return [];
  const name = DialectNames[dialect];
  if (!text.trim()) {
    return [
      {
        id: `${dialect}:help`,
        label: `Type a ${name} command`,
        detail:
          dialect === 'pymol'
            ? 'show sticks, resn HEM · color red, chain A · help'
            : 'mol addrep top · mol color Chain · help',
        icon: 'terminal',
        run: () => true,
      },
    ];
  }
  const env = area.env();
  const compiled = area.preview(text, dialect);
  const errors = compiled.diagnostics.filter((d) => d.severity === 'error');
  const c = complete(text, dialect, env);
  const completions: CommandMenuItem[] = c.items
    .filter((i) => i.insert !== text.slice(c.from, c.to))
    .slice(0, MaxCompletions)
    .map((i) => ({
      id: `${dialect}:complete:${i.label}`,
      label: i.label,
      ...(i.detail ? { detail: i.detail } : {}),
      run: () => ({ text: `${prefix}${applyCompletion(text, c, i)}` }),
    }));
  if (errors.length) {
    const first: CommandMenuItem = {
      id: `${dialect}:error`,
      label: formatDiagnostic(text, { ...errors[0], hint: undefined } as (typeof errors)[number]),
      detail: errors[0].hint ?? '',
      icon: 'circle-alert',
      // A line that does not read stays in the menu to be fixed.
      run: () => true,
    };
    // While the line is being typed (the error is about its end), the completions come first.
    const atEnd = errors[0].at.to >= text.trimEnd().length;
    return atEnd && completions.length ? [...completions, first] : [first, ...completions];
  }
  const lone = compiled.commands.length === 1 ? compiled.commands[0] : undefined;
  const label = compiled.commands.map(describeCommand).filter(Boolean).join('; ') || 'Nothing to run';
  const warnings = compiled.diagnostics.filter((d) => d.severity !== 'error').map((d) => d.message);
  const runRow: CommandMenuItem = {
    id: `${dialect}:run`,
    label,
    detail: warnings.join(' ') || `${name} · Enter runs it`,
    icon: 'terminal',
    run: () => {
      if (lone?.kind === 'history') void (lone.action === 'undo' ? api.history.undo() : api.history.redo());
      else if (!hasErrors(compiled.diagnostics)) void runLine(api, dialect, text, post);
    },
  };
  return [runRow, ...completions];
}

/** The prefixes: `pymol `, `vmd `, and `!` for the default dialect. */
export function commandPrefixes({ api, preferences, post }: CommandPrefixContext): CommandMenuPrefix[] {
  const word = (dialect: Dialect): CommandMenuPrefix => ({
    char: `${dialect} `,
    title: `${DialectNames[dialect]} command`,
    placeholder: `A ${DialectNames[dialect]} command`,
    items: (text) => prefixRows(api, `${dialect} `, dialect, text, post),
  });
  return [
    word('pymol'),
    word('vmd'),
    {
      char: '!',
      title: 'Command in the default dialect',
      placeholder: 'A command in your default dialect (PyMOL or VMD)',
      items: (text) => {
        const viewer = api.script;
        if (!viewer) return [];
        const dialect = defaultDialect(preferences, viewer);
        // The choice is the person's preference, made here (no setting elsewhere to look for).
        const choose = (d: Dialect, label: string): CommandMenuItem => ({
          id: `command-language:dialect:${d}`,
          label,
          detail: preferences ? 'Kept with your preferences' : 'For this session',
          icon: 'terminal',
          run: () => {
            setDialect(preferences, viewer, d);
            return { text: '!' };
          },
        });
        if (!dialect) {
          return Dialects.map((d) => choose(d, `Read ! as ${DialectNames[d]}`));
        }
        const rows = prefixRows(api, '!', dialect, text, post);
        if (text.trim()) return rows;
        const other = Dialects.find((d) => d !== dialect)!;
        return [...rows, choose(other, `Read ! as ${DialectNames[other]} instead of ${DialectNames[dialect]}`)];
      },
    },
  ];
}

/** Commands over the area's operations. */
export function commandLanguageCommands(api: ViewerApi): CommandDefinition<any>[] {
  const area = api.script;
  if (!area) return [];
  return [
    operationCommand(api, area.run, {
      category: 'Commands',
      icon: 'terminal',
      alias: ['pymol script', 'vmd script', 'run commands'],
      action: 'Run',
      order: ['dialect', 'text'],
      values: { files: [] },
    }),
    operationCommand(api, area.names, { category: 'Commands', icon: 'terminal', alias: ['named selections'] }),
  ];
}

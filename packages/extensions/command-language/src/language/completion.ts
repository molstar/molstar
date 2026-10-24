/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { all as PymolSymbols } from '@molstar/query-language/transpilers/pymol/symbols';
import { all as VmdSymbols } from '@molstar/query-language/transpilers/vmd/symbols';
import { PymolColorNames, pymolColor, VmdColorIds, vmdColor } from './colors.js';
import type { ScriptEnv } from './env.js';
import type { Dialect } from './ir.js';
import { PymolReps } from './pymol/args.js';
import { PymolCommands, PymolCommandsByName } from './pymol/commands.js';
import { PymolSettings } from './pymol/settings.js';
import { splitPymol, unquote } from './pymol/tokenizer.js';
import type { ArgKind, CommandHelp, ParamSpec } from './table.js';
import { VmdCommands, VmdGroups, VmdCommandsByName } from './vmd/commands.js';
import { VmdColorMethodNames, VmdStyleNames } from './vmd/styles.js';
import { literal, splitTcl } from './vmd/tcl.js';

// Completion at the cursor (_spec §8), from the tables: the command word, then by the argument's kind (rep names,
// colors, settings, enums, files, molecules), and in a selection the dialect's keywords, properties and operators
// (query-language's lists), the session's names and objects, and the structures' chains and ligands.

export interface CompletionItem {
  readonly label: string;
  /** What replaces the text from `from` to `to`. */
  readonly insert: string;
  readonly detail?: string;
  readonly kind: 'command' | ArgKind | 'keyword';
  /** A color's RGB (`0xRRGGBB`), for a swatch. */
  readonly rgb?: number;
}

/** Where the cursor is in a command: its help and which of its arguments is being typed (signature help). */
export interface Signature {
  readonly name: string;
  readonly usage: string;
  readonly summary: string;
  readonly params: readonly ParamSpec[];
  /** The argument being typed, an index into `params`; -1: none (past them, or the name). */
  readonly active: number;
}

export interface Completion {
  readonly from: number;
  readonly to: number;
  readonly items: readonly CompletionItem[];
  /** The command around the cursor, once its name is typed. */
  readonly signature?: Signature;
}

const signatureOf = (c: CommandHelp & { readonly params: readonly ParamSpec[] }, active: number): Signature => ({
  name: c.name,
  usage: c.usage,
  summary: c.summary,
  params: c.params,
  active,
});

const Limit = 40;

function match(word: string, candidates: readonly CompletionItem[]): CompletionItem[] {
  const w = word.toLowerCase();
  const starts = candidates.filter((c) => c.label.toLowerCase().startsWith(w));
  const inside = w
    ? candidates.filter((c) => !c.label.toLowerCase().startsWith(w) && c.label.toLowerCase().includes(w))
    : [];
  const seen = new Set<string>();
  return [...starts, ...inside].filter((c) => !seen.has(c.label) && seen.add(c.label)).slice(0, Limit);
}

const items = (labels: Iterable<string>, kind: CompletionItem['kind'], detail?: string): CompletionItem[] =>
  [...labels].map((label) => ({ label, insert: label, kind, ...(detail ? { detail } : {}) }));

/** In a selection: what fits after the word before (`chain` → chains, `resn` → ligands), else the language's words. */
function selectionItems(dialect: Dialect, before: string, env: ScriptEnv): CompletionItem[] {
  const prev = before.trim().split(/\s+/).pop()?.toLowerCase() ?? '';
  const chains = [...new Set(env.objects.flatMap((o) => o.chains ?? []))];
  const ligands = [...new Set(env.objects.flatMap((o) => o.ligands ?? []))];
  if (prev === 'chain' || prev === 'c.' || prev === 'segi') return items(chains, 'text', 'chain');
  if (prev === 'resn' || prev === 'r.' || prev === 'resname') return items(ligands, 'text', 'ligand');
  const symbols = dialect === 'pymol' ? PymolSymbols : VmdSymbols;
  return [
    ...(dialect === 'pymol' ? items(env.names.keys(), 'name', 'selection') : []),
    ...(dialect === 'pymol'
      ? items(
          env.objects.map((o) => o.name),
          'object',
          'object',
        )
      : []),
    ...items(symbols.Keywords, 'keyword', 'keyword'),
    ...items(symbols.Properties, 'keyword', 'property'),
    ...items(symbols.Operators, 'keyword', 'operator'),
  ];
}

function kindItems(dialect: Dialect, param: ParamSpec, env: ScriptEnv, before: string): CompletionItem[] {
  switch (param.kind) {
    case 'selection':
      return selectionItems(dialect, before, env);
    case 'rep':
      return dialect === 'pymol'
        ? items(Object.keys(PymolReps), 'rep', 'representation')
        : items(VmdStyleNames, 'rep', 'drawing style');
    case 'color':
      if (dialect === 'pymol') {
        return [
          { label: 'atomic', insert: 'atomic', kind: 'color', detail: 'by element' },
          ...PymolColorNames.map(
            (label): CompletionItem => ({ label, insert: label, kind: 'color', rgb: pymolColor(label)! }),
          ),
        ];
      }
      return [
        ...(param.values ?? VmdColorMethodNames).map(
          (label): CompletionItem => ({ label, insert: label, kind: 'color', detail: 'coloring method' }),
        ),
        ...(param.values
          ? []
          : Object.entries(VmdColorIds).map(
              ([label, id]): CompletionItem => ({
                label,
                insert: label,
                kind: 'color',
                detail: `ColorID ${id}`,
                rgb: vmdColor(label)!,
              }),
            )),
      ];
    case 'setting':
      return Object.entries(PymolSettings).map(([label, s]) => ({
        label,
        insert: label,
        kind: 'setting',
        detail: s.summary,
      }));
    case 'object':
      return [
        ...items(
          env.objects.map((o) => o.name),
          'object',
          'object',
        ),
        ...items(env.names.keys(), 'name', 'selection'),
      ];
    case 'molecule':
      return [
        { label: 'top', insert: 'top', kind: 'molecule', detail: env.top ?? '' },
        ...env.objects.map(
          (o): CompletionItem => ({
            label: String(o.molecule),
            insert: String(o.molecule),
            kind: 'molecule',
            detail: o.name,
          }),
        ),
        { label: 'all', insert: 'all', kind: 'molecule' },
      ];
    case 'file':
      return items(env.files, 'file');
    default:
      return param.values ? items(param.values, param.kind, param.name) : [];
  }
}

function pymolCompletion(text: string, env: ScriptEnv): Completion {
  const { calls } = splitPymol(text);
  const call = calls[calls.length - 1];
  // At a command's start: the commands.
  if (!call || (call.args.length === 0 && !call.open && call.at.to === text.length)) {
    const word = call?.name ?? '';
    const from = call ? call.nameAt.from : text.length;
    const commands = PymolCommands.map(
      (c): CompletionItem => ({ label: c.name, insert: c.name, kind: 'command', detail: c.summary }),
    );
    const named = call && PymolCommandsByName.get(call.name);
    return {
      from,
      to: text.length,
      items: match(word, commands),
      ...(named && named.name === word ? { signature: signatureOf(named, -1) } : {}),
    };
  }
  const command = PymolCommandsByName.get(call.name) ?? PymolCommandsByName.get(call.name.toLowerCase());
  if (!command) return { from: text.length, to: text.length, items: [] };
  const last = call.args[call.args.length - 1];
  // A new argument after a comma (or the first after the name); else the last one goes on, a space inside it included.
  const fresh = !last || /,\s*$/.test(text);
  const positional = call.args.filter((a) => !a.key).length;
  const param =
    !fresh && last.key
      ? command.params.find((p) => p.name === last.key)
      : command.params[fresh ? positional : positional - 1];
  const signature = signatureOf(command, param ? command.params.indexOf(param) : -1);
  if (!param) return { from: text.length, to: text.length, items: [], signature };
  const current = fresh ? '' : text.slice(last.at.from);
  // In a selection, the last word only; elsewhere the whole argument.
  const word = param.kind === 'selection' ? (current.match(/[\w.'+-]*$/)?.[0] ?? '') : unquote(current);
  const before = param.kind === 'selection' ? current.slice(0, current.length - word.length) : '';
  // A whole argument with more after it (`show sticks`, then the selection): its comma comes with it.
  const index = command.params.indexOf(param);
  const more = param.kind !== 'selection' && param.kind !== 'text' && index < command.params.length - 1;
  const found = match(word, kindItems('pymol', param, env, before));
  return {
    from: text.length - word.length,
    to: text.length,
    items: more ? found.map((i) => ({ ...i, insert: `${i.insert}, ` })) : found,
    signature,
  };
}

function vmdCompletion(text: string, env: ScriptEnv): Completion {
  const { commands } = splitTcl(text);
  const cmd = commands[commands.length - 1];
  const atEnd = !!cmd && cmd.at.to === text.trimEnd().length && !/[;\n]\s*$/.test(text);
  const words = atEnd ? cmd.words.map((w) => literal(w) ?? '') : [];
  const open = atEnd && cmd.open;
  const firsts = [...new Set(VmdCommands.map((c) => c.name.split(' ')[0]))];
  // Typing the first word.
  if (!words.length || (words.length === 1 && !open)) {
    const word = words[0] ?? '';
    const list = firsts.map((f): CompletionItem => {
      const own = VmdCommandsByName.get(f);
      return { label: f, insert: f, kind: 'command', ...(own ? { detail: own.summary } : {}) };
    });
    return { from: text.length - word.length, to: text.length, items: match(word, list) };
  }
  const group = words[0].toLowerCase();
  // Typing a subcommand.
  if (VmdGroups.has(group) && (words.length === 1 || (words.length === 2 && !open))) {
    const word = words.length === 2 ? words[1] : '';
    const subs = VmdCommands.filter((c) => c.name.startsWith(`${group} `)).map(
      (c): CompletionItem => ({
        label: c.name.split(' ')[1],
        insert: c.name.split(' ')[1],
        kind: 'command',
        detail: c.summary,
      }),
    );
    return { from: text.length - word.length, to: text.length, items: match(word, subs) };
  }
  const key = VmdGroups.has(group) ? `${group} ${words[1].toLowerCase()}` : group;
  const command = VmdCommandsByName.get(key);
  if (!command) return { from: text.length, to: text.length, items: [] };
  const args = words.slice(VmdGroups.has(group) ? 2 : 1);
  const index = open ? args.length : args.length - 1;
  const param =
    command.params[index] ??
    (command.params[command.params.length - 1]?.kind === 'selection'
      ? command.params[command.params.length - 1]
      : undefined);
  const signature = signatureOf(command, param ? command.params.indexOf(param) : -1);
  if (!param) return { from: text.length, to: text.length, items: [], signature };
  const current = open ? '' : (args[args.length - 1] ?? '');
  const word = param.kind === 'selection' ? (current.match(/[\w.'+-]*$/)?.[0] ?? '') : current;
  const before = param.kind === 'selection' ? current.slice(0, current.length - word.length) : '';
  return {
    from: text.length - word.length,
    to: text.length,
    items: match(word, kindItems('vmd', param, env, before)),
    signature,
  };
}

/** What may come at the end of `text` (a completion edits the text before the cursor). */
export function complete(text: string, dialect: Dialect, env: ScriptEnv): Completion {
  return dialect === 'pymol' ? pymolCompletion(text, env) : vmdCompletion(text, env);
}

/** The text with a completion applied (a space after a word that ends an argument). */
export function applyCompletion(text: string, c: Completion, item: CompletionItem): string {
  const space = item.kind === 'command' || item.kind === 'molecule' ? ' ' : '';
  return text.slice(0, c.from) + item.insert + space + text.slice(c.to);
}

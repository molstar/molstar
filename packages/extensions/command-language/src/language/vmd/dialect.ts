/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { fail, nearest, type Span } from '../diagnostics.js';
import type { DialectImpl, Unit } from '../dialect.js';
import type { Compiled, Compiler } from '../table.js';
import { vmdSelection } from '../selection.js';
import { moleculeOf, VmdCommands, VmdCommandsByName, VmdGroups, VmdWords } from './commands.js';
import { literal, splitTcl, type TclCommand, type TclWord } from './tcl.js';
import { VmdUnsupportedByName } from './unsupported.js';

// VMD as a dialect (_spec §3.2): Tcl words, a closed subset of Tcl evaluated while compiling. `set` and `unset` keep
// variables in the session; `$name` reads one; `[…]` runs `atomselect` or `molinfo top` for its value; an atomselect
// handle's methods (`$sel num`, `$sel text`, `$sel delete`) are commands of their own. The rest is VMD's commands.

interface VmdUnit extends Unit {
  readonly command: TclCommand;
  /** The whole submission: the words' places (and a substitution's script) are offsets into it. */
  readonly text: string;
}

/** A word's value: its text with variables and command substitutions in. */
function wordValue(w: TclWord, text: string, c: Compiler): string {
  let out = '';
  for (const part of w.parts) {
    if (part.kind === 'text') out += part.text;
    else if (part.kind === 'var') {
      const v = c.state.session.vars.get(part.name);
      if (v === undefined) {
        fail('no-name', `There is no variable '${part.name}'.`, part.at, 'Set it first: set name value.');
      }
      out += v;
    } else {
      const { commands, diagnostics } = splitTcl(text, part.from, part.at.to - 1);
      if (diagnostics.length) fail('syntax', diagnostics[0].message, diagnostics[0].at);
      if (commands.length !== 1) {
        fail('unsupported', 'A command substitution runs one command.', part.at);
      }
      out += substitute(commands[0], text, c, part.at);
    }
  }
  return out;
}

/** The value of `[…]`: an atomselect handle, or the top molecule's id. */
function substitute(cmd: TclCommand, text: string, c: Compiler, at: Span): string {
  const values = cmd.words.map((w) => wordValue(w, text, c));
  const head = values[0]?.toLowerCase();
  if (head === 'atomselect')
    return atomselect(
      values.slice(1),
      cmd.words.slice(1).map((w) => w.at),
      c,
      at,
    );
  if (head === 'molinfo' && values[1] === 'top') {
    const top = c.state.molecule('top');
    if (!top) fail('no-name', 'There is no top molecule: nothing is loaded.', at);
    return String(top.molecule);
  }
  return fail(
    'unsupported',
    `[${values[0] ?? ''} …] is not supported: substitutions run atomselect and molinfo top.`,
    at,
    'Write the value itself.',
  );
}

/** `atomselect molid "selection"`: a handle (atomselect0, …) for the selection on that molecule. */
function atomselect(args: readonly string[], ats: readonly Span[], c: Compiler, at: Span): string {
  if (args[0] === 'macro')
    fail('unsupported', 'atomselect macro is not supported yet.', ats[0] ?? at, 'Write the selection out.');
  if (args.length < 2) fail('argument', 'atomselect needs a molecule and a selection.', at, 'atomselect top "protein"');
  const object = moleculeOf(c, args[0], ats[0] ?? at);
  let text = args[1];
  if (args.length > 2) {
    if (args[2] === 'frame') {
      fail(
        'unsupported-argument',
        'atomselect … frame is not supported: selections are on the frame shown.',
        ats[2] ?? at,
      );
    }
    text = args.slice(1).join(' ');
  }
  const sel = vmdSelection(text, object.name, ats[1] ?? at);
  const name = `atomselect${c.state.session.nextHandle++}`;
  c.state.session.handles.set(name, { object: object.name, text, sel });
  return name;
}

/** A handle's method: `$sel num`, `$sel text`, `$sel molid`, `$sel get …`, `$sel delete`. */
function handleMethod(handleName: string, w: VmdWords, c: Compiler): Compiled {
  const handle = c.state.session.handles.get(handleName)!;
  const method = (w.opt(0) ?? '').toLowerCase();
  switch (method) {
    case 'num':
      return [{ kind: 'read', read: { kind: 'count', sel: handle.sel } }];
    case 'text':
      return [{ kind: 'print', text: handle.text }];
    case 'molid':
      return [{ kind: 'print', text: String(c.state.object(handle.object)?.molecule ?? '') }];
    case 'get':
      return [
        { kind: 'read', read: { kind: 'values', sel: handle.sel, properties: w.rest(1).split(/\s+/).filter(Boolean) } },
      ];
    case 'delete':
      c.state.session.handles.delete(handleName);
      return [];
    case 'set':
    case 'move':
    case 'moveby':
    case 'moveto':
    case 'writepdb':
      return fail('unsupported', `$sel ${method} changes or writes coordinates, which is not supported.`, w.atOf(0));
    default:
      return fail(
        'argument',
        `An atom selection has no method '${method}'.`,
        w.atOf(0),
        'Methods: num, text, molid, get, delete.',
      );
  }
}

export const VmdDialect: DialectImpl = {
  id: 'vmd',
  commands: VmdCommands,
  split(text) {
    const { commands, diagnostics } = splitTcl(text);
    const units: VmdUnit[] = commands.map((command) => ({
      name: literal(command.words[0]) ?? '',
      at: command.at,
      source: command.source,
      command,
      text,
    }));
    return { units, diagnostics };
  },
  includeOf(unit) {
    const { command } = unit as VmdUnit;
    // `source file.tcl` and `play file.vmd`: another script given with this one (a literal name).
    const head = literal(command.words[0])?.toLowerCase();
    const file = command.words[1] ? literal(command.words[1]) : undefined;
    if ((head === 'source' || head === 'play') && file && command.words.length === 2)
      return { name: file, at: command.at };
    return undefined;
  },
  compile(unit, c) {
    const { command, text: full } = unit as VmdUnit;
    const first = command.words[0];
    const head = literal(first)?.toLowerCase();
    // `set name value` and `unset name`: the session's variables (the value substituted now).
    if (head === 'set' || head === 'unset') {
      const name = command.words[1] ? literal(command.words[1]) : undefined;
      if (!name) fail('argument', `${head} needs a variable name.`, command.at);
      if (head === 'unset') {
        c.state.session.vars.delete(name);
        return [];
      }
      if (command.words.length < 3) {
        const v = c.state.session.vars.get(name);
        if (v === undefined) fail('no-name', `There is no variable '${name}'.`, command.words[1].at);
        return [{ kind: 'print', text: v }];
      }
      c.state.session.vars.set(
        name,
        command.words
          .slice(2)
          .map((w) => wordValue(w, full, c))
          .join(' '),
      );
      return [];
    }
    if (head === 'help') {
      const topic = command.words
        .slice(1)
        .map((w) => literal(w) ?? '')
        .join(' ')
        .trim()
        .toLowerCase();
      return [{ kind: 'print', text: helpText(topic, command.at) }];
    }
    if (head === 'source' || head === 'play') {
      // Not an include (`includeOf`): the file's name is not written out, or there is more.
      fail(
        'unsupported',
        `${head} runs a script given with this one, by its name as written: ${head} file.tcl.`,
        command.at,
        'Open the scripts together.',
      );
    }
    if (head === 'atomselect') {
      const values = command.words.slice(1).map((w) => wordValue(w, full, c));
      return [
        {
          kind: 'print',
          text: atomselect(
            values,
            command.words.slice(1).map((w) => w.at),
            c,
            command.at,
          ),
        },
      ];
    }
    const values = command.words.map((w) => wordValue(w, full, c));
    const ats = command.words.map((w) => w.at);
    if (c.state.session.handles.has(values[0])) {
      return handleMethod(values[0], new VmdWords(values.slice(1), ats.slice(1), command.at, values[0]), c);
    }
    const name0 = values[0].toLowerCase();
    if (VmdGroups.has(name0) && values.length > 1) {
      const key = `${name0} ${values[1].toLowerCase()}`;
      const sub = VmdCommandsByName.get(key);
      if (sub) return sub.compile(new VmdWords(values.slice(2), ats.slice(2), command.at, sub.name), c);
      const unsupported = VmdUnsupportedByName.get(key);
      if (unsupported)
        fail('unsupported', `VMD’s ${key}: ${unsupported.why}`, { from: ats[0].from, to: ats[1].to }, unsupported.hint);
      const subs = VmdCommands.filter((x) => x.name.startsWith(`${name0} `)).map((x) => x.name.split(' ')[1]);
      const near = nearest(values[1], subs);
      return fail(
        'unknown-command',
        `VMD’s ${name0} has no subcommand '${values[1]}'.`,
        ats[1],
        near ? `Did you mean ${name0} ${near}?` : `Subcommands: ${subs.join(', ')}.`,
      );
    }
    const command0 = VmdCommandsByName.get(name0);
    if (command0) return command0.compile(new VmdWords(values.slice(1), ats.slice(1), command.at, command0.name), c);
    const unsupported = VmdUnsupportedByName.get(name0);
    if (unsupported) fail('unsupported', `${name0}: ${unsupported.why}`, ats[0], unsupported.hint);
    if (VmdGroups.has(name0)) {
      fail(
        'argument',
        `${name0} needs a subcommand.`,
        ats[0],
        `Subcommands: ${VmdCommands.filter((x) => x.name.startsWith(`${name0} `))
          .map((x) => x.name.split(' ')[1])
          .join(', ')}.`,
      );
    }
    const near = nearest(name0, [...new Set(VmdCommands.map((x) => x.name.split(' ')[0]))]);
    return fail(
      'unknown-command',
      `VMD has no command '${values[0]}'.`,
      ats[0],
      near ? `Did you mean ${near}?` : '`help` lists the commands.',
    );
  },
};

function helpText(topic: string, at: Span): string {
  if (!topic) {
    const groups = [...new Set(VmdCommands.map((x) => x.name.split(' ')[0]))];
    return `VMD commands: ${groups.join(', ')}, set, atomselect. \`help mol\` lists mol’s; \`help mol addrep\` says more.`;
  }
  const exact = VmdCommandsByName.get(topic);
  if (exact) return `${exact.usage}: ${exact.summary}${exact.examples.length ? ` Example: ${exact.examples[0]}` : ''}`;
  const subs = VmdCommands.filter((x) => x.name.startsWith(`${topic} `));
  if (subs.length) return `${topic}: ${subs.map((x) => x.name.split(' ')[1]).join(', ')}.`;
  const u = VmdUnsupportedByName.get(topic);
  if (u) return `${topic}: ${u.why}${u.hint ? ` ${u.hint}` : ''}`;
  return fail('unknown-command', `VMD has no command '${topic}'.`, at, '`help` lists the commands.');
}

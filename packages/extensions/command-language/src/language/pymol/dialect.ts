/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { fail, nearest } from '../diagnostics.js';
import type { DialectImpl, Unit } from '../dialect.js';
import { PymolCommands, PymolCommandsByName, compilePymolCall } from './commands.js';
import { splitPymol, type PymolCall, unquote } from './tokenizer.js';
import { PymolUnsupportedByName } from './unsupported.js';

// PyMOL as a dialect (_spec §3.1): its lines split into calls, each looked up in the table, else in the list of
// commands Mol* does not do, else unknown.

interface PymolUnit extends Unit {
  readonly call: PymolCall;
}

export const PymolDialect: DialectImpl = {
  id: 'pymol',
  commands: PymolCommands,
  split(text) {
    const { calls, diagnostics } = splitPymol(text);
    const units: PymolUnit[] = calls.map((call) => ({ name: call.name, at: call.at, source: call.source, call }));
    return { units, diagnostics };
  },
  includeOf(unit) {
    const { call } = unit as PymolUnit;
    // `@file.pml`, and `load file.pml` (PyMOL runs a .pml it loads).
    if (call.name.startsWith('@')) return { name: call.name.slice(1), at: call.at };
    const first = call.args[0];
    if (call.name === 'load' && first && !first.key && /\.pml$/i.test(unquote(first.text))) {
      return { name: unquote(first.text), at: call.at };
    }
    return undefined;
  },
  compile(unit, c) {
    const { call } = unit as PymolUnit;
    const name = call.name;
    if (name === 'help') {
      const topic = call.args[0] ? unquote(call.args[0].text) : '';
      return [{ kind: 'print', text: helpText(topic, call.nameAt) }];
    }
    const command = PymolCommandsByName.get(name) ?? PymolCommandsByName.get(name.toLowerCase());
    if (command) return compilePymolCall(call, command, c);
    const head = name.split('.')[0].toLowerCase();
    const unsupported = PymolUnsupportedByName.get(name.toLowerCase()) ?? PymolUnsupportedByName.get(head);
    if (unsupported) {
      return fail('unsupported', `PyMOL’s ${unsupported.name}: ${unsupported.why}`, call.nameAt, unsupported.hint);
    }
    const near = nearest(name, PymolCommandsByName.keys());
    return fail(
      'unknown-command',
      `PyMOL has no command '${name}'.`,
      call.nameAt,
      near ? `Did you mean ${near}?` : '`help` lists the commands.',
    );
  },
};

function helpText(topic: string, at: { from: number; to: number }): string {
  if (!topic) {
    return `PyMOL commands: ${PymolCommands.map((c) => c.name).join(', ')}. \`help command\` says more.`;
  }
  const c = PymolCommandsByName.get(topic.toLowerCase());
  if (c) return `${c.usage}: ${c.summary}${c.examples.length ? ` Example: ${c.examples[0]}` : ''}`;
  const u = PymolUnsupportedByName.get(topic.toLowerCase());
  if (u) return `${u.name}: ${u.why}${u.hint ? ` ${u.hint}` : ''}`;
  return fail('unknown-command', `PyMOL has no command '${topic}'.`, at, '`help` lists the commands.');
}

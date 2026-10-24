/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { DialectImpl } from './dialect.js';
import { type Diagnostic, error, positionOf, ScriptError, type Span } from './diagnostics.js';
import { type AtomSelectHandle, CompileState, createScriptSession, type ScriptEnv, type ScriptSession } from './env.js';
import type { Command, Dialect } from './ir.js';
import { DialectNames } from './ir.js';
import { PymolDialect } from './pymol/dialect.js';
import { VmdDialect } from './vmd/dialect.js';

// A submission, start to end (_spec §5): split it, compile every command against a working copy of the env (each
// command's names visible to the next), and run nothing if anything did not read; then run the commands in order on
// the host, stopping at the first failure. The host runs them as one step, so a failure is all or nothing (_spec §6).

export const DialectImpls: Readonly<Record<Dialect, DialectImpl>> = { pymol: PymolDialect, vmd: VmdDialect };

/** What the host does with one command (`ScriptHost.execute`). */
export type ExecResult =
  | {
      readonly ok: true;
      /** One line: what it did ("Showed sticks for 87 atoms in 1cbs."). */
      readonly did?: string;
      /** What a read printed. */
      readonly output?: string;
      readonly warnings?: readonly string[];
    }
  | { readonly ok: false; readonly message: string; readonly hint?: string };

/** Where commands run (_spec §5): the Prism viewer's API, or a test's fake. */
export interface ScriptHost {
  /** What names mean now (read once a submission). */
  env(): ScriptEnv;
  execute(command: Command, signal?: AbortSignal): Promise<ExecResult>;
}

export interface ScriptOutput {
  readonly at: Span;
  readonly source: string;
  readonly kind: 'did' | 'value' | 'warning' | 'note';
  readonly text: string;
}

export interface CompiledScript {
  readonly dialect: Dialect;
  readonly commands: readonly Command[];
  /** Errors, warnings and notes; any error means nothing may run. */
  readonly diagnostics: readonly Diagnostic[];
  /** The session as it will be once these run (Tcl variables, atomselect handles). */
  readonly session: ScriptSession;
}

export interface ScriptResult {
  readonly ok: boolean;
  readonly outputs: readonly ScriptOutput[];
  readonly diagnostics: readonly Diagnostic[];
  /** The commands that ran (all of them when `ok`). */
  readonly ran: number;
  readonly commands: readonly Command[];
}

function cloneSession(s: ScriptSession): ScriptSession {
  return { vars: new Map(s.vars), handles: new Map<string, AtomSelectHandle>(s.handles), nextHandle: s.nextHandle };
}

/** A name another dialect has: the hint that says which prefix runs it. */
function otherDialectHint(dialect: Dialect, name: string): string | undefined {
  for (const other of Object.values(DialectImpls)) {
    if (other.id === dialect) continue;
    const known = other.commands.some((c) => c.name === name || c.name.split(' ')[0] === name);
    if (known) return `${name} is ${DialectNames[other.id]}’s: type ${other.id} ${name} …`;
  }
  return undefined;
}

/** Splits and compiles a submission without running it: what the menu previews and `runScript` runs. */
export function compileScript(
  text: string,
  dialect: Dialect,
  env: ScriptEnv,
  session: ScriptSession = createScriptSession(),
): CompiledScript {
  const impl = DialectImpls[dialect];
  const working = cloneSession(session);
  const state = new CompileState(env, working);
  const diagnostics: Diagnostic[] = [];
  const commands: Command[] = [];
  compileText(impl, text, state, commands, diagnostics, undefined, []);
  const history = commands.filter((c) => c.kind === 'history');
  if (history.length && commands.length > 1) {
    diagnostics.push(
      error('syntax', `${history[0].source} runs alone: it would undo this script’s own steps.`, history[0].at),
    );
  }
  return { dialect, commands, diagnostics, session: working };
}

/** Scripts include scripts this deep at most (a script that includes itself stops here). */
const MaxIncludeDepth = 8;

/** Where an included script's commands are: the include's place in the submission, and the file's name. */
interface Including {
  readonly at: Span;
  readonly name: string;
  readonly text: string;
}

/**
 * Compiles a text's commands into `commands` (an included script's too, `@file.pml`, `source file.tcl`: its commands
 * take the include's place, their text named by the file, and their messages say the file's line).
 */
function compileText(
  impl: DialectImpl,
  text: string,
  state: CompileState,
  commands: Command[],
  diagnostics: Diagnostic[],
  including: Including | undefined,
  stack: readonly string[],
): void {
  // A diagnostic of an included file: at the include, its message naming the file's line.
  const place = (d: Diagnostic): Diagnostic =>
    including
      ? {
          ...d,
          at: including.at,
          message: `${including.name} line ${positionOf(including.text, d.at.from).line}: ${d.message}`,
        }
      : d;
  const { units, diagnostics: split } = impl.split(text);
  for (const d of split) diagnostics.push(place(d));
  for (const unit of units) {
    try {
      const include = impl.includeOf?.(unit);
      if (include) {
        const found = state.script(include.name);
        if (!found) {
          throw new ScriptError(
            error(
              'no-name',
              `There is no script '${include.name}' with this one.`,
              include.at,
              'Open the scripts together (drop them all, with the files they load).',
            ),
          );
        }
        if (stack.includes(found.name))
          throw new ScriptError(error('syntax', `${found.name} includes itself.`, include.at));
        if (stack.length >= MaxIncludeDepth) {
          throw new ScriptError(
            error('syntax', `Scripts include scripts at most ${MaxIncludeDepth} deep.`, include.at),
          );
        }
        const at = including?.at ?? unit.at;
        compileText(impl, found.text, state, commands, diagnostics, { at, name: found.name, text: found.text }, [
          ...stack,
          found.name,
        ]);
        continue;
      }
      const bodies = impl.compile(unit, {
        dialect: impl.id,
        state,
        at: unit.at,
        source: unit.source,
        report: (d) => diagnostics.push(place(d)),
      });
      for (const body of bodies) {
        const cmd = {
          ...body,
          at: including?.at ?? unit.at,
          source: including ? `${including.name}: ${unit.source}` : unit.source,
          ...(including ? { file: { name: including.name, line: positionOf(including.text, unit.at.from).line } } : {}),
        } as Command;
        state.apply(cmd);
        commands.push(cmd);
      }
    } catch (e) {
      if (!(e instanceof ScriptError)) throw e;
      const d = e.diagnostic;
      const other = d.code === 'unknown-command' ? otherDialectHint(impl.id, unit.name) : undefined;
      diagnostics.push(place(other ? { ...d, hint: other } : d));
    }
  }
}

export const hasErrors = (ds: readonly Diagnostic[]) => ds.some((d) => d.severity === 'error');

/**
 * Runs a submission on a host (_spec §5): compiled whole first (any error: nothing runs), then each command in order;
 * the first failure stops it. The session (Tcl variables, handles) takes the script's changes when it ran whole.
 */
export async function runScript(
  host: ScriptHost,
  text: string,
  dialect: Dialect,
  options: { readonly session?: ScriptSession; readonly signal?: AbortSignal } = {},
): Promise<ScriptResult> {
  const session = options.session ?? createScriptSession();
  const compiled = compileScript(text, dialect, host.env(), session);
  const diagnostics = [...compiled.diagnostics];
  const outputs: ScriptOutput[] = diagnostics
    .filter((d) => d.severity !== 'error')
    .map((d) => ({
      at: d.at,
      source: text.slice(d.at.from, d.at.to),
      kind: d.severity as 'warning' | 'note',
      text: d.message,
    }));
  if (hasErrors(diagnostics)) return { ok: false, outputs, diagnostics, ran: 0, commands: compiled.commands };
  let ran = 0;
  for (const cmd of compiled.commands) {
    if (options.signal?.aborted) {
      diagnostics.push(error('runtime', 'Stopped.', cmd.at));
      return { ok: false, outputs, diagnostics, ran, commands: compiled.commands };
    }
    const r =
      cmd.kind === 'print' ? ({ ok: true, output: cmd.text } as const) : await host.execute(cmd, options.signal);
    if (!r.ok) {
      // In an included script: its file and line, as its compile errors say them.
      const where = cmd.file ? `${cmd.file.name} line ${cmd.file.line}: ` : '';
      diagnostics.push(error('runtime', `${where}${r.message}`, cmd.at, r.hint));
      return { ok: false, outputs, diagnostics, ran, commands: compiled.commands };
    }
    ran++;
    if (r.did) outputs.push({ at: cmd.at, source: cmd.source, kind: 'did', text: r.did });
    if (r.output) outputs.push({ at: cmd.at, source: cmd.source, kind: 'value', text: r.output });
    for (const w of r.warnings ?? []) outputs.push({ at: cmd.at, source: cmd.source, kind: 'warning', text: w });
  }
  session.vars.clear();
  for (const [k, v] of compiled.session.vars) session.vars.set(k, v);
  session.handles.clear();
  for (const [k, v] of compiled.session.handles) session.handles.set(k, v);
  session.nextHandle = compiled.session.nextHandle;
  return { ok: true, outputs, diagnostics, ran, commands: compiled.commands };
}

/**
 * The scripts among several that no other one includes: what to run when they are opened together (the others run
 * where they are included).
 */
export function rootScripts(scripts: ReadonlyMap<string, string>): string[] {
  const included = new Set<string>();
  for (const [name, text] of scripts) {
    const dialect = dialectOfFile(name);
    if (!dialect) continue;
    const impl = DialectImpls[dialect];
    for (const unit of impl.split(text).units) {
      const inc = impl.includeOf?.(unit);
      if (!inc) continue;
      const base = inc.name.split(/[\\/]/).pop();
      for (const k of scripts.keys()) if (k === inc.name || k.split(/[\\/]/).pop() === base) included.add(k);
    }
  }
  return [...scripts.keys()].filter((k) => !included.has(k));
}

/** The dialect a line names by its first word (`pymol show …`, `vmd mol …`), and the rest. */
export function dialectPrefix(
  text: string,
): { readonly dialect: Dialect; readonly rest: string; readonly offset: number } | undefined {
  const m = text.match(/^\s*(pymol|vmd)(\s+|$)/i);
  if (!m) return undefined;
  return { dialect: m[1].toLowerCase() as Dialect, rest: text.slice(m[0].length), offset: m[0].length };
}

/** The dialect a script file is in, by its extension. */
export function dialectOfFile(name: string): Dialect | undefined {
  const ext = name.toLowerCase().split('.').pop();
  if (ext === 'pml') return 'pymol';
  if (ext === 'tcl' || ext === 'vmd') return 'vmd';
  return undefined;
}

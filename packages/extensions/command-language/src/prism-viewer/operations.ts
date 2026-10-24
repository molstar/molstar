/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ApiArea, type ApiAreaContext, ApiResult } from '@molstar/prism-viewer/api/area';
import { formatDiagnostic } from '../language/diagnostics.js';
import { createScriptSession, EmptyEnv, type ScriptEnv, type ScriptSession } from '../language/env.js';
import { type Dialect, DialectNames } from '../language/ir.js';
import {
  type CompiledScript,
  compileScript,
  dialectOfFile,
  hasErrors,
  runScript,
  type ScriptOutput,
} from '../language/script.js';
import { PrismScriptHost } from './host.js';
import { envOf, StateHolder } from './state.js';

declare module '@molstar/prism-viewer/api' {
  // The `script` area: the extension's own (api.md §2, extensions add areas).
  interface ExtensionApis {
    script: ScriptApi;
  }
}

// The `script` area (_spec §5.1): `script.run` runs a submission of PyMOL or VMD commands as one step through the
// viewer's operations (its calls join its step: Undo takes it back whole, a failure rolls it back), and `script.names`
// lists what names mean. The area keeps the session state (custom state) and a console session (Tcl variables).

export const DialectParam = () =>
  PD.Select<Dialect>('pymol', [
    ['pymol', DialectNames.pymol],
    ['vmd', DialectNames.vmd],
  ]);

/** What `script.run` gives back: what each command printed and said. */
export interface ScriptRunValue {
  readonly outputs: readonly ScriptOutput[];
  readonly ran: number;
}

/** The addresses a submission's loads reach: none for entries (the known services), a URL's own for a URL. */
export function loadAddresses(text: string, dialect: Dialect): string[] {
  // A fresh session, as the agent's run has (Tcl variables are only those the text sets). A command that does not
  // compile here (it names an object) adds no address; the run, compiled against the scene, refuses a load to any
  // address not among these.
  const { commands } = compileScript(text, dialect, EmptyEnv);
  return commands.flatMap((c) => (c.kind === 'load' && 'url' in c.from ? [c.from.url] : []));
}

/** The texts of the scripts among files (by extension), by name. */
export async function scriptTexts(files: readonly File[]): Promise<Map<string, string>> {
  const out = new Map<string, string>();
  for (const f of files) if (dialectOfFile(f.name)) out.set(f.name, await f.text());
  return out;
}

export class ScriptApi extends ApiArea {
  readonly id = 'script';
  readonly #state: StateHolder;
  /** Tcl variables and atomselect handles, across submissions (not state: not undone, not saved, as in VMD). */
  readonly session: ScriptSession = createScriptSession();

  constructor(ctx: ApiAreaContext) {
    super(ctx);
    this.#state = new StateHolder(ctx.plugin);
  }

  /** What names mean now: for the menu's preview and completion (no run). */
  env(files: readonly string[] = []): ScriptEnv {
    return envOf(this.ctx, this.#state.value, files);
  }

  /**
   * A submission compiled against the scene now, without running it: `files`, the names of those given with it;
   * `scripts`, the texts of the scripts among them.
   */
  preview(
    text: string,
    dialect: Dialect,
    files: readonly string[] = [],
    scripts: ReadonlyMap<string, string> = new Map(),
  ): CompiledScript {
    return compileScript(text, dialect, { ...this.env(files), scripts }, this.session);
  }

  readonly run = this.operation({
    id: 'script.run',
    title: 'Run PyMOL or VMD commands',
    params: {
      dialect: DialectParam(),
      text: PD.Text('', { description: 'The commands: one a line, or separated by ; (help lists them)' }),
      files: PD.Value<readonly File[]>([], {
        isHidden: true,
        description: 'Files the commands load by name (load model.pdb, mol new model.pdb)',
      }),
    },
    // What the agent's policy reads (agent.md §4): the addresses its loads reach (entries: the known services; URLs:
    // themselves). The operations it runs inside join its step and are not asked again.
    approval: 'network',
    addresses: ({ dialect, text }) => loadAddresses(text, dialect),
    run: async ({ dialect, text, files }, ctx): Promise<ApiResult<ScriptRunValue>> => {
      if (!text.trim()) return ApiResult.fail('invalid', 'There are no commands.', 'Type one: show sticks, ligand.');
      // The agent's run (agent.md §4): its own session each call, so a Tcl variable set before (by the person or an
      // earlier call) cannot name a URL or a molecule its approval did not see; the loads it will make are exactly the
      // addresses approved; what would ask the person on its own (removing structures, saving a file) is refused, as
      // nothing inside asks.
      const agent = ctx.origin.kind === 'agent';
      const session = agent ? createScriptSession() : this.session;
      if (agent) {
        const compiled = compileScript(text, dialect, this.env((files ?? []).map((f) => f.name)), session);
        const approved = new Set(loadAddresses(text, dialect));
        const reach = compiled.commands.flatMap((c) => (c.kind === 'load' && 'url' in c.from ? [c.from.url] : []));
        const unapproved = reach.find((u) => !approved.has(u));
        if (unapproved) {
          return ApiResult.fail('blocked', `${unapproved} was not among the addresses approved for this run.`);
        }
        const asks = compiled.commands.find((c) => c.kind === 'object.remove' || c.kind === 'export.image');
        if (asks) {
          return ApiResult.fail(
            'blocked',
            `${asks.source} is not run for the assistant: it would ask the person first.`,
            asks.kind === 'object.remove' ? 'Use structure.remove, which asks.' : 'Use export.image, which asks.',
          );
        }
      }
      const byName = new Map((files ?? []).map((f) => [f.name, f]));
      const host = new PrismScriptHost(ctx, this.#state, byName, await scriptTexts(files ?? []));
      const preview = compileScript(text, dialect, host.env(), session);
      if (preview.commands.length === 1 && preview.commands[0].kind === 'history' && !hasErrors(preview.diagnostics)) {
        return ApiResult.fail(
          'invalid',
          `${preview.commands[0].source} is the viewer's own Undo and Redo, not a command of a step.`,
          'Run history.undo or history.redo, or press ⌘Z.',
        );
      }
      host.keepObjects();
      const r = await runScript(host, text, dialect, { session, signal: ctx.signal });
      const failed = r.diagnostics.find((d) => d.severity === 'error');
      if (!r.ok && failed) {
        const others = r.diagnostics.filter((d) => d.severity === 'error').length - 1;
        return ApiResult.fail(
          failed.code === 'runtime' ? 'failed' : 'invalid',
          `${formatDiagnostic(text, { ...failed, hint: undefined } as typeof failed)}${others > 0 ? ` (and ${others} more)` : ''}`,
          failed.hint,
        );
      }
      const did = r.outputs.filter((o) => o.kind === 'did' || o.kind === 'value').map((o) => o.text);
      const warnings = r.outputs.filter((o) => o.kind === 'warning' || o.kind === 'note').map((o) => o.text);
      const value: ScriptRunValue = { outputs: r.outputs, ran: r.ran };
      const changes = r.commands.some((c) => c.kind !== 'read' && c.kind !== 'print');
      const line = did.join(' ') || `Ran ${r.ran} command${r.ran === 1 ? '' : 's'}.`;
      if (!changes) return ApiResult.already(line, { value });
      return ApiResult.ok(line, {}, { value, ...(warnings.length ? { warnings } : {}) });
    },
  });

  readonly names = this.operation({
    id: 'script.names',
    title: 'List the objects and named selections of the commands',
    params: {},
    read: true,
    run: async () => {
      const env = this.env();
      const objects = env.objects.map((o) => `${o.molecule} ${o.name}`);
      const names = [...env.names.entries()].map(([k, s]) => `${k}: ${s.text}`);
      return ApiResult.ok(
        [
          `Objects (VMD molecule ids): ${objects.join(', ') || 'none'}${env.top ? `; top: ${env.top}` : ''}.`,
          `Named selections: ${names.join(', ') || 'none'}.`,
        ].join('\n'),
        {},
        {
          value: {
            objects: env.objects.map(({ name, molecule }) => ({ name, molecule })),
            names: Object.fromEntries(env.names),
          },
        },
      );
    },
  });

  override dispose(): void {
    this.#state.dispose();
  }
}

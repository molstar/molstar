/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ApiArea, type ApiAreaContext, ApiResult } from '@molstar/prism-viewer/api/area';
import { formatDiagnostic } from '../language/diagnostics.js';
import { createScriptSession, type ScriptEnv, type ScriptSession } from '../language/env.js';
import { type Dialect, DialectNames } from '../language/ir.js';
import { type CompiledScript, compileScript, hasErrors, runScript, type ScriptOutput } from '../language/script.js';
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

  /** A submission compiled against the scene now, without running it. */
  preview(text: string, dialect: Dialect): CompiledScript {
    return compileScript(text, dialect, this.env(), this.session);
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
    run: async ({ dialect, text, files }, ctx): Promise<ApiResult<ScriptRunValue>> => {
      if (!text.trim()) return ApiResult.fail('invalid', 'There are no commands.', 'Type one: show sticks, ligand.');
      const byName = new Map((files ?? []).map((f) => [f.name, f]));
      const host = new PrismScriptHost(ctx, this.#state, byName);
      const preview = compileScript(text, dialect, host.env(), this.session);
      if (preview.commands.length === 1 && preview.commands[0].kind === 'history' && !hasErrors(preview.diagnostics)) {
        return ApiResult.fail(
          'invalid',
          `${preview.commands[0].source} is the viewer's own Undo and Redo, not a command of a step.`,
          'Run history.undo or history.redo, or press ⌘Z.',
        );
      }
      host.keepObjects();
      const r = await runScript(host, text, dialect, { session: this.session, signal: ctx.signal });
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

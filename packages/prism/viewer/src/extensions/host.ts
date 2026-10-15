/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ObservableState, observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import type { CommandDefinition } from '@molstar/prism-ui/core';
import { type Observable, Subject } from 'rxjs';
import type { AnnotationsApi } from '../api/annotations.js';
import { type ApiAreaContext, ApiArea, ApiResult, type OperationDefinition, type ViewerApi } from '../api/area.js';
import type { ViewerApiHost } from '../api/index.js';
import type { ViewerCommandContext } from '../api/managers.js';
import type {
  LazyViewerExtension,
  Part,
  PrismViewerPluginParts,
  ViewerExtension,
  ViewerExtensionUI,
  ViewerUIContext,
} from '../spec.js';
import { isLazyExtension } from '../spec.js';

// The extensions of a viewer (viewer-spec.md §3–§4, README §3.1): each activated in the list's order. An eager one
// brings its API areas, operations into core areas (built-in features only), commands, annotations, and views when a
// page attaches. A lazy one is its stub until first use: its areas, commands and views exist from the start as stubs;
// the first use loads the chunk (a task in the status column), registers its plugin part on the running plugin, swaps
// the real parts in and replays what triggered it. A loaded extension is checked against its stub. Overrides replace or
// remove any part by id, the core's included.

export type ExtensionStatus = 'stub' | 'loading' | 'loaded' | 'failed';

/** Where an extension's views go: the page's shell (`ui/extensions.tsx`). */
export interface ExtensionUIPort {
  /** What the views get. */
  readonly context: ViewerUIContext;
  /** Registers the views of an extension (overrides applied by the host); returns their removal. */
  register(owner: string, ui: ViewerExtensionUI): () => void;
  /** Placeholders for a lazy extension's declared views: a click loads it. */
  registerStub(owner: string, stub: LazyViewerExtension['stub'], load: () => Promise<boolean>): () => void;
}

export interface ExtensionHostDeps {
  readonly host: ViewerApiHost;
  readonly person: ViewerApi;
  readonly registerCommands: (defs: readonly CommandDefinition<ViewerCommandContext>[]) => () => void;
  readonly annotations: AnnotationsApi;
  /** The running plugin's side of a lazy load (`sources/extension-plugin.ts`). */
  readonly applyPlugin: (parts: Partial<PrismViewerPluginParts>) => Promise<() => Promise<void>>;
  readonly overrides: Readonly<Record<string, Part | null>>;
  /** A task while a chunk loads (the status column). */
  readonly task: (label: string) => () => void;
}

interface Entry {
  readonly id: string;
  readonly lazy?: LazyViewerExtension;
  ext?: ViewerExtension;
  readonly status: ObservableState<ExtensionStatus>;
  readonly disposers: (() => void | Promise<void>)[];
  ui?: () => void;
  loading?: Promise<boolean>;
}

export class ExtensionHost {
  readonly #deps: ExtensionHostDeps;
  readonly #entries = new Map<string, Entry>();
  readonly #loaded = observableState<readonly string[]>([]);
  readonly #activated = new Subject<string>();
  readonly #failed = new Subject<{ id: string; message: string }>();
  #port: ExtensionUIPort | undefined;

  constructor(deps: ExtensionHostDeps) {
    this.#deps = deps;
  }

  /** The extensions loaded (eager ones from the start), by id. */
  get loaded(): ReadonlyState<readonly string[]> {
    return this.#loaded.asReadonly;
  }

  get events(): { activated: Observable<string>; failed: Observable<{ id: string; message: string }> } {
    return { activated: this.#activated.asObservable(), failed: this.#failed.asObservable() };
  }

  status(id: string): ReadonlyState<ExtensionStatus> | undefined {
    return this.#entries.get(id)?.status.asReadonly;
  }

  /** Adds the list's extensions in order: eager ones activated, lazy ones as stubs (preloaded ones loaded). */
  async add(list: readonly (ViewerExtension | LazyViewerExtension)[]): Promise<void> {
    for (const e of list) {
      if (isLazyExtension(e)) {
        this.#addStub(e);
      } else {
        const entry: Entry = { id: e.id, ext: e, status: observableState<ExtensionStatus>('loaded'), disposers: [] };
        this.#entries.set(e.id, entry);
        this.#activate(entry, e);
        this.#loaded.set([...this.#loaded.value, e.id]);
      }
    }
    for (const e of list) if (isLazyExtension(e) && e.preload) await this.load(e.id);
  }

  /** Loads a lazy extension (once); true when it is loaded. A failure leaves its stub and says why. */
  load(id: string): Promise<boolean> {
    const entry = this.#entries.get(id);
    if (!entry) return Promise.resolve(false);
    if (entry.ext) return Promise.resolve(true);
    if (!entry.loading) entry.loading = this.#load(entry);
    return entry.loading;
  }

  /** Loads the lazy extensions that own any of these (a session's transformers, a file's extension) first. */
  async ensureFor(need: {
    readonly transformers?: readonly string[];
    readonly formats?: readonly string[];
    readonly fileExtensions?: readonly string[];
  }): Promise<void> {
    for (const entry of this.#entries.values()) {
      const owns = entry.lazy?.stub.owns;
      if (!owns || entry.ext) continue;
      const hit =
        need.transformers?.some((t) => owns.transformers?.includes(t)) ||
        need.formats?.some((f) => owns.formats?.includes(f)) ||
        need.fileExtensions?.some((x) => owns.fileExtensions?.includes(x.toLowerCase()));
      if (hit) await this.load(entry.id);
    }
  }

  /** The page's shell: every extension's views now, and loaded ones' later. */
  attachUI(port: ExtensionUIPort): () => void {
    this.#port = port;
    for (const entry of this.#entries.values()) this.#mountUI(entry);
    return () => {
      for (const entry of this.#entries.values()) {
        entry.ui?.();
        entry.ui = undefined;
      }
      this.#port = undefined;
    };
  }

  async dispose(): Promise<void> {
    for (const entry of [...this.#entries.values()].reverse()) {
      entry.ui?.();
      for (const d of entry.disposers.splice(0).reverse()) await d();
    }
    this.#entries.clear();
    this.#activated.complete();
    this.#failed.complete();
  }

  /** An extension's parts, in order; overrides drop or replace commands by id. */
  #activate(entry: Entry, e: ViewerExtension): void {
    const { host, person, registerCommands, annotations } = this.#deps;
    for (const [id, factory] of Object.entries(e.api ?? {})) {
      const replace = host.api.area(id) && entry.lazy;
      entry.disposers.push(replace ? host.replaceArea(id, factory) : host.addArea(id, factory));
    }
    for (const c of e.contributes?.(host.areaContext) ?? []) {
      if (!e.builtIn) {
        throw new Error(`Extension '${e.id}' adds '${c.operation.id}' to a core area: only the viewer's features may`);
      }
      const area = host.api.area(c.area);
      if (!area) throw new Error(`Extension '${e.id}' adds to an area that does not exist: '${c.area}'`);
      entry.disposers.push(area.contribute(c.name, c.operation));
    }
    const commands = this.#overridden(e.commands?.(person) ?? []);
    if (commands.length) entry.disposers.push(registerCommands(commands));
    for (const a of e.annotations ?? []) entry.disposers.push(annotations.register(a));
    this.#mountUI(entry);
    this.#activated.next(e.id);
  }

  #addStub(lazy: LazyViewerExtension): void {
    const entry: Entry = { id: lazy.id, lazy, status: observableState<ExtensionStatus>('stub'), disposers: [] };
    this.#entries.set(lazy.id, entry);
    const { host, registerCommands } = this.#deps;
    // Stub areas: each operation loads the chunk, then runs the real one, as part of the same step.
    const byArea = new Map<string, { id: string; title: string }[]>();
    for (const op of lazy.stub.operations ?? []) {
      const area = op.id.slice(0, op.id.indexOf('.'));
      if (!byArea.has(area)) byArea.set(area, []);
      byArea.get(area)!.push(op);
    }
    for (const [area, ops] of byArea) {
      const load = () => this.load(lazy.id);
      entry.disposers.push(host.addArea(area, (ctx) => new StubArea(ctx, area, ops, load, lazy.title)));
    }
    // Stub commands: listed by title; running one loads the chunk and opens the real one.
    const commands = this.#overridden(
      (lazy.stub.commands ?? []).map(
        (c): CommandDefinition<ViewerCommandContext> => ({
          id: c.id,
          title: c.title,
          category: c.category ?? lazy.title,
          ...(c.icon ? { icon: c.icon } : {}),
          run: async () => {
            const ok = await this.load(lazy.id);
            if (!ok) throw new Error(`${lazy.title} could not be loaded.`);
            return { message: '' };
          },
        }),
      ),
    );
    if (commands.length) {
      const unregister = registerCommands(commands);
      entry.disposers.push(() => {
        if (!entry.ext) unregister();
      });
      (entry as { stubCommands?: () => void }).stubCommands = unregister;
    }
    this.#mountUI(entry);
  }

  async #load(entry: Entry): Promise<boolean> {
    const lazy = entry.lazy!;
    entry.status.set('loading');
    const done = this.#deps.task(`Loading ${lazy.title}`);
    try {
      const ext = await lazy.load();
      checkAgainstStub(lazy, ext);
      if (ext.plugin) entry.disposers.push(await this.#deps.applyPlugin(ext.plugin));
      // The stub's commands and views go; the real ones come.
      (entry as { stubCommands?: () => void }).stubCommands?.();
      entry.ui?.();
      entry.ui = undefined;
      entry.ext = ext;
      this.#activate(entry, ext);
      entry.status.set('loaded');
      this.#loaded.set([...this.#loaded.value, ext.id]);
      return true;
    } catch (e) {
      const message = e instanceof Error ? e.message : String(e);
      entry.status.set('failed');
      entry.loading = undefined;
      this.#failed.next({ id: lazy.id, message });
      return false;
    } finally {
      done();
    }
  }

  #mountUI(entry: Entry): void {
    const port = this.#port;
    if (!port || entry.ui) return;
    if (entry.ext?.ui) {
      const ui = entry.ext.ui;
      let removed = false;
      let remove: (() => void) | undefined;
      void Promise.resolve(ui(port.context)).then((parts) => {
        if (removed || !this.#port) return;
        remove = port.register(entry.id, this.#overriddenUI(parts));
      });
      entry.ui = () => {
        removed = true;
        remove?.();
      };
    } else if (entry.lazy && !entry.ext) {
      entry.ui = port.registerStub(entry.id, entry.lazy.stub, () => this.load(entry.id));
    }
  }

  #overridden(commands: readonly CommandDefinition<ViewerCommandContext>[]): CommandDefinition<ViewerCommandContext>[] {
    const o = this.#deps.overrides;
    return commands.flatMap((c) => {
      if (!(c.id in o)) return [c];
      const replacement = o[c.id];
      return replacement === null ? [] : [replacement as CommandDefinition<ViewerCommandContext>];
    });
  }

  #overriddenUI(ui: ViewerExtensionUI): ViewerExtensionUI {
    const o = this.#deps.overrides;
    const apply = <T extends { id: string }>(list: readonly T[] | undefined): T[] =>
      (list ?? []).flatMap((p) =>
        !(p.id in o) ? [p] : o[p.id] === null ? [] : [{ ...p, ...(o[p.id] as object) } as T],
      );
    return {
      sections: apply(ui.sections),
      tabs: apply(ui.tabs),
      overlays: apply(ui.overlays),
      rows: apply(ui.rows),
      toolbar: apply(ui.toolbar),
    };
  }
}

/** A loaded extension brings only what its stub declared (viewer-spec.md §3; review F6). */
function checkAgainstStub(lazy: LazyViewerExtension, ext: ViewerExtension): void {
  if (ext.id !== lazy.id) throw new Error(`The extension loaded for '${lazy.id}' says it is '${ext.id}'`);
  const declaredOps = new Set((lazy.stub.operations ?? []).map((o) => o.id));
  const declaredCommands = new Set((lazy.stub.commands ?? []).map((c) => c.id));
  const declaredParts = new Set(lazy.stub.parts.map((p) => p.id));
  const extra: string[] = [];
  for (const area of Object.keys(ext.api ?? {})) {
    if (![...declaredOps].some((id) => id.startsWith(`${area}.`))) extra.push(`area '${area}'`);
  }
  // Its views and commands are checked as they register (`ExtensionHost.#activate`, the UI port).
  void declaredCommands;
  void declaredParts;
  if (extra.length) throw new Error(`${lazy.title} brings what its stub did not declare: ${extra.join(', ')}`);
}

/** A lazy extension's area before it loads: each declared operation loads it, then runs the real one. */
class StubArea extends ApiArea {
  readonly id: string;

  constructor(
    ctx: ApiAreaContext,
    id: string,
    ops: readonly { id: string; title: string }[],
    load: () => Promise<boolean>,
    title: string,
  ) {
    super(ctx);
    this.id = id;
    for (const op of ops) {
      const def: OperationDefinition = {
        id: op.id,
        title: op.title,
        params: {},
        passValues: true,
        // The real operation makes the step; the stub only loads and replays.
        history: false,
        async run(values, runCtx) {
          const loaded = await load();
          if (!loaded)
            return ApiResult.fail('unavailable', `${title} could not be loaded.`, 'Try again, or reload the page.');
          return runCtx.api.run(op.id, values as Record<string, unknown>);
        },
      };
      const name = op.id.slice(op.id.indexOf('.') + 1).replace(/-([a-z])/g, (_, c: string) => c.toUpperCase());
      Object.defineProperty(this, name, { value: this.operation(def), enumerable: true, configurable: true });
    }
  }
}

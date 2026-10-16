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
import {
  type ApiAreaContext,
  ApiArea,
  ApiResult,
  memberName,
  type OperationDefinition,
  type ViewerApi,
} from '../api/area.js';
import type { ViewerApiHost } from '../api/index.js';
import type { ViewerCommandContext } from '../api/managers.js';
import type {
  LazyOperationStub,
  LazyViewerExtension,
  PrismViewerPluginParts,
  SequenceBehavior,
  ViewerExtension,
  ViewerExtensionUI,
  ViewerUIContext,
} from '../spec.js';
import { isLazyExtension } from '../spec.js';
import type { Overrides } from './overrides.js';

// The extensions of a viewer (viewer-spec.md §3–§4, README §3.1): each activated in the list's order. An eager one
// brings its API areas, operations into core areas (built-in features only), commands, annotations, and views when a
// page attaches. A lazy one is its stub until first use: its operations (in its own areas, or a built-in feature's in a
// core area), commands and views exist from the start as stubs; the first use loads the chunk (a task in the status
// column), checks it against the stub, registers its plugin part on the running plugin, swaps the real parts in and
// replays what triggered it. Activation is all or nothing: a failure undoes what it did and puts the stub back.
// Overrides (`overrides.ts`) apply to every part, the core's included.

export type ExtensionStatus = 'stub' | 'loading' | 'loaded' | 'failed';

/** Where an extension's views go: the page's shell (`ui/extensions.tsx`). */
export interface ExtensionUIPort {
  /** What the views get. */
  readonly context: ViewerUIContext;
  /** Registers the views of an extension (overrides applied by the host); returns their removal. */
  register(owner: string, ui: ViewerExtensionUI): () => void;
  /** Placeholders for a lazy extension's declared views (overrides applied by the host): a click loads it. */
  registerStub(owner: string, stub: LazyViewerExtension['stub'], load: () => Promise<boolean>): () => void;
}

export interface ExtensionHostDeps {
  readonly host: ViewerApiHost;
  readonly person: ViewerApi;
  /** Registers commands, overrides applied. */
  readonly registerCommands: (defs: readonly CommandDefinition<ViewerCommandContext>[]) => () => void;
  /** Opens a command as the command menu does (its form, or runs it). */
  readonly openCommand: (id: string) => void;
  /** The `annotations` area, when the spec has its feature (viewer-spec.md §3.1); read as each extension activates. */
  readonly annotations: () => AnnotationsApi | undefined;
  /** Attaches a sequence behavior to the sequence viewer; absent: the viewer has no sequence (sequence.md §3). */
  readonly addSequenceBehavior?: (behavior: SequenceBehavior) => () => void;
  /** The running plugin's side of a lazy load (`sources/extension-plugin.ts`). */
  readonly applyPlugin: (parts: Partial<PrismViewerPluginParts>) => Promise<() => Promise<void>>;
  readonly overrides: Overrides;
  /** A task while a chunk loads (the status column). */
  readonly task: (label: string) => () => void;
}

type Disposer = () => void | Promise<void>;

interface Entry {
  readonly id: string;
  readonly lazy?: LazyViewerExtension;
  ext?: ViewerExtension;
  readonly status: ObservableState<ExtensionStatus>;
  /** The active extension's parts. */
  readonly disposers: Disposer[];
  /** The stub's parts, until the extension loads. */
  readonly stub: Disposer[];
  ui?: () => void;
  loading?: Promise<boolean>;
}

/** The definitions of stub operations: their runs are the real ones' replays (the status column shows the result). */
const StubDefinitions = new WeakSet<OperationDefinition<any, any>>();

export function isStubDefinition(def: OperationDefinition<any, any>): boolean {
  return StubDefinitions.has(def);
}

export class ExtensionHost {
  readonly #deps: ExtensionHostDeps;
  readonly #entries = new Map<string, Entry>();
  readonly #loaded = observableState<readonly string[]>([]);
  readonly #activated = new Subject<string>();
  readonly #failed = new Subject<{ id: string; message: string }>();
  #port: ExtensionUIPort | undefined;
  /** Who brought each area an extension added (operations.md §2: an extension names its own areas, once). */
  readonly #areaOwners = new Map<string, string>();
  /** Set by `dispose`: a load in flight then stops at its next step and undoes what it did. */
  #disposed = false;

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
        const entry = this.#entry(e.id, 'stub', e);
        this.#installStub(entry);
        this.#mountUI(entry);
      } else {
        const entry = this.#entry(e.id, 'loaded');
        checkContributes(e);
        entry.disposers.push(...(await this.#activate(e)));
        entry.ext = e;
        this.#mountUI(entry);
        this.#loaded.set([...this.#loaded.value, e.id]);
        this.#activated.next(e.id);
      }
    }
    for (const e of list) if (isLazyExtension(e) && e.preload) await this.load(e.id);
  }

  /**
   * Adds an extension to the running viewer (viewer-spec.md §3): an eager one is activated, its plugin part on the
   * running plugin; a lazy one gets its stub (and loads when asked, or now when `preload`). The extension must be new to
   * the viewer and its areas' names free. Returns its removal.
   */
  async mount(ext: ViewerExtension | LazyViewerExtension): Promise<() => Promise<void>> {
    if (this.#disposed) throw new Error('The viewer is gone.');
    if (this.#entries.has(ext.id)) throw new Error(`Extension '${ext.id}' is in the viewer already`);
    if (isLazyExtension(ext)) {
      await this.add([ext]);
    } else {
      checkContributes(ext);
      const parts = await this.#activate(ext, undefined, true);
      const entry = this.#entry(ext.id, 'loaded');
      entry.disposers.push(...parts);
      entry.ext = ext;
      this.#mountUI(entry);
      this.#loaded.set([...this.#loaded.value, ext.id]);
      this.#activated.next(ext.id);
    }
    return () => this.#unmount(ext.id);
  }

  async #unmount(id: string): Promise<void> {
    const entry = this.#entries.get(id);
    if (!entry) return;
    await entry.loading?.catch(() => false);
    this.#entries.delete(id);
    entry.ui?.();
    entry.ui = undefined;
    await runAll(entry.disposers);
    await runAll(entry.stub);
    this.#loaded.set(this.#loaded.value.filter((x) => x !== id));
  }

  /** An area an extension brings: a name another has fails, naming who has it. */
  #addArea(owner: string, id: string, factory: Parameters<ViewerApiHost['addArea']>[1]): Disposer {
    const { host } = this.#deps;
    if (host.api.area(id)) {
      const by = this.#areaOwners.get(id);
      throw new Error(
        `Extension '${owner}' brings the API area '${id}', which ${by ? `extension '${by}'` : 'the viewer'} has already: an extension names its own areas`,
      );
    }
    const remove = host.addArea(id, factory);
    this.#areaOwners.set(id, owner);
    return () => {
      if (this.#areaOwners.get(id) === owner) this.#areaOwners.delete(id);
      remove();
    };
  }

  /** Loads a lazy extension (once); true when it is loaded. A failure leaves its stub and says why. */
  load(id: string): Promise<boolean> {
    const entry = this.#entries.get(id);
    if (!entry || this.#disposed) return Promise.resolve(false);
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
    const lower = (xs: readonly string[] | undefined) => xs?.map((x) => x.toLowerCase().replace(/^\./, ''));
    const files = lower(need.fileExtensions);
    for (const entry of [...this.#entries.values()]) {
      const owns = entry.lazy?.stub.owns;
      if (!owns || entry.ext) continue;
      const owned = lower(owns.fileExtensions);
      const hit =
        need.transformers?.some((t) => owns.transformers?.includes(t)) ||
        need.formats?.some((f) => owns.formats?.includes(f)) ||
        files?.some((x) => owned?.includes(x));
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
    this.#disposed = true;
    for (const entry of [...this.#entries.values()].reverse()) {
      entry.ui?.();
      await runAll(entry.disposers);
      await runAll(entry.stub);
    }
    this.#entries.clear();
    this.#activated.complete();
    this.#failed.complete();
  }

  #entry(id: string, status: ExtensionStatus, lazy?: LazyViewerExtension): Entry {
    const entry: Entry = { id, ...(lazy ? { lazy } : {}), status: observableState(status), disposers: [], stub: [] };
    this.#entries.set(id, entry);
    return entry;
  }

  /**
   * An extension's parts, in order, all or nothing: on a failure, what was done is undone and the error rethrown. The
   * disposers are returned to the caller, which keeps them once the activation is accepted.
   */
  async #activate(e: ViewerExtension, lazy?: LazyViewerExtension, running = !!lazy): Promise<Disposer[]> {
    const { host, person, registerCommands, annotations } = this.#deps;
    const done: Disposer[] = [];
    try {
      // An eager extension's plugin part is in the composition; a lazy or mounted one's goes on the running plugin.
      if (running && e.plugin) done.push(await this.#deps.applyPlugin(e.plugin));
      for (const [id, factory] of Object.entries(e.api ?? {})) done.push(this.#addArea(e.id, id, factory));
      for (const c of e.contributes?.(host.areaContext) ?? []) {
        const area = host.api.area(c.area);
        if (!area) throw new Error(`Extension '${e.id}' adds to an area that does not exist: '${c.area}'`);
        if (c.name !== memberName(c.operation.id)) {
          throw new Error(
            `Extension '${e.id}' adds '${c.operation.id}' as '${c.name}': its name is '${memberName(c.operation.id)}'`,
          );
        }
        done.push(area.contribute(c.name, c.operation));
      }
      // Its commands are made once its operations exist (they read their params).
      const commands = e.commands?.(person) ?? [];
      if (lazy) checkCommands(lazy, commands);
      if (commands.length) done.push(registerCommands(commands));
      // Its annotations, where the viewer has the `annotations` feature to offer them; without it they are not shown.
      const annotationsArea = e.annotations?.length ? annotations() : undefined;
      for (const a of annotationsArea ? (e.annotations ?? []) : []) done.push(annotationsArea!.register(a));
      // Its sequence behaviors, overrides applied (sequence.md §3), each `<id>.sequence.<name>`.
      for (const b of this.#deps.overrides.parts(e.sequence)) {
        if (!b.id.startsWith(`${e.id}.sequence.`)) {
          throw new Error(
            `Extension '${e.id}' brings a sequence behavior named '${b.id}': not '${e.id}.sequence.<name>'`,
          );
        }
        if (this.#deps.addSequenceBehavior) done.push(this.#deps.addSequenceBehavior(b));
      }
      return done;
    } catch (error) {
      await runAll(done);
      throw error;
    }
  }

  /** A lazy extension's stub: its operations, commands (its views come with the UI). */
  #installStub(entry: Entry): void {
    const lazy = entry.lazy!;
    const { host, registerCommands } = this.#deps;
    const load = () => this.load(lazy.id);
    // Stub operations: each loads the chunk, then runs the real one, as part of the same step. In the extension's own
    // area a stub area; in an area that exists (a built-in feature's in a core area) contributed to it.
    const byArea = new Map<string, LazyOperationStub[]>();
    for (const op of lazy.stub.operations ?? []) {
      const area = op.id.slice(0, op.id.indexOf('.'));
      if (!byArea.has(area)) byArea.set(area, []);
      byArea.get(area)!.push(op);
    }
    for (const [area, ops] of byArea) {
      const existing = host.api.area(area);
      if (existing) {
        if (!lazy.builtIn) {
          const by = this.#areaOwners.get(area);
          throw new Error(
            by
              ? `Extension '${lazy.id}' adds '${ops[0].id}' to the area '${area}', which extension '${by}' has: an extension names its own areas`
              : `Extension '${lazy.id}' adds '${ops[0].id}' to a core area: only the viewer's features may`,
          );
        }
        for (const op of ops) entry.stub.push(existing.contribute(memberName(op.id), stubDefinition(op, load, lazy)));
      } else {
        entry.stub.push(this.#addArea(lazy.id, area, (ctx) => new StubArea(ctx, area, ops, load, lazy)));
      }
    }
    // Stub commands: listed by title; running one loads the chunk, then opens the real one.
    const commands = (lazy.stub.commands ?? []).map(
      (c): CommandDefinition<ViewerCommandContext> => ({
        id: c.id,
        title: c.title,
        category: c.category ?? lazy.title,
        ...(c.icon ? { icon: c.icon } : {}),
        run: async () => {
          if (!(await this.load(lazy.id))) throw new Error(`${lazy.title} could not be loaded.`);
          // Opened once this run has settled: the menu takes one run at a time.
          setTimeout(() => this.#deps.openCommand(c.id), 0);
          return { message: `Loaded ${lazy.title}.` };
        },
      }),
    );
    if (commands.length) entry.stub.push(registerCommands(commands));
  }

  async #load(entry: Entry): Promise<boolean> {
    const lazy = entry.lazy!;
    entry.status.set('loading');
    const done = this.#deps.task(`Loading ${lazy.title}`);
    try {
      const ext = await lazy.load();
      // The viewer went while the chunk loaded: nothing of it goes on the disposed viewer.
      if (this.#disposed) return false;
      // Checked before anything changes (its commands once its operations exist, in `#activate`).
      checkAgainstStub(lazy, ext);
      checkContributes(ext);
      // The stub's parts go and the real ones come; a failure puts the stub back.
      entry.ui?.();
      entry.ui = undefined;
      await runAll(entry.stub);
      if (this.#disposed) return false;
      let parts: Disposer[];
      try {
        parts = await this.#activate(ext, lazy);
        // Disposed while it activated: its parts are undone here, as `dispose` no longer holds this entry.
        if (this.#disposed) {
          await runAll(parts);
          return false;
        }
        const definitions = this.#deps.host.runner.definitions.value;
        const missing = (lazy.stub.operations ?? []).filter((o) => !definitions.has(o.id)).map((o) => o.id);
        if (missing.length) {
          await runAll(parts);
          throw new Error(`${lazy.title} does not bring what its stub declared: ${missing.join(', ')}`);
        }
        // The policy declared is the policy it has: tools and commands read it before the load.
        const differing = (lazy.stub.operations ?? [])
          .filter((o) => !sameFlags(o, definitions.get(o.id)!))
          .map((o) => o.id);
        if (differing.length) {
          await runAll(parts);
          throw new Error(`${lazy.title}'s stub declares other flags than it has: ${differing.join(', ')}`);
        }
      } catch (error) {
        if (this.#disposed) return false;
        this.#installStub(entry);
        this.#mountUI(entry);
        throw error;
      }
      entry.disposers.push(...parts);
      entry.ext = ext;
      this.#mountUI(entry);
      entry.status.set('loaded');
      this.#loaded.set([...this.#loaded.value, ext.id]);
      this.#activated.next(ext.id);
      return true;
    } catch (e) {
      entry.status.set('failed');
      entry.loading = undefined;
      this.#failed.next({ id: lazy.id, message: e instanceof Error ? e.message : String(e) });
      return false;
    } finally {
      done();
    }
  }

  #mountUI(entry: Entry): void {
    const port = this.#port;
    if (!port || entry.ui) return;
    const { overrides } = this.#deps;
    if (entry.ext?.ui) {
      const ui = entry.ext.ui;
      let removed = false;
      let remove: (() => void) | undefined;
      Promise.resolve()
        .then(() => ui(port.context))
        .then((views) => {
          if (removed || !this.#port) return;
          const { checked, problems } = checkViews(entry, views);
          if (problems.length) this.#failed.next({ id: entry.id, message: problems.join('; ') });
          remove = port.register(entry.id, {
            sections: overrides.parts(checked.sections),
            tabs: overrides.parts(checked.tabs),
            overlays: overrides.parts(checked.overlays),
            pages: overrides.parts(checked.pages),
            rows: overrides.parts(checked.rows),
            toolbar: overrides.parts(checked.toolbar),
          });
        })
        .catch((e) => {
          if (!removed) this.#failed.next({ id: entry.id, message: `Its views did not load: ${messageOf(e)}` });
        });
      entry.ui = () => {
        removed = true;
        remove?.();
      };
    } else if (entry.lazy && !entry.ext) {
      const stub = { ...entry.lazy.stub, parts: entry.lazy.stub.parts.filter((p) => !overrides.removes(p.id)) };
      entry.ui = port.registerStub(entry.id, stub, () => this.load(entry.id));
    }
  }
}

/** Whether a stub's declared flags are the definition's (read, history, undoable, approval, addresses, modules). */
function sameFlags(stub: LazyOperationStub, def: OperationDefinition<any, any>): boolean {
  return (
    !!stub.read === !!def.read &&
    (stub.history === false) === (def.history === false) &&
    (stub.undoable === false) === (def.undoable === false) &&
    !!stub.addresses === !!def.addresses &&
    stub.approval === def.approval &&
    [...(stub.modules ?? [])].sort().join() === [...(def.modules ?? [])].sort().join()
  );
}

/** Disposers run in reverse, emptied. */
async function runAll(disposers: Disposer[]): Promise<void> {
  for (const d of disposers.splice(0).reverse()) await d();
}

const messageOf = (e: unknown) => (e instanceof Error ? e.message : String(e));

/** Operations into core areas: built-in features only, under the member name their id gives. */
function checkContributes(e: ViewerExtension): void {
  if (!e.contributes) return;
  if (!e.builtIn) throw new Error(`Extension '${e.id}' adds to a core area: only the viewer's features may`);
}

/** A loaded extension brings only what its stub declared, and every operation it declared (viewer-spec.md §3). */
function checkAgainstStub(lazy: LazyViewerExtension, ext: ViewerExtension): void {
  if (ext.id !== lazy.id) throw new Error(`The extension loaded for '${lazy.id}' says it is '${ext.id}'`);
  if (!!ext.builtIn !== !!lazy.builtIn) throw new Error(`${lazy.title}: its stub and the extension differ in builtIn`);
  const declaredOps = (lazy.stub.operations ?? []).map((o) => o.id);
  const extra = Object.keys(ext.api ?? {})
    .filter((area) => !declaredOps.some((id) => id.startsWith(`${area}.`)))
    .map((area) => `area '${area}'`);
  if (extra.length) throw new Error(`${lazy.title} brings what its stub did not declare: ${extra.join(', ')}`);
}

/** Its commands, as its stub declared them. Its views are checked when they register (`checkViews`). */
function checkCommands(lazy: LazyViewerExtension, commands: readonly CommandDefinition<ViewerCommandContext>[]): void {
  const declared = new Set((lazy.stub.commands ?? []).map((c) => c.id));
  const extra = commands.filter((c) => !declared.has(c.id)).map((c) => `command '${c.id}'`);
  if (extra.length) throw new Error(`${lazy.title} brings what its stub did not declare: ${extra.join(', ')}`);
}

const ViewKinds = {
  sections: 'section',
  tabs: 'tab',
  overlays: 'overlay',
  pages: 'page',
  rows: 'row',
  toolbar: 'toolbar',
} as const;

/**
 * An extension's views: each id `<owner>.<kind>.<name>` with its own id and kind, and, for a lazy one, declared by its
 * stub. The others are not shown, and said.
 */
function checkViews(entry: Entry, views: ViewerExtensionUI): { checked: ViewerExtensionUI; problems: string[] } {
  const declared = entry.lazy ? new Set(entry.lazy.stub.parts.map((p) => p.id)) : undefined;
  const problems: string[] = [];
  const checked: Record<string, unknown[]> = {};
  for (const [key, kind] of Object.entries(ViewKinds) as [keyof typeof ViewKinds, string][]) {
    checked[key] = (views[key] ?? []).filter((p) => {
      const [owner, k] = p.id.split('.');
      if (owner !== entry.id || k !== kind) {
        problems.push(`its ${kind} '${p.id}' is not named '${entry.id}.${kind}.<name>'`);
        return false;
      }
      if (declared && !declared.has(p.id)) {
        problems.push(`its ${kind} '${p.id}' is not declared by its stub`);
        return false;
      }
      return true;
    });
  }
  return { checked: checked as ViewerExtensionUI, problems };
}

/** A stub operation: loads the extension, then runs the real one with the same values (one step, the real one's). */
function stubDefinition(
  op: LazyOperationStub,
  load: () => Promise<boolean>,
  lazy: LazyViewerExtension,
): OperationDefinition {
  const def: OperationDefinition = {
    id: op.id,
    title: op.title,
    params: {},
    passValues: true,
    // Its flags are the declared ones (the agent's policy reads them before the load); the real operation, replayed
    // inside it, joins its step.
    ...(op.read ? { read: true as const } : {}),
    ...(op.history === false ? { history: false as const } : {}),
    ...(op.undoable === false ? { undoable: false as const } : {}),
    ...(op.approval ? { approval: op.approval } : {}),
    ...(op.addresses ? { addresses: op.addresses } : {}),
    ...(op.modules ? { modules: op.modules } : {}),
    async run(values, ctx) {
      if (!(await load())) {
        return ApiResult.fail('unavailable', `${lazy.title} could not be loaded.`, 'Try again, or reload the page.');
      }
      return ctx.api.run(op.id, values as Record<string, unknown>);
    },
  };
  StubDefinitions.add(def);
  return def;
}

/** A lazy extension's area before it loads: each declared operation loads it, then runs the real one. */
class StubArea extends ApiArea {
  readonly id: string;

  constructor(
    ctx: ApiAreaContext,
    id: string,
    ops: readonly LazyOperationStub[],
    load: () => Promise<boolean>,
    lazy: LazyViewerExtension,
  ) {
    super(ctx);
    this.id = id;
    for (const op of ops) {
      const fn = this.operation(stubDefinition(op, load, lazy));
      Object.defineProperty(this, memberName(op.id), { value: fn, enumerable: true, configurable: true });
    }
  }
}

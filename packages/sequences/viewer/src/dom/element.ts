/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/// <reference path="./css.d.ts" />
// The vanilla entry (seqspec:interaction#the-dom-element): `mountSeqViewer` places the chrome's composition as React's
// `SeqViewerLayout` does (and returns its `ui`, for parts a host mounts elsewhere), and the `<seq-viewer>` element mounts
// it in a shadow root, the skin as adopted stylesheets.
// Properties, not attributes; each member changes only what it names: `chrome`, `ui` and `slots` remount the UI on the
// same context, `spec` and `ports` remake the context the element made, `document` and `data` load once per (context,
// object). A context the element made outlives a move in the DOM (its disposal waits a task, which a reconnect
// cancels); a host's outlives the element.

import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import type { Loci } from '@molstar/seq-core/model';
import { mountSection, windowScheduler } from '../canvas/index.js';
import { isDomPart, labelText, placedParts, resolveSlots } from '../chrome.js';
import {
  createSeqViewer,
  DefaultSeqViewerSpec,
  type LegendMode,
  type LociSize,
  lociSize,
  type SeqDataInput,
  SeqViewerConfig,
  type SeqViewerContext,
  type SeqViewerPorts,
  type SeqViewerSpec,
  type Unsubscribable,
} from '../index.js';
import skin from '../skin.css?inline';
import {
  cx,
  type DomPart,
  type DomSlotName,
  type DomUI,
  type DomUISpec,
  Hint,
  h,
  later,
  Menu,
  off,
  type PartProps,
  ScrollArea,
} from './core.js';
import { withDownloadExport } from './export.js';
import { HoverCardPart } from './hover-card.js';
import * as P from './parts.js';
import { ViewSurface } from './section.js';

const DEFAULTS: DomUI['slots'] = {
  ViewTabs: P.ViewTabs,
  SectionHeader: P.SectionHeader,
  SectionSwitches: P.SectionSwitches,
  TrackActions: P.TrackActions,
  HiddenChips: P.HiddenChips,
  RowControl: null,
  Tooltip: P.Tooltip,
  HoverCard: HoverCardPart,
  Legend: P.Legend,
  SnapshotPlayer: P.SnapshotPlayer,
  Description: P.Description,
  Search: P.Search,
  Export: P.Export,
  Issues: P.IssueBadge,
  Busy: P.Busy,
  EmptyState: P.EmptyState,
  ErrorState: P.ErrorState,
  Hint,
  ScrollArea,
  Menu,
};

/**
 * A UI for parts to share (React's provider): the slots over the defaults (`resolveSlots`), the composition the layout
 * places (`placedParts`, by `spec.chrome`; the context's config is never written, seqspec:interaction#slots-and-chrome),
 * and the surfaces map `ViewSurface` fills (a host's own, kept across remounts, or a new one).
 */
export function createDomUI(
  context: SeqViewerContext,
  spec: DomUISpec = {},
  surfaces: Map<string, HTMLElement> = new Map(),
): DomUI {
  return {
    slots: resolveSlots(context, DEFAULTS, spec.slots, isDomPart),
    placed: placedParts(context, spec.chrome, spec.slots, spec.legend),
    className: spec.className ?? '',
    legend: spec.legend,
    surface: spec.surface ?? mountSection,
    search: spec.search,
    surfaces,
  };
}

/**
 * The viewer's root in `container` for a composition: `.seqv-root` (the box the canvas overlays are placed in, focusable
 * by script: where focus lands when nothing of its own remains, seqspec:interaction#accessibility), what `compose`
 * mounts into it, then a visually hidden polite live region with the keyboard cursor's label, never the pointer's.
 * React's `SeqViewerRoot`.
 */
export function mountSeqViewerRoot(
  container: HTMLElement | ShadowRoot,
  context: SeqViewerContext,
  ui: DomUI,
  compose: (root: HTMLElement) => readonly { dispose(): void }[],
): { readonly root: HTMLElement; dispose(): void } {
  const root = h('div', { class: cx('seqv-root', ui.className), tabindex: -1 });
  container.append(root);
  const parts = compose(root);
  const live = h('span', { 'aria-live': 'polite', class: 'seqv-sr-only' });
  root.append(live);
  const sub = context.labels.current.subscribe((l) => {
    if (context.hover.live.current.owner === 'keyboard') live.textContent = labelText(l);
  });
  return {
    root,
    dispose() {
      off([sub, ...parts]);
      root.remove();
    },
  };
}

/**
 * The default composition (`ChromeCompositions`, by `spec.chrome`): a bar (views, search, busy, export, issues), the
 * player, the description, the states, the sections, the legend, and the canvas overlays. Returns its `ui`: a part the
 * host mounts elsewhere (its slot "none" here) shares its slots and surfaces with `{ ui }`.
 */
export function mountSeqViewer(
  container: HTMLElement | ShadowRoot,
  context: SeqViewerContext,
  spec: DomUISpec = {},
  surfaces?: Map<string, HTMLElement>,
): { readonly ui: DomUI; dispose(): void } {
  const ui = createDomUI(context, spec, surfaces);
  const mounted = mountSeqViewerRoot(container, context, ui, (root) => {
    const bar = h('div', { class: 'seqv-bar' });
    root.append(bar);
    const parts: { dispose(): void }[] = [];
    const part = (name: DomSlotName | 'ViewSurface', host = root) => {
      const Part = (name === 'ViewSurface' ? ViewSurface : ui.placed.has(name) ? ui.slots[name] : null) as DomPart<
        PartProps & { ui: DomUI }
      > | null;
      if (Part) parts.push(Part(host, context, { ui }));
    };
    part('ViewTabs', bar);
    bar.append(h('span', { class: 'seqv-spacer' }));
    for (const p of ['Search', 'Busy', 'Export', 'Issues'] as const) part(p, bar);
    for (const p of [
      'SnapshotPlayer',
      'Description',
      'EmptyState',
      'ErrorState',
      'ViewSurface',
      'Legend',
      'Tooltip',
      'HoverCard',
    ] as const)
      part(p);
    return parts;
  });
  return { ui, dispose: mounted.dispose };
}

/** The skin in the shadow root; its `--seq-*` defaults on the document at zero specificity, so `data-theme` on the
 * page's root and a host's tokens on any ancestor (or the element) reach the parts, as with the React UI. */
let sheets: { shadow: CSSStyleSheet; tokens: CSSStyleSheet } | null = null;
function adoptSkin(doc: Document, shadow: ShadowRoot) {
  if (!sheets) {
    const all = new CSSStyleSheet();
    all.replaceSync(`:host { display: block; }\n${skin}`);
    const tokens = new CSSStyleSheet();
    for (const r of Array.from(all.cssRules))
      if (r.cssText.includes(':root')) tokens.insertRule(r.cssText, tokens.cssRules.length);
    sheets = { shadow: all, tokens };
  }
  const { tokens, shadow: s } = sheets;
  if (!doc.adoptedStyleSheets.includes(tokens)) doc.adoptedStyleSheets = [...doc.adoptedStyleSheets, tokens];
  if (!shadow.adoptedStyleSheets.includes(s)) shadow.adoptedStyleSheets = [s];
}

/**
 * Loci as the edge speaks them (seqspec:interaction#focus): `Coordinates` with the descriptor of each system, by key.
 */
function described(ctx: SeqViewerContext, loci: readonly Loci[]) {
  const keys: string[] = [];
  const runs: number[] = [];
  for (const l of loci) {
    const c = ctx.systems.coordinatesOf(l);
    for (const k of c.systems) if (!keys.includes(k)) keys.push(k);
    for (let i = 0; i < c.runs.length; i += 3)
      runs.push(keys.indexOf(c.systems[c.runs[i]]), c.runs[i + 1], c.runs[i + 2]);
  }
  const systems: Record<string, ReturnType<SeqViewerContext['systems']['descriptor']> & object> = {};
  for (const k of keys) {
    const d = ctx.systems.descriptor(k);
    if (d) systems[k] = d;
  }
  return {
    coordinates: Object.freeze({ systems: Object.freeze(keys), runs: Int32Array.from(runs) }),
    systems: Object.freeze(systems),
  };
}

/**
 * What `seq-hover` carries, frozen: the hovered element's numbers, copied when the event is built; `loci`,
 * `coordinates` and `systems` computed on their first read (a host that reads none allocates nothing for them) from
 * the hover then, and kept. Read them in the listener: once the hover has moved on they are null and `{}`.
 */
export interface SeqHoverDetail {
  readonly owner: string;
  readonly section: string | null;
  readonly track: string | null;
  /** `HitRegion`. */
  readonly region: number;
  readonly block: number;
  readonly index: number;
  /** What is lit, closed (`hover.toLoci()`); null: nothing is hovered. */
  readonly loci: Loci | null;
  /** Its size (`hover.size()`: items, covered positions), known before `loci` resolves: a large one can be left alone. */
  readonly size: LociSize | null;
  readonly coordinates: ReturnType<typeof described>['coordinates'] | null;
  readonly systems: ReturnType<typeof described>['systems'];
}
/** What `seq-selection` carries, frozen: the selection value and its loci as coordinates with their systems. */
export interface SeqSelectionDetail {
  readonly selection: SeqViewerContext['selection']['current']['value'];
  readonly origin: SeqViewerContext['selection']['current']['value']['origin'];
  readonly coordinates: ReturnType<typeof described>['coordinates'];
  readonly systems: ReturnType<typeof described>['systems'];
  /** The selection's size (items, covered positions), as `lociSize` reads its parts. */
  readonly size: LociSize;
}

const NO_SYSTEMS: SeqHoverDetail['systems'] = Object.freeze({});
/** The hover at `version`; its loci and their coordinates are a reader's cost (an event a frame, read or not). */
function hoverDetail(ctx: SeqViewerContext, version: number): SeqHoverDetail {
  const h = ctx.hover.live.current;
  const view = ctx.graph.assembled.value?.view;
  let lit: { loci: Loci | null; d: ReturnType<typeof described> | null } | null = null;
  let sized: { size: LociSize | null } | null = null;
  const read = () => {
    if (lit) return lit;
    const loci =
      h.owner === '' || ctx.hover.live.version.value !== version ? null : (ctx.hover.toLoci() as Loci | null);
    lit = { loci, d: loci ? described(ctx, [loci]) : null };
    return lit;
  };
  return Object.freeze({
    owner: h.owner,
    section: ctx.layout.sectionRef(h.section),
    track: view && h.track >= 0 ? view.set.strings.get(view.resolved.tracks.ref[h.track]) : null,
    region: h.region,
    block: h.block,
    index: h.index,
    get loci() {
      return read().loci;
    },
    get size() {
      sized ??= { size: h.owner === '' || ctx.hover.live.version.value !== version ? null : ctx.hover.size() };
      return sized.size;
    },
    get coordinates() {
      return read().d?.coordinates ?? null;
    },
    get systems() {
      return read().d?.systems ?? NO_SYSTEMS;
    },
  });
}
function selectionDetail(ctx: SeqViewerContext): SeqSelectionDetail {
  const s = ctx.selection.current.value;
  const parts = [s.seq, s.aligned, s.columns, s.rows, s.items, s.boundaries].filter((l) => l !== null) as Loci[];
  const d = described(ctx, parts);
  const size = { items: 0, length: 0 };
  for (const l of parts) {
    const z = lociSize(l);
    size.items += z?.items ?? 0;
    size.length += z?.length ?? 0;
  }
  return Object.freeze({ selection: s, origin: s.origin, ...d, size: Object.freeze(size) });
}

const PROPS = ['context', 'spec', 'ports', 'ui', 'slots', 'chrome', 'legend', 'document', 'data'] as const;
const LEGEND_MODES: readonly string[] = ['collapsed', 'expanded', 'hidden'];

/**
 * `<seq-viewer>`: set `context` (outlives the element) or `spec` (a context made on connect), then `document` or
 * `data`. Re-dispatches `seq-context`, `seq-selection`, `seq-hover` and `seq-snapshot`, composed and bubbling, each
 * with a frozen detail. Light-DOM slots: `header` above the viewer, `toolbar` in its bar. `chrome` (property or
 * attribute, "full" or "clean") and `legend` (property or attribute, "collapsed", "expanded" or "hidden";
 * seqspec:interaction#the-legend) are this UI's options; the context's config is never written.
 */
export class SeqViewerElement extends HTMLElement {
  static readonly observedAttributes = ['chrome', 'legend'];
  #host: SeqViewerContext | null = null;
  #made: SeqViewerContext | null = null;
  #spec: SeqViewerSpec | null = null;
  #ports: Partial<SeqViewerPorts> | null = null;
  #ui: DomUISpec | null = null;
  #slots: DomUISpec['slots'] = {};
  #chrome: 'full' | 'clean' | null = null;
  #legend: LegendMode | null = null;
  #document: SeqViewSpecDoc | null = null;
  #data: SeqDataInput | null = null;
  /** What was last loaded into each context: setting the same object again does nothing, nor does a context set again. */
  #loaded = new WeakMap<SeqViewerContext, { document: object | null; data: object | null }>();
  #told: SeqViewerContext | null = null;
  /** The UI mounted (`ui` null: it failed, its error shown in its place). */
  #mounted: { readonly ui: DomUI | null; dispose(): void } | null = null;
  /** One surfaces map across remounts: a host's part on `mountedUI` reads the surfaces the next mount fills. */
  readonly #surfaces = new Map<string, HTMLElement>();
  #subs: Unsubscribable[] = [];
  #leaving = 0;

  constructor() {
    super();
    this.attachShadow({ mode: 'open' });
    // Properties a host set before the element was defined shadow the accessors: pass them through.
    for (const k of PROPS) {
      if (!Object.hasOwn(this, k)) continue;
      const v = (this as Record<string, unknown>)[k];
      delete (this as Record<string, unknown>)[k];
      (this as Record<string, unknown>)[k] = v;
    }
  }

  /** The host's context, else the one made from `spec` while connected. */
  get context(): SeqViewerContext | null {
    return this.#host ?? this.#made;
  }
  /** The UI is remounted on it; a context the element made is disposed a task later; `document` and `data` load into it. */
  set context(c: SeqViewerContext | null) {
    if (c === this.#host) return;
    this.#host = c;
    if (c) this.#retire();
    if (this.#mounted) this.#remountUI();
    else if (c) this.#announce(c);
  }
  get spec() {
    return this.#spec;
  }
  /** Read when the element makes its context: a later spec remakes it (a host's context stays). */
  set spec(s: SeqViewerSpec | null) {
    this.#spec = s;
    this.#respec();
  }
  get ports() {
    return this.#ports;
  }
  set ports(p: Partial<SeqViewerPorts> | null) {
    this.#ports = p;
    this.#respec();
  }
  get ui() {
    return this.#ui;
  }
  set ui(u: DomUISpec | null) {
    this.#ui = u;
    this.#remountUI();
  }
  get slots(): DomUISpec['slots'] {
    return this.#slots;
  }
  set slots(s: DomUISpec['slots']) {
    this.#slots = s ?? {};
    this.#remountUI();
  }
  get chrome(): 'full' | 'clean' {
    return this.#chrome ?? this.#ui?.chrome ?? this.context?.config.get(SeqViewerConfig.UI.Chrome) ?? 'clean';
  }
  set chrome(c: 'full' | 'clean' | null) {
    this.#chrome = c === 'clean' || c === 'full' ? c : null;
    this.#remountUI();
  }
  /**
   * The mounted UI (its slots and surfaces): parts a host mounts outside the element share it; null while disconnected
   * or when the mount failed. Replaced on each remount (`chrome`, `ui`, `slots`, `legend`, `spec`, `ports`), which
   * `seq-ui` announces; its surfaces map stays the element's.
   */
  get mountedUI(): DomUI | null {
    return this.#mounted?.ui ?? null;
  }
  get legend(): LegendMode {
    return this.#legend ?? this.#ui?.legend ?? this.context?.config.get(SeqViewerConfig.UI.Legend) ?? 'collapsed';
  }
  set legend(m: LegendMode | null) {
    this.#legend = m !== null && LEGEND_MODES.includes(m) ? m : null;
    this.#remountUI();
  }
  attributeChangedCallback(name: string, _old: string | null, value: string | null) {
    if (name === 'legend') this.legend = value as LegendMode | null;
    else this.chrome = value as 'full' | 'clean' | null;
  }
  get document() {
    return this.#document;
  }
  set document(d: SeqViewSpecDoc | null) {
    this.#document = d;
    this.#load();
  }
  get data() {
    return this.#data;
  }
  set data(d: SeqDataInput | null) {
    this.#data = d;
    this.#load();
  }

  connectedCallback() {
    window.clearTimeout(this.#leaving);
    if (!this.#mounted) this.#mount();
  }

  /** The UI goes at once; a context the element made stays for a task, so a move in the DOM keeps it. */
  disconnectedCallback() {
    this.#unmount();
    const made = this.#made;
    if (!made) return;
    this.#leaving = window.setTimeout(() => {
      if (this.#made === made) this.#made = null;
      made.dispose();
    }, 0);
  }

  #make(): SeqViewerContext {
    const spec = withDownloadExport(this.#spec ?? DefaultSeqViewerSpec());
    this.#made = createSeqViewer(spec, { scheduler: windowScheduler(), ...this.#ports });
    return this.#made;
  }

  /** Disposes the context the element made, a task later (nothing it was showing is cut off mid-frame). */
  #retire() {
    window.clearTimeout(this.#leaving);
    const old = this.#made;
    this.#made = null;
    if (old) window.setTimeout(() => old.dispose(), 0);
  }

  #mount() {
    const shadow = this.shadowRoot as ShadowRoot;
    adoptSkin(this.ownerDocument, shadow);
    shadow.append(h('slot', { name: 'header' }));
    let ctx: SeqViewerContext;
    try {
      ctx = this.#host ?? this.#made ?? this.#make();
      const ui = this.#ui;
      const spec = { ...ui, slots: { ...ui?.slots, ...this.#slots } };
      this.#mounted = mountSeqViewer(
        shadow,
        ctx,
        { ...spec, chrome: this.#chrome ?? ui?.chrome, legend: this.#legend ?? ui?.legend },
        this.#surfaces,
      );
    } catch (e) {
      return this.#fail(shadow, e);
    }
    // The bar takes its look while the host's toolbar has something in it (a slot's content is not the bar's child).
    const toolbar = h('slot', { name: 'toolbar' });
    const bar = shadow.querySelector('.seqv-bar');
    bar?.querySelector(':scope > .seqv-spacer')?.after(toolbar);
    toolbar.addEventListener('slotchange', () =>
      bar?.toggleAttribute('data-seqv-filled', toolbar.assignedElements().length > 0),
    );
    this.#subs = this.#relay(ctx);
    this.#announce(ctx);
    this.dispatchEvent(
      new CustomEvent('seq-ui', { detail: Object.freeze({ ui: this.mountedUI }), bubbles: true, composed: true }),
    );
    this.#load();
  }

  /**
   * A mount that threw (a spec whose registration conflicts, a part's error): never an empty element. The error is shown
   * in the viewer's place and dispatched as `seq-error`; setting `spec`, `ports`, `ui` or `context` tries again.
   */
  #fail(shadow: ShadowRoot, e: unknown) {
    const message = e instanceof Error ? e.message : String(e);
    const root = h(
      'div',
      { class: 'seqv-root', tabindex: -1 },
      h('div', { role: 'alert', class: 'seqv-state seqv-error' }, message),
    );
    shadow.append(root);
    this.#mounted = { ui: null, dispose: () => root.remove() };
    this.dispatchEvent(
      new CustomEvent('seq-error', { detail: Object.freeze({ error: e }), bubbles: true, composed: true }),
    );
  }

  #unmount() {
    off([...this.#subs, this.#mounted]);
    [this.#subs, this.#mounted] = [[], null];
    this.shadowRoot?.replaceChildren();
  }

  /** `chrome`, `ui`, `slots`: the UI alone, on the same context (its document, selection, viewports and history stay). */
  #remountUI() {
    if (!this.#mounted) return;
    this.#unmount();
    this.#mount();
  }

  /** `spec`, `ports`: a new context for the element's own, loaded again; a host's context ignores them. */
  #respec() {
    if (this.#host) return;
    if (!this.#mounted) return this.#retire();
    this.#unmount();
    this.#retire();
    this.#mount();
  }

  #announce(ctx: SeqViewerContext) {
    if (this.#told === ctx) return;
    this.#told = ctx;
    this.dispatchEvent(
      new CustomEvent('seq-context', { detail: Object.freeze({ context: ctx }), bubbles: true, composed: true }),
    );
  }

  /** `document` and `data` go in with `setDocument` / `setData`, each once per (context, object); a failure is the context's issue. */
  #load() {
    const ctx = this.context;
    if (!ctx || !this.#mounted) return;
    const done = this.#loaded.get(ctx) ?? { document: null, data: null };
    this.#loaded.set(ctx, done);
    if (this.#data && done.data !== this.#data) {
      done.data = this.#data;
      ctx.setData(this.#data).catch(() => {});
    }
    if (this.#document && done.document !== this.#document) {
      done.document = this.#document;
      ctx.setDocument(this.#document).catch(() => {});
    }
  }

  /** The context's changes a host commonly wants, as DOM events (not their replay on subscribe). */
  #relay(ctx: SeqViewerContext): Unsubscribable[] {
    const emit = (type: string, detail: unknown) =>
      this.dispatchEvent(new CustomEvent(type, { detail, bubbles: true, composed: true }));
    // The hover moves at pointer rate: at most one event a frame, after a change (`frame` is pending until it runs), the
    // detail built then for the element hovered by then.
    let frame = 0;
    const send = () => {
      frame = 0;
      emit('seq-hover', hoverDetail(ctx, ctx.hover.live.version.value));
    };
    const hover = later(ctx.hover.live.version, () => {
      frame ||= ctx.scheduler.requestFrame(send);
    });
    return [
      later(ctx.selection.current, () => emit('seq-selection', selectionDetail(ctx))),
      hover,
      { unsubscribe: () => void (frame && ctx.scheduler.cancelFrame(frame)) },
      later(ctx.snapshots.current, (ref) => emit('seq-snapshot', Object.freeze({ ref }))),
    ];
  }
}

/** Registers the element once (`seq-viewer` by default; another tag gets a subclass). */
export function defineSeqViewer(tag = 'seq-viewer'): void {
  if (customElements.get(tag)) return;
  customElements.define(tag, tag === 'seq-viewer' ? SeqViewerElement : class extends SeqViewerElement {});
}

declare global {
  interface HTMLElementTagNameMap {
    'seq-viewer': SeqViewerElement;
  }
}

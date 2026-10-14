/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The provider, the root, the default composition of the parts, and `SeqViewer` (seqspec:interaction#the-default-ui): a
// host that keeps the context in its own model renders the parts under `SeqViewerProvider`, in `SeqViewerLayout` (the
// chrome's composition) or in its own arrangement around `SeqViewerRoot` and `ViewSurface`; `SeqViewer` makes a context
// per mount for everyone else.

import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { type ComponentType, type ReactNode, useEffect, useMemo, useRef, useState } from 'react';
import { mountSection, windowScheduler } from '../canvas/index.js';
import { isReactPart, labelText, placedParts, resolveSlots } from '../chrome.js';
import { withDownloadExport } from '../dom/export.js';
import {
  createSeqViewer,
  DefaultSeqViewerSpec,
  type SeqDataInput,
  type SeqViewerContext,
  type SeqViewerPorts,
  type SeqViewerSpec,
} from '../index.js';
import {
  type PartProps,
  type ResolvedUISpec,
  ScopeContext,
  type SeqViewerUISpec,
  type SlotName,
  useScope,
} from './hooks.js';
import { HoverCard } from './hover-card.js';
import * as P from './parts.js';
import { ViewSurface } from './section.js';

const DEFAULTS: ResolvedUISpec['components'] = {
  ViewTabs: P.ViewTabs,
  SectionHeader: P.SectionHeader,
  SectionSwitches: P.SectionSwitches,
  TrackActions: P.TrackActions,
  HiddenChips: P.HiddenChips,
  RowControl: null,
  Tooltip: P.Tooltip,
  HoverCard,
  Legend: P.Legend,
  SnapshotPlayer: P.SnapshotPlayer,
  Description: P.Description,
  Search: P.Search,
  Export: P.Export,
  Issues: P.IssueBadge,
  Busy: P.Busy,
  EmptyState: P.EmptyState,
  ErrorState: P.ErrorState,
  Hint: P.Hint,
  ScrollArea: P.ScrollArea,
  Menu: P.Menu,
};
/**
 * Slots over the defaults (`resolveSlots`), and what the layout places: this UI's chrome, else the viewer's
 * (`placedParts`, seqspec:interaction#slots-and-chrome).
 */
function resolve(context: SeqViewerContext, ui: SeqViewerUISpec | undefined): ResolvedUISpec {
  return {
    components: resolveSlots(context, DEFAULTS, ui?.components, isReactPart),
    placed: placedParts(context, ui?.chrome, ui?.components, ui?.legend),
    className: ui?.className ?? '',
    legend: ui?.legend,
    surface: ui?.surface ?? mountSection,
    search: ui?.search,
  };
}

/** `ui` two levels deep: `className`, `chrome`, `legend`, `surface` and each `components` entry by identity. */
function sameUI(a: SeqViewerUISpec | undefined, b: SeqViewerUISpec | undefined): boolean {
  if (a === b) return true;
  if (
    a?.className !== b?.className ||
    a?.chrome !== b?.chrome ||
    a?.legend !== b?.legend ||
    a?.surface !== b?.surface ||
    a?.search !== b?.search
  )
    return false;
  const [x, y] = [a?.components ?? {}, b?.components ?? {}] as Record<string, unknown>[];
  const keys = new Set([...Object.keys(x), ...Object.keys(y)]);
  for (const k of keys) if (x[k] !== y[k]) return false;
  return true;
}

/**
 * `ui` may be an inline object: while it is equal in two levels to the last one (`sameUI`) the scope is the same
 * object, so what reads it sees no change and every mounted surface and its keyboard focus stay. The surfaces map is
 * the provider instance's own, never remade with `ui`. A `surface` that is a new function each render does remount them.
 */
export function SeqViewerProvider(props: { context: SeqViewerContext; ui?: SeqViewerUISpec; children: ReactNode }) {
  const { context } = props;
  const kept = useRef(props.ui);
  if (!sameUI(kept.current, props.ui)) kept.current = props.ui;
  const ui = kept.current;
  const [surfaces] = useState(() => new Map<string, HTMLElement>());
  const scope = useMemo(() => ({ context, ui: resolve(context, ui), surfaces }), [context, ui, surfaces]);
  return <ScopeContext.Provider value={scope}>{props.children}</ScopeContext.Provider>;
}

/**
 * A visually hidden polite live region: the keyboard cursor's label, never the pointer's
 * (seqspec:interaction#accessibility).
 */
function Announcer({ context: ctx }: PartProps) {
  const el = useRef<HTMLSpanElement>(null);
  useEffect(() => {
    const sub = ctx.labels.current.subscribe((l) => {
      if (el.current && ctx.hover.live.current.owner === 'keyboard') el.current.textContent = labelText(l);
    });
    return () => sub.unsubscribe();
  }, [ctx]);
  return <span ref={el} aria-live="polite" className="seqv-sr-only" />;
}

/**
 * The viewer's root under the provider: the box its canvas overlays (the tooltip, the hover card) are placed in, focusable
 * by script (where focus lands when nothing of its own remains), with the keyboard cursor's live region. `children` is the
 * composition: `SeqViewerLayout`'s, or a host's own (`ViewSurface`, the overlays, any parts).
 */
export function SeqViewerRoot({ children }: { children?: ReactNode }) {
  const { context, ui } = useScope();
  return (
    <div className={P.cx('seqv-root', ui.className)} tabIndex={-1}>
      {children}
      <Announcer context={context} />
    </div>
  );
}

/**
 * The default composition (`ChromeCompositions`, by this UI's `chrome`): a bar (views, search, busy, export, issues),
 * the player, the description, the states, the sections (their headers and track menus), the legend, and the canvas
 * overlays. A part a host places elsewhere is set to "none" here.
 */
export function SeqViewerLayout() {
  const { context, ui } = useScope();
  const part = (name: SlotName) => {
    const Part = (ui.placed.has(name) ? ui.components[name] : null) as ComponentType<PartProps> | null;
    return Part && <Part context={context} />;
  };
  return (
    <SeqViewerRoot>
      <div className="seqv-bar">
        {part('ViewTabs')}
        <span className="seqv-spacer" />
        {part('Search')}
        {part('Busy')}
        {part('Export')}
        {part('Issues')}
      </div>
      {part('SnapshotPlayer')}
      {part('Description')}
      {part('EmptyState')}
      {part('ErrorState')}
      <ViewSurface context={context} />
      {part('Legend')}
      {part('Tooltip')}
      {part('HoverCard')}
    </SeqViewerRoot>
  );
}

export interface SeqViewerProps {
  readonly spec?: SeqViewerSpec;
  readonly ui?: SeqViewerUISpec;
  readonly document?: SeqViewSpecDoc;
  readonly data?: SeqDataInput;
  readonly ports?: Partial<SeqViewerPorts>;
  onReady?(context: SeqViewerContext): void;
}
/**
 * A context per mount, made from the props of that mount and disposed in its cleanup (StrictMode's second mount gets its
 * own); `spec` and `ports` are read at mount. `document` and `data` go in by identity, one effect each: a new `document`
 * object with `setDocument`, a new `data` object with `setData`, and neither applies the other again (a new `data` never
 * reverts edits committed to the document). Like the element, it registers the browser's `view.exportSvg` offering a
 * download for the Export part, unless the spec registers its own.
 */
export function SeqViewer(props: SeqViewerProps) {
  const [context, setContext] = useState<SeqViewerContext | null>(null);
  const initial = useRef(props);
  initial.current = props;
  useEffect(() => {
    const p = initial.current;
    const ctx = createSeqViewer(withDownloadExport(p.spec ?? DefaultSeqViewerSpec()), {
      scheduler: windowScheduler(),
      ...p.ports,
    });
    setContext(ctx);
    p.onReady?.(ctx);
    return () => {
      setContext(null);
      ctx.dispose();
    };
  }, []);
  const { document: doc, data } = props;
  // `data` first, as the element loads them.
  useEffect(() => {
    if (context && data) context.setData(data).catch(() => {});
  }, [context, data]);
  useEffect(() => {
    if (context && doc) context.setDocument(doc).catch(() => {});
  }, [context, doc]);
  if (!context) return <div className={P.cx('seqv-root', props.ui?.className)} />;
  return (
    <SeqViewerProvider context={context} ui={props.ui}>
      <SeqViewerLayout />
    </SeqViewerProvider>
  );
}

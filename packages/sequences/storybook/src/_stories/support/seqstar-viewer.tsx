/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer in the Seq* stories: a context per case on the fixtures, shown by the default React UI
// (`@molstar/seq-viewer/react`) with the Storybook's own hint and scroll area in its slots, as a host maps them. Pan, zoom, hover and selection are the default behaviors over the engine's DOM input.

import type { DocNode, Loaders, SeqViewSpecDoc } from '@molstar/seq-core/format';
import {
  CleanParts,
  createSeqViewer,
  DefaultPresets,
  DefaultRegistry,
  DefaultSeqViewerSpec,
  type LegendMode,
  type PresetEntry,
  type RenderProps,
  SeqViewerConfig,
  type SeqViewerConfigItem,
  type SeqViewerContext,
  type SeqViewerRegistryEntry,
  type TrackActionList,
} from '@molstar/seq-viewer';
import { mountSection, type SectionSurface, type SurfaceOptions, windowScheduler } from '@molstar/seq-viewer/canvas';
import { type DomUISpec, defineSeqViewer, type SeqHoverDetail } from '@molstar/seq-viewer/dom';
import { InteractionsViewerRegistry } from '@molstar/seq-viewer/interactions/viewer';
import {
  type HintProps,
  type ScrollAreaProps,
  SeqViewerLayout,
  SeqViewerProvider,
  type SeqViewerUISpec,
} from '@molstar/seq-viewer/react';
import { fixtureLoaders, MolstarFixture } from '../../../../core/_test/format/load-fixtures.js';
import '@molstar/seq-viewer/skin.css';
import { type CSSProperties, type ReactNode, useEffect, useMemo, useRef, useState } from 'react';
import { Fold, ScrollArea, TextAction, Tooltip } from '../ui/index.js';
import { caseSource } from './seqstar-cases.js';
import f from './seqstar-format.module.css';
import { PatchBanner } from './seqstar-patch-banner.js';
import { PrismCanvas } from './seqstar-prism.js';
import s from './seqstar-viewer.module.css';

/** The viewer's chrome (`UI.Chrome`). */
export type Chrome = 'full' | 'clean';
/** The chrome every Seq* story starts from: the canvas, hover, selection, tooltip and navigator; a story brings parts back. */
export const STORY_CHROME: Chrome = 'clean';
/** The parts `clean` leaves out (`CleanParts`): `ViewTabs`, `SectionHeader`, `TrackActions`, `Legend`, `Export`, `Issues`. */
export type CleanPart = (typeof CleanParts)[number];
/** Every part `clean` leaves out, for a story about the viewer's whole UI. */
export const ALL_PARTS: readonly CleanPart[] = CleanParts;
/** A UI's slots that bring the parts back under `clean` chrome (`"default"`: the part's default); undefined for none. */
export function partSlots(parts: readonly CleanPart[] | undefined): Partial<Record<CleanPart, 'default'>> | undefined {
  if (!parts || parts.length === 0) return undefined;
  return Object.fromEntries(parts.map((p) => [p, 'default' as const]));
}

/** A compiled case's parts under `clean`: the view tabs when it has more than one view (else nothing to switch), and `extra`. */
export function caseParts(c: { readonly views: readonly unknown[] }, extra: readonly CleanPart[] = []): CleanPart[] {
  return [...(c.views.length > 1 ? (['ViewTabs'] as const) : []), ...extra];
}

/** The story's choices over the built-in presets. */
export interface ViewerOptions {
  /** The linear presets' navigator (default: theirs, `rulerRail`); another name is a registered one (`registry`). */
  readonly navigator?: NavigatorChoice | (string & {});
  /** The `msa` preset's vertical navigation (default: the map beside a scrollbar) and row zoom handles (off). */
  readonly verticalNav?: RenderProps['verticalNav'];
  readonly verticalZoom?: boolean;
  /** Caps each box section's height (CSS px; default 560): a lower cap makes rows overflow and scroll. */
  readonly maxHeight?: number;
  /** The `rows` layers' color theme (default: the document's). */
  readonly rowsColor?: RowsColor;
  /** The view shown first (default: the document's). */
  readonly view?: string;
  /** Render props over every preset's (`controlWidth` for row controls). */
  readonly render?: Partial<RenderProps>;
  /** The default UI drawn by React (`@molstar/seq-viewer/react`) or by the `<seq-viewer>` element (`/dom`). */
  readonly ui?: UiChoice;
  /** The viewer's chrome (`UI.Chrome`): `clean` (the stories' default, `STORY_CHROME`) or every control. */
  readonly chrome?: Chrome;
  /** The parts `clean` leaves out that this story needs (a story about legends asks for `Legend`), brought back by their slots. */
  readonly parts?: readonly CleanPart[];
  /** The legend (`UI.Legend`): collapsed (the default), expanded or hidden. */
  readonly legend?: LegendMode;
  /** The track menu's items (`UI.TrackActions`; default: the built-ins). */
  readonly trackActions?: TrackActionList;
  /** Registry entries after the story's own (a host's actions, such as `browserExportAction`); keep the array stable. */
  readonly registry?: readonly SeqViewerRegistryEntry[];
  /** The loaders for the document's sources (default: the format's fixtures); keep the object stable. */
  readonly loaders?: Loaders;
  /**
   * False: the fixture structure entries (hand-written outputs, load-fixtures.ts) stay out, for a story that registers
   * the Mol* adapter's format half itself (`registry`), whose formats and services have the same names.
   */
  readonly fixtureStructures?: false;
}
export type UiChoice = 'React' | 'DOM';

/** Color themes to compare on alignment rows the document's, or one set on every `rows` layer. */
export const ROWS_COLORS = {
  document: null,
  clustal_x: { theme: 'clustal_x' },
  'letter (Zappo)': { theme: 'letter' },
  'letter (Taylor)': { theme: 'letter', params: { colors: 'Taylor' } },
  identity_to_row: { theme: 'identity_to_row' },
  by_unit: { theme: 'by_unit' },
} as const;
export type RowsColor = keyof typeof ROWS_COLORS;
/** The document with every `rows` layer's color replaced (a copy; edge-rate). */
export function withRowsColor(doc: SeqViewSpecDoc, color: RowsColor | undefined): SeqViewSpecDoc {
  const theme = ROWS_COLORS[color ?? 'document'];
  if (!theme) return doc;
  const map = (n: DocNode): DocNode => ({
    ...n,
    params: n.kind === 'layer' && n.params?.type === 'rows' ? { ...n.params, color: theme } : n.params,
    children: n.children?.map(map),
  });
  return { ...doc, root: map(doc.root) };
}

/** The built-in navigators and `band` with `axis: ["labels"]`: the window drags to pan, its ends zoom. */
export const NAVIGATORS = {
  /** The ruler's bottom border, the window a thicker segment; a section without a ruler gets the `rail`. */
  rulerRail: 'rulerRail',
  /** A slim row of its own: a hairline per block, the window a segment. */
  rail: 'rail',
  /** The axis in miniature (the section's tracks but letters, ruler and rows), outside the window faded. */
  strip: 'strip',
  /** The tracks in miniature in lanes, the window an outlined box. */
  band: 'band',
  /** The labeled block strip: block names in boxes. */
  'band-labels': { height: 36, thickness: 28, lanes: true, box: true, fade: 0, axis: ['labels'] },
  none: 'none',
} as const;
export type NavigatorChoice = keyof typeof NAVIGATORS;

/** T1's presets with the story's navigator on the linear ones, its vertical navigation on `msa`, its render props. */
export function presets(o: ViewerOptions, custom?: PresetEntry['navigator']): PresetEntry[] {
  const navigator =
    custom ?? (o.navigator ? ((NAVIGATORS as Record<string, string | object>)[o.navigator] ?? o.navigator) : undefined);
  return DefaultPresets.map((p) => {
    if (p.geometry !== 'linear' || p.name === 'strip')
      return o.render ? { ...p, render: { ...p.render, ...o.render } } : p;
    const msa = p.name === 'msa' ? { verticalNav: o.verticalNav ?? 'both', verticalZoom: o.verticalZoom ?? false } : {};
    return { ...p, render: { ...p.render, ...msa, ...o.render }, navigator: navigator ?? p.navigator };
  });
}

/** Dev only: every mounted viewer and its surfaces, for benchmarks from the console (`__seqViewers`). */
export interface Mounted {
  readonly ctx: SeqViewerContext;
  readonly surfaces: Map<string, SectionSurface>;
  /** The parts this mount brings back under `clean` chrome (`ViewerOptions.parts`). */
  readonly parts?: readonly CleanPart[];
}
const mounted = new Set<Mounted>();
(globalThis as { __seqViewers?: Set<Mounted> }).__seqViewers = mounted;

/** A context over `doc` on the fixtures, made per mount; null until made. */
export function useCaseContext(
  doc: SeqViewSpecDoc,
  o: ViewerOptions = {},
  navigator?: PresetEntry['navigator'],
): Mounted | null {
  const [m, setM] = useState<Mounted | null>(null);
  const { rowsColor, view, chrome = STORY_CHROME, legend, trackActions, registry: extraRegistry, loaders } = o;
  const fixtures = o.fixtureStructures !== false;
  const key = JSON.stringify([o.navigator, o.verticalNav, o.verticalZoom, o.render, navigator]);
  const partsKey = (o.parts ?? []).join();
  // biome-ignore lint/correctness/useExhaustiveDependencies: by content (`key`): callers pass options inline.
  const presetList = useMemo(() => presets(o, navigator), [key]);
  useEffect(() => {
    const spec = DefaultSeqViewerSpec();
    const registry = [
      { ...DefaultRegistry, presets: presetList },
      ...(fixtures ? [{ svs: MolstarFixture }] : []),
      // The interactions extension, both halves in one entry case 6's fingerprint as KLIFS draws it.
      InteractionsViewerRegistry,
      ...(extraRegistry ?? []),
    ];
    const config: [SeqViewerConfigItem, unknown][] = [[SeqViewerConfig.UI.Chrome, chrome]];
    if (trackActions) config.push([SeqViewerConfig.UI.TrackActions, trackActions]);
    if (legend) config.push([SeqViewerConfig.UI.Legend, legend]);
    const ctx = createSeqViewer(
      { ...spec, registry, config },
      {
        loaders: loaders ?? fixtureLoaders(),
        scheduler: windowScheduler(),
      },
    );
    const next: Mounted = {
      ctx,
      surfaces: new Map(),
      parts: partsKey ? (partsKey.split(',') as CleanPart[]) : undefined,
    };
    mounted.add(next);
    void ctx.setDocument(withRowsColor(doc, rowsColor)).then(() => !ctx.disposed && view && ctx.views.setCurrent(view));
    setM(next);
    return () => {
      mounted.delete(next);
      setM(null);
      ctx.dispose();
    };
  }, [doc, rowsColor, view, presetList, chrome, partsKey, legend, trackActions, extraRegistry, loaders, fixtures]);
  return m;
}

/** The Storybook's tooltip as the viewer's `Hint` (the trigger a span, so any child works). */
function StoryHint({ label, children }: HintProps) {
  return (
    <Tooltip content={label}>
      <span className={s.hint}>{children}</span>
    </Tooltip>
  );
}
function StoryScrollArea({ children, viewportRef, className }: ScrollAreaProps) {
  return (
    <ScrollArea className={className} viewportRef={viewportRef as never}>
      {children}
    </ScrollArea>
  );
}

/** The slots a host maps, here the Storybook's hint and scroll area; the surfaces recorded for `__seqViewers`. */
export function useStoryUi(
  m: Mounted | null,
  extra?: Partial<SeqViewerUISpec>,
  surfaceOptions?: SurfaceOptions,
): SeqViewerUISpec {
  return useMemo(
    () => ({
      className: s.viewer,
      ...extra,
      components: { Hint: StoryHint, ScrollArea: StoryScrollArea, ...partSlots(m?.parts), ...extra?.components },
      surface: (ctx, section, host) => {
        const surface = mountSection(ctx, section, host, surfaceOptions) as SectionSurface;
        m?.surfaces.set(section, surface);
        return {
          element: surface.element,
          dispose: () => {
            m?.surfaces.delete(section);
            surface.dispose();
          },
        };
      },
    }),
    [m, extra, surfaceOptions],
  );
}

/**
 * The vanilla UI (`@molstar/seq-viewer/dom`): a `<seq-viewer>` on the story's context (the host's: it outlives the
 * element), themed from outside (the story frame's `--seq-*` mapping inherits into its shadow root; the panel class on
 * the element itself), with the events it re-dispatches shown under it.
 */
function DomViewer({ m, ui, slots }: { m: Mounted; ui: SeqViewerUISpec; slots?: DomUISpec['slots'] }) {
  const host = useRef<HTMLDivElement>(null);
  const log = useRef<HTMLDivElement>(null);
  useEffect(() => {
    defineSeqViewer();
    const el = document.createElement('seq-viewer');
    el.className = s.viewer;
    el.ui = { surface: ui.surface, slots };
    el.context = m.ctx;
    const last: Record<string, string> = {};
    const show = (type: string, text: string) => {
      last[type] = text;
      if (log.current)
        log.current.textContent = `Events: ${Object.entries(last)
          .map(([k, v]) => `${k} ${v}`)
          .join(' · ')}`;
    };
    const events: [string, (e: Event) => void][] = [
      [
        'seq-hover',
        (e) => {
          const d = (e as CustomEvent<SeqHoverDetail>).detail;
          show('seq-hover', `(${d.owner}, ${d.section ?? '–'}, ${d.track ?? '–'})`);
        },
      ],
      [
        'seq-selection',
        (e) => show('seq-selection', `(origin ${(e as CustomEvent<{ origin: string }>).detail.origin})`),
      ],
      ['seq-snapshot', (e) => show('seq-snapshot', (e as CustomEvent<{ ref: string | null }>).detail.ref ?? 'none')],
    ];
    const box = host.current as HTMLDivElement;
    // Listened for on the story's own element: the events bubble out of the shadow root.
    for (const [type, fn] of events) box.addEventListener(type, fn);
    box.append(el);
    return () => {
      for (const [type, fn] of events) box.removeEventListener(type, fn);
      el.remove();
    };
  }, [m, ui, slots]);
  return (
    <div className={f.stack}>
      <div ref={host} />
      <div ref={log} className={f.code}>
        Events: hover, select or play a snapshot to see the element's events.
      </div>
    </div>
  );
}

/** The viewer on a canvas in the UI theme (the toolbar's): the default of every card. */
const uiCanvas = (viewer: ReactNode) => <PrismCanvas>{viewer}</PrismCanvas>;

/** A case's viewer in a card: every part of the default UI over a context on the fixtures, by React or the element. */
export function ViewerCard({
  doc,
  label,
  options = {},
  ui,
  domSlots,
  canvas,
  children,
}: {
  doc: SeqViewSpecDoc;
  label: ReactNode;
  options?: ViewerOptions;
  ui?: Partial<SeqViewerUISpec>;
  /** The DOM variant's slots (`ui.components` are React's). */
  domSlots?: DomUISpec['slots'];
  /**
   * Wraps the viewer alone (not the card's title or children): its canvas. Default: a `PrismCanvas` in the UI theme; a
   * story comparing themes passes one pinned to a theme (`<PrismCanvas theme="dark">`).
   */
  canvas?: (viewer: ReactNode) => ReactNode;
  children?: ReactNode;
}) {
  const m = useCaseContext(doc, options);
  const spec = useStoryUi(m, ui);
  const style = { '--seq-section-max-height': `${options.maxHeight ?? 560}px` } as CSSProperties;
  const dom = options.ui === 'DOM';
  const slots = useMemo(() => ({ ...partSlots(m?.parts), ...domSlots }), [m, domSlots]);
  return (
    <section className={f.card} style={style}>
      <h2 className={f.heading}>
        {label}
        {dom && <span className={f.code}> · {'<seq-viewer>'}</span>}
      </h2>
      {m ? (
        <SeqViewerProvider context={m.ctx} ui={spec}>
          {(canvas ?? uiCanvas)(dom ? <DomViewer m={m} ui={spec} slots={slots} /> : <SeqViewerLayout />)}
          {children}
          <PatchBanner />
        </SeqViewerProvider>
      ) : (
        <div className={f.status}>Loading…</div>
      )}
    </section>
  );
}

/**
 * The case's builder function, the source the story runs (cases.ts, by its region markers), behind a toggle; or a
 * story's own `source` (its code and what it uses).
 */
export function BuilderCode({ name, source }: { name?: string; source?: { code: string; uses: string } }) {
  const [open, setOpen] = useState(false);
  const { code, helpers, fixtures } = useMemo(
    () => (source ? { code: source.code, helpers: [], fixtures: [] } : caseSource(name ?? '')),
    [name, source],
  );
  const uses = [
    source?.uses,
    helpers.length > 0 && `from cases.ts: ${helpers.join(', ')}`,
    fixtures.length > 0 && `from the fixtures (load-fixtures.ts): ${fixtures.join(', ')}`,
  ].filter(Boolean);
  return (
    <div className={f.stack}>
      <div className={f.row}>
        <TextAction tone="quiet" onClick={() => setOpen((o) => !o)}>
          {open ? 'Hide the builder code' : 'Builder code'}
        </TextAction>
      </div>
      <Fold open={open}>
        <div className={f.stack}>
          <div className={f.code}>
            {uses.length > 0 ? `Uses ${uses.join('; ')}.` : 'Uses only the builder (@molstar/seq-core/builder).'}
          </div>
          <ScrollArea className={f.box} orientation="both">
            <pre className={f.pre}>{open ? code : ''}</pre>
          </ScrollArea>
        </div>
      </Fold>
    </div>
  );
}

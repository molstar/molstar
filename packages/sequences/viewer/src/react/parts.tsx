/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The parts of seqspec:interaction#the-default-ui that are DOM chrome, and the default primitives (`Hint`,
// `ScrollArea`, `Menu`). Each part receives the context and reads the smallest state it renders; it holds no state
// beyond leaf state (a menu open). Each is also a component a host places anywhere (seqspec:interaction#slots-and-chrome):
// under the provider it takes its slots and primitives, without one the defaults. The structure, class names, text and `data-seqv-key`s are the DOM parts'
// (src/dom/parts.ts), from the same helpers (src/chrome.ts, src/dom/shared.ts); parity is tested.

import {
  type ButtonHTMLAttributes,
  cloneElement,
  type ReactElement,
  type RefObject,
  useCallback,
  useEffect,
  useId,
  useLayoutEffect,
  useRef,
  useState,
} from 'react';
import { createPortal, flushSync } from 'react-dom';
import {
  type ColorsDrill,
  dataSpan,
  editRange,
  type IssueGroup,
  issueBadge,
  issueList,
  type LegendGroup,
  type LegendRange,
  legendExpanded,
  legendGroups,
  legendShown,
  legendShownGroups,
  legendTitle,
  plainDescription,
  rampEnd,
  resetRange,
  rootErrors,
  type SectionHighlightEntry,
  type SectionSwitch,
  sameMenu,
  schemeItems,
  sectionChrome,
  sectionChromeStates,
  sectionHighlightEntries,
  trackMenuModel,
} from '../chrome.js';
import { cx, h, landFocus, later } from '../dom/core.js';
import {
  busyShown,
  CROSS,
  handOff,
  hintBehavior,
  legendGradient,
  legendGrid,
  legendSizeStyle,
  PALETTE,
  PIN,
  type Popover,
  placeTooltip,
  popover,
  sizeCaption,
  tipSpot,
} from '../dom/shared.js';
import {
  type LegendRamp,
  type LegendSize,
  type SearchFn,
  type SeqIssue,
  SeqViewerConfig,
  type SeqViewerContext,
  type Unsubscribable,
} from '../index.js';
import { sameData } from '../same.js';

import { CLEAR_KEPT_HINT, KEEP_HINT, SearchBox } from '../search-box.js';
import {
  type HintProps,
  type MenuProps,
  type PartProps,
  type ScrollAreaProps,
  type SectionProps,
  type TrackProps,
  useDerived,
  useLanding,
  useScopeMaybe,
  useSeqState,
} from './hooks.js';

export { cx };

/**
 * A callback ref that runs `attach`'s returned cleanup when React detaches the element. React 19 would take a cleanup
 * returned from the ref, React 18 ignores it and calls the ref with `null` instead, so this keeps the cleanup and runs
 * it on that `null` (both detach before the node leaves the document, in the commit of the removal).
 */
function useAttach<T extends HTMLElement>(attach: (el: T) => (() => void) | undefined) {
  const off = useRef<(() => void) | undefined>(undefined);
  return useCallback(
    (el: T | null) => {
      const prev = off.current;
      off.current = undefined;
      prev?.();
      if (el) off.current = attach(el);
    },
    [attach],
  );
}

/**
 * A callback ref that lands focus (seqspec:interaction#accessibility, `landFocus`) when React removes the element with
 * focus inside it: for the root of a part that comes and goes (`ref` of a node that may not be rendered); the landing
 * runs before the removal.
 */
function useLand(near?: (el: HTMLElement) => Node) {
  return useAttach<HTMLElement>(useCallback((el) => () => landFocus(el, near?.(el) ?? el), [near]));
}

/**
 * An accessible hint on hover and on keyboard focus, never in the pointer's way (`hintBehavior`); the control keeps its
 * own accessible name.
 */
export function Hint({ label, children }: HintProps) {
  const bubble = useRef<HTMLSpanElement>(null);
  const behave = useAttach<HTMLSpanElement>(
    useCallback((el) => (bubble.current ? hintBehavior(el, bubble.current) : undefined), []),
  );
  return (
    <span ref={behave} className="seqv-hint">
      {children}
      <span ref={bubble} className="seqv-hint-bubble" aria-hidden="true">
        {label}
      </span>
    </span>
  );
}

/** A plain region: the default UI scrolls nothing in the DOM (a host maps this to its own scroll area). */
export function ScrollArea({ children, viewportRef, className }: ScrollAreaProps) {
  return (
    <div ref={viewportRef as RefObject<HTMLDivElement>} className={cx('seqv-scroll', className)}>
      {children}
    </div>
  );
}

/**
 * `popover` (dom/shared.ts: the trigger, Escape, Tab, a press outside, focus leaving, the arrows with `roving`, and
 * focus back to the trigger) over a trigger React renders and a list it fills but does not own: the popover puts the
 * list into the tree and out of it, so the list is a detached node and its content a portal into it. `root` holds the
 * trigger as its first child. `open` is the owner's: the popover asks `onOpenChange` and follows `open`
 * (seqspec:interaction#slots-and-chrome). Each commit checks the trigger: one React replaced (another element type)
 * gets a new popover, so its clicks and ARIA state are never a gone node's. A StrictMode remount disposes the popover
 * and the next commit's check makes it again.
 */
function usePopover(
  root: RefObject<HTMLElement | null>,
  list: HTMLElement,
  o: { roving: boolean; top?: boolean; open: boolean; onOpenChange(open: boolean): void },
) {
  const ask = useRef(o.onOpenChange);
  const pop = useRef<{ p: Popover; trigger: Element } | null>(null);
  useLayoutEffect(() => {
    ask.current = o.onOpenChange;
    const r = root.current;
    const trigger = r?.firstElementChild;
    if (pop.current && pop.current.trigger !== trigger) {
      pop.current.p.dispose();
      pop.current = null;
    }
    if (!pop.current && r && trigger instanceof HTMLElement)
      pop.current = {
        p: popover(r, trigger, list, { roving: o.roving, top: o.top, ask: (open) => ask.current(open) }),
        trigger,
      };
    if (o.open) pop.current?.p.open();
    else pop.current?.p.close();
  });
  useLayoutEffect(
    () => () => {
      pop.current?.p.dispose();
      pop.current = null;
    },
    [],
  );
}

/**
 * The trigger, then a list of buttons while open (a detached list the popover moves in; see `usePopover`): opening
 * focuses the first item, the arrows move between items, Escape, Tab, a choice and a press outside ask `onOpenChange`,
 * and closing with focus inside returns it to the trigger.
 */
export function Menu({ trigger, items, open, onOpenChange }: MenuProps) {
  const root = useRef<HTMLSpanElement>(null);
  const id = useId();
  const [list] = useState(() => h('span', { id: `seqv-menu-${id}`, role: 'menu', class: 'seqv-menu-list' }));
  useLanding(root);
  usePopover(root, list, { roving: true, top: true, open, onOpenChange });
  const t = trigger as ReactElement<ButtonHTMLAttributes<HTMLButtonElement>>;
  return (
    <span ref={root} className="seqv-menu">
      {cloneElement(t, { 'aria-haspopup': 'menu' })}
      {createPortal(
        items.map((i) => (
          <button
            key={i.key}
            type="button"
            {...(i.checked === undefined ? { role: 'menuitem' } : { role: 'menuitemradio', 'aria-checked': i.checked })}
            className={cx('seqv-menu-item', i.danger && 'seqv-danger')}
            data-seqv-key={`item:${i.key}`}
            tabIndex={-1}
            onClick={() => {
              onOpenChange(false);
              i.run();
            }}
          >
            {i.label}
          </button>
        )),
        list,
      )}
    </span>
  );
}

const usePrimitives = () => {
  // Tolerant of no provider: every part also stands alone, with the default primitives.
  const c = useScopeMaybe()?.ui.components;
  return { Hint: c?.Hint ?? Hint, Menu: c?.Menu ?? Menu, ScrollArea: c?.ScrollArea ?? ScrollArea };
};

export function ViewTabs({ context: ctx }: PartProps) {
  const views = useSeqState(ctx.views.list);
  const current = useSeqState(ctx.views.current);
  const land = useLand();
  if (views.length < 2) return null;
  const go = (at: number, el: HTMLElement) => {
    const v = views[(at + views.length) % views.length];
    ctx.actions.run('view.switch', v.ref);
    (el.parentElement?.children[views.indexOf(v)] as HTMLElement | undefined)?.focus();
  };
  return (
    <div ref={land} role="tablist" aria-label="Views" className="seqv-tabs">
      {views.map((v, i) => (
        <button
          key={v.ref}
          type="button"
          role="tab"
          aria-selected={v.ref === current}
          tabIndex={v.ref === current ? 0 : -1}
          className="seqv-tab"
          data-seqv-key={`tab:${v.ref}`}
          data-seqv-group="tabs"
          onClick={() => ctx.actions.run('view.switch', v.ref)}
          onKeyDown={(e) => {
            if (e.key === 'ArrowRight' || e.key === 'ArrowLeft')
              go(i + (e.key === 'ArrowRight' ? 1 : -1), e.currentTarget);
          }}
        >
          <span className="seqv-tab-label">{v.label}</span>
        </button>
      ))}
    </div>
  );
}

/** Two or three choices as pressed buttons; more as a menu (its open state is the switch's own). */
function Switch({ sw, run }: { sw: SectionSwitch; run: (value: string | null) => void }) {
  const { Menu } = usePrimitives();
  const [open, setOpen] = useState(false);
  const shown = sw.options.find((o) => o.value === sw.value)?.label ?? sw.value ?? 'default';
  if (sw.options.length > 3)
    return (
      <Menu
        trigger={
          <button
            type="button"
            className="seqv-switch-trigger"
            aria-label={`${sw.name}: ${shown}`}
            data-seqv-key={`switch:${sw.action}`}
          >
            <span className="seqv-quiet">{sw.name}</span> {shown}
          </button>
        }
        items={sw.options.map((o) => ({
          key: String(o.value),
          label: o.label,
          checked: o.value === sw.value,
          run: () => run(o.value),
        }))}
        open={open}
        onOpenChange={setOpen}
      />
    );
  return (
    <fieldset aria-label={sw.name} className="seqv-switch">
      <span className="seqv-quiet">{sw.name}</span>
      {sw.options.map((o) => (
        <button
          key={String(o.value)}
          type="button"
          aria-pressed={o.value === sw.value}
          className="seqv-choice"
          data-seqv-key={`choice:${sw.action}:${o.value}`}
          data-seqv-group={`choices:${sw.action}`}
          onClick={() => run(o.value)}
        >
          {o.label}
        </button>
      ))}
    </fieldset>
  );
}

/**
 * A section's label, its filter chip, then its `HiddenChips` and `SectionSwitches` (the provider's slots, else the
 * defaults); nothing when it would be empty.
 */
export function SectionHeader({ context: ctx, section }: SectionProps) {
  const c = useScopeMaybe()?.ui.components;
  const [Chips, Switches] = c ? [c.HiddenChips, c.SectionSwitches] : [HiddenChips, SectionSwitches];
  // What the row shows of its own; the switches are their part's (a change of one renders only that part).
  // The hidden tracks' chips sit in this row, left of the switches: a header with only those is still a header.
  const info = useDerived(
    sectionChromeStates(ctx),
    () => {
      const i = sectionChrome(ctx, section);
      return {
        label: i.label,
        filter: i.filter,
        switches: i.switches.length > 0,
        hidden: !!Chips && i.hidden.length > 0,
      };
    },
    { deps: [ctx, section, Chips] },
  );
  const hidden = info.hidden;
  const land = useLand();
  if (!info.label && !info.filter && !(Switches && info.switches) && !hidden) return null;
  return (
    <div ref={land} className="seqv-section-header">
      {info.label && <span className="seqv-section-label">{info.label}</span>}
      {info.filter && (
        <button
          type="button"
          className="seqv-chip"
          aria-label={`Show all positions (now only ${info.filter})`}
          data-seqv-key="filter"
          onClick={() => ctx.actions.run('section.clearFilter', section)}
        >
          Only {info.filter} <span aria-hidden="true">✕</span>
        </button>
      )}
      <span className="seqv-spacer" />
      {Chips && <Chips context={ctx} section={section} />}
      {Switches && <Switches context={ctx} section={section} />}
    </div>
  );
}

/** A section's switches, each running its `section.*` action: in its header, or wherever a host places them. */
export function SectionSwitches({ context: ctx, section }: SectionProps) {
  const switches = useDerived(sectionChromeStates(ctx), () => sectionChrome(ctx, section).switches, {
    deps: [ctx, section],
  });
  const land = useLand();
  if (switches.length === 0) return null;
  return (
    <div ref={land} className="seqv-switches">
      {switches.map((sw) => (
        <Switch key={sw.action} sw={sw} run={(value) => ctx.actions.run(sw.action, { section, value })} />
      ))}
    </div>
  );
}

/**
 * The menu at a track header: filter view, focus, hide, remove. It holds its menu's `open` (leaf state) and opens it on
 * `events.menu` for its track (Shift+F10 or ContextMenu on the keyboard cursor, seqspec:viewer#behaviors): a host's
 * `Menu` gets it as its `open`.
 */
export function TrackActions({ context: ctx, section, track, row }: TrackProps) {
  const { Menu, Hint } = usePrimitives();
  const [open, setOpen] = useState(false);
  const [layer, setLayer] = useState<string | null>(null);
  // Closed, the menu opens at its top again, not in its "Colors" submenu.
  const onOpenChange = (o: boolean) => {
    setOpen(o);
    if (!o) setLayer(null);
  };
  // The "Colors" submenu: an item chosen closed the menu; the close is committed first, so the menu opens again with
  // the other list and focus on its first item, as the DOM UI's does.
  const drill: ColorsDrill = {
    layer,
    open: (l) => {
      flushSync(() => onOpenChange(false));
      setLayer(l);
      setOpen(true);
    },
  };
  // "Retry" comes with a failed region: the menu is read again on the streams' version, but only while it is open (the
  // third state is the second again until then). A new view: which layers have schemes, and the one picked.
  const menu = useDerived(
    [ctx.document.shown, ctx.ui.expanded, ctx.graph.assembled, open ? ctx.streams.version : ctx.ui.expanded],
    () => trackMenuModel(ctx, section, track, drill, row?.() ?? -1),
    // Opening reads the menu again for the line it opens on (a split row's own items), as the DOM part does.
    { compare: sameMenu, deps: [ctx, section, track, row, layer, open] },
  );
  const has = menu.items.length > 0;
  useEffect(() => {
    const sub = ctx.events.menu.subscribe((e) => {
      if (has && e.section === section && e.track === track) setOpen(true);
    });
    return () => sub.unsubscribe();
  }, [ctx, section, track, has]);
  if (!has) return null;
  const label = `Actions for ${menu.label}`;
  return (
    <Hint label={label}>
      <Menu
        trigger={
          <button type="button" className="seqv-icon" aria-label={label} data-seqv-key="track-actions">
            ⋯
          </button>
        }
        items={menu.items}
        open={open}
        onOpenChange={onOpenChange}
      />
    </Hint>
  );
}

/** One chip: focus on it when it goes lands on the next chip, else the previous, else the surface. */
function Chip({ children, ...p }: ButtonHTMLAttributes<HTMLButtonElement> & { 'data-seqv-key': string }) {
  const land = useLand();
  return (
    <button ref={land} type="button" className="seqv-chip" data-seqv-group="chips" {...p}>
      {children}
    </button>
  );
}
export function HiddenChips({ context: ctx, section }: SectionProps) {
  const hidden = useDerived(sectionChromeStates(ctx), () => sectionChrome(ctx, section).hidden, {
    deps: [ctx, section],
  });
  if (hidden.length === 0) return null;
  return (
    <fieldset aria-label="Hidden tracks" className="seqv-chips">
      {hidden.map((t) => (
        <Chip key={t.ref} data-seqv-key={`chip:${t.ref}`} onClick={() => ctx.actions.run('track.show', t.ref)}>
          Show {t.label}
        </Chip>
      ))}
    </fieldset>
  );
}

/**
 * The hover's label, placed at its owner's anchor (`placeTooltip`): the pointer's, or the keyboard cursor's cell, which
 * it follows while the viewport moves under it. It moves by subscriptions that write styles; it renders again only when
 * the label changes.
 */
/** No UI around a canvas overlay: it is never placed (it renders nothing). */
export const NO_SURFACES: ReadonlyMap<string, HTMLElement> = new Map();
export function Tooltip({ context: ctx }: PartProps) {
  // A canvas overlay, placed at its UI's surfaces: nothing without a provider (as the DOM's without a UI).
  const scope = useScopeMaybe();
  const surfaces = (scope?.surfaces ?? NO_SURFACES) as Map<string, HTMLElement>;
  // The same words again (the hover refreshed under a new document) render nothing: equal by content.
  const label = useDerived([ctx.labels.current], () => ctx.labels.current.value, {
    compare: sameData,
    deps: [ctx],
  });
  const el = useRef<HTMLDivElement>(null);
  const spot = useRef(tipSpot());
  const shown = label !== null;
  // The tip is a new node when a label comes after none: its subscriptions are made again.
  useLayoutEffect(() => {
    const tip = el.current;
    if (!tip) return;
    const at = spot.current;
    [at.x, at.y, at.shown] = [Number.NaN, Number.NaN, true];
    const place = () => placeTooltip(ctx, surfaces, tip, at);
    let follow: Unsubscribable | undefined;
    const subs = [
      later(ctx.pointer.version, place),
      // The replay is the first run: it follows the cursor's section if there is one, and places.
      ctx.ui.keyCursor.subscribe((c) => {
        follow?.unsubscribe();
        const live = c && ctx.viewports.get(c.section)?.live.version;
        follow = live ? later(live, place) : undefined;
        place();
      }),
    ];
    return () => {
      for (const s of subs) s.unsubscribe();
      follow?.unsubscribe();
    };
  }, [ctx, surfaces, shown]);
  // New words change its size: placed again once they are in.
  useLayoutEffect(() => {
    if (label && el.current) placeTooltip(ctx, surfaces, el.current, spot.current);
  }, [ctx, surfaces, label]);
  if (!label || !scope) return null;
  return (
    <div ref={el} role="tooltip" className="seqv-tooltip">
      {label.title && <div className="seqv-tooltip-title">{label.title}</div>}
      {label.lines.map((l) => (
        <div key={`${l.key}:${l.value}`} className="seqv-tooltip-line">
          <span className="seqv-quiet">{l.key}</span> {l.value}
        </div>
      ))}
      {label.legend && (
        <div className="seqv-tooltip-legend">
          {label.legend.items.map((i) => (
            <span key={`${i.label}:${i.color}`} className="seqv-legend-item">
              <span className="seqv-swatch" style={{ background: i.color ?? 'var(--seq-ink)' }} />
              {i.label}
            </span>
          ))}
          {label.legend.ramps.map((r) => (
            <Ramp key={`${r.lo}:${r.hi}:${r.caption}`} ramp={r} />
          ))}
          {label.legend.more > 0 && <span className="seqv-quiet">+{label.legend.more} more</span>}
        </div>
      )}
    </div>
  );
}

/**
 * On the canvas's grid: each group's name in the label column, its items from the data area's left (`legendGrid`); a
 * group past `UI.LegendItems` ends in "+n more" (`ui.legendOpen`), as does the legend past `UI.LegendGroups` groups; a
 * continuous palette is a ramp with its ends and caption, a sample size a dot. A group with a color scheme picker has
 * its title as the picker's trigger (seqspec:interaction#color-schemes-at-runtime).
 */
export function Legend({ context: ctx }: PartProps) {
  const mode = useScopeMaybe()?.ui.legend;
  const { groups, kept, expanded } = useDerived(
    [ctx.graph.assembled, ctx.ui.legendOpen, ctx.ui.legendExpanded],
    () => ({
      groups: legendGroups(ctx),
      kept: sectionHighlightEntries(ctx),
      // The open groups are in the value, so a toggle renders (the groups alone are one memoized list per view; the
      // groups, the kept entries and the open set are each the same object until they change).
      open: ctx.ui.legendOpen.value,
      expanded: legendExpanded(ctx, mode),
    }),
    { deps: [ctx, mode] },
  );
  return groups.length === 0 && kept.length === 0 ? null : (
    <LegendBox context={ctx} groups={groups} kept={kept} expanded={expanded} />
  );
}
function LegendBox({
  context: ctx,
  groups,
  kept,
  expanded,
}: PartProps & { groups: readonly LegendGroup[]; kept: readonly SectionHighlightEntry[]; expanded: boolean }) {
  const mode = useScopeMaybe()?.ui.legend;
  const root = useRef<HTMLDetailsElement>(null);
  const first = groups[0]?.section ?? kept[0]?.section;
  useLanding(root);
  useLayoutEffect(
    () => (first !== undefined && root.current ? legendGrid(ctx, first, root.current) : undefined),
    [ctx, first],
  );
  const shown = legendShownGroups(ctx, groups);
  return (
    <details
      ref={root}
      className="seqv-legend"
      open={expanded}
      // The person's toggle is UI state (`ui.legendExpanded`), as `view.legend` is.
      onToggle={(e) => {
        const o = e.currentTarget.open;
        if (o !== legendExpanded(ctx, mode)) ctx.ui.setLegendExpanded(o);
      }}
    >
      <summary data-seqv-key="legend">
        <span className="seqv-chevron" aria-hidden="true" />
        Legend
      </summary>
      {shown.groups.map((g) => (
        <LegendRow key={g.key} context={ctx} group={g} />
      ))}
      {shown.toggle && <LegendRest context={ctx} section={groups[0].section} name={shown.key} more={shown.more} />}
      {kept.length > 0 && <KeptHighlights context={ctx} entries={kept} />}
    </details>
  );
}
/**
 * "Highlights": each entry of the sections' `highlights` (a kept result set, bookmarks), its style's swatch, its label,
 * and a ✕ that removes it (`highlight.remove`, one patch: Undo brings it back).
 */
function KeptHighlights({ context: ctx, entries }: PartProps & { entries: readonly SectionHighlightEntry[] }) {
  const { Hint } = usePrimitives();
  const el = useRef<HTMLDivElement>(null);
  const section = entries[0].section;
  useLayoutEffect(() => (el.current ? legendGrid(ctx, section, el.current) : undefined), [ctx, section]);
  return (
    <div ref={el} className="seqv-legend-group" data-section={section}>
      <span className="seqv-legend-track">
        <span className="seqv-legend-name">Highlights</span>
      </span>
      <span className="seqv-legend-items">
        {entries.map((e) => (
          <span key={e.key} className="seqv-legend-item">
            <span className="seqv-swatch seqv-swatch-set" data-style={e.style} aria-hidden="true" />
            <span>{e.label}</span>
            <Hint label={`Remove ${e.label} from the highlights`}>
              <button
                type="button"
                className="seqv-legend-remove"
                aria-label={`Remove ${e.label}`}
                data-seqv-key={`hl:${e.key}`}
                onClick={() => ctx.actions.run('highlight.remove', { section: e.section, data: e.data }, 'pointer')}
              >
                ✕
              </button>
            </Hint>
          </span>
        ))}
      </span>
    </div>
  );
}
function More({ context: ctx, name, more }: PartProps & { name: string; more: number }) {
  return (
    <button
      type="button"
      className="seqv-legend-more"
      aria-expanded={more === 0}
      data-seqv-key={`more:${name}`}
      onClick={() => ctx.ui.toggleLegend(name)}
    >
      {more > 0 ? `+${more} more` : 'fewer'}
    </button>
  );
}
/** The legend's own "+n more" (groups past `UI.LegendGroups`), on the first group's grid. */
function LegendRest({
  context: ctx,
  section,
  name,
  more,
}: PartProps & { section: string; name: string; more: number }) {
  const el = useRef<HTMLDivElement>(null);
  useLayoutEffect(() => (el.current ? legendGrid(ctx, section, el.current) : undefined), [ctx, section]);
  return (
    <div ref={el} className="seqv-legend-group" data-section={section}>
      <span className="seqv-legend-track" />
      <span className="seqv-legend-items">
        <More context={ctx} name={name} more={more} />
      </span>
    </div>
  );
}
/** The legend title's palette icon (`PALETTE`), an image labeled with the scheme shown. */
function PaletteIcon({ label }: { label: string }) {
  return (
    <svg className="seqv-legend-icon" role="img" aria-label={label} viewBox="0 0 16 16">
      <path d={PALETTE.outline} fill="none" stroke="currentColor" strokeWidth={1.3} />
      {PALETTE.dots.map(([cx, cy]) => (
        <circle key={`${cx}:${cy}`} cx={cx} cy={cy} r={1.1} fill="currentColor" />
      ))}
    </svg>
  );
}
/**
 * A group's title: its name, or with a picker the trigger of its scheme menu (the name, dotted underline, then a palette
 * icon on hover or focus); in a hint with the full name only while the name is cut (`legendGrid`'s `onCut`).
 */
function LegendTitle({ context: ctx, group: g, cut }: PartProps & { group: LegendGroup; cut: boolean }) {
  const { Menu, Hint } = usePrimitives();
  const [open, setOpen] = useState(false);
  const t = legendTitle(g);
  const body =
    g.colors && t.colors ? (
      <Menu
        trigger={
          <button type="button" className="seqv-legend-title" data-seqv-key={`colors:${g.key}`}>
            <span className="seqv-legend-name">{t.name}</span>
            <PaletteIcon label={t.colors} />
          </button>
        }
        items={schemeItems(ctx, g.colors.layers, g.colors)}
        open={open}
        onOpenChange={setOpen}
      />
    ) : (
      <span className="seqv-legend-name">{t.name}</span>
    );
  return <span className="seqv-legend-track">{cut ? <Hint label={t.name}>{body}</Hint> : body}</span>;
}
function LegendRow({ context: ctx, group: g }: PartProps & { group: LegendGroup }) {
  const el = useRef<HTMLDivElement>(null);
  const [cut, setCut] = useState(false);
  useLayoutEffect(() => (el.current ? legendGrid(ctx, g.section, el.current, setCut) : undefined), [ctx, g.section]);
  const { items, more } = legendShown(ctx, g);
  const expanded = ctx.ui.legendOpen.value.has(g.key);
  const toggle = more > 0 || (expanded && g.items.length > ctx.config.get(SeqViewerConfig.UI.LegendItems));
  return (
    <div ref={el} className="seqv-legend-group" data-section={g.section}>
      <LegendTitle context={ctx} group={g} cut={cut} />
      <span className="seqv-legend-items">
        {items.map((i) => (
          <span key={`${i.label}:${i.color}`} className="seqv-legend-item">
            <span className="seqv-swatch" style={{ background: i.color ?? 'var(--seq-ink)' }} />
            {i.label}
          </span>
        ))}
        {g.ramps.map((r) =>
          g.range ? (
            <RangeRamp key={`${r.lo}:${r.hi}:${r.caption}`} context={ctx} group={g} ramp={r} />
          ) : (
            <Ramp key={`${r.lo}:${r.hi}:${r.caption}`} ramp={r} />
          ),
        )}
        {g.sizes.map((s, i) => (
          <SizeSample key={`${s.channel}:${s.field}:${s.label}`} size={s} caption={sizeCaption(g.sizes, i)} />
        ))}
        {toggle && <More context={ctx} name={g.key} more={more} />}
      </span>
    </div>
  );
}
/** A ramp with its range editable (the DOM part's `rangeRamp`): fields for the ends, the data's span marked, a reset. */
function RangeRamp({ context: ctx, group: g, ramp: r }: PartProps & { group: LegendGroup; ramp: LegendRamp }) {
  const range = g.range as LegendRange;
  const span = dataSpan(range);
  const field = (end: 'lo' | 'hi') => (
    <input
      key={`${end}:${range[end]}`}
      className="seqv-range-input"
      type="text"
      inputMode="decimal"
      aria-label={`${end === 'lo' ? 'Lowest' : 'Highest'} value of the color scale`}
      data-seqv-key={`range-${end}:${g.key}`}
      defaultValue={rampEnd(range[end])}
      onBlur={(e) => {
        if (e.currentTarget.value !== rampEnd(range[end]) && !editRange(ctx, range, end, e.currentTarget.value))
          e.currentTarget.value = rampEnd(range[end]);
      }}
      onKeyDown={(e) => e.key === 'Enter' && e.currentTarget.blur()}
    />
  );
  return (
    <span className="seqv-ramp seqv-ramp-edit">
      {field('lo')}
      <span className="seqv-ramp-track">
        <span className="seqv-ramp-bar" aria-hidden="true" style={{ background: legendGradient(r) }} />
        {span && (
          <span
            className="seqv-ramp-data"
            aria-hidden="true"
            style={{ left: `${span[0] * 100}%`, width: `${(span[1] - span[0]) * 100}%` }}
          />
        )}
      </span>
      {field('hi')}
      {range.edited && (
        <button
          type="button"
          className="seqv-legend-more"
          aria-label="Reset the color scale"
          data-seqv-key={`range-reset:${g.key}`}
          onClick={() => resetRange(ctx, range)}
        >
          reset
        </button>
      )}
      {r.caption && <span className="seqv-quiet">{r.caption}</span>}
    </span>
  );
}
function Ramp({ ramp: r }: { ramp: LegendRamp }) {
  return (
    <span className="seqv-ramp">
      <span className="seqv-quiet">{rampEnd(r.lo)}</span>
      <span className="seqv-ramp-bar" aria-hidden="true" style={{ background: legendGradient(r) }} />
      <span className="seqv-quiet">{rampEnd(r.hi)}</span>
      {r.caption && <span className="seqv-quiet">{r.caption}</span>}
    </span>
  );
}
function SizeSample({ size: s, caption }: { size: LegendSize; caption: string | null }) {
  return (
    <>
      {caption !== null && <span className="seqv-quiet">{caption}</span>}
      <span className="seqv-legend-item">
        <span className="seqv-size" aria-hidden="true" style={legendSizeStyle(s)} />
        {s.label}
      </span>
    </>
  );
}

export function SnapshotPlayer({ context: ctx }: PartProps) {
  const { Hint } = usePrimitives();
  const list = useSeqState(ctx.snapshots.list);
  const current = useSeqState(ctx.snapshots.current);
  const playing = useSeqState(ctx.snapshots.playing);
  const land = useLand();
  if (list.length === 0) return null;
  const at = list.findIndex((s) => s.ref === current);
  /** An unusable control stays focusable: `aria-disabled`, never `disabled`, and does nothing when pressed. */
  const icon = (
    key: string,
    label: string,
    text: string,
    onClick: () => void,
    disabled?: boolean,
    pressed?: boolean,
  ) => (
    <Hint label={label}>
      <button
        type="button"
        className="seqv-icon"
        aria-label={label}
        aria-pressed={pressed}
        aria-disabled={disabled}
        data-seqv-key={key}
        onClick={() => disabled || onClick()}
      >
        {text}
      </button>
    </Hint>
  );
  return (
    <fieldset ref={land} aria-label="Snapshots" className="seqv-player">
      {icon('prev', 'Previous snapshot', '‹', () => void ctx.snapshots.next(-1), at <= 0)}
      {icon(
        'play',
        playing ? 'Pause' : 'Play',
        playing ? '❚❚' : '▶',
        () => (playing ? ctx.snapshots.stop() : ctx.snapshots.play()),
        undefined,
        playing,
      )}
      {icon('next', 'Next snapshot', '›', () => void ctx.snapshots.next(1), at >= list.length - 1)}
      <span className="seqv-player-title">{at >= 0 ? (list[at].title ?? list[at].ref) : 'Start'}</span>
      <ol className="seqv-snapshots">
        {list.map((s, i) => (
          <li key={s.ref}>
            <button
              type="button"
              className="seqv-dot"
              aria-current={s.ref === current ? 'step' : undefined}
              aria-label={`${i + 1}. ${s.title ?? s.ref}${s.valid ? '' : ' (no longer applies)'}`}
              aria-disabled={!s.valid}
              data-invalid={s.valid ? undefined : ''}
              data-seqv-key={`dot:${s.ref}`}
              data-seqv-group="dots"
              onClick={() => s.valid && void ctx.snapshots.apply(s.ref)}
            >
              {i + 1}
            </button>
          </li>
        ))}
      </ol>
      <span className="seqv-quiet">
        {at + 1} / {list.length}
      </span>
    </fieldset>
  );
}

/**
 * A `plaintext` description as paragraphs, text only. The viewer has no markdown: a `markdown` description is shown by
 * a host's part through the `Description` slot, with its own renderer and `@molstar/seq-core/links`.
 */
export function Description({ context: ctx }: PartProps) {
  const paragraphs = useDerived(
    [ctx.document.shown, ctx.document.base, ctx.snapshots.current],
    () => plainDescription(ctx),
    { deps: [ctx] },
  );
  if (!paragraphs?.length) return null;
  return (
    <div className="seqv-description">
      {paragraphs.map((p, i) => (
        <p key={i}>{p}</p>
      ))}
    </div>
  );
}

export interface SearchProps extends PartProps {
  /** An application's search; else the UI spec's. With neither, the part draws nothing. */
  readonly search?: SearchFn;
}
/**
 * A search box and the current results: `2 / 14`, previous, next, save. In the bar it needs the UI spec's `search`;
 * anywhere else, `<Search context={ctx} search={fn} />` (inside a provider or not; give `fn` a stable identity, a new
 * one remakes the box). The status sits in a slot of a fixed width inside the field, so no state moves the field or the
 * buttons; Save is always there, hidden while it cannot run.
 */
export function Search({ context: ctx, search: given }: SearchProps) {
  const scope = useScopeMaybe();
  const search = given ?? scope?.ui.search;
  const [box, setBox] = useState<SearchBox | null>(null);
  useEffect(() => {
    if (!search) return;
    const b = new SearchBox(ctx, search);
    setBox(b);
    return () => {
      setBox(null);
      b.dispose();
    };
  }, [ctx, search]);
  return box && <SearchBar context={ctx} box={box} />;
}
function SearchBar({ context: ctx, box }: { context: SeqViewerContext; box: SearchBox }) {
  const { Hint } = usePrimitives();
  const v = useSeqState(box.view);
  useSeqState(ctx.document.shown);
  const land = useLand();
  const input = useRef<HTMLInputElement>(null);
  const save = useRef<HTMLButtonElement>(null);
  const clear = useRef<HTMLButtonElement>(null);
  const idle = !box.canSave();
  // Keep or Clear kept hiding while it has focus hands focus to the field (`handOff`), as the DOM part does.
  useLayoutEffect(() => {
    handOff(save.current, idle, input.current);
    handOff(clear.current, v.kept === 0, input.current);
  });
  const icon = (key: string, label: string, text: string, dir: 1 | -1) => (
    <Hint label={label}>
      <button
        type="button"
        className="seqv-icon"
        aria-label={label}
        aria-disabled={v.count === 0}
        data-seqv-key={key}
        onClick={() => v.count === 0 || box.step(dir)}
      >
        {text}
      </button>
    </Hint>
  );
  return (
    <div ref={land} className="seqv-search">
      <div className="seqv-search-field">
        <input
          ref={input}
          type="search"
          className="seqv-search-input"
          aria-label="Search"
          placeholder="Search"
          data-seqv-key="input"
          value={v.text}
          onChange={(e) => box.input(e.target.value)}
          onKeyDown={(e) => {
            if (e.key !== 'Enter') return;
            e.preventDefault();
            box.step(e.shiftKey ? -1 : 1);
          }}
        />
        <span className="seqv-search-status seqv-quiet" aria-live="polite">
          {box.status(v)}
        </span>
      </div>
      {icon('prev', 'Previous result', '‹', -1)}
      {icon('next', 'Next result', '›', 1)}
      <Hint label={KEEP_HINT}>
        <button
          type="button"
          className="seqv-switch-trigger seqv-keep"
          data-seqv-key="save"
          ref={save}
          data-idle={idle ? '' : undefined}
          onClick={() => box.save()}
        >
          <svg className="seqv-pin" aria-hidden="true" viewBox="0 0 16 16">
            <path d={PIN.outline} fill="none" stroke="currentColor" strokeWidth={1.3} strokeLinejoin="round" />
            <path d={PIN.needle} fill="none" stroke="currentColor" strokeWidth={1.3} strokeLinejoin="round" />
          </svg>
          Keep
        </button>
      </Hint>
      <Hint label={CLEAR_KEPT_HINT}>
        <button
          type="button"
          className="seqv-icon"
          aria-label="Clear kept"
          data-seqv-key="clear-kept"
          ref={clear}
          data-idle={v.kept > 0 ? undefined : ''}
          onClick={() => box.clearKept()}
        >
          <svg className="seqv-pin" aria-hidden="true" viewBox="0 0 16 16">
            <path d={CROSS} fill="none" stroke="currentColor" strokeWidth={1.4} strokeLinecap="round" />
          </svg>
        </button>
      </Hint>
    </div>
  );
}

const EXPORT = 'view.exportSvg';
/** "Export as SVG" when the context registers `view.exportSvg` (`browserExportAction`): a menu in the bar. */
export function Export({ context: ctx }: PartProps) {
  const { Menu } = usePrimitives();
  const surfaces = useScopeMaybe()?.surfaces;
  // Read again when the view changes, as the DOM part does.
  useSeqState(ctx.views.current);
  const [open, setOpen] = useState(false);
  const land = useLand();
  const entry = ctx.registry.get('actions', EXPORT) as { title?: string } | undefined;
  if (!entry) return null;
  const width = () => Math.round(surfaces?.values().next().value?.clientWidth || 1000);
  return (
    <span ref={land} className="seqv-export">
      <Menu
        trigger={
          <button type="button" className="seqv-switch-trigger" data-seqv-key="export">
            Export
          </button>
        }
        items={[
          {
            key: 'svg',
            label: entry.title ?? 'Export as SVG',
            run: () => void ctx.actions.run(EXPORT, { width: width() }),
          },
        ]}
        open={open}
        onOpenChange={setOpen}
      />
    </span>
  );
}

/**
 * A badge with the count; the list on demand, the most severe first and grouped by ref, notes last under their heading
 * (`issueBadge`, `issueList`), opened with the same `popover` as the menus. The `Issues` slot's default, named as the
 * DOM UI's.
 */
export function IssueBadge({ context: ctx }: PartProps) {
  const issues = useSeqState(ctx.issues.current);
  const badge = issueBadge(issues);
  return badge ? <IssueBox issues={issues} badge={badge} /> : null;
}
function IssueGroups({ groups }: { groups: readonly IssueGroup[] }) {
  return groups.map((g, k) => (
    <div key={k} className="seqv-issue-group">
      <div className="seqv-quiet">{g.ref ?? 'Viewer'}</div>
      {g.issues.map((i, n) => (
        <div key={n} className="seqv-issue" data-severity={i.severity}>
          <code>{i.code}</code> {i.message}
        </div>
      ))}
    </div>
  ));
}
function IssueBox({
  issues,
  badge,
}: {
  issues: readonly SeqIssue[];
  badge: NonNullable<ReturnType<typeof issueBadge>>;
}) {
  const { ScrollArea } = usePrimitives();
  const [open, setOpen] = useState(false);
  const root = useRef<HTMLSpanElement>(null);
  const id = useId();
  const [list] = useState(() => h('div', { id: `seqv-issues-${id}`, class: 'seqv-issue-list' }));
  useLanding(root);
  usePopover(root, list, { roving: false, open, onOpenChange: setOpen });
  const shown = open ? issueList(issues) : null;
  return (
    <span ref={root} className="seqv-issues">
      <button type="button" className="seqv-badge" data-seqv-key="issues" data-severity={badge.severity}>
        {badge.text}
      </button>
      {shown &&
        createPortal(
          <ScrollArea className="seqv-issue-scroll">
            <IssueGroups groups={shown.problems} />
            {shown.notes.length > 0 && (
              <div className="seqv-issue-notes">
                <div className="seqv-issue-heading">Notes</div>
                <IssueGroups groups={shown.notes} />
              </div>
            )}
            {shown.more > 0 && <div className="seqv-quiet">and {shown.more} more</div>}
          </ScrollArea>,
          list,
        )}
    </span>
  );
}

/** A quiet indicator while tasks run or regions in view are in flight, held briefly after the last region (`busyShown`). */
export function Busy({ context: ctx }: PartProps) {
  const [busy, setBusy] = useState(false);
  useEffect(() => {
    const b = busyShown(ctx, setBusy);
    return () => b.dispose();
  }, [ctx]);
  return busy ? <span role="status" aria-label="Loading" className="seqv-busy" /> : null;
}

export function EmptyState({ context: ctx }: PartProps) {
  const empty = useDerived([ctx.document.shown], () => !ctx.document.shown.value, { deps: [ctx] });
  return empty ? <div className="seqv-state">No sequences to show.</div> : null;
}

export function ErrorState({ context: ctx }: PartProps) {
  const errors = useDerived([ctx.graph.statuses, ctx.graph.assembled], () => rootErrors(ctx), { deps: [ctx] });
  if (errors.length === 0) return null;
  return (
    <div role="alert" className="seqv-state seqv-error">
      {errors.map((e) => (
        <div key={e}>{e}</div>
      ))}
    </div>
  );
}

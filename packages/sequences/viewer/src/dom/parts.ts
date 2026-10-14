/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The parts of seqspec:interaction#the-default-ui that are DOM chrome, as `DomPart`s
// (seqspec:interaction#the-dom-element): the same structure, class names and text as the React parts, from the same
// helpers (src/chrome.ts). Each also stands alone: a host mounts it into an element of its own, with a mounted UI's
// `ui` (its slots and surfaces) or none (the default primitives). Each subscribes to the smallest state it shows and writes only its own nodes; leaf state (a
// menu or the issue list open) stays in the part, never in the context. A node a part rebuilds goes through
// `spot(host).set`, which keeps keyboard focus by `data-seqv-key` (seqspec:interaction#accessibility).

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
  type MenuItem,
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
  type TrackMenuModel,
  trackMenuModel,
} from "../chrome.ts";
import {
  type HoverLabel,
  type LegendItem,
  type LegendRamp,
  type LegendSize,
  type SearchFn,
  SeqViewerConfig,
  type SeqViewerContext,
  type Unsubscribable,
} from "../index.ts";
import { sameData } from "../same.ts";
import { CLEAR_KEPT_HINT, KEEP_HINT, SearchBox } from "../search-box.ts";
import {
  button,
  type DomHandle,
  type DomPart,
  type DomUI,
  derived,
  type HintProps,
  h,
  later,
  lose,
  type MenuProps,
  off,
  ownedMenu,
  type PartProps,
  primitives,
  type SectionProps,
  setAttr,
  setText,
  spot,
  type TrackProps,
  watch,
} from "./core.ts";
import {
  busyShown,
  CROSS,
  handOff,
  legendGradient,
  legendGrid,
  legendSizeStyle,
  PALETTE,
  PIN,
  placeTooltip,
  popover,
  sizeCaption,
  tipSpot,
} from "./shared.ts";

type Handles = { dispose(): void }[];

/** A legend item: its swatch and label (the legend's, a track header's tooltip). */
function legendItem(i: LegendItem) {
  const swatch = h("span", { class: "seqv-swatch" });
  swatch.style.background = i.color ?? "var(--seq-ink)";
  return h("span", { class: "seqv-legend-item" }, swatch, i.label);
}
/** A continuous palette: its two ends around the gradient bar, and its caption. */
function legendRamp(r: LegendRamp) {
  const bar = h("span", { class: "seqv-ramp-bar", "aria-hidden": "true" });
  bar.style.background = legendGradient(r);
  const end = (n: number) => h("span", { class: "seqv-quiet" }, rampEnd(n));
  return h(
    "span",
    { class: "seqv-ramp" },
    end(r.lo),
    bar,
    end(r.hi),
    r.caption && h("span", { class: "seqv-quiet" }, r.caption),
  );
}

/**
 * A continuous palette with its range editable: the ends as fields around the gradient bar, the data's own span marked
 * under it (a range wider than the data shows where the data sits), and a reset once the range is the person's.
 */
function rangeRamp(ctx: SeqViewerContext, g: LegendGroup, r: LegendRamp) {
  const range = g.range as LegendRange;
  const field = (end: "lo" | "hi") => {
    const label = `${end === "lo" ? "Lowest" : "Highest"} value of the color scale`;
    const input = h("input", {
      class: "seqv-range-input",
      type: "text",
      inputmode: "decimal",
      "aria-label": label,
      "data-seqv-key": `range-${end}:${g.key}`,
    });
    input.value = rampEnd(range[end]);
    input.addEventListener("change", () => {
      if (!editRange(ctx, range, end, input.value)) input.value = rampEnd(range[end]);
    });
    return input;
  };
  const bar = h("span", { class: "seqv-ramp-bar", "aria-hidden": "true" });
  bar.style.background = legendGradient(r);
  const span = dataSpan(range);
  const mark = span && h("span", { class: "seqv-ramp-data", "aria-hidden": "true" });
  if (mark && span) {
    mark.style.left = `${span[0] * 100}%`;
    mark.style.width = `${(span[1] - span[0]) * 100}%`;
  }
  const reset = range.edited
    ? button(
        { class: "seqv-legend-more", "data-seqv-key": `range-reset:${g.key}`, "aria-label": "Reset the color scale" },
        () => resetRange(ctx, range),
        "reset",
      )
    : null;
  return h(
    "span",
    { class: "seqv-ramp seqv-ramp-edit" },
    field("lo"),
    h("span", { class: "seqv-ramp-track" }, bar, mark),
    field("hi"),
    reset,
    r.caption && h("span", { class: "seqv-quiet" }, r.caption),
  );
}

const SVG = "http://www.w3.org/2000/svg";
const svgEl = <K extends keyof SVGElementTagNameMap>(tag: K, attrs: Record<string, string | number>) => {
  const e = document.createElementNS(SVG, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, String(v));
  return e;
};
/** The `Search` part's "Keep" pin (`PIN`), decorative: the button says "Keep". */
function pinIcon(): SVGSVGElement {
  const svg = svgEl("svg", { class: "seqv-pin", "aria-hidden": "true", viewBox: "0 0 16 16" });
  const stroke = { fill: "none", stroke: "currentColor", "stroke-width": 1.3, "stroke-linejoin": "round" };
  svg.append(svgEl("path", { d: PIN.outline, ...stroke }), svgEl("path", { d: PIN.needle, ...stroke }));
  return svg;
}
/** The legend title's palette icon (`PALETTE`), an image labeled with the scheme shown. */
function paletteIcon(label: string): SVGSVGElement {
  const el = svgEl;
  const svg = el("svg", { class: "seqv-legend-icon", role: "img", "aria-label": label, viewBox: "0 0 16 16" });
  svg.append(
    el("path", { d: PALETTE.outline, fill: "none", stroke: "currentColor", "stroke-width": 1.3 }),
    ...PALETTE.dots.map(([cx, cy]) => el("circle", { cx, cy, r: 1.1, fill: "currentColor" })),
  );
  return svg;
}

export const ViewTabs: DomPart<PartProps> = (host, ctx) => {
  const at = spot(host);
  let tabs: HTMLButtonElement[] = [];
  const select = () => {
    const current = ctx.views.current.value;
    ctx.views.list.value.forEach((v, i) => {
      setAttr(tabs[i], "aria-selected", v.ref === current);
      setAttr(tabs[i], "tabindex", v.ref === current ? 0 : -1);
    });
  };
  const subs = [
    ctx.views.list.subscribe((views) => {
      const go = (k: number) => {
        const v = (k + views.length) % views.length;
        ctx.actions.run("view.switch", views[v].ref);
        tabs[v].focus();
      };
      tabs = views.map((v, i) => {
        const b = button(
          { role: "tab", class: "seqv-tab", "data-seqv-key": `tab:${v.ref}`, "data-seqv-group": "tabs" },
          () => ctx.actions.run("view.switch", v.ref),
          h("span", { class: "seqv-tab-label" }, v.label),
        );
        b.addEventListener("keydown", (e) => {
          if (e.key === "ArrowRight" || e.key === "ArrowLeft") go(i + (e.key === "ArrowRight" ? 1 : -1));
        });
        return b;
      });
      at.set(
        views.length < 2 ? null : h("div", { role: "tablist", "aria-label": "Views", class: "seqv-tabs" }, ...tabs),
      );
      select();
    }),
    // `select` ran in the list's replay.
    later(ctx.views.current, select),
  ];
  return { dispose: () => off([...subs, at]) };
};

/** Two or three choices as pressed buttons; more as a menu (open state kept by the switch, in its closure). */
function addSwitch(
  root: HTMLElement,
  ctx: SeqViewerContext,
  ui: DomUI | undefined,
  sw: SectionSwitch,
  run: (v: string | null) => void,
) {
  const shown = sw.options.find((o) => o.value === sw.value)?.label ?? sw.value ?? "default";
  const name = h("span", { class: "seqv-quiet" }, sw.name);
  if (sw.options.length > 3) {
    const trigger = button(
      { class: "seqv-switch-trigger", "aria-label": `${sw.name}: ${shown}`, "data-seqv-key": `switch:${sw.action}` },
      () => {},
      name,
      ` ${shown}`,
    );
    const items = sw.options.map((o) => ({
      key: String(o.value),
      label: o.label,
      checked: o.value === sw.value,
      run: () => run(o.value),
    }));
    return ownedMenu(primitives(ui).Menu, root, ctx, trigger, items);
  }
  const choice = (o: SectionSwitch["options"][number]) =>
    button(
      {
        "aria-pressed": o.value === sw.value,
        class: "seqv-choice",
        "data-seqv-key": `choice:${sw.action}:${o.value}`,
        "data-seqv-group": `choices:${sw.action}`,
      },
      () => run(o.value),
      o.label,
    );
  root.append(h("fieldset", { "aria-label": sw.name, class: "seqv-switch" }, name, ...sw.options.map(choice)));
  return null;
}

/**
 * A section's label, its filter chip, then its `HiddenChips` and `SectionSwitches` (the UI's slots, else the defaults);
 * nothing when it would be empty. The hidden tracks' chips sit in this row, left of the switches: a header with only
 * those is still a header.
 */
export const SectionHeader: DomPart<SectionProps> = (host, ctx, { ui, section }) => {
  const at = spot(host);
  const [Chips, Switches] = ui ? [ui.slots.HiddenChips, ui.slots.SectionSwitches] : [HiddenChips, SectionSwitches];
  let inner: Handles = [];
  const subs = derived(
    sectionChromeStates(ctx),
    () => {
      const info = sectionChrome(ctx, section);
      const hidden = !!Chips && info.hidden.length > 0;
      return { label: info.label, filter: info.filter, switches: !!Switches && info.switches.length > 0, hidden };
    },
    (info) => {
      const old = inner;
      inner = [];
      if (!info.label && !info.filter && !info.switches && !info.hidden) {
        at.set(null);
        return off(old);
      }
      const clear = () => ctx.actions.run("section.clearFilter", section);
      const chip =
        info.filter &&
        button(
          {
            class: "seqv-chip",
            "aria-label": `Show all positions (now only ${info.filter})`,
            "data-seqv-key": "filter",
          },
          clear,
          `Only ${info.filter} `,
          h("span", { "aria-hidden": "true" }, "✕"),
        );
      const label = info.label && h("span", { class: "seqv-section-label" }, info.label);
      const root = h("div", { class: "seqv-section-header" }, label, chip, h("span", { class: "seqv-spacer" }));
      if (Chips) inner.push(Chips(root, ctx, { ui, section }));
      if (Switches) inner.push(Switches(root, ctx, { ui, section }));
      // The old header goes after the swap: a menu of it that was open must not take focus back.
      at.set(root);
      off(old);
    },
  );
  return { dispose: () => off([...inner, ...subs, at]) };
};

/** A section's switches, each running its `section.*` action: in its header, or wherever a host mounts them. */
export const SectionSwitches: DomPart<SectionProps> = (host, ctx, { ui, section }) => {
  const at = spot(host);
  let menus: Handles = [];
  const subs = derived(
    sectionChromeStates(ctx),
    () => sectionChrome(ctx, section).switches,
    (switches) => {
      const old = menus;
      menus = [];
      const root = h("div", { class: "seqv-switches" });
      for (const sw of switches) {
        const m = addSwitch(root, ctx, ui, sw, (value) => ctx.actions.run(sw.action, { section, value }));
        if (m) menus.push(m);
      }
      at.set(switches.length === 0 ? null : root);
      off(old);
    },
  );
  return { dispose: () => off([...menus, ...subs, at]) };
};

/**
 * The menu at a track header: filter view, focus, hide, remove. Alone in the box `ViewSurface` places at the header,
 * so it appends there. It holds its menu's `open` (leaf state) and opens it on `events.menu` for its track (Shift+F10
 * or ContextMenu on the keyboard cursor, seqspec:viewer#behaviors): a host's `Menu` gets it through `update`.
 */
export const TrackActions: DomPart<TrackProps> = (host, ctx, { ui, section, track, row }) => {
  let [open, items] = [false, [] as readonly MenuItem[]];
  let hint: DomHandle<HintProps> | null = null;
  let menu: DomHandle<MenuProps> | null = null;
  let trigger: HTMLButtonElement | null = null;
  const props = (): MenuProps => ({ trigger: trigger as HTMLButtonElement, items, open, onOpenChange });
  const onOpenChange = (o: boolean) => {
    open = o;
    // Opened: built again for the line it opens on (a split row's own items).
    if (o) sync();
    // Closed, the menu opens at its top again, not in its "Colors" submenu.
    if (!o && drill.layer !== null) {
      drill.layer = null;
      sync();
    }
    menu?.update?.(props());
  };
  // The "Colors" submenu: an item chosen closed the menu; it opens again with the other list (focus on its first item).
  const drill: { layer: string | null } & ColorsDrill = {
    layer: null,
    open(layer) {
      drill.layer = layer;
      sync();
      onOpenChange(true);
    },
  };
  const { Menu, Hint } = primitives(ui);
  const children = (x: HTMLElement) => {
    menu = Menu(x, ctx, props());
  };
  let last: TrackMenuModel | undefined;
  const apply = (m: TrackMenuModel) => {
    items = m.items;
    if (items.length === 0) {
      // Nothing applies any more (the track was hidden by its own menu): focus on the trigger lands on the surface.
      lose(host, host, () => off([menu, hint]));
      [menu, hint, trigger, open] = [null, null, null, false];
      return;
    }
    const label = `Actions for ${m.label}`;
    if (trigger) setAttr(trigger, "aria-label", label);
    else trigger = button({ class: "seqv-icon", "aria-label": label, "data-seqv-key": "track-actions" }, () => {}, "⋯");
    if (hint) hint.update?.({ label, children });
    else hint = Hint(host, ctx, { label, children });
    menu?.update?.(props());
  };
  // Written when what it shows changed (`sameMenu`).
  const sync = () => {
    const m = trackMenuModel(ctx, section, track, drill, row?.() ?? -1);
    if (last && sameMenu(last, m)) return;
    last = m;
    apply(m);
  };
  const subs = [
    // A new view: which layers have schemes, and the one picked (here, in the legend or by command).
    ...watch([ctx.document.shown, ctx.ui.expanded, ctx.graph.assembled], sync),
    // "Retry" comes with a failed region: the menu is read again on the streams' version, but only while it is open.
    later(ctx.streams.version, () => open && sync()),
    ctx.events.menu.subscribe((e) => {
      if (menu && e.section === section && e.track === track) onOpenChange(true);
    }),
  ];
  return { dispose: () => off([menu, hint, ...subs]) };
};

export const HiddenChips: DomPart<SectionProps> = (host, ctx, { section }) => {
  const at = spot(host);
  const subs = derived(
    sectionChromeStates(ctx),
    () => sectionChrome(ctx, section).hidden,
    (hidden) => {
      const chip = (t: (typeof hidden)[number]) =>
        button(
          { class: "seqv-chip", "data-seqv-key": `chip:${t.ref}`, "data-seqv-group": "chips" },
          () => ctx.actions.run("track.show", t.ref),
          `Show ${t.label}`,
        );
      const root = h("fieldset", { "aria-label": "Hidden tracks", class: "seqv-chips" }, ...hidden.map(chip));
      at.set(hidden.length === 0 ? null : root);
    },
  );
  return { dispose: () => off([...subs, at]) };
};

/**
 * The hover's label, placed at its owner's anchor (`placeTooltip`): the pointer's, or the keyboard cursor's cell
 *, which it follows while the viewport moves under it. Its content is rewritten only when the label changes.
 */
export const Tooltip: DomPart<PartProps> = (host, ctx, { ui }) => {
  // A canvas overlay, placed at its UI's surfaces: nothing without a UI (as React's without a provider).
  if (!ui) return { dispose() {} };
  const surfaces = ui.surfaces;
  const at = spot(host);
  const last = tipSpot();
  let tip: HTMLElement | null = null;
  let shown: HoverLabel | null = null;
  let follow: Unsubscribable | undefined;
  const place = () => tip && placeTooltip(ctx, surfaces, tip, last);
  const subs = [
    ctx.labels.current.subscribe((label) => {
      if (!label) {
        [tip, shown] = [null, null];
        return at.set(null);
      }
      if (!tip) {
        tip = h("div", { role: "tooltip", class: "seqv-tooltip" });
        [last.x, last.y, last.shown] = [Number.NaN, Number.NaN, true];
      }
      // The same label again (the hover refreshed under a new document) writes nothing; a legend's swatches or ramp
      // that change while it shows are written, as React's Tooltip compares (`sameData`).
      if (!shown || !sameData(shown, label)) {
        const line = (l: (typeof label.lines)[number]) =>
          h("div", { class: "seqv-tooltip-line" }, h("span", { class: "seqv-quiet" }, l.key), ` ${l.value}`);
        const title = label.title && h("div", { class: "seqv-tooltip-title" }, label.title);
        const lg = label.legend;
        const legend =
          lg &&
          h(
            "div",
            { class: "seqv-tooltip-legend" },
            ...lg.items.map(legendItem),
            ...lg.ramps.map(legendRamp),
            lg.more > 0 && h("span", { class: "seqv-quiet" }, `+${lg.more} more`),
          );
        tip.replaceChildren(...[title, ...label.lines.map(line), legend].filter((n) => !!n));
        shown = label;
      }
      at.set(tip);
      place();
    }),
    later(ctx.pointer.version, place),
    // The replay is the first run: it follows the cursor's section if there is one, and places.
    ctx.ui.keyCursor.subscribe((c) => {
      follow?.unsubscribe();
      const live = c && ctx.viewports.get(c.section)?.live.version;
      follow = live ? later(live, place) : undefined;
      place();
    }),
  ];
  return { dispose: () => off([...subs, { dispose: () => follow?.unsubscribe() }, at]) };
};

/**
 * On the canvas's grid: each group's name in the label column, its items from the data area's left (`legendGrid`);
 * a group past `UI.LegendItems` ends in "+n more" (`ui.legendOpen`), as does the legend past `UI.LegendGroups` groups;
 * a continuous palette is a ramp with its ends and caption, a sample size a dot. A group with a color scheme picker has
 * its title as the picker's trigger (seqspec:interaction#color-schemes-at-runtime). The groups are rebuilt on a change;
 * focus stays on the same group's control.
 */
export const Legend: DomPart<PartProps> = (host, ctx, { ui }) => {
  const at = spot(host);
  const { Menu, Hint } = primitives(ui);
  const grids: (() => void)[] = [];
  let parts: Handles = [];
  /**
   * The title: its name, or with a picker the trigger of the group's scheme menu (the name, dotted underline, then a
   * palette icon on hover or focus); in a hint with the full name only while the name is cut (`legendGrid`'s `onCut`).
   * Its handles go into `own`.
   */
  const title = (g: LegendGroup, cut: boolean, own: Handles) => {
    const cell = h("span", { class: "seqv-legend-track" });
    const t = legendTitle(g);
    const name = h("span", { class: "seqv-legend-name" }, t.name);
    const body = (x: HTMLElement) => {
      if (!g.colors || !t.colors) return x.append(name);
      const trigger = button(
        { class: "seqv-legend-title", "data-seqv-key": `colors:${g.key}` },
        () => {},
        name,
        paletteIcon(t.colors),
      );
      own.push(ownedMenu(Menu, x, ctx, trigger, schemeItems(ctx, g.colors.layers, g.colors)));
    };
    if (cut) own.push(Hint(cell, ctx, { label: t.name, children: body }));
    else body(cell);
    return cell;
  };
  const item = legendItem;
  const ramp = legendRamp;
  const size = (s: LegendSize, i: number, all: readonly LegendSize[]) => {
    const swatch = h("span", { class: "seqv-size", "aria-hidden": "true" });
    Object.assign(swatch.style, legendSizeStyle(s));
    const caption = sizeCaption(all, i);
    const item = h("span", { class: "seqv-legend-item" }, swatch, s.label);
    return caption === null ? item : [h("span", { class: "seqv-quiet" }, caption), item];
  };
  const more = (key: string, n: number, expanded: boolean) =>
    button(
      { class: "seqv-legend-more", "aria-expanded": expanded, "data-seqv-key": `more:${key}` },
      () => ctx.ui.toggleLegend(key),
      n > 0 ? `+${n} more` : "fewer",
    );
  const row = (g: LegendGroup) => {
    const { items, more: n } = legendShown(ctx, g);
    const expanded = ctx.ui.legendOpen.value.has(g.key);
    const toggle =
      (n > 0 || (expanded && g.items.length > ctx.config.get(SeqViewerConfig.UI.LegendItems))) &&
      more(g.key, n, expanded);
    let own: Handles = [];
    let cell = title(g, false, own);
    parts.push({ dispose: () => off(own) });
    const el = h(
      "div",
      { class: "seqv-legend-group", "data-section": g.section },
      cell,
      h(
        "span",
        { class: "seqv-legend-items" },
        ...items.map(item),
        ...g.ramps.map((r) => (g.range ? rangeRamp(ctx, g, r) : ramp(r))),
        ...g.sizes.flatMap(size),
        toggle,
      ),
    );
    // Cut, or whole again: the title is made again with its hint, or without.
    const recut = (cut: boolean) => {
      const old = own;
      own = [];
      const next = title(g, cut, own);
      el.replaceChild(next, cell);
      cell = next;
      off(old);
    };
    grids.push(legendGrid(ctx, g.section, el, recut));
    return el;
  };
  /**
   * "Highlights": each entry of the sections' `highlights` (a kept result set, bookmarks), its style's swatch, its label,
   * and a ✕ that removes it (`highlight.remove`, one patch: Undo brings it back).
   */
  const highlights = (entries: readonly SectionHighlightEntry[]) => {
    const item = (e: SectionHighlightEntry) => {
      const swatch = h("span", { class: "seqv-swatch seqv-swatch-set", "data-style": e.style, "aria-hidden": "true" });
      const x = button(
        { class: "seqv-legend-remove", "aria-label": `Remove ${e.label}`, "data-seqv-key": `hl:${e.key}` },
        () => ctx.actions.run("highlight.remove", { section: e.section, data: e.data }, "pointer"),
        "✕",
      );
      const el = h("span", { class: "seqv-legend-item" }, swatch, h("span", null, e.label));
      const hint = Hint(el, ctx, {
        label: `Remove ${e.label} from the highlights`,
        children: (c) => c.append(x),
      });
      parts.push(hint);
      return el;
    };
    const el = h(
      "div",
      { class: "seqv-legend-group", "data-section": entries[0].section },
      h("span", { class: "seqv-legend-track" }, h("span", { class: "seqv-legend-name" }, "Highlights")),
      h("span", { class: "seqv-legend-items" }, ...entries.map(item)),
    );
    grids.push(legendGrid(ctx, entries[0].section, el));
    return el;
  };
  /** The legend's own "+n more" (groups past `UI.LegendGroups`), on the first group's grid. */
  const rest = (all: readonly LegendGroup[]) => {
    const s = legendShownGroups(ctx, all);
    if (!s.toggle) return null;
    const el = h(
      "div",
      { class: "seqv-legend-group", "data-section": all[0].section },
      h("span", { class: "seqv-legend-track" }),
      h("span", { class: "seqv-legend-items" }, more(s.key, s.more, s.more === 0)),
    );
    grids.push(legendGrid(ctx, all[0].section, el));
    return el;
  };
  const subs = derived(
    [ctx.graph.assembled, ctx.ui.legendOpen, ctx.ui.legendExpanded],
    // The open groups are in the value, so a toggle writes (the groups alone are one memoized list per view; the
    // groups, the kept entries and the open set are each the same object until they change).
    () => ({
      groups: legendGroups(ctx),
      kept: sectionHighlightEntries(ctx),
      open: ctx.ui.legendOpen.value,
      expanded: legendExpanded(ctx, ui?.legend),
    }),
    ({ groups, kept, expanded }) => {
      for (const d of grids.splice(0)) d();
      const old = parts;
      parts = [];
      const chevron = h("span", { class: "seqv-chevron", "aria-hidden": "true" });
      const root = h(
        "details",
        { class: "seqv-legend", open: expanded },
        h("summary", { "data-seqv-key": "legend" }, chevron, "Legend"),
        ...legendShownGroups(ctx, groups).groups.map(row),
        groups.length > 0 && rest(groups),
        kept.length > 0 && highlights(kept),
      );
      // The person's toggle is UI state (`ui.legendExpanded`), as `view.legend` is.
      root.addEventListener("toggle", () => {
        if (root.open !== legendExpanded(ctx, ui?.legend)) ctx.ui.setLegendExpanded(root.open);
      });
      const first = groups[0]?.section ?? kept[0]?.section;
      if (first !== undefined) grids.push(legendGrid(ctx, first, root));
      at.set(groups.length === 0 && kept.length === 0 ? null : root);
      off(old);
    },
  );
  return {
    dispose: () => {
      for (const d of grids.splice(0)) d();
      off([...parts, ...subs, at]);
    },
  };
};

export const SnapshotPlayer: DomPart<PartProps> = (host, ctx, { ui }) => {
  const at = spot(host);
  const { Hint } = primitives(ui);
  let hints: DomHandle<HintProps>[] = [];
  let shown: {
    prev: HTMLButtonElement;
    play: HTMLButtonElement;
    next: HTMLButtonElement;
    title: HTMLElement;
    dots: HTMLButtonElement[];
    pos: HTMLElement;
  } | null = null;
  const S = ctx.snapshots;
  /** An unusable control stays focusable: `aria-disabled`, never `disabled`, and does nothing when pressed. */
  const usable = (act: () => void) => (e: MouseEvent) =>
    (e.currentTarget as HTMLElement).getAttribute("aria-disabled") === "true" || act();
  const icon = (root: HTMLElement, key: string, label: string, text: string, onClick: () => void) => {
    const b = button({ class: "seqv-icon", "aria-label": label, "data-seqv-key": key }, usable(onClick), text);
    hints.push(Hint(root, ctx, { label, children: (x) => x.append(b) }));
    return b;
  };
  const current = () => {
    if (!shown) return;
    const list = S.list.value;
    const i = list.findIndex((s) => s.ref === S.current.value);
    setAttr(shown.prev, "aria-disabled", i <= 0);
    setAttr(shown.next, "aria-disabled", i >= list.length - 1);
    setText(shown.title, i >= 0 ? (list[i].title ?? list[i].ref) : "Start");
    for (let k = 0; k < shown.dots.length; k++) setAttr(shown.dots[k], "aria-current", k === i ? "step" : null);
    setText(shown.pos, `${i + 1} / ${list.length}`);
  };
  const playing = () => {
    if (!shown) return;
    const p = S.playing.value;
    setAttr(shown.play, "aria-label", p ? "Pause" : "Play");
    setAttr(shown.play, "aria-pressed", p);
    setText(shown.play, p ? "❚❚" : "▶");
    hints[1].update?.({ label: p ? "Pause" : "Play", children: () => {} });
  };
  const build = (list: readonly { ref: string; title: string | null; valid: boolean }[]) => {
    const old = hints;
    hints = [];
    shown = null;
    if (list.length === 0) {
      at.set(null);
      return off(old);
    }
    const root = h("fieldset", { "aria-label": "Snapshots", class: "seqv-player" });
    const prev = icon(root, "prev", "Previous snapshot", "‹", () => void S.next(-1));
    const play = icon(root, "play", "Play", "▶", () => (S.playing.value ? S.stop() : S.play()));
    const next = icon(root, "next", "Next snapshot", "›", () => void S.next(1));
    const dots = list.map((s, i) =>
      button(
        {
          class: "seqv-dot",
          "aria-label": `${i + 1}. ${s.title ?? s.ref}${s.valid ? "" : " (no longer applies)"}`,
          "aria-disabled": !s.valid,
          "data-invalid": s.valid ? null : "",
          "data-seqv-key": `dot:${s.ref}`,
          "data-seqv-group": "dots",
        },
        usable(() => void S.apply(s.ref)),
        String(i + 1),
      ),
    );
    shown = {
      prev,
      play,
      next,
      dots,
      title: h("span", { class: "seqv-player-title" }),
      pos: h("span", { class: "seqv-quiet" }),
    };
    root.append(shown.title, h("ol", { class: "seqv-snapshots" }, ...dots.map((d) => h("li", null, d))), shown.pos);
    current();
    playing();
    at.set(root);
    off(old);
  };
  const subs = [
    ...derived(
      [S.list],
      () => S.list.value.map((s) => ({ ref: s.ref, title: s.title, valid: s.valid })),
      build,
      sameData,
    ),
    // `build` ran both in the list's replay.
    later(S.current, current),
    later(S.playing, playing),
  ];
  return { dispose: () => off([...hints, ...subs, at]) };
};

/**
 * A `plaintext` description as paragraphs, text only. The viewer has no markdown: a `markdown` description is shown by
 * a host's part through the `Description` slot, with its own renderer and `@molstar/seq-core/links`.
 */
export const Description: DomPart<PartProps> = (host, ctx) => {
  const at = spot(host);
  const subs = derived(
    [ctx.document.shown, ctx.document.base, ctx.snapshots.current],
    () => plainDescription(ctx),
    (paragraphs) =>
      at.set(
        !paragraphs?.length ? null : h("div", { class: "seqv-description" }, ...paragraphs.map((p) => h("p", null, p))),
      ),
  );
  return { dispose: () => off([...subs, at]) };
};

const EXPORT = "view.exportSvg";
/** What the `Search` part takes: in the bar, the UI (`ui.search` and its `Hint`); anywhere else, just `search`. */
export interface SearchProps {
  readonly ui?: DomUI;
  /** An application's search; else the UI spec's. With neither, the part draws nothing. */
  readonly search?: SearchFn;
}
/**
 * A search box and the current results: `2 / 14`, previous, next, save. In the bar it needs the UI spec's `search`;
 * anywhere else, `Search(host, ctx, { search })` (with the skin's styles). The status sits in a slot of a fixed width
 * inside the field, so no state moves the field or the buttons; Save is always there, hidden while it cannot run.
 */
export const Search: DomPart<SearchProps> = (host, ctx, props) => {
  const search = props.search ?? props.ui?.search;
  if (!search) return { dispose() {} };
  const box = new SearchBox(ctx, search);
  const { Hint } = primitives(props.ui);
  const attrs = { type: "search", class: "seqv-search-input", "aria-label": "Search", placeholder: "Search" };
  const input = h("input", { ...attrs, "data-seqv-key": "input" });
  input.addEventListener("input", () => box.input(input.value));
  input.addEventListener("keydown", (e) => {
    if (e.key !== "Enter") return;
    e.preventDefault();
    box.step(e.shiftKey ? -1 : 1);
  });
  const status = h("span", { class: "seqv-search-status seqv-quiet", "aria-live": "polite" });
  const root = h("div", { class: "seqv-search" }, h("div", { class: "seqv-search-field" }, input, status));
  const hints: Handles = [];
  const icon = (key: string, label: string, text: string, dir: 1 | -1) => {
    const b = button(
      { class: "seqv-icon", "aria-label": label, "data-seqv-key": key },
      () => b.getAttribute("aria-disabled") === "true" || box.step(dir),
      text,
    );
    hints.push(Hint(root, ctx, { label, children: (x) => x.append(b) }));
    return b;
  };
  const steps = [icon("prev", "Previous result", "‹", -1), icon("next", "Next result", "›", 1)];
  const save = button(
    { class: "seqv-switch-trigger seqv-keep", "data-seqv-key": "save" },
    () => box.save(),
    pinIcon(),
    "Keep",
  );
  hints.push(Hint(root, ctx, { label: KEEP_HINT, children: (x) => x.append(save) }));
  const cross = svgEl("svg", { class: "seqv-pin", "aria-hidden": "true", viewBox: "0 0 16 16" });
  cross.append(
    svgEl("path", { d: CROSS, fill: "none", stroke: "currentColor", "stroke-width": 1.4, "stroke-linecap": "round" }),
  );
  const clear = button(
    { class: "seqv-icon", "aria-label": "Clear kept", "data-seqv-key": "clear-kept" },
    () => box.clearKept(),
    cross,
  );
  hints.push(Hint(root, ctx, { label: CLEAR_KEPT_HINT, children: (x) => x.append(clear) }));
  host.append(root);
  const show = () => {
    const v = box.view.value;
    if (input.value !== v.text) input.value = v.text;
    setText(status, box.status(v));
    for (const b of steps) setAttr(b, "aria-disabled", v.count === 0);
    handOff(save, !box.canSave(), input);
    handOff(clear, v.kept === 0, input);
    setAttr(save, "data-idle", !box.canSave());
    setAttr(clear, "data-idle", v.kept === 0);
  };
  const subs = watch([box.view, ctx.document.shown], show);
  return {
    dispose() {
      off([...subs, ...hints]);
      box.dispose();
      root.remove();
    },
  };
};

/** "Export as SVG" when the context registers `view.exportSvg` (`browserExportAction`): a menu in the bar. */
export const Export: DomPart<PartProps> = (host, ctx, { ui }) => {
  const at = spot(host);
  const root = h("span", { class: "seqv-export" });
  const trigger = button({ class: "seqv-switch-trigger", "data-seqv-key": "export" }, () => {}, "Export");
  const width = () => Math.round(ui?.surfaces.values().next().value?.clientWidth || 1000);
  const title = (ctx.registry.get("actions", EXPORT) as { title?: string } | undefined)?.title ?? "Export as SVG";
  const menu = ownedMenu(primitives(ui).Menu, root, ctx, trigger, [
    { key: "svg", label: title, run: () => void ctx.actions.run(EXPORT, { width: width() }) },
  ]);
  const sub = ctx.views.current.subscribe(() => at.set(ctx.registry.get("actions", EXPORT) ? root : null));
  return { dispose: () => off([sub, menu, at]) };
};

let lists = 0;
/**
 * A badge with the count; the list on demand, the most severe first and grouped by ref, notes last under their heading
 * (`issueBadge`, `issueList`). The `Issues` slot's default.
 */
export const IssueBadge: DomPart<PartProps> = (host, ctx, { ui }) => {
  const at = spot(host);
  let scroll: Handles[number] | null = null;
  const badge = button({ class: "seqv-badge", "data-seqv-key": "issues" }, () => {});
  const list = h("div", { id: `seqv-issues-${++lists}`, class: "seqv-issue-list" });
  const root = h("span", { class: "seqv-issues" }, badge);
  const group = (g: IssueGroup) =>
    h(
      "div",
      { class: "seqv-issue-group" },
      h("div", { class: "seqv-quiet" }, g.ref ?? "Viewer"),
      ...g.issues.map((i) =>
        h("div", { class: "seqv-issue", "data-severity": i.severity }, h("code", null, i.code), ` ${i.message}`),
      ),
    );
  const fill = (v: HTMLElement) => {
    const { problems, notes, more } = issueList(ctx.issues.current.value);
    v.append(...problems.map(group));
    if (notes.length > 0)
      v.append(
        h(
          "div",
          { class: "seqv-issue-notes" },
          h("div", { class: "seqv-issue-heading" }, "Notes"),
          ...notes.map(group),
        ),
      );
    if (more > 0) v.append(h("div", { class: "seqv-quiet" }, `and ${more} more`));
  };
  const content = () => {
    off([scroll]);
    list.replaceChildren();
    scroll = primitives(ui).ScrollArea(list, ctx, { className: "seqv-issue-scroll", children: fill });
  };
  const pop = popover(root, badge, list, {
    roving: false,
    ask(open) {
      if (!open) return pop.close();
      content();
      pop.open();
    },
  });
  const show = () => {
    const b = issueBadge(ctx.issues.current.value);
    if (!b) {
      pop.close();
      return at.set(null);
    }
    setAttr(badge, "data-severity", b.severity);
    setText(badge, b.text);
    if (pop.isOpen) content();
    at.set(root);
  };
  const subs = [ctx.issues.current.subscribe(show)];
  return {
    dispose() {
      pop.dispose();
      off([scroll, ...subs, at]);
    },
  };
};

/** A quiet indicator while tasks run or regions in view are in flight. */
export const Busy: DomPart<PartProps> = (host, ctx) => {
  const at = spot(host);
  const el = h("span", { role: "status", "aria-label": "Loading", class: "seqv-busy" });
  const shown = busyShown(ctx, (busy) => at.set(busy ? el : null));
  return { dispose: () => off([shown, at]) };
};

export const EmptyState: DomPart<PartProps> = (host, ctx) => {
  const at = spot(host);
  const el = h("div", { class: "seqv-state" }, "No sequences to show.");
  const subs = derived(
    [ctx.document.shown],
    () => !ctx.document.shown.value,
    (empty) => at.set(empty ? el : null),
  );
  return { dispose: () => off([...subs, at]) };
};

export const ErrorState: DomPart<PartProps> = (host, ctx) => {
  const at = spot(host);
  const subs = derived(
    [ctx.graph.statuses, ctx.graph.assembled],
    () => rootErrors(ctx),
    (errors) => {
      const root = h("div", { role: "alert", class: "seqv-state seqv-error" }, ...errors.map((e) => h("div", null, e)));
      at.set(errors.length === 0 ? null : root);
    },
  );
  return { dispose: () => off([...subs, at]) };
};

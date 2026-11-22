/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Value } from '@molstar/core/repr/tree';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import { type Box, HeaderStateNames, type HiddenTrack } from '@molstar/seqx-core/render';
import type { IconName } from '@molstar/ui/components/icon/icon';
import type { HeadersBehavior } from './behaviors/headers.js';
import type { CardProviderEntry } from './card-provider.js';
import { FilteredHeaderEntry } from './filter-chip.js';
import type { SeqViewer } from './viewer.js';

// A track's and a section's header entries (_spec/README.md "Behaviors", Headers): registry entries of kind
// `header-action` naming a command, made into slots, chips and menu items once a layout from DOM-free facts. The
// built-ins are entries like a host's (`BuiltInHeaderEntries`), in the default spec's registry. Data and types only:
// the `headers` behavior reads them.

/** A layer of a header's track, as the entries read it. */
export interface HeaderLayerFacts {
  /** Its address (a layer with no id: its key, `<track>:<n>`, which no command takes). */
  readonly address: string;
  readonly type: string;
  readonly label: string;
  readonly hidden: boolean;
  /** A rows layer ordered by a tree it has. */
  readonly tree?: boolean;
  /** The rows its user filter leaves out (0: none). */
  readonly filtered: number;
  /** Its color as written (absent: its type's default). */
  readonly color?: Value;
  /** The schemes that fit it (what "Color ▸" offers): none for a layer with no color channel. */
  readonly schemes: readonly { readonly name: string; readonly label: string }[];
}

/**
 * What a header's entries read (_spec/README.md "Behaviors", Headers: `HeaderFacts`): plain data from the resolved view
 * and the last layout, made once a layout a shown header (and when a menu opens). Pure functions of it decide whether
 * an entry shows, its label, its args and a state's text.
 */
export interface HeaderFacts {
  /** Whose header: a track's, a section's, or an element's menu. */
  readonly on: 'track' | 'section' | 'element';
  /** The track's address (a section's for its own header). */
  readonly address: string;
  readonly label: string;
  /** The section's address and the system of its axis; whether it is linear (no ring). */
  readonly section: string;
  readonly system: string;
  readonly linear: boolean;
  /** The section's tracks' ids in spec order, the hidden ones too (`track.move`'s siblings). */
  readonly tracks: readonly string[];
  /** The track's layers (none on a section's header). */
  readonly layers: readonly HeaderLayerFacts[];
  /** The section's tracks the spec hides, in spec order. */
  readonly hiddenTracks: readonly HiddenTrack[];
  /** The element a menu was opened on (`on: 'element'`): its layer and key. */
  readonly element?: { readonly layer: string; readonly key: string | number };
  /**
   * The block a section's menu was opened on (its label, its band: genome coordinates): its address, label, its
   * system's name and whether it reads against the section.
   */
  readonly block?: {
    readonly address: string;
    readonly label: string;
    readonly name: string;
    readonly reversed: boolean;
  };
  /** The track has a rows layer with a tree (`order: 'tree'`): whether the reader opened it, and whether the frame lacks the room (the tree is collapsed to its toggle, or kept open past it). */
  readonly tree?: { readonly open: boolean; readonly lacks: boolean };
  /** A section of assemblies (its `include` can be set), and whether it has one. */
  readonly assemblies?: boolean;
  readonly included?: boolean;
}

/** An item of an entry's submenu: a label and the command it runs with its args (a layer's checkbox, a scheme). */
export interface HeaderSubEntry {
  readonly label: string;
  readonly command: string;
  readonly args: Readonly<Record<string, unknown>>;
  /** A checked item (a shown layer), or the chosen one of a choice. */
  readonly checked?: boolean;
}

/**
 * A header entry (_spec/README.md "Behaviors", Headers): an action (a slot when it has an icon, and a menu item) or a
 * state (a chip), on a track's header, a section's or an element's menu; namespaced (`seq:hide`, `myapp:color-by`).
 */
export interface HeaderEntry {
  /** Namespaced: two entries under one name are a registration error. */
  readonly name: string;
  readonly on: 'track' | 'section' | 'element';
  /** `action` (default): a slot and a menu item; `state`: a chip ("12 filtered"). */
  readonly kind?: 'action' | 'state';
  /** `@molstar/ui`'s; an action without one is in the menu only (the command's icon when absent). */
  readonly icon?: IconName;
  /**
   * The slot's icon as SVG markup in a 24-unit box (Lucide's: `@molstar/ui`'s per-icon exports, `lucideEyeOff`), read
   * into paths when the entry is laid out; absent: a built-in's icon by `icon`, else the label's first letter.
   */
  readonly iconMarkup?: string;
  /** Its label (the command's title when absent). */
  readonly label?: string | ((f: HeaderFacts) => string);
  /** Shown at all (absent: always). The command's own `when` says whether it runs now. */
  readonly when?: (f: HeaderFacts) => boolean;
  /** A command the viewer has (a state may have none: it says, it does not run). */
  readonly command?: string;
  /** The command's context (default `{ address: f.address }`). */
  readonly args?: (f: HeaderFacts) => Readonly<Record<string, unknown>>;
  /** A state's chip text ("12 filtered"); none: no chip. */
  readonly text?: (f: HeaderFacts) => string | undefined;
  /** Slots by order, then name. */
  readonly order?: number;
  /** false: the menu only (default: true when it has an icon). */
  readonly slot?: boolean;
  /**
   * A submenu of other commands or args (a layer's checkboxes, a chip's list of what it hides): each runs its own
   * command. Its slot opens the menu at it; a chip's secondary click lists them before the entry itself ("Show all").
   * Without it, a command whose one param is a select makes the submenu: its options, each the entry's args with the
   * choice (the command owns its choices: the menu, a palette and an agent read the same).
   */
  readonly items?: (f: HeaderFacts) => readonly HeaderSubEntry[];
}

/** A registry entry of the viewer: a header entry (`SeqViewerSpec.registry`, beside the core's kinds). */
export interface HeaderActionEntry {
  readonly kind: 'header-action';
  readonly entry: HeaderEntry;
}

/** What a viewer's registry lists: the core's entries and the viewer's own kinds (header entries, card providers). */
export type SeqViewerRegistryEntry = SeqRegistryEntry | HeaderActionEntry | CardProviderEntry;

/** What a header (or a menu) is about: its section's address, the track's, the element a menu was opened on. */
export interface HeaderTarget {
  readonly section: string;
  readonly track?: string;
  readonly element?: { readonly layer: string; readonly key: string | number };
  /** The block under the point (its index in the section's blocks), on a section of blocks. */
  readonly block?: number;
}

/**
 * An item of a header's menu or parts (_spec/README.md "Behaviors", Headers: one item list for both): DOM-free data. A
 * choice runs its `command` with `args` through `viewer.execute` (`HeadersBehavior.runItem`).
 */
export interface HeaderMenuItem {
  /** Its entry's name. */
  readonly name: string;
  readonly kind: 'action' | 'state';
  readonly label: string;
  readonly icon?: IconName;
  /** A state's chip text. */
  readonly text?: string;
  /** Its command's key, as text (`⇧F10`). */
  readonly key?: string;
  /** Its command can run now; else `why` says why not. */
  readonly enabled: boolean;
  readonly why?: string;
  readonly checked?: boolean;
  /** Whether it is drawn in a slot (an action with an icon, in the first `header.slots`). */
  readonly slot?: boolean;
  readonly command?: string;
  readonly args?: Readonly<Record<string, unknown>>;
  /** A submenu. */
  readonly items?: readonly HeaderMenuItem[];
}

/** A header's menu: its target, and its groups (the element's entries, the track's, the section's), each non-empty. */
export interface HeaderMenu {
  readonly target: HeaderTarget;
  readonly groups: readonly (readonly HeaderMenuItem[])[];
}

/** A part's box and its entry's name. */
export interface HeaderPartBox {
  readonly name: string;
  readonly box: Box;
}

/**
 * A header's geometry (_spec/README.md "Behaviors", Headers: `viewer.headers`), in the canvas's CSS px as the hits are:
 * its target, its line's box, its label's, each slot's (the "⋯" `seq:menu` among them) and each chip's with its
 * entry's name, its track's well, and whether it is on the page.
 */
export interface HeaderGeometry {
  readonly target: HeaderTarget;
  readonly line: Box;
  readonly label: Box;
  readonly slots: readonly HeaderPartBox[];
  readonly states: readonly (HeaderPartBox & { readonly text: string; readonly mark: boolean })[];
  readonly well?: Box;
  /** Where the track's label text ends (its x), and the width of an ellipsis in its font: for the actions' backdrop. */
  readonly labelEnd?: number;
  readonly ellipsis?: number;
  readonly shown: boolean;
}

/** A header's parts for a host's DOM (`headerParts`): its geometry, whether its actions show now, and its items. */
export interface HeaderPartsHeader {
  readonly geometry: HeaderGeometry;
  readonly actions: boolean;
  readonly items: readonly HeaderMenuItem[];
}

/** What a parts port is given (`HeadersBehavior.parts`): every shown header's parts. */
export interface HeaderPartsState {
  readonly headers: readonly HeaderPartsHeader[];
}

// ---- the built-ins ----

/** The layer "Color ▸" acts on: the first with more than one fitting scheme. */
const colorLayer = (f: HeaderFacts) => f.layers.find((l) => l.schemes.length > 1 && !l.hidden);
/** A track's place among its siblings, by id, and the id the move puts it after (null: first). */
const moveAfter = (f: HeaderFacts, delta: -1 | 1): string | null | undefined => {
  const id = f.address.slice(f.address.lastIndexOf('/') + 1);
  const i = f.tracks.indexOf(id);
  if (i < 0 || i + delta < 0 || i + delta >= f.tracks.length) return undefined;
  return delta > 0 ? f.tracks[i + 1] : i >= 2 ? f.tracks[i - 2] : null;
};
const count = (n: number, one: string, many: string) => (n === 1 ? `1 ${one}` : `${n} ${many}`);

/** The built-in entries' names. */
export const HeaderEntryNames = {
  hide: 'seq:hide',
  color: 'seq:color',
  layers: 'seq:layers',
  moveUp: 'seq:move-up',
  moveDown: 'seq:move-down',
  filtered: HeaderStateNames.filtered,
  hiddenLayers: HeaderStateNames.hiddenLayers,
  hiddenTracks: 'seq:hidden-tracks',
  reverse: 'seq:reverse',
  collapse: 'seq:collapse',
  reverseBlock: 'seq:reverse-block',
  onlyBlock: 'seq:only-block',
  allBlocks: 'seq:all-blocks',
  tree: 'seq:tree',
} as const;

/**
 * The built-in header entries (_spec/README.md "Behaviors", Headers), registry entries like a host's: in
 * `DefaultSeqViewerSpec().registry`; a host leaves one out of the registry it spreads to remove it.
 */
export const BuiltInHeaderEntries: readonly HeaderActionEntry[] = /* @__PURE__ */ (
  [
    { name: HeaderEntryNames.hide, on: 'track', icon: 'eye-off', label: 'Hide track', command: 'track.hide', order: 1 },
    {
      name: HeaderEntryNames.color,
      on: 'track',
      label: 'Color',
      command: 'layer.channel',
      // the menu only ("Color ▸"): one place for it, never also a slot (a host that wants a slot registers its own entry)
      slot: false,
      order: 2,
      when: (f) => colorLayer(f) !== undefined,
      // its submenu: the command's choices (`layer.channel`'s select param, the layer's fitting schemes)
      args: (f) => ({ address: colorLayer(f)!.address, channel: 'color' }),
    },
    {
      name: HeaderEntryNames.layers,
      on: 'track',
      label: 'Layers',
      command: 'layer.hide',
      slot: false,
      order: 3,
      when: (f) => f.layers.filter((l) => !l.address.includes(':')).length > 1,
      args: (f) => ({ address: f.layers[0].address }),
      items: (f) =>
        f.layers
          .filter((l) => !l.address.includes(':'))
          .map((l) => ({
            label: l.label,
            command: l.hidden ? 'layer.show' : 'layer.hide',
            args: { address: l.address },
            checked: !l.hidden,
          })),
    },
    {
      name: HeaderEntryNames.moveUp,
      on: 'track',
      label: 'Move up',
      command: 'track.move',
      slot: false,
      order: 4,
      when: (f) => moveAfter(f, -1) !== undefined,
      args: (f) => ({ address: f.address, after: moveAfter(f, -1) }),
    },
    {
      name: HeaderEntryNames.moveDown,
      on: 'track',
      label: 'Move down',
      command: 'track.move',
      slot: false,
      order: 5,
      when: (f) => moveAfter(f, 1) !== undefined,
      args: (f) => ({ address: f.address, after: moveAfter(f, 1) }),
    },
    {
      // a narrow frame's tree: shown whole from its toggle, folded back from the toggle or here (`rows.tree`, Shift+T)
      name: HeaderEntryNames.tree,
      on: 'track',
      label: (f) => (f.tree?.open ? 'Fold tree' : 'Show tree'),
      command: 'rows.tree',
      slot: false,
      order: 6,
      when: (f) => f.tree?.lacks === true,
      args: (f) => ({ track: f.address }),
    },
    FilteredHeaderEntry,
    {
      name: HeaderEntryNames.hiddenLayers,
      on: 'track',
      kind: 'state',
      label: 'Show every layer',
      command: 'layer.show',
      when: (f) => f.layers.some((l) => l.hidden),
      args: (f) => ({ addresses: f.layers.filter((l) => l.hidden).map((l) => l.address) }),
      text: (f) => `${count(f.layers.filter((l) => l.hidden).length, 'layer', 'layers')} hidden`,
      items: (f) =>
        f.layers
          .filter((l) => l.hidden)
          .map((l) => ({ label: `Show ${l.label}`, command: 'layer.show', args: { address: l.address } })),
    },
    {
      // menu only: nothing in the view says a track is hidden (the application lists them, `mountHiddenTracks`)
      name: HeaderEntryNames.hiddenTracks,
      on: 'section',
      order: 5.5,
      label: 'Show hidden',
      command: 'track.show',
      slot: false,
      when: (f) => f.hiddenTracks.length > 0,
      args: (f) => ({ addresses: f.hiddenTracks.map((t) => t.address) }),
      items: (f) => [
        ...f.hiddenTracks.map((t) => ({ label: t.label, command: 'track.show', args: { address: t.address } })),
        ...(f.hiddenTracks.length > 1
          ? [{ label: 'Show all', command: 'track.show', args: { addresses: f.hiddenTracks.map((t) => t.address) } }]
          : []),
      ],
    },
    {
      name: HeaderEntryNames.reverse,
      on: 'section',
      order: 6,
      label: 'Reverse',
      command: 'section.reverse',
      slot: false,
      when: (f) => f.linear,
    },
    {
      name: HeaderEntryNames.collapse,
      on: 'section',
      order: 7,
      label: 'Collapse or expand hidden stretches',
      command: 'section.collapse',
      slot: false,
      when: (f) => f.linear,
    },
    // a block's (genome coordinates: its label's menu, Shift+F10 on it): turned over, shown alone, every block again
    {
      name: HeaderEntryNames.reverseBlock,
      on: 'section',
      order: 8,
      label: (f) => (f.block!.reversed ? `Turn ${f.block!.label} back` : `Reverse ${f.block!.label}`),
      command: 'block.reverse',
      slot: false,
      when: (f) => f.block !== undefined,
      args: (f) => ({ address: f.block!.address }),
    },
    {
      name: HeaderEntryNames.onlyBlock,
      on: 'section',
      order: 9,
      label: (f) => `Show only ${f.block!.name}`,
      command: 'section.include',
      slot: false,
      when: (f) => f.block !== undefined && !!f.assemblies,
      args: (f) => ({ address: f.section, include: [f.block!.name] }),
    },
    {
      name: HeaderEntryNames.allBlocks,
      on: 'section',
      order: 10,
      label: 'Show all',
      command: 'section.include',
      slot: false,
      when: (f) => !!f.included,
      args: (f) => ({ address: f.section, include: null }),
    },
  ] satisfies HeaderEntry[]
).map((entry) => ({ kind: 'header-action', entry }));

// ---- the ports (the `headers` behavior's params) ----

/** Where a menu is shown: a point in the canvas (CSS px from its top left). */
export interface HeaderMenuPlace {
  readonly x: number;
  readonly y: number;
}

/** The menu's port (_spec/README.md "Ports", Header menu): shows a menu at a place, one at a time. */
export interface HeaderMenuPort {
  show(menu: HeaderMenu, place: HeaderMenuPlace, behavior: HeadersBehavior): void;
  hide(): void;
  /** Its DOM's chunk asked for and not landed (the default's at its first open): `settle` waits for it. */
  loaded?(): Promise<void> | undefined;
  dispose?(): void;
  /** false: the host owns it and the viewer never disposes it. */
  readonly disposeWithViewer?: boolean;
}

/** Made once a viewer: with the canvas's element when mounted, else none (headless). */
export type HeaderMenuPortFactory = (container: HTMLElement | undefined) => HeaderMenuPort;

/** A host's DOM parts over the headers' boxes (_spec/README.md "Ports", Header parts). */
export interface HeaderPartsPort {
  update(parts: HeaderPartsState): void;
  /**
   * Resolves once its DOM is in place (the default's lazy chunk, asked for here): the canvas draws the parts until
   * then, and nothing in their boxes after (absent: at once). `settle` waits for it.
   */
  loaded?(): Promise<void>;
  /** The keyboard cursor came onto a header: its parts take the focus. */
  focus?(target: HeaderTarget): void;
  dispose?(): void;
  readonly disposeWithViewer?: boolean;
}

/** Made at the mount, with the canvas's element and the behavior its presses run through. */
export type HeaderPartsPortFactory = (container: HTMLElement, behavior: HeadersBehavior) => HeaderPartsPort;

export interface HeadersParams {
  /** The menu's adapter (the default spec gives `@molstar/ui`'s context menu; none: the menu's state moves alone). */
  readonly menu?: HeaderMenuPortFactory;
  /**
   * `canvas` (default): the parts drawn; a factory: a host's DOM parts over the laid-out boxes (the slots, the "⋯" and
   * the chips; the canvas draws none of them once the port's `loaded` resolved, and an export draws the states).
   */
  readonly parts?: 'canvas' | HeaderPartsPortFactory;
  /** What a secondary click opens: `everywhere` (default), `header` (the label column), `none` (the browser's). */
  readonly contextClick?: 'everywhere' | 'header' | 'none';
  /** A long press opens the menu after this many ms (500: touch's context menu). */
  readonly longPressMs?: number;
  /**
   * A submenu's flyout opens after the pointer rests on its row this many ms (300: long enough that a pointer passing
   * over the row on its way down the menu opens nothing, short enough to feel like hover; a click, →, Enter and Space
   * open it at once); the same rest on another row of the menu closes it, so a diagonal move into it keeps it. 0 turns
   * resting off.
   */
  readonly submenuRestMs?: number;
  /** A count chip's click: `all` (default) shows every hidden thing in one patch; `menu` lists them. */
  readonly hiddenChip?: 'all' | 'menu';
}

/** The menu's data for a target (DOM-free; _spec/README.md "Behaviors", Headers): undefined without the behavior. */
export function headerMenu(viewer: SeqViewer, target: HeaderTarget): HeaderMenu | undefined {
  return (viewer.behaviors.get('headers') as HeadersBehavior | undefined)?.menuOf(target);
}

/** Every shown header's parts (DOM-free): what a parts port is given. */
export function headerParts(viewer: SeqViewer): HeaderPartsState | undefined {
  return (viewer.behaviors.get('headers') as HeadersBehavior | undefined)?.partsState();
}

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { PartialCanvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import type { PluginConfigItem } from '@molstar/plugin/config';
import type { PluginRegistryEntry, PluginSpec } from '@molstar/plugin/spec';
import type {
  CommandDefinition,
  DockId,
  IconName,
  KeyChord,
  ShellOptions,
  ShellViewContext,
  UITheme,
  ViewRenderer,
} from '@molstar/prism-ui';
import type { ApiArea, ApiAreaContext, RunOrigin, ViewerApi } from './api/area.js';
import type { PrismViewerOptionValues } from './options.js';
import type { SceneNodeKind } from './model/scene.js';

// What a Prism viewer is made from (viewer-spec.md): plain data, explicit ordered lists, nothing resolved by name, a
// fresh default from a function. Types only: no DOM, no value import of the plugin.

export interface PrismViewerSpec {
  /** The engine's base, before the extensions' parts (§2.1). */
  readonly plugin: PrismViewerPluginParts;
  /** One PD definition read by the JS options, the URL and Settings (§2.2). */
  readonly options?: PartialOptions;
  /** In order: the viewer's built-in features first, then Mol*'s and the app's (§3). */
  readonly extensions?: readonly (ViewerExtension | LazyViewerExtension)[];
  /** Replace or remove any part by id (§4). */
  readonly overrides?: Readonly<Record<PartId, Part | null>>;
  /** The one sequence viewer, or none (§2.3). */
  readonly sequence?: PrismViewerSequenceSpec | false;
  /** The tools and skills, and the assistant that runs them (§2.4). */
  readonly agent?: PrismViewerAgentSpec | false;
  /** Layout, theme, keybindings; ignored headless (§2.5). */
  readonly ui?: PrismViewerUISpec;
  /** Light state in localStorage and sessionStorage, kept records (§2.6). */
  readonly persist?: PrismViewerPersist;
  /** The app's last word on the composed plugin spec (§2.1). */
  readonly finalize?: (spec: PluginSpec) => PluginSpec;
}

/** Options by group, each group partial. */
export type PartialOptions = { readonly [G in keyof PrismViewerOptionValues]?: Partial<PrismViewerOptionValues[G]> };

/** §2.1: the engine's base. */
export interface PrismViewerPluginParts {
  readonly registry: readonly PluginRegistryEntry[];
  readonly behaviors: readonly PluginSpec.Behavior[];
  readonly config?: readonly (readonly [PluginConfigItem, unknown])[];
  readonly canvas3d?: PartialCanvas3DProps;
}

/** §2.3 */
export interface PrismViewerSequenceSpec {
  /** The viewer's sequence behaviors (follow, secondary structure); extensions add theirs. */
  readonly behaviors: readonly SequenceBehavior[];
  /** The initial style; Tracks by default. */
  readonly style?: 'tracks' | 'text';
}

/**
 * A sequence behavior: what keeps part of the sequence's document in step with the plugin's state (sequence.md §3).
 * Its shape is settled with the sequence on the API (plan phase 2.4).
 */
export interface SequenceBehavior {
  readonly id: PartId;
  readonly title: string;
}

/** §2.4 */
export interface PrismViewerAgentSpec {
  /** Tool ids offered from the first turn (the viewer's default core). */
  readonly core?: readonly string[];
  /** The viewer's tools (agent.md §3); extensions add theirs. */
  readonly tools?: readonly ViewerTool[];
  /** The viewer's skills; extensions add theirs. */
  readonly skills?: readonly ViewerSkill[];
  /** The UI: `() => import('@molstar/prism-assistant/viewer')`. */
  readonly assistant?: () => Promise<unknown>;
}

/** A tool over the viewer's API (agent.md §3); typed with `@molstar/ai-tools`' registrations in plan phase 4. */
export interface ViewerTool {
  readonly id: string;
}

/** A skill: guidance and tools, disclosed by name (agent.md §3). */
export interface ViewerSkill {
  readonly id: string;
  readonly title: string;
  readonly summary: string;
}

/** §2.5 */
export interface PrismViewerUISpec {
  readonly layout?: Omit<ShellOptions, 'buttons' | 'viewState' | 'navigation'> & {
    readonly tabs?: Readonly<Record<string, { readonly dock?: DockId; readonly order?: number }>>;
    /** The tab open at start; none: no panel open. */
    readonly start?: string;
  };
  readonly theme?: UITheme;
  readonly keybindings?: readonly { readonly command: string; readonly key: KeyChord | null }[];
}

/** §2.6 */
export interface PrismViewerPersist {
  /** The namespace of every stored key ('prism-viewer'); one per live viewer on a page. */
  readonly key: string;
  /** localStorage: layout, theme, settings (small, lasting). */
  readonly preferences?: PrismStore;
  /** sessionStorage: this tab's light restore after a reload (journal.md §6). */
  readonly restore?: PrismStore;
  /** IndexedDB: a session record the person keeps (journal.md §1). */
  readonly rows?: unknown;
}

export interface PrismStore {
  get(key: string): unknown;
  set(key: string, value: unknown): void;
  remove(key: string): void;
}

// --- Extensions (§3) ---------------------------------------------------------------------------------------------

/** `<owner>.<kind>.<name>`: `core.toolbar.xr`, `snapshots.tab.snapshots`, `geo-export.section.export` (§4). */
export type PartId = string;

export type PartKind =
  | 'tab'
  | 'section'
  | 'row'
  | 'toolbar'
  | 'overlay'
  | 'command'
  | 'tool'
  | 'skill'
  | 'sequence'
  | 'api'
  | 'settings';

/** The places an extension's views go: closed and documented; the viewer decides where things go. */
export type ViewerSlot =
  | 'home.import'
  | 'export'
  | 'selection'
  | 'canvas-settings'
  | 'animation'
  | 'settings'
  | `inspector:${SceneNodeKind}`;

export interface ViewerExtension {
  /** 'geo-export'; every part id starts with it. */
  readonly id: string;
  readonly title: string;
  /** The engine side: registry entries, behaviors, config. */
  readonly plugin?: Partial<PrismViewerPluginParts>;
  /** Its operations, under its id (api.md §6): an area bound to an origin, made for each origin as the core's are. */
  readonly api?: (ctx: ApiAreaContext, origin: RunOrigin, group?: object) => ApiArea;
  readonly commands?: (api: ViewerApi) => readonly CommandDefinition<any>[];
  readonly tools?: readonly ViewerTool[];
  /** Never core tools: skills are disclosed. */
  readonly skills?: readonly ViewerSkill[];
  readonly sequence?: readonly SequenceBehavior[];
  /** Its Settings group. */
  readonly settings?: PD.Params;
  /** Views in the viewer's slots; a function, so the views are a chunk of their own. */
  readonly ui?: () => Promise<ViewerExtensionUI> | ViewerExtensionUI;
}

export interface ViewerExtensionUI {
  readonly sections?: readonly {
    readonly id: PartId;
    readonly slot: ViewerSlot;
    readonly title: string;
    readonly icon?: IconName;
    readonly order?: number;
    readonly keywords?: readonly string[];
    readonly render: ViewRenderer<ShellViewContext>;
  }[];
  readonly tabs?: readonly {
    readonly id: PartId;
    readonly title: string;
    readonly icon: IconName;
    readonly dock?: DockId;
    readonly order?: number;
    readonly render: ViewRenderer<ShellViewContext>;
  }[];
  /** Cards over the canvas, opened by a toolbar button or a command. */
  readonly overlays?: readonly {
    readonly id: PartId;
    readonly title: string;
    readonly render: ViewRenderer<ShellViewContext>;
  }[];
  /** Structure rows. */
  readonly rows?: readonly {
    readonly id: PartId;
    readonly when: (node: unknown) => boolean;
    readonly render: ViewRenderer<ShellViewContext>;
  }[];
  /** The canvas stack only: a command by id, or a card to open. */
  readonly toolbar?: readonly {
    readonly id: PartId;
    readonly icon: IconName;
    readonly title: string;
    readonly command?: string;
    readonly overlay?: PartId;
    readonly group?: string;
    readonly order?: number;
    readonly place?: 'rail' | 'canvas';
  }[];
}

export interface LazyViewerExtension {
  readonly id: string;
  readonly title: string;
  /** Every part it will bring, and what it answers for. */
  readonly stub: LazyViewerExtensionStub;
  readonly load: () => Promise<ViewerExtension>;
  /** Loaded during creation, before `ready`. */
  readonly preload?: boolean;
}

export interface LazyViewerExtensionStub {
  /** Checked by overrides; a loaded extension that brings a part not listed here fails its load. */
  readonly parts: readonly {
    readonly id: PartId;
    readonly kind: PartKind;
    readonly slot?: ViewerSlot;
    readonly title?: string;
  }[];
  /** Calling one of these loads the chunk and replays the call. */
  readonly operations?: readonly { readonly id: string; readonly title: string }[];
  readonly commands?: readonly {
    readonly id: string;
    readonly title: string;
    readonly category?: string;
    readonly icon?: IconName;
  }[];
  /** In the catalog; disclosing one loads the chunk. */
  readonly skills?: readonly ViewerSkill[];
  /** A session or file that needs one of these loads it first. */
  readonly owns?: {
    readonly transformers?: readonly string[];
    readonly formats?: readonly string[];
    readonly fileExtensions?: readonly string[];
    readonly representations?: readonly string[];
    readonly themes?: readonly string[];
  };
}

/** What an override puts in a part's place: another of its kind. */
export type Part = unknown;

export function isLazyExtension(e: ViewerExtension | LazyViewerExtension): e is LazyViewerExtension {
  return 'stub' in e && typeof (e as LazyViewerExtension).load === 'function';
}

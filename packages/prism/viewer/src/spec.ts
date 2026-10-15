/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ModelSource } from '@molstar/ai-session';
import type { CommandMenuAssistant } from '@molstar/prism-ui';
import type { Child } from '@molstar/prism-ui/dom/jsx-types';
import type { Observable } from 'rxjs';
import type { createViewerAgent, ViewerAgent } from './agent/session.js';
import type { PrismViewer } from './viewer.js';
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
import type { AnnotationProvider } from './api/annotations.js';
import type { ApiArea, ApiAreaContext, ApiModuleKey, OperationDefinition, ViewerApi } from './api/area.js';
import type { PrismViewerOptionValues } from './options.js';
import type { SceneNode, SceneNodeKind } from './model/scene.js';
import type { SequenceBehavior } from './sequence/behaviors.js';
import type { SkillDeclaration } from './agent/toolkit.js';
import type { OperationToolSpec } from './agent/tools.js';

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
  /**
   * Seq* and the views over it, its own chunk (`() => import('@molstar/prism-viewer/sequence/runtime')`): loaded when a
   * structure is first shown, or an operation or a saved state needs it; without a sequence in the spec, none of it is
   * bundled.
   */
  readonly runtime: () => Promise<typeof import('./sequence/runtime.js')>;
  /** The viewer's sequence behaviors (follow, secondary structure); extensions add theirs. */
  readonly behaviors: readonly SequenceBehavior[];
  /** The initial style; Tracks by default. */
  readonly style?: 'tracks' | 'text';
}

/** A sequence behavior: what keeps part of the sequence's document in step with what it shows (sequence.md §3). */
export type {
  SequenceBehavior,
  SequenceBehaviorContext,
  SequenceOffer,
  SequenceShown,
  SequenceTrackMenuItem,
} from './sequence/behaviors.js';

/** §2.4 */
export interface PrismViewerAgentSpec {
  /** The tools offered from the first turn, over operations (default the viewer's `CoreTools`, agent.md §3). */
  readonly core?: readonly ViewerToolSpec[];
  /** The viewer's skills (default its `Skills`); extensions add theirs. */
  readonly skills?: readonly ViewerSkill[];
  /** The UI: `() => import('@molstar/prism-assistant/viewer')` (its `createAssistant`). */
  /**
   * Its module brings the agent too (`createViewerAgent`, from `@molstar/prism-viewer/agent/session`): the viewer itself
   * never imports the agent's code (the AI SDK, zod), so a page without an assistant, the single file, is without it.
   */
  readonly assistant?: () => Promise<{
    readonly createAssistant: AssistantFactory;
    readonly createViewerAgent: typeof createViewerAgent;
  }>;
}

/**
 * What the assistant gets (§2.4): the viewer, and the agent's factory, which it calls with its models (the keys are the
 * person's, kept by the assistant: agent.md §4).
 */
export interface AssistantHost {
  readonly viewer: PrismViewer;
  readonly createAgent: (models: ModelSource) => ViewerAgent;
}

/** The assistant's parts in the page: its views, the command menu's Ask row, its turn card under the status line. */
export interface AssistantUI {
  readonly ui: ViewerExtensionUI;
  readonly commandMenu: CommandMenuAssistant;
  readonly agentView: () => Child;
  /** Its commands (the Inspector's: Assistant trace, tools, limits), registered once it has loaded. */
  readonly commands?: readonly CommandDefinition<any>[];
  dispose(): void;
}

export type AssistantFactory = (host: AssistantHost) => Promise<AssistantUI> | AssistantUI;

/** A tool over an operation, in the model's words (agent.md §2). */
export type ViewerToolSpec = OperationToolSpec;

/**
 * A skill: guidance and tools over operations, disclosed by name (agent.md §3). An extension's skills are over its own
 * operations (its `feature` is the extension); a lazy one's are in the catalog and load it when disclosed.
 */
export type ViewerSkill = SkillDeclaration;

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
  | 'page'
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
  /** Besides an image, what renders the view (an animation): the Image section's switch; its subject an `ExportFrame`. */
  | 'export.picture'
  | 'selection'
  | 'canvas-settings'
  | 'animation'
  | 'settings'
  | `inspector:${SceneNodeKind}`;

/** The Image section's shared frame, an `export.picture` view's subject: the resolution and the axes the person set. */
export interface ExportFrame {
  /** `viewport`, or the frame's pixels (a preset's or a custom size). */
  get(): ExportFrameValues;
}

export interface ExportFrameValues {
  readonly size: 'viewport' | 'custom';
  readonly width: number;
  readonly height: number;
  readonly axes: boolean;
  /** The crop, as fractions of the frame (the preview's box, auto or the person's). */
  readonly box: { readonly x: number; readonly y: number; readonly width: number; readonly height: number };
}

export interface ViewerExtension {
  /** 'geo-export'; every part id starts with it. */
  readonly id: string;
  readonly title: string;
  /** The engine side: registry entries, behaviors, config. */
  readonly plugin?: Partial<PrismViewerPluginParts>;
  /**
   * Its areas, by area id (api.md §6): each made once. Their operation ids start with the area's id; an extension's
   * areas are its own (`geo-export`), a built-in feature's may be named for what they do (`measure`, `snapshot`).
   */
  readonly api?: Readonly<Record<string, (ctx: ApiAreaContext) => ApiArea>>;
  /**
   * Operations it adds to core areas (the `loaders` feature's `load.entry`, `canvas-settings`' `canvas.lighting`):
   * built-in features only (`builtIn`), so the core's areas do not grow with Mol*'s extensions (operations.md §2).
   */
  readonly contributes?: (ctx: ApiAreaContext) => readonly {
    readonly area: string;
    readonly name: string;
    readonly operation: OperationDefinition<any, any>;
  }[];
  /** One of the viewer's own features (viewer-spec.md §3.1). */
  readonly builtIn?: true;
  readonly commands?: (api: ViewerApi) => readonly CommandDefinition<any>[];
  /** Its tools, in skills: never core tools, they are disclosed (agent.md §3). */
  readonly skills?: readonly ViewerSkill[];
  readonly sequence?: readonly SequenceBehavior[];
  /** What it can show of a structure (`annotations.apply`; README decision 27). */
  readonly annotations?: readonly AnnotationProvider[];
  /** Its Settings group. */
  readonly settings?: PD.Params;
  /** Views in the viewer's slots; a function, so the views are a chunk of their own. */
  readonly ui?: (ctx: ViewerUIContext) => Promise<ViewerExtensionUI> | ViewerExtensionUI;
}

/** What an extension's views get: the viewer's API as the person, its sources and managers (no plugin). */
export interface ViewerUIContext {
  readonly api: ViewerApi;
  readonly sources: import('./sources/viewer-sources.js').ViewerSources;
  readonly managers: import('./api/managers.js').ViewerManagers;
}

export interface ViewerExtensionUI {
  readonly sections?: readonly {
    readonly id: PartId;
    readonly slot: ViewerSlot;
    readonly title: string;
    readonly icon?: IconName;
    readonly order?: number;
    readonly keywords?: readonly string[];
    /** Open until the person folds it, in a slot whose sections start folded (Export's Image). */
    readonly open?: boolean;
    readonly render: ViewRenderer<ShellViewContext>;
  }[];
  readonly tabs?: readonly {
    readonly id: PartId;
    readonly title: string;
    readonly icon: IconName;
    readonly dock?: DockId;
    readonly order?: number;
    /** On the rail: among the tabs (default), or at its foot with Settings (the assistant). */
    readonly place?: 'strip' | 'foot';
    readonly render: ViewRenderer<ShellViewContext>;
  }[];
  /** Cards over the canvas, opened by a toolbar button or a command. */
  readonly overlays?: readonly {
    readonly id: PartId;
    readonly title: string;
    readonly render: ViewRenderer<ShellViewContext>;
  }[];
  /**
   * Pages over the whole app, for what needs more room than a panel or a card (the assistant's Inspector): opened by a
   * command or a link (`navigation.setPage({ id, focus })`), closed with Esc or their ✕.
   */
  readonly pages?: readonly {
    readonly id: PartId;
    readonly title: string;
    readonly icon?: IconName;
    /** What goes in the header after the title (the page's own tabs). */
    readonly header?: ViewRenderer<ShellViewContext>;
    readonly render: ViewRenderer<ShellViewContext>;
  }[];
  /**
   * Rows in the scene tree, under the nodes `when` holds for (a structure's assembly symmetry, its membrane): the view
   * gets the node as its `subject`; `SceneRow` draws it as the tree's own rows.
   */
  readonly rows?: readonly {
    readonly id: PartId;
    readonly when: (node: SceneNode) => boolean;
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
    /** Shown while true (always, when not given). */
    readonly visible?: Observable<boolean>;
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
  /** One of the viewer's own features: its stub may declare operations in core areas (as its `contributes`). */
  readonly builtIn?: true;
}

export interface LazyViewerExtensionStub {
  /** Checked by overrides; a loaded extension that brings a part not listed here fails its load. */
  readonly parts: readonly {
    readonly id: PartId;
    readonly kind: PartKind;
    readonly slot?: ViewerSlot;
    readonly title?: string;
    /** As the real part has them, so its placeholder stands where it will be. */
    readonly icon?: IconName;
    readonly order?: number;
    readonly place?: 'strip' | 'foot';
  }[];
  /**
   * Calling one of these loads the chunk and replays the call. Their flags are declared too, so tools and commands see
   * the policy before the load (checked against the real ones when it loads).
   */
  readonly operations?: readonly LazyOperationStub[];
  readonly commands?: readonly {
    readonly id: string;
    readonly title: string;
    readonly category?: string;
    readonly icon?: IconName;
  }[];
  /** In the catalog; disclosing one loads the chunk (its tools are over the stub's declared operations). */
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

/** A lazy extension's operation before it loads: its id, title and the flags its definition will have. */
export interface LazyOperationStub {
  readonly id: string;
  readonly title: string;
  readonly read?: true;
  /** A view change (no history step), as the real operation says. */
  readonly history?: false;
  readonly undoable?: false;
  readonly approval?: 'network' | 'upload' | 'destructive';
  /** What a call reaches (`OperationDefinition.addresses`): plain code, no import of the extension's chunk. */
  addresses?(values: Readonly<Record<string, any>>): readonly string[];
  readonly modules?: readonly ApiModuleKey[];
}

/** What an override puts in a part's place: another of its kind. */
export type Part = unknown;

export function isLazyExtension(e: ViewerExtension | LazyViewerExtension): e is LazyViewerExtension {
  return 'stub' in e && typeof (e as LazyViewerExtension).load === 'function';
}

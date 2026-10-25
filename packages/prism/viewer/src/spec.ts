/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { QueryLanguageEntry } from './targets/languages.js';
import type { CanvasMenuItem, ShellModel } from '@molstar/prism-ui/shell/shell';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ModelSource } from '@molstar/ai';
import type { CommandMenuAssistant } from '@molstar/prism-ui';
import type { CommandMenuPrefix } from '@molstar/prism-ui/commands/types';
import type { Child } from '@molstar/prism-ui/dom/jsx-types';
import type { Observable } from 'rxjs';
import type { DomainRules } from './agent/policies.js';
import type { createViewerAgent, ViewerAgent } from './agent/session.js';
import type { PrismViewer } from './viewer.js';
import type { PartialCanvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import type { PluginConfigItem } from '@molstar/plugin/config';
import type { PluginContext } from '@molstar/plugin/context';
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
import type { SequenceTrackProvider } from './sequence/tracks.js';
import type { AgentListChange, OperationToolSpec, SkillDeclaration } from './agent/skills.js';

// What a Prism viewer is made from (viewer-spec.md): plain data, explicit ordered lists, nothing resolved by name, a
// fresh default from a function. Types only: no DOM, no value import of the plugin.

export interface PrismViewerSpec {
  /** The engine's base, before the extensions' parts (§2.1). */
  readonly plugin: PrismViewerPluginParts;
  /** One PD definition read by the JS options, the URL and Settings (§2.2). */
  readonly options?: PartialOptions;
  /** In order: the viewer's built-in features first, then Mol*'s and the app's (§3). */
  readonly extensions?: readonly (ViewerExtension | LazyViewerExtension)[];
  /**
   * The query languages target text may be written in (targets.md §1), PyMOL first if listed (the default for text):
   * only these transpilers are bundled. None: keywords, address lists and MolQL selectors only (§2.7).
   */
  readonly languages?: readonly QueryLanguageEntry[];
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
  /** The lists the data palettes take (§2.8); `options.representation.defaultPalette` picks one. */
  readonly palettes?: PrismViewerPalettes;
  /** The app's last word on the composed plugin spec (§2.1). */
  readonly finalize?: (spec: PluginSpec) => PluginSpec;
}

/** The options a spec gives: any group, any value in it, and any value of a nested group (`ui.buttons`), as `resolveOptions` merges them. */
export type PartialOptions = {
  readonly [G in keyof PrismViewerOptionValues]?: {
    readonly [K in keyof PrismViewerOptionValues[G]]?: PrismViewerOptionValues[G][K] extends Readonly<
      Record<string, boolean>
    >
      ? Partial<PrismViewerOptionValues[G][K]>
      : PrismViewerOptionValues[G][K];
  };
};

/** §2.1: the engine's base. */
export interface PrismViewerPluginParts {
  readonly registry: readonly PluginRegistryEntry[];
  readonly behaviors: readonly PluginSpec.Behavior[];
  readonly config?: readonly (readonly [PluginConfigItem, unknown])[];
  readonly canvas3d?: PartialCanvas3DProps;
}

/**
 * §2.8: the lists the data palettes take (accessibility.md), the categorical themes' colors (chains, entities, …), in
 * the order categories take them: each color `0xRRGGBB` or `'#rrggbb'`. `options.representation.defaultPalette` picks among them
 * (`prism` by default; `classic` is Mol*'s own many-distinct).
 */
export interface PrismViewerPalettes {
  /** Prism's own; Many-Distinct (Prism) when not given (Mol*'s green and orange first). */
  readonly prism?: readonly (number | string)[];
  /** The colorblind-safe one; Okabe–Ito when not given. */
  readonly colorblind?: readonly (number | string)[];
}

/** §2.3 */
export interface PrismViewerSequenceSpec {
  /**
   * Seq* and the views over it, its own chunk (`() => import('@molstar/prism-viewer/sequence/runtime')`): loaded when a
   * structure is first shown, or an operation or a saved state needs it; without a sequence in the spec, none of it is
   * bundled.
   */
  readonly runtime: () => Promise<typeof import('./sequence/runtime.js')>;
  /** The viewer's sequence behaviors; extensions add theirs. */
  readonly behaviors: readonly SequenceBehavior[];
  /**
   * The viewer's track providers (sequence.md §3: ruler, letters, secondary structure); extensions add theirs. Absent:
   * `DefaultSequenceTracks`.
   */
  readonly tracks?: readonly SequenceTrackProvider[];
  /** The initial style; Tracks by default. */
  readonly style?: 'tracks' | 'text';
  /**
   * Each structure's Sequence node in the structure tree (sequence.md §8.1); false: none (the band, its menu, the
   * commands and the API stay; Edit in structure tree goes). True by default.
   */
  readonly tree?: boolean;
}

/** A sequence behavior: what keeps part of the sequence's document in step with what it shows (sequence.md §3). */
export type {
  SequenceBehavior,
  SequenceBehaviorContext,
  SequenceOffer,
  SequenceShown,
  SequenceTrackMenuItem,
} from './sequence/behaviors.js';
/** A track provider: a track a structure's sequence can have, which the viewer puts in and takes out (sequence.md §3). */
export type { SequenceTrackContext, SequenceTrackProvider } from './sequence/tracks.js';

/** §2.4 */
export interface PrismViewerAgentSpec {
  /**
   * The tools offered from the first turn, over operations: the viewer's `CoreTools` with these added or removed by id
   * (`{ add: [tool], remove: ['load_url'] }`; agent.md §3). The declarations are plain data
   * (`@molstar/prism-viewer/agent/skills`), so a spec that changes them does not bring the toolkit.
   */
  readonly core?: AgentListChange<ViewerToolSpec>;
  /** The viewer's skills: its `Skills` with these added (an id of a default's replaces it) or removed; extensions add theirs. */
  readonly skills?: AgentListChange<ViewerSkill>;
  /**
   * The UI, loaded lazily: `assistant({ provider, model })` from `@molstar/prism-assistant/spec` (or
   * `() => import('@molstar/prism-assistant/viewer')`). Its module brings the agent too (`createViewerAgent`, from
   * `@molstar/prism-viewer/agent/session`): the viewer itself never imports the agent's code (the AI SDK, zod), so a
   * page without an assistant, the single file, is without it.
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
  /** The viewer's agent over these models; `rules`: the domains the person allowed (agent.md §4, the assistant's). */
  readonly createAgent: (models: ModelSource, options?: { readonly rules?: DomainRules }) => ViewerAgent;
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
  /** localStorage: the theme and the layout (small, lasting; journal.md §6). */
  readonly preferences?: PrismStore;
  /**
   * Keep the layout in `preferences` (docks, sizes, hidden tabs, what is open in them) and restore it at start
   * (journal.md §6). On by default, except for an embedded layout (`ui.layout.mode: 'embedded'` or `compactRail`): that
   * one is the host page's, not the person's.
   */
  readonly layout?: boolean;
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
  | 'selection-action'
  | 'node-action'
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
  /** Rows of the canvas's playback card, with the trajectory's (the snapshots', their transition's). */
  | 'playback'
  | `inspector:${SceneNodeKind}`;

/** The Image section's shared frame, an `export.picture` view's subject: the resolution and the axes the person set. */
export interface ExportFrame {
  /** `viewport`, or the frame's pixels (a preset's or a custom size). */
  get(): ExportFrameValues;
}

export interface ExportFrameValues {
  /** The viewport's, or its pixels (api.md §1: a choice with its settings is a union). */
  readonly size:
    | { readonly name: 'viewport' }
    | { readonly name: 'custom'; readonly params: { readonly width: number; readonly height: number } };
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
   * built-in features only (`builtIn`), so the core's areas do not grow with Mol*'s extensions (api.md §2).
   */
  readonly contributes?: (ctx: ApiAreaContext) => readonly {
    readonly area: string;
    readonly name: string;
    readonly operation: OperationDefinition<any, any>;
  }[];
  /** One of the viewer's own features (viewer-spec.md §3.1). */
  readonly builtIn?: true;
  /** The extensions it needs in the spec before it (by id: the \`annotations\` feature needs \`structure\`). */
  readonly requires?: readonly string[];
  readonly commands?: (api: ViewerApi) => readonly CommandDefinition<any>[];
  /**
   * Prefixes it adds to the command menu (`pymol `, ui's commands.md §5): typed at the start, they turn the menu into
   * the prefix's list. Made once the extension's operations exist; an eager extension's only (a lazy stub declares
   * none: the rows are made as the person types).
   */
  readonly commandPrefixes?: (ctx: CommandPrefixContext) => readonly CommandMenuPrefix[];
  /** Its tools, in skills: never core tools, they are disclosed (agent.md §3). */
  readonly skills?: readonly ViewerSkill[];
  /**
   * What it adds to the sequence (sequence.md §3): track providers (the way to add sequence data: a track the person
   * can add, remove and hide per structure) and behaviors (rules that follow the scene), ids `<id>.sequence.<name>`.
   */
  readonly sequence?: {
    readonly tracks?: readonly SequenceTrackProvider[];
    readonly behaviors?: readonly SequenceBehavior[];
  };
  /** What it can show of a structure (`annotations.apply`; README decision 27). */
  readonly annotations?: readonly AnnotationProvider[];
  /** Its Settings group. */
  readonly settings?: PD.Params;
  /** Views in the viewer's slots; a function, so the views are a chunk of their own. */
  readonly ui?: (ctx: ViewerUIContext) => Promise<ViewerExtensionUI> | ViewerExtensionUI;
  /**
   * A built-in feature's view models (\`builtIn\` only; viewer-spec.md §3.1): made when the viewer is, headless too, and
   * read as \`viewer.views\` (the Structure panel's \`scene\`, Home's \`homeActions\`).
   */
  readonly views?: (ctx: BuiltInViewsContext) => Partial<BuiltInViews>;
  /**
   * A built-in feature's part of the page's shell (\`builtIn\` only): its tabs, sections and commands as the core's own,
   * with the core's ids (\`home\`, \`structure\`, \`selection\`, \`help\`; \`core.tab.*\` as parts). Returns their removal.
   */
  readonly shell?: (ctx: BuiltInShellContext) => () => void;
}

/** The view models the built-in features make (\`ViewerExtension.views\`), read as \`viewer.views\`. */
export interface BuiltInViews {
  /** The scene tree and inspector (the Structure panel). */
  readonly scene: import('./model/scene.js').SceneModel;
  /** The representation presets, from the registry (the Structure panel's Preset menu). */
  readonly presets: import('./model/scene.js').PresetsModel;
  /** The look's scene defaults and the overrides kept (the Structure panel's sliders; api.md §2, "Look"). */
  readonly look: import('./model/look.js').SceneLookModel;
  /** Quick styles (the Structure panel's Styles). */
  readonly quickStyles: import('./model/scene.js').QuickStylesModel;
  /** Home's loading actions, from the registry. */
  readonly homeActions: import('@molstar/prism-ui/core').ReadonlyState<readonly import('./model/home.js').HomeAction[]>;
}

/** What a built-in feature's view models are made with: the viewer's scope (they go with it), its plugin and API. */
export interface BuiltInViewsContext {
  readonly scope: import('./sources/state.js').SourceScope;
  readonly plugin: PluginContext;
  /** The API as the person. */
  readonly person: ViewerApi;
  /** Opens a command as the command menu does (its form, or runs it). */
  readonly openCommand: (id: string, values?: Readonly<Record<string, unknown>>) => void;
  /** Each structure's sequence for the structure tree (sequence.md §8.1); none: the spec has no sequence. */
  readonly sequence?: import('./model/scene.js').SceneSequence;
}

/** What a built-in feature's part of the shell gets: the shell, the viewer's core, and the shell's own services. */
export interface BuiltInShellContext {
  readonly shell: ShellModel;
  readonly core: import('./viewer.js').ViewerCore;
  /** The API as the person. */
  readonly person: ViewerApi;
  /** The core's parts after \`overrides\` (viewer-spec.md §4): by kind, with their short ids (\`tab\`, \`section\`). */
  readonly own: <T extends { id: string }>(kind: string, parts: readonly T[]) => T[];
  /** Opens the command menu, at a group. */
  readonly openCommandMenu: (group?: string) => void;
  /** Opens files as a drop does (a session file replaces what is loaded). */
  readonly openFiles: (files: readonly File[]) => void;
}

/** What an extension's views get: the viewer's API as the person, its sources and managers (no plugin). */
export interface ViewerUIContext {
  readonly api: ViewerApi;
  readonly sources: import('./sources/viewer-sources.js').ViewerSources;
  readonly managers: import('./api/managers.js').ViewerManagers;
  /** The person's preferences (`persist.preferences`, journal.md §6): small and lasting; none, the session's only. */
  readonly preferences?: PrismStore;
}

/** What an extension's command-menu prefixes are made with: the API as the person, the person's preferences. */
export interface CommandPrefixContext {
  readonly api: ViewerApi;
  /** The person's preferences (`persist.preferences`); none, the session's only. */
  readonly preferences?: PrismStore;
  /**
   * Says a row's outcome on the status line (and the log), as a command's is: a row runs after the menu closes, so
   * this is where its result goes.
   */
  readonly post: (text: string, kind?: 'info' | 'success' | 'warning' | 'error') => void;
}

export interface ViewerExtensionUI {
  readonly sections?: readonly {
    readonly id: PartId;
    readonly slot: ViewerSlot;
    readonly title: string;
    /** Its header's words when they follow what it offers (Export's "Image & animation" once an extension adds one). */
    readonly heading?: (shell: ShellModel) => Observable<string>;
    readonly icon?: IconName;
    readonly order?: number;
    readonly keywords?: readonly string[];
    /** Open until the person folds it, in a slot whose sections start folded (Export's Image). */
    readonly open?: boolean;
    readonly render: ViewRenderer<ShellViewContext>;
  }[];
  readonly tabs?: readonly {
    readonly id: PartId;
    /** Its name: the rail's tooltip, and the panel's header unless `heading` says otherwise. */
    readonly title: string;
    /** The panel's header when the title is longer than it needs (the PyMOL / VMD console: "Console"). */
    readonly heading?: string;
    readonly icon: IconName;
    readonly dock?: DockId;
    readonly order?: number;
    /** On the rail: among the tabs (default), or at its foot with Settings (the assistant). */
    readonly place?: 'strip' | 'foot';
    /**
     * On the rail while true (always, when not given): a panel for what is loaded (a kinemage's), which an app with the
     * extension does not show until there is something in it. Not the person's hidden tabs: those they choose.
     */
    readonly visible?: Observable<boolean>;
    readonly render: ViewRenderer<ShellViewContext>;
  }[];
  /** Cards over the canvas, opened by a toolbar button or a command. */
  readonly overlays?: readonly {
    readonly id: PartId;
    readonly title: string;
    /** Stays open on a click outside it (Esc or its ✕ close it): what is read beside the scene (the Story). */
    readonly persistent?: boolean;
    /** Its view takes the card's whole height and scrolls its own parts. */
    readonly fill?: boolean;
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
    /**
     * The nodes it is for. A row whose `when` holds takes a slot in the tree (its rows are of one height), so its view
     * renders a row whenever `when` holds: one that renders nothing leaves a gap.
     */
    readonly when: (node: SceneNode) => boolean;
    /**
     * Told when `when` may answer otherwise without a change of the scene (a PAE fetched): the tree asks it again for
     * the nodes it shows. Not needed when `when` reads only the node and its structure.
     */
    readonly whenChanged?: import('@molstar/prism-ui/core').ReadonlyState<unknown>;
    readonly render: ViewRenderer<ShellViewContext>;
  }[];
  /**
   * Rows of the selection bar's menu (what to do with the selection: the clipping feature's Clip…), after the
   * viewer's own (Focus, Show in structure tree, Make a component) and before More in the Selection panel; by `order`.
   * Each runs a command by id (its form, as the command menu opens it), or its own `run`.
   */
  readonly selectionActions?: readonly {
    readonly id: PartId;
    readonly title: string;
    /** One line under the title: what it does. */
    readonly hint?: string;
    readonly icon?: IconName;
    readonly order?: number;
    readonly command?: string;
    readonly run?: () => void;
  }[];
  /**
   * Rows of the scene tree's node menus (and the inspector's, which shows the node's menu): the clipping feature's
   * Clip this… on a component or a representation. On the nodes of `kinds` that pass `when`, after the node's own rows
   * and before Remove, by `order`; `run` gets the node (its `id` is the state cell's ref for a component or a
   * representation).
   */
  readonly nodeActions?: readonly {
    readonly id: PartId;
    readonly title: string;
    readonly icon?: IconName;
    readonly order?: number;
    readonly kinds: readonly SceneNodeKind[];
    readonly when?: (node: SceneNode) => boolean;
    readonly run: (node: SceneNode) => void;
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
    /** Its own action (in place of a command or a card). */
    readonly run?: () => void;
    /** A face that changes with what it does (play, pause), in place of `icon` and `title`. */
    readonly face?: Observable<{ icon: IconName; title: string }>;
    /** Variants beside it on hover (canvas items): stop, close. */
    readonly menu?: readonly CanvasMenuItem[];
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
    /**
     * Files it reads that their extension does not tell apart (mmCIF holding particles, not a structure): the files
     * with one of these `extensions` are read, and `test` (plain code, no import of the extension's chunk) says
     */
    readonly sniff?: {
      readonly extensions: readonly string[];
      test(file: { readonly name: string; readonly ext: string }, content: string | Uint8Array): boolean;
    };
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

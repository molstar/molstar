/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map, type Observable, of } from 'rxjs';
import { formatChord, isMac } from '../commands/keys.js';
import type { CommandManager } from '../commands/manager.js';
import type { CommandDefinition } from '../commands/types.js';
import type { IconName } from '../components/icon/icon.js';
import { Dynamic, Show } from '../dom/control-flow.js';
import type { Child } from '../dom/jsx-types.js';
import { onCleanup } from '../dom/scope.js';
import { domView } from '../dom/view.js';
import { parseColor, themeForColor } from '../model/canvas-theme.js';
import { observableState } from '../reactive/observable-state.js';
import { FullscreenModel } from '../shell/fullscreen.js';
import { ShellModel, type ShellTab } from '../shell/shell.js';
import { StatusModel } from '../shell/status.js';
import { ThemeModel, type UITheme, UIThemes } from '../shell/theme.js';
import { CommandMenu, type CommandMenuAssistant, type CommandMenuPrefix } from '../views/command-menu/command-menu.js';
import { CanvasPlaceholder } from '../views/shell/canvas-placeholder.js';
import { type ModeChoice, ModeChooser } from '../views/shell/mode-chooser.js';
import { CommandManager as Commands } from '../commands/manager.js';
import { SettingsTab } from '../views/settings/settings-tab.js';
import { LogView } from '../views/shell/log-view.js';
import { ShellView } from '../views/shell/shell-view.js';
import { mockSettings } from './settings.js';

// The mock app, in two layers: the frame (the shell and its view, the status column, the theme, fullscreen, the
// command menu, Log and Settings, and the add-ons such as the agent) and the one mode the app runs: the viewer
// (@molstar/prism-viewer) or the story mode (@molstar/prism-stories). A mode is chosen before the app exists (the
// splash, `PrismAppStory`); switching later starts a new app in the other mode (design/modes.md).

/** What the frame gives its mode. */
export interface MockFrameContext {
  shell: ShellModel;
  status: StatusModel;
  loaded: boolean;
  /** The canvas background (a CSS color). */
  background?: string;
  openCommandMenu(group?: string): void;
  /** The app's modes (to offer the choice, as Home does while nothing is loaded) and the one running. */
  modes: readonly ModeChoice[];
  mode: string;
  /**
   * Starts a new app: in that mode at once (nothing loaded, nothing lost: Home's choice), or, without one, after asking
   * (the splash over the app, saying what goes).
   */
  switchMode(id?: string): void;
}

/**
 * What a mode brings to the frame: its commands and the canvas parts (all optional but the commands). The parts are
 * render functions: the frame renders them in their slots, in the slot's scope (a mode activates from a click, outside
 * any render).
 */
export interface MockModeParts {
  commands: CommandManager<any>;
  prefixes?: CommandMenuPrefix[];
  canvas?: () => Child;
  overlays?: () => Child;
  playback?: () => Child;
  context?: () => Child;
  canvasLead?: () => Child;
  canvasAside?: () => Child;
  empty?: Observable<boolean>;
  onDropFiles?: (files: readonly File[]) => void;
  dropFormats?: readonly string[];
  /** The tab shown when the mode starts. */
  start?: string;
  /** Called when the app goes. */
  dispose?(): void;
}

export interface MockMode {
  id: string;
  title: string;
  icon: IconName;
  /** What it is for, in a line (the mode chooser). */
  description?: string;
  activate(ctx: MockFrameContext): MockModeParts;
}

/**
 * What a module adds to the mock shell, the way an app would contribute it: tabs (one with an existing id replaces it,
 * as History replaces Log), the command menu's assistant, the canvas's agent slot, and the person's finished commands.
 */
export interface MockShellAddOnParts {
  tabs?: readonly ShellTab[];
  assistant?: CommandMenuAssistant;
  agent?: Child;
  commandDone?(message: string, commandId: string): void;
}

export type MockShellAddOn = (ctx: { shell: ShellModel; status: StatusModel; loaded: boolean }) => MockShellAddOnParts;

function withAddOnTabs(tabs: ShellTab[], extra: readonly ShellTab[] | undefined): ShellTab[] {
  if (!extra?.length) return tabs;
  const replaced = new Set(extra.map((t) => t.id));
  return [...tabs.filter((t) => !replaced.has(t.id)), ...extra];
}

/** The frame alone, with an empty canvas: for stories about the frame or what is added to it. */
export function mockBlankMode(): MockMode {
  return {
    id: 'blank',
    title: 'Blank',
    icon: 'square',
    activate: ({ background }) => ({
      commands: new Commands<object>({ context: observableState({}).asReadonly, formats: { get: () => undefined } }),
      canvas: () => <CanvasPlaceholder empty={false} background={background} />,
      empty: of(false),
    }),
  };
}

export interface ShellStoryProps {
  loaded: boolean;
  background?: string;
  /** Opens the Export card at a section (a mode registers the card). */
  openExportAt?: string;
  /** The right panel pinned beside the canvas. */
  pinned?: boolean;
  /** Opens with the command menu open. */
  commandMenu?: boolean;
  /** What an add-on brings (the agent module's mock: `@molstar/prism-agent/testing`). */
  addOn?: MockShellAddOn;
}

/** The app's side of switching (`PrismAppStory`): its modes, and starting a new app. */
export interface MockApp {
  modes: readonly ModeChoice[];
  switchMode(id?: string): void;
}

/**
 * The mock frame running one mode: Mol*'s own UI is `mockViewerMode` (@molstar/prism-viewer), the story mode
 * `mockStoryMode` (@molstar/prism-stories). Inside `PrismAppStory`, it can start a new app in another mode.
 */
export function ShellStory(props: ShellStoryProps & { mode: MockMode; app?: MockApp }) {
  const { mode, app } = props;
  const shell = new ShellModel({ docks: { right: { pinned: !!props.pinned } } });
  const status = new StatusModel();
  const commandMenuOpen = observableState(false);
  // The group the command menu opens in (Preset, Styles); none for ⌘K.
  const commandMenuGroup = observableState<string | undefined>(undefined);
  const openCommandMenu = (group?: string) => {
    commandMenuGroup.set(group);
    commandMenuOpen.set(true);
  };
  // The plugin root goes fullscreen (in the stories, the story page).
  const fullscreen = new FullscreenModel(() => document.documentElement);
  // The UI theme, starting from Storybook's (its toolbar sets the same attribute on the page).
  const theme = new ThemeModel((document.documentElement.dataset.mspTheme as UITheme | undefined) ?? 'light');
  theme.applyTo(document.documentElement);
  // What an add-on brings (the agent module's tabs, its Ask port and turn card): only when the story includes one.
  const addOn = props.addOn?.({ shell, status, loaded: props.loaded });
  const settings = mockSettings();

  // The frame's own tabs and buttons: in every mode.
  shell.registerTabs(
    withAddOnTabs(
      [
        {
          id: 'log',
          title: 'Log',
          icon: 'scroll-text',
          order: 2.5,
          // Not a working panel: with Settings at the rail's foot.
          place: 'foot',
          render: domView(() => <LogView status={status} />),
        },
        {
          id: 'settings',
          title: 'Settings',
          icon: 'settings',
          order: 3,
          place: 'foot',
          render: domView(() => <SettingsTab model={settings} />),
        },
      ],
      addOn?.tabs,
    ),
  );

  const parts = mode.activate({
    shell,
    status,
    loaded: props.loaded,
    background: props.background,
    openCommandMenu,
    modes: app?.modes ?? [mode],
    mode: mode.id,
    switchMode: (id) => app?.switchMode(id),
  });
  // The canvas tools act on what is shown: the mode says whether anything is.
  const empty = parts.empty ?? of(!props.loaded);
  const hasData = empty.pipe(map((e) => !e));
  shell.registerToolbar([
    // The UI theme, by Settings at the rail's foot: a click goes to the next one; the icon is the current one.
    {
      id: 'theme',
      icon: 'sun',
      title: 'Theme',
      face: theme.theme.pipe(
        map((t) => {
          const current = UIThemes.find((u) => u.id === t)!;
          return { icon: current.icon, title: `Theme: ${current.title}` };
        }),
      ),
      run: () => theme.next(),
    },
    {
      id: 'palette',
      icon: 'command',
      title: 'Command menu',
      caption: formatChord({ key: 'k', mod: true }),
      run: () => openCommandMenu(),
    },
    // Expand by context (FullscreenModel.mode): on a page, Expand fills it and Fullscreen is beside it; in an iframe
    // that allows it, the button is Fullscreen (expanding would only fill the frame); in one that does not, Expand fills
    // the frame and Fullscreen shows why it is not possible.
    fullscreen.mode === 'fullscreen'
      ? {
          id: 'expand',
          place: 'canvas',
          group: 'canvas',
          visible: hasData,
          icon: 'fullscreen' as const,
          title: 'Fullscreen',
          pressed: fullscreen.active,
          run: () => fullscreen.toggle(),
        }
      : {
          id: 'expand',
          place: 'canvas',
          group: 'canvas',
          visible: hasData,
          icon: 'maximize-2' as const,
          title: 'Expand',
          run: () => status.post('Expanded'),
          menu: [
            {
              id: 'fullscreen',
              label: 'Fullscreen',
              icon: 'fullscreen' as const,
              hint: 'The whole screen, not just the page',
              disabled: fullscreen.unavailableReason,
              pressed: fullscreen.active,
              run: () => fullscreen.toggle(),
            },
          ],
        },
  ]);

  // Switching modes, with more than one: a command that asks (the splash over the app, saying what goes).
  const switching: CommandDefinition<any>[] =
    app && app.modes.length > 1
      ? [
          {
            id: 'app.switch-mode',
            title: 'Switch mode…',
            alias: ['mode', 'viewer', 'story', 'start over'],
            category: 'View',
            icon: 'prism',
            run: () => {
              app.switchMode();
              return { message: '' };
            },
          },
        ]
      : [];
  if (switching.length) parts.commands.register(switching);
  // What commands did goes to the status line (and, with the agent, to History).
  const outcomes = parts.commands.outcomes.subscribe((o) => {
    if (o.kind === 'done' && o.result.message) status.post(o.result.message);
    if (o.kind === 'done') addOn?.commandDone?.(o.result.message ?? o.command.title, o.command.id);
    else if (o.kind === 'failed' || o.kind === 'blocked') status.post(o.message, 'warning');
  });

  // ⌘K / ⇧⌘P open the command menu; the mode's commands' own shortcuts run them.
  const onKey = (e: KeyboardEvent) => {
    const mod = isMac() ? e.metaKey : e.ctrlKey;
    if (mod && (e.key.toLowerCase() === 'k' || (e.shiftKey && e.key.toLowerCase() === 'p'))) {
      e.preventDefault();
      openCommandMenu();
    } else if (!commandMenuOpen.value && parts.commands.handleKey(e, isMac())) e.preventDefault();
  };
  window.addEventListener('keydown', onKey);
  if (parts.start) shell.show(parts.start);
  shell.mount();
  status.mount();
  if (props.openExportAt) shell.openOverlay('export', undefined, props.openExportAt);
  if (props.commandMenu) commandMenuOpen.set(true);
  onCleanup(() => {
    window.removeEventListener('keydown', onKey);
    outcomes.unsubscribe();
    parts.dispose?.();
    shell.unmount();
    status.unmount();
  });

  const background = props.background ?? '#ffffff';
  return (
    <ShellView
      shell={shell}
      status={status}
      canvas={parts.canvas?.() ?? null}
      overlays={parts.overlays?.() ?? null}
      playback={parts.playback?.() ?? null}
      canvasColor={background}
      canvasTheme={themeForColor(parseColor(background) ?? 0xffffff)}
      context={parts.context}
      commandMenu={
        <Show when={commandMenuOpen}>
          {() => (
            <CommandMenu
              manager={parts.commands}
              prefixes={parts.prefixes}
              group={commandMenuGroup.value}
              assistant={addOn?.assistant}
              onClose={() => {
                commandMenuOpen.set(false);
                commandMenuGroup.set(undefined);
              }}
            />
          )}
        </Show>
      }
      agent={addOn?.agent}
      canvasLead={parts.canvasLead?.() ?? null}
      canvasAside={parts.canvasAside?.() ?? null}
      onOpenLog={() => shell.show('log')}
      empty={empty}
      onDropFiles={parts.onDropFiles}
      dropFormats={parts.dropFormats}
    />
  );
}

/**
 * The app as it starts (design/modes.md): the splash asks which mode before any app exists; choosing one creates the
 * app in it. Switching later shows the splash over the app, saying what goes (what is loaded; stories are kept in the
 * browser); choosing creates a new app, Cancel goes back. With `start`, the app starts in that mode (no splash).
 */
export function PrismAppStory(props: ShellStoryProps & { modes: readonly MockMode[]; start?: string }) {
  const find = (id: string | undefined) => props.modes.find((m) => m.id === id);
  // Each start counted: the app is created again even in the same mode (states skip equal records).
  let starts = 0;
  const start = (mode: MockMode) => ({ mode, n: ++starts });
  const running = observableState<{ mode: MockMode; n: number } | undefined>(
    find(props.start) ? start(find(props.start)!) : undefined,
  );
  const asking = observableState(false);
  const app: MockApp = {
    modes: props.modes,
    switchMode: (id) => {
      const next = find(id);
      if (next) queueMicrotask(() => running.set(start(next)));
      else asking.set(true);
    },
  };
  const choose = (id: string) => {
    asking.set(false);
    const next = find(id);
    if (next && next.id !== running.value?.mode.id) queueMicrotask(() => running.set(start(next)));
  };
  return (
    <div style="position: relative; height: 100%">
      <Dynamic
        component={running.asReadonly.pipe(
          map((r) =>
            r
              ? () => <ShellStory {...props} mode={r.mode} app={app} />
              : () => (
                  <ModeChooser variant="splash" modes={props.modes} onChoose={(id) => running.set(start(find(id)!))} />
                ),
          ),
        )}
      />
      <Show when={asking}>
        {() => (
          <ModeChooser
            variant="welcome"
            modes={props.modes}
            active={running.asReadonly.pipe(map((r) => r?.mode.id ?? ''))}
            notice="Switching starts a new session: what is loaded now goes. Stories are kept in this browser."
            onChoose={choose}
            onCancel={() => asking.set(false)}
          />
        )}
      </Show>
    </div>
  );
}

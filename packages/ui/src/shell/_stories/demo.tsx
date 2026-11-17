/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState, State } from '@molstar/core/reactive';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ActionFooter } from '../../components/action-footer/action-footer.js';
import { Button } from '../../components/button/button.js';
import { Expandable, ExpandableList } from '../../components/expandable/expandable.js';
import { Icon } from '../../components/icon/icon.js';
import { GlassWordmark } from '../../components/logo/glass-wordmark.js';
import { ParamsForm } from '../../components/params/params-form.js';
import { SkeletonRows } from '../../components/skeleton/skeleton.js';
import { TextInput } from '../../components/text-input/text-input.js';
import type { ToolbarGroup } from '../../components/toolbar/toolbar.js';
import { type Child, domView, For, onCleanup, Show } from '../../dom/index.js';
import { BandPlaceholder, StagePlaceholder, structureLines } from '../../stage/placeholders.js';
import { Stage } from '../../stage/stage.js';
import { StageModel } from '../../stage/stage-model.js';
import { StageToolbar } from '../../stage/toolbar.js';
import { StatusChips, StatusLine } from '../../status/status-line.js';
import { StatusModel } from '../../status/status-model.js';
import { ParamsFormModel } from '../../forms/params-form.js';
import { PaneMap } from '../pane-map.js';
import { type PanePort, ShellModel, type ShellOptions, type ShellViewDef } from '../shell-model.js';
import { type ShellBusy, ShellView } from '../shell-view.js';
import { SectionList } from '../section-list.js';
import { ViewEmpty, ViewFooter, ViewHeading, ViewRow, ViewStack, ViewToolbar } from '../view-frame.js';

// Generic views for the shell's stories: what an app (or an extension) contributes, written only with the view
// primitives (shell foundations, styles §8): no padding, gaps or type of their own.

const form = <P extends PD.Params>(params: P) => (
  <ParamsForm form={new ParamsFormModel(params, PD.getDefaultValues(params))} />
);

const appearance = {
  quality: PD.Select('auto', [
    ['auto', 'Auto'],
    ['high', 'High'],
    ['low', 'Low'],
  ]),
  outline: PD.Boolean(true),
  opacity: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
};

const lighting = {
  ambient: PD.Numeric(0.4, { min: 0, max: 1, step: 0.01 }),
  key: PD.Numeric(0.6, { min: 0, max: 1, step: 0.01 }),
  shadows: PD.Boolean(true),
};

/** The demo's layers: a tree (level), each with a description; some long enough to be cut off in a panel. */
const layers: { name: string; description: string; level: number; icon: 'layers' | 'image' | 'tag' | 'ruler' }[] = [
  { name: 'Background', description: 'Gradient, 3 stops', level: 0, icon: 'image' },
  { name: 'Annotations of the selected region', description: '128 items, 4 hidden', level: 0, icon: 'tag' },
  { name: 'Labels', description: 'Names and identifiers, shown on hover', level: 1, icon: 'tag' },
  { name: 'Measurements', description: 'Distances and angles, 6 items', level: 1, icon: 'ruler' },
  { name: 'Selection', description: 'Highlighted with an outline', level: 0, icon: 'layers' },
  { name: 'Overlay with a very long name that does not fit', description: 'Drawn last', level: 0, icon: 'layers' },
];

function HomeView() {
  return (
    <ViewStack>
      <ViewHeading>Start</ViewHeading>
      <ViewEmpty action={<Button variant="primary">Open a file…</Button>}>
        Nothing is open yet. Open a file, or drop one anywhere.
      </ViewEmpty>
      <ViewHeading>Recent</ViewHeading>
      <div class="msui-list">
        {['session-03.json', 'notes.md', 'survey.csv'].map((name) => (
          <div class="msui-row">
            <span class="msui-slot">
              <Icon name="file" />
            </span>
            <span>{name}</span>
          </div>
        ))}
      </div>
    </ViewStack>
  );
}

function LayersView() {
  const filter = new State('');
  const shown = derive(filter, (f) => layers.filter((l) => l.name.toLowerCase().includes(f.toLowerCase())));
  return (
    <ViewStack>
      <ViewToolbar label="Layers">
        <TextInput value={filter} onInput={(v) => filter.set(v)} placeholder="Filter layers" />
        <Button variant="ghost" size="sm" icon={<Icon name="plus" />} title="Add a layer" />
      </ViewToolbar>
      <div class="msui-list" role="tree" aria-label="Layers">
        <For each={shown} by={(l) => l.name}>
          {(l) => (
            // Cut off, its name and description show whole in a card beside the view (the shell's side tooltip).
            <ViewRow
              role="treeitem"
              tabIndex={0}
              aria-level={String(l.level + 1)}
              level={l.level}
              icon={l.icon}
              label={l.name}
              description={l.description}
              trailing={<Button variant="ghost" size="xs" icon={<Icon name="eye" />} title={`Hide ${l.name}`} />}
            />
          )}
        </For>
      </div>
    </ViewStack>
  );
}

function InspectView() {
  return (
    <ExpandableList>
      <Expandable title="Appearance" open>
        {() => form(appearance)}
      </Expandable>
      <Expandable title="Lighting">{() => form(lighting)}</Expandable>
    </ExpandableList>
  );
}

/** A form and its footer: the footer after the stack (the body's own child), so it sits at the bottom. */
function SettingsView() {
  return [
    <ViewStack>
      {form({
        theme: PD.Select('auto', [
          ['auto', 'Auto'],
          ['light', 'Light'],
          ['dark', 'Dark'],
        ]),
        animations: PD.Boolean(true),
      })}
    </ViewStack>,
    <ViewFooter>
      <ActionFooter label="Apply" onApply={() => {}} onReset={() => {}} />
    </ViewFooter>,
  ];
}

function HelpView() {
  return (
    <ViewStack>
      <ViewHeading>Keys</ViewHeading>
      <div class="msui-list">
        {[
          ['F6', 'Next region'],
          ['Esc', 'Close the innermost'],
          ['Shift F10', 'A rail button’s menu'],
        ].map(([k, v]) => (
          <div class="msui-row msui-row-plain">
            <span>{v}</span>
            <kbd>{k}</kbd>
          </div>
        ))}
      </div>
    </ViewStack>
  );
}

/** The details of one item (an instance per item: two side by side to compare): its properties as rows. */
function DetailsView(props: { subject: ReadonlyState<unknown> }) {
  const name = derive(props.subject, (s) => String(s ?? 'Nothing selected'));
  const rows: [string, (n: string) => string][] = [
    ['Name', (n) => n],
    ['Kind', () => 'Layer'],
    ['Items', (n) => String(12 + n.length * 7)],
    ['Visible', () => 'Yes'],
    ['Updated', () => 'Today, 10:42'],
  ];
  return (
    <ViewStack>
      <ViewHeading>Properties</ViewHeading>
      <div class="msui-form">
        {rows.map(([label, value]) => (
          <>
            <span class="msui-label">{label}</span>
            <span class="msui-value">{derive(name, value)}</span>
          </>
        ))}
      </div>
    </ViewStack>
  );
}

function LoadingView() {
  return <SkeletonRows count={6} icon />;
}

/** The demo's views: a head view, strip views, foot views, a card with sections, a page, a modal, a plot by instance. */
export function demoViews(): ShellViewDef[] {
  return [
    {
      id: 'home',
      title: 'Home',
      // The Mol* wordmark as glass: it assembles each time Home opens.
      heading: () => <GlassWordmark />,
      icon: 'molstar',
      place: 'head',
      order: 0,
      compactBar: true,
      render: domView(HomeView),
    },
    { id: 'layers', title: 'Layers', icon: 'layers', order: 1, compactBar: true, render: domView(LayersView) },
    {
      id: 'inspect',
      title: 'Inspect',
      icon: 'sliders-horizontal',
      order: 2,
      compactBar: true,
      render: domView(InspectView),
    },
    { id: 'loading', title: 'Loading', icon: 'hourglass', order: 3, render: domView(LoadingView) },
    { id: 'settings', title: 'Settings', icon: 'settings', place: 'foot', order: 10, render: domView(SettingsView) },
    { id: 'help', title: 'Help', icon: 'circle-help', place: 'foot', order: 11, render: domView(HelpView) },
    {
      id: 'export',
      title: 'Export',
      icon: 'download',
      home: 'card',
      render: domView(({ ctx }) => <SectionList shell={ctx.shell} slot="export" ctx={ctx} />),
    },
    {
      id: 'details',
      title: 'Details',
      icon: 'info',
      home: 'float',
      instances: 'many',
      size: { width: 320, height: 240 },
      render: domView(({ ctx }) => <DetailsView subject={ctx.subject} />),
    },
    {
      id: 'report',
      title: 'Report',
      icon: 'scroll-text',
      home: 'page',
      render: domView(() => (
        <ViewStack>
          <ViewEmpty>A page: opaque over the main area, for what needs the room. Escape closes it.</ViewEmpty>
        </ViewStack>
      )),
    },
    {
      id: 'confirm',
      title: 'Discard changes?',
      icon: 'triangle-alert',
      home: 'modal',
      render: domView(({ ctx }) => [
        <ViewStack>
          <ViewEmpty>The layout and its windows go back to how they were.</ViewEmpty>
        </ViewStack>,
        <ViewFooter>
          <ActionFooter
            label="Discard"
            onApply={() => ctx.close()}
            secondary={{ label: 'Keep', onClick: () => ctx.close() }}
          />
        </ViewFooter>,
      ]),
    },
  ];
}

/** The main area's placeholder: what an app's view would be (a canvas), a quiet grid with its name. */
/**
 * The app's main area: a stage over a canvas stand-in, the status bottom-left, the stage's toolbar (Expand, Fullscreen,
 * Clear) bottom-right. Empty, it shows the glass mark; a drop "loads" it.
 */
/** Where the demo's stage toolbars go (shell §2.1): one a pane, one shared over the main area, or both split. */
export type DemoToolbars = 'pane' | 'shared' | 'mixed';

/** What a camera tool acts on (shell §2.1): a pane's stage, its pane (where its message goes), its name. */
interface CameraTarget {
  readonly stage: StageModel;
  /** Its scope in the status (its pane's id). */
  readonly pane: string;
  readonly title: string;
}

/** What the demo keeps a pane (app side, `PaneMap`): its stage and its name. */
interface DemoPane {
  readonly stage: StageModel;
  readonly title: string;
}

/**
 * The stage toolbar's camera tools, acting on `target()`: a pane's own stage (its toolbar), or the active pane's
 * (a shared toolbar: the app's `PaneMap` of the shell's active pane). The message goes to that pane's status line.
 */
function cameraGroup(target: () => CameraTarget | undefined, status: StatusModel, own?: StageModel): ToolbarGroup {
  // A pane's own toolbar greys out on its empty stage; a shared one says so when used (its pane changes under it).
  const disabled = own ? derive(own.empty, (e) => (e ? 'Nothing is loaded' : false)) : false;
  const act = (verb: string, fn?: (s: StageModel) => void) => () => {
    const t = target();
    if (!t) return;
    if (t.stage.empty.value) {
      status.post(`Nothing is loaded in ${t.title}`, 'warning', { scope: t.pane });
      return;
    }
    fn?.(t.stage);
    status.post(`${verb}: ${t.title}`, 'info', { scope: t.pane });
  };
  return {
    id: 'scene',
    label: 'Scene',
    items: [
      { id: 'reset', icon: 'refresh-cw', title: 'Reset the camera', disabled, run: act('Camera reset') },
      { id: 'clear', icon: 'trash-2', title: 'Clear', disabled, run: act('Cleared', (s) => s.setEmpty(true)) },
    ],
  };
}

/** The view tools: Expand, Fullscreen (the shell's). */
/**
 * The view tools (shell §3.9): Expand, Fullscreen. A pane's own toolbar expands its pane (`pane`); a toolbar over the
 * whole main area expands it all (no pane: the global expand).
 */
function viewGroup(shell: ShellModel, pane?: PanePort): ToolbarGroup {
  // A pane's own toolbar maximizes its pane (the layout's, through its port); one over the whole main area expands the
  // shell (the shell's chrome away). Both from one button would be the app's choice: call both.
  const item = pane
    ? {
        id: 'maximize',
        icon: 'maximize' as const,
        title: 'Maximize',
        pressed: derive(pane.status, (p) => p.maximized),
        run: () => pane.toggleMaximized(),
      }
    : {
        id: 'expand',
        icon: 'maximize-2' as const,
        title: 'Expand the main area',
        pressed: shell.expanded,
        run: () => shell.toggleExpanded(),
      };
  return {
    id: 'view',
    label: 'View',
    items: [
      item,
      {
        id: 'fullscreen',
        icon: 'fullscreen',
        title: 'Fullscreen',
        pressed: shell.fullscreen,
        disabled: derive(shell.fullscreenAvailable, (a) => (a === true ? false : a)),
        run: () => shell.requestFullscreen(),
      },
    ],
  };
}

/**
 * A stage (the app's content, or a stage view beside it): a canvas stand-in, its own toolbar as the placement says,
 * the status where the app put it (the content's stage, or over the main area).
 */
function DemoStage(props: {
  shell: ShellModel;
  stage: StageModel;
  status: StatusModel;
  name: string;
  /** Its pane's port (the layout's side, shell §2.1): its id (its messages' scope), status, insets, Expand. */
  port: PanePort;
  toolbars: DemoToolbars;
  /** Its status line's scope: its pane's messages, and with `rest` the rest (the one line that takes it). */
  rest?: boolean;
  top?: () => Child;
  colors?: readonly string[];
  /** The toolbar only while the pane is active, hovered or focused (its port's status: the stage decides). */
  activeOnly?: boolean;
}) {
  const { shell, stage, status } = props;
  const own = (): CameraTarget => ({ stage, pane: props.port.id, title: props.name });
  const groups: ToolbarGroup[] =
    props.toolbars === 'pane'
      ? [cameraGroup(own, status, stage), viewGroup(shell, props.port)]
      : props.toolbars === 'mixed'
        ? [cameraGroup(own, status, stage)]
        : [];
  return (
    <Stage
      model={stage}
      top={props.top?.()}
      canvas={
        <Show
          when={stage.empty}
          fallback={() => (
            <StagePlaceholder
              colors={props.colors}
              onHover={(i) => stage.setHover(i === undefined ? [] : structureLines(i))}
            />
          )}
        >
          {() => <StagePlaceholder empty />}
        </Show>
      }
      status={<StatusLine status={status} scope={{ id: props.port.id, rest: props.rest }} />}
      chips={props.rest ? <StatusChips status={status} /> : undefined}
      corner={groups.length ? <StageToolbar groups={groups} stage={stage} /> : undefined}
      toolbarShown={props.activeOnly ? derive(props.port.status, (p) => p.active || p.hovered || p.focused) : true}
      insets={props.port.insets}
    />
  );
}

/** Another look a stage view's instance: its chains' colors by its key. */
const PaneColors: readonly (readonly string[])[] = [
  ['#5fae8b', '#c98a4b', '#7b74c9'],
  ['#c56b8a', '#5b9bd6', '#d1b04a'],
  ['#8a8f99', '#d66b5b', '#4bb3c4'],
];

/**
 * The demo's stage views (shell §2.1): a stage beside the app's (`stage`, one a key, its own model), a sequence for the
 * top (overlaid or split by the shell's option), a panel for a pane beside the stage.
 */
function stageViews(
  shell: ShellModel,
  status: StatusModel,
  panes: PaneMap<DemoPane>,
  toolbars: DemoToolbars,
  activeOnly: boolean,
): ShellViewDef[] {
  let n = 0;
  return [
    {
      id: 'stage',
      title: 'Stage',
      icon: 'box',
      stage: true,
      home: 'stage',
      instances: 'many',
      fill: true,
      render: domView(({ ctx }) => {
        const model = new StageModel({ empty: false });
        onCleanup(() => model.dispose());
        const name = `Stage ${(ctx.key ?? '').toUpperCase()}`;
        // Its title to the layout (its port: tabs, a window's header); its stage to the app's map (a shared toolbar).
        onCleanup(ctx.pane.setTitle(name));
        onCleanup(panes.set(ctx.ref, { stage: model, title: name }));
        return (
          <DemoStage
            shell={shell}
            stage={model}
            status={status}
            name={name}
            port={ctx.pane}
            toolbars={toolbars}
            colors={PaneColors[n++ % PaneColors.length]}
            activeOnly={activeOnly}
          />
        );
      }),
    },
    {
      id: 'sequence',
      title: 'Sequence',
      icon: 'panel-top',
      home: 'top',
      fill: true,
      render: domView(() => <BandPlaceholder rows={2} length={160} selection={[11, 18]} />),
    },
    {
      id: 'panel',
      title: 'Panel',
      icon: 'panel-right',
      home: 'stage',
      stageSide: 'right',
      stageSize: 280,
      render: domView(() => (
        <ViewStack>
          <ViewEmpty>A view docked beside the stage.</ViewEmpty>
        </ViewStack>
      )),
    },
  ];
}

export interface DemoOptions extends ShellOptions {
  /** Called with the shell once its views are registered (to put it in a state). */
  setup?(shell: ShellModel): void;
  busy?: ReadonlyState<ShellBusy | undefined>;
  width?: string;
  height?: string;
  /** The stage starts empty (the glass mark) until files are dropped. */
  empty?: boolean;
  /** Where the stage toolbars go (default one a pane). */
  toolbars?: DemoToolbars;
  /** A pane's toolbar only on the active or hovered stage. */
  activeOnly?: boolean;
  /** The Mol* star that shows a hidden rail (default on). */
  railToggle?: boolean;
  /** What the drop overlay says. */
  dropLabel?: string;
  /** What can be dropped (a hint under the label). */
  dropFormats?: string;
  /** How far Expand reaches. */
  fill?: 'shell' | 'page';
}

/** A shell with the demo's views, sections and rail items, in a box of the given size. */
export function DemoShell(props: { options?: DemoOptions; children?: Child }) {
  const o = props.options ?? {};
  const shell = new ShellModel(o);
  const toolbars = o.toolbars ?? 'pane';
  const stage = new StageModel({ empty: o.empty ?? false });
  const status = new StatusModel();
  onCleanup(() => {
    stage.dispose();
    status.retire();
  });
  // What the demo keeps a pane (app side): the shell says which pane is active, this says what it holds.
  const panes = new PaneMap<DemoPane>(shell.activePane);
  shell.registerViews([...demoViews(), ...stageViews(shell, status, panes, toolbars, !!o.activeOnly)]);
  shell.registerSections([
    {
      id: 'image',
      slot: 'export',
      title: 'Image',
      order: 1,
      render: domView(() =>
        form({ width: PD.Numeric(1920, { min: 16, max: 8192, step: 1 }), transparent: PD.Boolean(false) }),
      ),
    },
    {
      id: 'data',
      slot: 'export',
      title: 'Data',
      order: 2,
      render: domView(() =>
        form({
          format: PD.Select('json', [
            ['json', 'JSON'],
            ['csv', 'CSV'],
          ]),
        }),
      ),
    },
  ]);
  // The theme: a choice, from the rail item's menu (on a press); auto follows the system.
  const theme = new State<'auto' | 'light' | 'dark' | 'blue' | 'contrast'>('auto');
  const themes = [
    ['light', 'Light', 'sun'],
    ['dark', 'Dark', 'moon'],
    ['blue', 'Blue', 'palette'],
    ['contrast', 'High contrast', 'eye'],
    ['auto', 'Auto (system)', 'sun-moon'],
  ] as const;
  shell.registerRailItems([
    {
      id: 'theme',
      icon: 'sun-moon',
      title: 'Theme',
      // In the rail's bottom section (no group: the foot's tools).
      order: 1,
      menu: () =>
        themes.map(([id, label, icon]) => ({
          id,
          label,
          icon,
          checked: theme.value === id,
          run: () => theme.set(id),
        })),
    },
    // The command menu: in the corner, under the rail's Hide button (shell §3.2).
    {
      id: 'commands',
      icon: 'command',
      title: 'Commands',
      caption: '⌘K',
      corner: true,
      run: () => shell.show('confirm'),
    },
  ]);
  o.setup?.(shell);
  // A shared toolbar acts on the active pane: the app's entry for it (`PaneMap`), the shell saying only which pane.
  const activeTarget = (): CameraTarget | undefined => {
    const a = panes.active.value;
    return a && { stage: a.value.stage, pane: a.pane, title: a.value.title };
  };
  // The demo has no loader: a drop fills the stage and names the files in the status.
  const onDropFiles = (files: File[]) => {
    stage.setEmpty(false);
    status.post(
      `Opened ${files.length === 1 ? files[0].name : `${files.length} files: ${files.map((f) => f.name).join(', ')}`}`,
      'success',
    );
  };
  return (
    <div
      data-msui-theme={derive(theme, (t) => (t === 'auto' ? undefined : t))}
      style={`width: ${o.width ?? '100%'}; height: ${o.height ?? '100vh'}; display: flex; background: var(--msui-bg); color: var(--msui-text);`}
    >
      <ShellView
        shell={shell}
        main={({ top, pane }) => {
          // The app's content as a pane: its title to the layout, its stage to the app's map.
          onCleanup(pane.setTitle('Stage A'));
          onCleanup(panes.set('content', { stage, title: 'Stage A' }));
          return (
            <DemoStage
              shell={shell}
              stage={stage}
              status={status}
              name="Stage A"
              toolbars={toolbars}
              port={pane}
              rest={toolbars === 'pane'}
              top={top}
              activeOnly={o.activeOnly}
            />
          );
        }}
        over={
          toolbars === 'pane'
            ? undefined
            : {
                status: () => (
                  <>
                    <StatusLine status={status} scope={{ rest: true }} />
                    <StatusChips status={status} />
                  </>
                ),
                corner: () => (
                  <StageToolbar
                    groups={[...(toolbars === 'shared' ? [cameraGroup(activeTarget, status)] : []), viewGroup(shell)]}
                    stage={stage}
                  />
                ),
              }
        }
        busy={o.busy}
        onDropFiles={onDropFiles}
        dropLabel={o.dropLabel}
        dropFormats={o.dropFormats}
        fill={o.fill}
        railToggle={o.railToggle}
      />
    </div>
  );
}

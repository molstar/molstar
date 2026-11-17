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
import { type Child, domView, For } from '../../dom/index.js';
import { ParamsFormModel } from '../../forms/params-form.js';
import { ShellModel, type ShellOptions, type ShellViewDef } from '../shell-model.js';
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
export function MainPlaceholder() {
  return (
    <div
      style="flex: 1; display: grid; place-items: center; color: var(--msui-text-subtle); background-color: var(--msui-surface-sunken); background-image: radial-gradient(circle, var(--msui-line-strong) 1px, transparent 1.5px); background-size: 16px 16px; background-position: 8px 8px;"
      aria-label="The app's main view"
    >
      <span>Main area</span>
    </div>
  );
}

export interface DemoOptions extends ShellOptions {
  /** Called with the shell once its views are registered (to put it in a state). */
  setup?(shell: ShellModel): void;
  busy?: ReadonlyState<ShellBusy | undefined>;
  width?: string;
  height?: string;
  status?: boolean;
  /** The Mol* star that shows a hidden rail (default on). */
  railToggle?: boolean;
}

/** A shell with the demo's views, sections and rail items, in a box of the given size. */
export function DemoShell(props: { options?: DemoOptions; children?: Child }) {
  const o = props.options ?? {};
  const shell = new ShellModel(o);
  shell.registerViews(demoViews());
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
  const theme = new State(false);
  shell.registerRailItems([
    {
      id: 'theme',
      icon: 'sun-moon',
      title: 'Theme',
      group: 'app',
      pressed: theme.asReadonly,
      run: () => theme.set(!theme.value),
    },
    { id: 'commands', icon: 'command', title: 'Commands', caption: '⌘K', run: () => shell.show('confirm') },
  ]);
  o.setup?.(shell);
  return (
    <div style={`width: ${o.width ?? '100%'}; height: ${o.height ?? '100vh'}; display: flex;`}>
      <ShellView
        shell={shell}
        main={() => <MainPlaceholder />}
        status={o.status ? () => <span class="msui-label">Ready</span> : undefined}
        busy={o.busy}
        onDropFiles={() => {}}
        railToggle={o.railToggle}
      />
    </div>
  );
}

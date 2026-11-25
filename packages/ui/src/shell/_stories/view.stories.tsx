/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, State } from '@molstar/core/reactive';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { Heading, Note, page, story } from '../../_stories/support.js';
import { ActionFooter } from '../../components/action-footer/action-footer.js';
import { Button } from '../../components/button/button.js';
import { Expandable, ExpandableList } from '../../components/expandable/expandable.js';
import { Icon } from '../../components/icon/icon.js';
import { ParamsForm } from '../../components/params/params-form.js';
import { SkeletonRows } from '../../components/skeleton/skeleton.js';
import { TextInput } from '../../components/text-input/text-input.js';
import { installTooltips } from '../../components/tooltip/tooltip.js';
import { type Child, domView, For, onCleanup } from '../../dom/index.js';
import { ParamsFormModel } from '../../forms/params-form.js';
import { ShellModel, type ShellViewDef, type ViewPlacement } from '../shell-model.js';
import {
  ViewEmpty,
  ViewError,
  ViewFooter,
  ViewFrame,
  ViewHeading,
  ViewRow,
  ViewStack,
  ViewToolbar,
} from '../view-frame.js';
import { viewMenu } from '../view-menu.js';

/**
 * The view frame on its own (shell foundations, styles §8): the header and the body every view is shown in, at a
 * docked group's width and as a window. The views are written as an extension would write them: only the view
 * primitives (`ViewStack`, `ViewToolbar`, `ViewEmpty`, `ViewError`, `ViewFooter`, `ViewHeading`), the components and
 * the row classes, with no padding, gaps or type of their own. The boxes around them are the story's.
 */
const meta: Meta = { title: 'Shell/View' };
export default meta;

// ---------------------------------------------------------------------------------------------------------------------
// The views (content only: what an extension writes)

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

const exportImage = {
  format: PD.Select('png', [
    ['png', 'PNG'],
    ['jpeg', 'JPEG'],
    ['webp', 'WebP'],
  ]),
  width: PD.Numeric(1920, { min: 16, max: 8192, step: 1 }),
  height: PD.Numeric(1080, { min: 16, max: 8192, step: 1 }),
  transparent: PD.Boolean(false),
  axes: PD.Boolean(false),
  scale: PD.Numeric(1, { min: 0.25, max: 4, step: 0.25 }),
  quality: PD.Numeric(0.9, { min: 0, max: 1, step: 0.01 }),
  label: PD.Text('figure-1'),
};

const items = ['Background', 'Grid', 'Annotations', 'Selection', 'Measurements', 'Labels', 'Overlay'];

function Paragraph() {
  return (
    <ViewStack>
      <ViewHeading>About</ViewHeading>
      <ViewEmpty>
        A view with one line of content: the header above it, the body padded a unit at top and bottom.
      </ViewEmpty>
    </ViewStack>
  );
}

function SectionsView() {
  return (
    <ExpandableList>
      <Expandable title="Appearance" open>
        {() => form(appearance)}
      </Expandable>
      <Expandable title="Lighting" open>
        {() => form(lighting)}
      </Expandable>
      <Expandable title="Advanced">{() => form({ fog: PD.Boolean(false), depth: PD.Numeric(3) })}</Expandable>
    </ExpandableList>
  );
}

function FormView() {
  return <ViewStack>{form(appearance)}</ViewStack>;
}

function ListView() {
  const filter = new State('');
  const shown = derive(filter, (f) => items.filter((n) => n.toLowerCase().includes(f.toLowerCase())));
  return (
    <ViewStack>
      <ViewToolbar label="Layers">
        <TextInput value={filter} onInput={(v) => filter.set(v)} placeholder="Filter" aria-label="Filter" />
        <Button variant="ghost" size="sm" icon={<Icon name="plus" />} title="Add a layer" />
      </ViewToolbar>
      <div class="msui-list" role="list">
        <For each={shown} by={(n) => n}>
          {(name) => (
            <div class="msui-row" role="listitem">
              <span class="msui-slot">
                <Icon name="layers" />
              </span>
              <span>{name}</span>
              <Button variant="ghost" size="xs" icon={<Icon name="eye" />} title={`Hide ${name}`} />
            </div>
          )}
        </For>
      </div>
    </ViewStack>
  );
}

function EmptyView() {
  return (
    <ViewStack>
      <ViewEmpty action={<Button variant="primary">Open a file…</Button>}>
        No layers yet. Open a file, or drop one anywhere, and its layers show here.
      </ViewEmpty>
    </ViewStack>
  );
}

function LoadingView() {
  return (
    <ViewStack>
      <ViewHeading>Layers</ViewHeading>
      <SkeletonRows count={6} icon />
    </ViewStack>
  );
}

function ErrorView() {
  return (
    <ViewStack>
      <ViewError action={<Button size="xs">Retry</Button>}>Could not load the layers.</ViewError>
      <ViewHeading>Kept from before</ViewHeading>
      <div class="msui-list">
        {items.slice(0, 3).map((name) => (
          <div class="msui-row">
            <span class="msui-slot">
              <Icon name="layers" />
            </span>
            <span>{name}</span>
          </div>
        ))}
      </div>
    </ViewStack>
  );
}

/** A form longer than the body, its Download kept at the bottom while it scrolls: the footer after the content, not in its stack. */
function FooterView() {
  return [
    <ViewStack>
      <ExpandableList>
        <Expandable title="Image" open>
          {() => form(exportImage)}
        </Expandable>
        <Expandable title="Lighting" open>
          {() => form(lighting)}
        </Expandable>
      </ExpandableList>
    </ViewStack>,
    <ViewFooter>
      <ActionFooter label="Download" onApply={() => {}} onReset={() => {}} info="1920 × 1080, about 2 MB" />
    </ViewFooter>,
  ];
}

/** Everything a view is composed of, in its order: a toolbar, a list, a section with a form, a footer. */
function ComposedView() {
  return [
    <ViewStack>
      <ViewToolbar label="Measurements">
        <TextInput value="" placeholder="Search" aria-label="Search" />
        <Button variant="ghost" size="sm" icon={<Icon name="plus" />} title="Add" />
      </ViewToolbar>
      <div class="msui-list">
        {[
          ['Distance', '12.40 Å'],
          ['Angle', '104.5°'],
        ].map(([name, value]) => (
          <div class="msui-row">
            <span class="msui-slot">
              <Icon name="ruler" />
            </span>
            <span>{name}</span>
            <span class="msui-num">{value}</span>
          </div>
        ))}
      </div>
      <ExpandableList>
        <Expandable title="Options" open>
          {() => form(appearance)}
        </Expandable>
      </ExpandableList>
    </ViewStack>,
    <ViewFooter>
      <ActionFooter label="Apply" onApply={() => {}} onReset={() => {}} />
    </ViewFooter>,
  ];
}

// ---------------------------------------------------------------------------------------------------------------------
// The story's stage: a frame in a box of a docked group's width, or of a window's

const headerActions = domView(() => (
  <>
    <Button variant="ghost" size="sm" icon={<Icon name="refresh-cw" />} title="Refresh" />
    <Button variant="ghost" size="sm" icon={<Icon name="plus" />} title="Add" />
  </>
));

interface FrameOptions {
  title?: string;
  heading?: string;
  actions?: boolean;
  menu?: boolean;
  close?: boolean;
  placement?: ViewPlacement;
  width?: number;
  height?: number;
  fill?: boolean;
}

let next = 0;

/** A view in its frame, in a box: a docked group's surface, or a window's raised one with its edge and shadow. */
function Framed(props: { view: () => Child; options?: FrameOptions }) {
  const o = props.options ?? {};
  const placement = o.placement ?? 'docked';
  const shell = new ShellModel();
  const id = `story-view-${next++}`;
  const def: ShellViewDef = {
    id,
    title: o.title ?? 'Layers',
    heading: o.heading,
    icon: 'layers',
    fill: o.fill,
    actions: o.actions ? headerActions : undefined,
    render: domView(props.view),
  };
  shell.registerViews([def]);
  const width = o.width ?? (placement === 'window' ? 360 : 340);
  const height = o.height ?? (placement === 'window' ? 320 : 420);
  const box =
    placement === 'window'
      ? 'background: var(--msui-surface); border-radius: var(--msui-radius-lg); box-shadow: var(--msui-edge), var(--msui-shadow-float); overflow: hidden;'
      : 'background: var(--msui-surface); box-shadow: -1px 0 0 var(--msui-line);';
  return (
    <div
      ref={(el: HTMLElement) => onCleanup(installTooltips(el))}
      data-tip-edge=""
      style={`display: flex; flex-direction: column; width: ${width}px; height: ${height}px; ${box}`}
    >
      <ViewFrame
        shell={shell}
        ref={id}
        def={def}
        placement={placement}
        menu={o.menu ? () => viewMenu(shell, id) : undefined}
        onClose={o.close ? () => {} : undefined}
      />
    </div>
  );
}

/** A caption over a framed view. */
function Labelled(props: { label: string; children: Child }) {
  return (
    <div style="display: grid; gap: var(--msui-space-2); align-content: start;">
      <Heading>{props.label}</Heading>
      {props.children}
    </div>
  );
}

const row = 'display: flex; flex-wrap: wrap; gap: var(--msui-space-5); align-items: start;';

const stage = (content: () => Child): StoryObj => ({
  render: () => story(() => <div style={page.replace('max-width: 960px;', '')}>{content()}</div>),
});

// ---------------------------------------------------------------------------------------------------------------------
// Stories

/**
 * The header's variants: the title alone; the view's own actions after it (ghost `sm` icon buttons); the layout menu
 * (⋯: Float, Hide, Reset layout); ✕ (a card, a window); a long title that fades; a heading that differs from the title
 * (the rail's tooltip and a tab say the title, the header the heading).
 */
export const Headers = stage(() => (
  <div style={row}>
    <Labelled label="Title only">
      <Framed view={Paragraph} options={{ height: 120 }} />
    </Labelled>
    <Labelled label="With actions">
      <Framed view={Paragraph} options={{ height: 120, actions: true }} />
    </Labelled>
    <Labelled label="With the menu">
      <Framed view={Paragraph} options={{ height: 120, menu: true }} />
    </Labelled>
    <Labelled label="With ✕">
      <Framed view={Paragraph} options={{ height: 120, close: true }} />
    </Labelled>
    <Labelled label="Everything">
      <Framed view={Paragraph} options={{ height: 120, actions: true, menu: true, close: true }} />
    </Labelled>
    <Labelled label="A long title fades">
      <Framed
        view={Paragraph}
        options={{ height: 120, title: 'Measurements, labels and annotations of the selection', menu: true }}
      />
    </Labelled>
    <Labelled label="A heading other than the title">
      <Framed
        view={Paragraph}
        options={{ height: 120, title: 'PyMOL / VMD console', heading: 'Console', menu: true }}
      />
    </Labelled>
  </div>
));

/** Sections: an `ExpandableList` of forms, docked and in a window. */
export const Sections = stage(() => (
  <div style={row}>
    <Labelled label="Docked">
      <Framed view={SectionsView} options={{ title: 'Inspect', menu: true }} />
    </Labelled>
    <Labelled label="Window">
      <Framed view={SectionsView} options={{ title: 'Inspect', placement: 'window', menu: true, close: true }} />
    </Labelled>
  </div>
));

/** A form: labels on the inset, values from the label column. */
export const Form = stage(() => (
  <div style={row}>
    <Labelled label="Docked">
      <Framed view={FormView} options={{ title: 'Appearance', height: 200 }} />
    </Labelled>
    <Labelled label="Window">
      <Framed view={FormView} options={{ title: 'Appearance', placement: 'window', height: 200, close: true }} />
    </Labelled>
  </div>
));

/** A list with a toolbar above it (a filter and Add): the toolbar and the rows' text on the inset. */
export const ListWithToolbar = stage(() => (
  <div style={row}>
    <Labelled label="Docked">
      <Framed view={ListView} options={{ actions: true, menu: true }} />
    </Labelled>
    <Labelled label="Window">
      <Framed view={ListView} options={{ placement: 'window', menu: true, close: true }} />
    </Labelled>
  </div>
));

const treeRows: { name: string; description: string; level: number; icon: 'layers' | 'tag' | 'ruler' }[] = [
  { name: 'Annotations of the selected region', description: '128 items, 4 hidden', level: 0, icon: 'tag' },
  { name: 'Labels', description: 'Names and identifiers, shown on hover', level: 1, icon: 'tag' },
  { name: 'Measurements', description: 'Distances and angles', level: 1, icon: 'ruler' },
  { name: 'Overlay with a very long name that does not fit', description: 'Drawn last', level: 0, icon: 'layers' },
  { name: 'Grid', description: 'Every 10 units', level: 0, icon: 'layers' },
];

function RowsView() {
  return (
    <ViewStack>
      <div class="msui-list" role="tree" aria-label="Layers">
        {treeRows.map((r) => (
          <ViewRow
            role="treeitem"
            tabIndex={0}
            aria-level={String(r.level + 1)}
            level={r.level}
            icon={r.icon}
            label={r.name}
            description={r.description}
            trailing={<Button variant="ghost" size="xs" icon={<Icon name="eye" />} title={`Hide ${r.name}`} />}
          />
        ))}
      </div>
    </ViewStack>
  );
}

/**
 * Rows (`ViewRow`): a slot, the label, a muted description, what trails; a level indents a slot (a tree). Cut off, a
 * row's name and description fade, and a card beside the view shows them whole: hover a row's text, or Tab into the
 * list and move with the arrows. The card goes outside the view, on the side with room (in the shell: toward the main
 * area).
 */
export const Rows = stage(() => (
  // One above the other: each view's cards go to its open side (beside each other, one's cards would cover the other).
  <div style="display: grid; gap: var(--msui-space-5); justify-items: start;">
    <Labelled label="Wide (420)">
      <Framed view={RowsView} options={{ width: 420, height: 240 }} />
    </Labelled>
    <Labelled label="Narrow (280)">
      <Framed view={RowsView} options={{ width: 280, height: 240 }} />
    </Labelled>
  </div>
));

/** Empty, loading and error: what goes here and how (one line, an action), a skeleton in place, an error line. */
export const EmptyLoadingError = stage(() => (
  <div style={row}>
    <Labelled label="Empty">
      <Framed view={EmptyView} options={{ height: 240 }} />
    </Labelled>
    <Labelled label="Loading">
      <Framed view={LoadingView} options={{ height: 240 }} />
    </Labelled>
    <Labelled label="Error">
      <Framed view={ErrorView} options={{ height: 240 }} />
    </Labelled>
  </div>
));

/**
 * A footer of actions under the body's scroll, at the view's bottom, docked and in a window: wherever the view writes
 * its `ViewFooter` (after its `ViewStack` or inside it), it moves into the body's footer slot.
 */
export const Footer = stage(() => (
  <div style={row}>
    <Labelled label="Docked">
      <Framed view={FooterView} options={{ title: 'Export', height: 360, menu: true }} />
    </Labelled>
    <Labelled label="Window">
      <Framed view={FooterView} options={{ title: 'Export', placement: 'window', menu: true, close: true }} />
    </Labelled>
    <Labelled label="Short content: the footer at the bottom">
      <Framed
        view={() => [
          <ViewStack>{form({ animations: PD.Boolean(true) })}</ViewStack>,
          <ViewFooter>
            <ActionFooter label="Apply" onApply={() => {}} />
          </ViewFooter>,
        ]}
        options={{ title: 'Settings', height: 240 }}
      />
    </Labelled>
    <Labelled label="The footer written inside a ViewStack: still at the bottom">
      <Framed
        view={() => (
          <ViewStack>
            {form({ animations: PD.Boolean(true) })}
            <ViewFooter>
              <ActionFooter label="Apply" onApply={() => {}} />
            </ViewFooter>
          </ViewStack>
        )}
        options={{ title: 'Settings', height: 240 }}
      />
    </Labelled>
  </div>
));

/** The same view at its narrowest (280) and wide (560): labels never cut, names fade. */
export const NarrowAndWide = stage(() => (
  <div style={row}>
    <Labelled label="Narrow (280)">
      <Framed view={ComposedView} options={{ title: 'Measurements', width: 280, menu: true, actions: true }} />
    </Labelled>
    <Labelled label="Docked (340)">
      <Framed view={ComposedView} options={{ title: 'Measurements', menu: true, actions: true }} />
    </Labelled>
    <Labelled label="Wide (560)">
      <Framed view={ComposedView} options={{ title: 'Measurements', width: 560, menu: true, actions: true }} />
    </Labelled>
  </div>
));

/** The guides: the inset text column and the trailing column, from the frame's edges. */
const Guides: readonly [label: string, side: 'left' | 'right', at: string, color: string][] = [
  ['text column (inset)', 'left', 'var(--msui-inset)', 'var(--msui-accent)'],
  ['text after an icon slot', 'left', 'calc(var(--msui-inset) + var(--msui-indent))', 'var(--msui-info)'],
  ['trailing column (inset from the right)', 'right', 'var(--msui-inset)', 'var(--msui-danger)'],
];

function GuideLines(props: { show: State<boolean> }) {
  return (
    <div style="position: absolute; inset: 0; pointer-events: none;">
      {Guides.map(([label, side, at, color]) => (
        <div
          title={label}
          style={{
            position: 'absolute',
            top: '0',
            bottom: '0',
            [side]: at,
            width: '1px',
            background: color,
            opacity: '0.8',
            display: derive(props.show, (s) => (s ? 'block' : 'none')),
          }}
        />
      ))}
    </div>
  );
}

/**
 * Columns (shell foundations, as `Foundations/Alignment` does for panels): the inset text column and the trailing
 * column drawn over views. The header's title, a section's header, a form's labels, the toolbar and the rows' text
 * start on the inset; the header's actions, a section's chevron and a row's trailing button end on the trailing
 * column.
 */
export const Columns = stage(() => {
  const show = new State(true);
  const guided = (content: Child) => (
    <div style="position: relative;">
      {content}
      <GuideLines show={show} />
    </div>
  );
  return (
    <>
      <label style="display: flex; gap: var(--msui-space-2); align-items: center;">
        <input type="checkbox" checked={show} onChange={(e) => show.set(e.currentTarget.checked)} />
        Column guides
      </label>
      <div style={row}>
        <Labelled label="Docked">
          {guided(<Framed view={ComposedView} options={{ title: 'Measurements', menu: true, actions: true }} />)}
        </Labelled>
        <Labelled label="Window">
          {guided(
            <Framed
              view={SectionsView}
              options={{ title: 'Inspect', placement: 'window', height: 420, menu: true, close: true }}
            />,
          )}
        </Labelled>
        <dl style="display: grid; gap: var(--msui-space-2); margin: 0; align-content: start;">
          {Guides.map(([label, , , color]) => (
            <div style="display: flex; gap: var(--msui-space-2); align-items: center;">
              <span style={{ width: '12px', height: '2px', background: color }} />
              <span>{label}</span>
            </div>
          ))}
        </dl>
      </div>
      <Note>
        A view never sets its own padding: the body pads the surface's padding, and each primitive pads the rest of the
        inset.
      </Note>
    </>
  );
});

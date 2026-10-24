/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { derive } from '@molstar/prism-ui/reactive/state';
import type { SequenceChainsEntity, SequenceChainsPreset, SequenceInspectorModel } from '../../model/scene.js';
import { classes as quick } from './_generated/representation-quick.module.css.js';
import { classes } from './_generated/sequence-inspector.module.css.js';

// A structure's sequence in its Sequence node's inspector (sequence.md §8.1), drawn as a representation's quick
// settings are (representation-quick.tsx): what the band shows and Add track (the command menu at its group), then each
// setting a row, its name and value over its choices: the style as two icons in one strip, the letter colors as tiles
// painted with their colors, the chains as All / None / Polymers in a strip and a line per entity (a box for its chains,
// checked, empty or mixed; a dot for its kind; its name, fading where cut; how many are shown), an entity of several
// chains opening to its chains as toggles that wrap. Hundreds of entities: the first eight, then Show all; an entity's
// chains only once it is opened.

const PageSize = 8;

const Kinds: Readonly<Record<SequenceChainsEntity['kind'], { label: string; class: string }>> = {
  polymer: { label: 'Polymer', class: classes.polymer },
  ligand: { label: 'Ligand', class: classes.ligand },
  branched: { label: 'Carbohydrate', class: classes.branched },
};

const Styles = [
  { id: 'tracks', label: 'Tracks', icon: 'chart-gantt' },
  { id: 'text', label: 'Text', icon: 'align-left' },
] as const;

const Presets: readonly { id: SequenceChainsPreset; label: string; tip: string }[] = [
  { id: 'all', label: 'All', tip: 'Every chain' },
  { id: 'none', label: 'None', tip: 'No chain, then show them one by one' },
  { id: 'polymer', label: 'Polymers', tip: 'The polymer chains, no ligands' },
];

export function SequenceInspector(props: { model: SequenceInspectorModel }) {
  const { model } = props;
  const chainsValue = derive([model.shown, model.entities], (s, es) => {
    const total = es.reduce((n, e) => n + e.chains.length, 0);
    return s.size === total ? 'all' : `${s.size} of ${total}`;
  });
  return (
    <section class={`${quick.quick} ${classes.root}`} aria-label="Sequence">
      <div class={classes.status}>
        <span class={classes.statusText}>{model.status}</span>
        <Show when={model.canAdd}>
          {() => (
            <Button size="xs" variant="ghost" icon={<Icon name="plus" />} onClick={() => model.addTrack()}>
              Add track
            </Button>
          )}
        </Show>
      </div>
      {/* The band's style and letters side by side: two short settings, one row. */}
      <div class={classes.pair}>
        <div class={quick.setting}>
          <div class={quick.head}>
            <span class={quick.name}>Style</span>
            <span class={quick.value}>{derive(model.style, (s) => (s === 'text' ? 'Text' : 'Tracks'))}</span>
          </div>
          <div class={`${quick.segments} ${classes.styles}`} role="radiogroup" aria-label="Style">
            {Styles.map((s) => (
              <button
                type="button"
                role="radio"
                class={`${quick.segment} ${classes.segment}`}
                aria-checked={derive(model.style, (v) => v === s.id)}
                aria-label={s.label}
                data-tip={s.label}
                onClick={() => model.setStyle(s.id)}
              >
                <Icon name={s.icon} size={16} />
              </button>
            ))}
          </div>
        </div>
        <div class={quick.setting}>
          <div class={quick.head}>
            <span class={quick.name}>Letters</span>
            <span class={quick.value}>
              {derive([model.colors, model.schemes], (c, list) => list.find((x) => x.id === c)?.label ?? '')}
            </span>
          </div>
          <div class={`${quick.tiles} ${classes.tiles}`} role="radiogroup" aria-label="Letter colors">
            <For each={model.schemes} by={(x) => x.id ?? 'none'}>
              {(x) => (
                <button
                  type="button"
                  role="radio"
                  class={quick.tile}
                  classList={{ [quick.unpainted]: !x.swatch }}
                  style={{ background: x.swatch ?? '' }}
                  aria-checked={derive(model.colors, (c) => c === x.id)}
                  aria-label={x.label}
                  data-tip={x.label}
                  onClick={() => model.setColors(x.id)}
                />
              )}
            </For>
          </div>
        </div>
      </div>
      <div class={quick.setting}>
        <div class={classes.chainsHead}>
          <div class={quick.head}>
            <span class={quick.name}>Chains</span>
            <span class={quick.value}>{chainsValue}</span>
          </div>
          <div class={`${quick.segments} ${classes.presets}`} role="radiogroup" aria-label="Chains">
            {Presets.map((p) => (
              <button
                type="button"
                role="radio"
                class={`${quick.segment} ${classes.segment} ${classes.word}`}
                aria-checked={derive(model.preset, (v) => v === p.id)}
                data-tip={p.tip}
                onClick={() => model.apply(p.id)}
              >
                {p.label}
              </button>
            ))}
          </div>
        </div>
        <Entities model={model} />
        <button
          type="button"
          role="checkbox"
          class={classes.line}
          aria-checked={derive(model.each, (e) => (e ? 'true' : 'false'))}
          data-tip="Chains of one entity each a block of their own; merged, one block says ×N"
          onClick={() => model.setEach(!model.each.value)}
        >
          <Box checked={derive(model.each, (e) => (e ? 'true' : 'false'))} />
          <span class={classes.lineText}>Every copy apart</span>
        </button>
      </div>
    </section>
  );
}

/** A checkbox's box: checked, empty or mixed. */
function Box(props: { checked: ReadonlyState<'true' | 'false' | 'mixed'> }) {
  return (
    <span class={classes.box} data-checked={props.checked} aria-hidden="true">
      <Show when={derive(props.checked, (c) => c === 'true')}>{() => <Icon name="check" size={12} />}</Show>
      <Show when={derive(props.checked, (c) => c === 'mixed')}>{() => <Icon name="minus" size={12} />}</Show>
    </span>
  );
}

function Entities(props: { model: SequenceInspectorModel }) {
  const { model } = props;
  const all = observableState(false);
  const listed = derive([all.asReadonly, model.entities], (a, es) =>
    a || es.length <= PageSize + 1 ? es : es.slice(0, PageSize),
  );
  return (
    <div class={classes.entities}>
      <div class={classes.list} role="list" aria-label="Entities">
        <For each={listed} by={(e) => e.id}>
          {(e) => <EntityRow model={model} entity={e} />}
        </For>
      </div>
      <Show when={derive([all.asReadonly, model.entities], (a, es) => !a && es.length > PageSize + 1)}>
        {() => (
          <button type="button" class={classes.more} onClick={() => all.set(true)}>
            {derive(model.entities, (es) => `Show all ${es.length}`)}
            <span class={classes.moreCount}>{derive(model.entities, (es) => `${es.length - PageSize} more`)}</span>
          </button>
        )}
      </Show>
    </div>
  );
}

function EntityRow(props: { model: SequenceInspectorModel; entity: SequenceChainsEntity }) {
  const { model, entity } = props;
  const kind = Kinds[entity.kind];
  const several = entity.chains.length > 1;
  const open = observableState(false);
  const count = derive(model.shown, (s) => entity.chains.filter((c) => s.has(c.name)).length);
  const checked = derive(count, (n): 'true' | 'false' | 'mixed' =>
    n === 0 ? 'false' : n === entity.chains.length ? 'true' : 'mixed',
  );
  const meta = derive(count, (n) =>
    several
      ? `${kind.label} · ${n} of ${entity.chains.length} chains`
      : `${kind.label} · chain ${entity.chains[0].label}`,
  );
  return (
    <div role="listitem" class={classes.entity}>
      <div class={classes.row}>
        <button
          type="button"
          role="checkbox"
          class={classes.boxButton}
          aria-checked={checked}
          aria-label={`Show ${entity.name}`}
          onClick={() => model.toggleEntity(entity.id)}
        >
          <Box checked={checked} />
        </button>
        <span class={`${classes.dot} ${kind.class}`} data-tip={kind.label} aria-label={kind.label} />
        <button
          type="button"
          class={classes.text}
          aria-expanded={several ? open.asReadonly : undefined}
          onClick={() => (several ? open.set(!open.value) : model.toggleEntity(entity.id))}
        >
          <span class={classes.name} data-tip-overflow="">
            {entity.name}
          </span>
          <span class={classes.meta}>{meta}</span>
        </button>
        {several ? (
          <span class={classes.chevron} classList={{ [classes.chevronOpen]: open.asReadonly }}>
            <Icon name="chevron-right" size={14} />
          </span>
        ) : null}
      </div>
      <Show when={open.asReadonly}>
        {() => (
          <div class={classes.chips} role="group" aria-label={`Chains of ${entity.name}`}>
            {entity.chains.map((c) => (
              <button
                type="button"
                class={classes.chip}
                aria-pressed={derive(model.shown, (s) => s.has(c.name))}
                data-tip={c.label}
                onClick={() => model.toggleChain(c.name)}
              >
                {c.name}
              </button>
            ))}
          </div>
        )}
      </Show>
    </div>
  );
}

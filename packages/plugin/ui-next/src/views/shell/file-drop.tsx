/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map } from 'rxjs';
import { Icon } from '../../components/icon/icon.js';
import { onCleanup } from '../../dom/scope.js';
import { observableState } from '../../reactive/observable-state.js';
import { classes } from './_generated/file-drop.module.css.js';

/**
 * Files dragged anywhere over the plugin: the whole plugin becomes the drop target (Mol*'s drag overlay; the spec turns
 * it off, `components.disableDragOverlay`). Shown only for drags that carry files; dropped files go to `onDrop` (the
 * plugin's drag-and-drop manager, which picks the format).
 */
export function FileDropOverlay(props: {
  /** Where drags are watched for: the plugin's root. */
  target: HTMLElement;
  onDrop(files: readonly File[]): void;
  /** File formats it reads, as a hint ("mmCIF, PDB, CCP4, …"). */
  formats?: readonly string[];
}) {
  const shown = observableState(false);
  const onEnter = (e: DragEvent) => {
    if (hasFiles(e)) shown.set(true);
  };
  props.target.addEventListener('dragenter', onEnter);
  onCleanup(() => props.target.removeEventListener('dragenter', onEnter));
  // Over the overlay: copy, not the browser opening the file; leaving it or dropping hides it.
  const over = (e: DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    if (e.dataTransfer) e.dataTransfer.dropEffect = 'copy';
  };
  return (
    <div
      class={classes.overlay}
      hidden={shown.asReadonly.pipe(map((s) => !s))}
      onDragEnter={over}
      onDragOver={over}
      onDragLeave={() => shown.set(false)}
      onDrop={(e) => {
        over(e);
        shown.set(false);
        const files = [...(e.dataTransfer?.files ?? [])];
        if (files.length) props.onDrop(files);
      }}
    >
      <div class={classes.target}>
        <Icon name="folder-open" size={24} />
        <span class={classes.title}>Drop to open</span>
        {props.formats?.length ? <span class={classes.hint}>{props.formats.join(', ')}</span> : null}
      </div>
    </div>
  );
}

function hasFiles(e: DragEvent): boolean {
  const types = e.dataTransfer?.types;
  if (!types) return false;
  for (let i = 0; i < types.length; i++) if (types[i] === 'Files') return true;
  return false;
}

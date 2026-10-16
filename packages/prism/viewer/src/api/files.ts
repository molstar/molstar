/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Asset } from '@molstar/core/util/assets';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

// File params as the API takes them: DOM `File`s (what a drop, a picker or a script has). PD's own file params are
// typed as Mol*'s assets, so an operation holds plain values; its command shows a file picker (`filePickers`), whose
// assets `filesOf` reads too.

type Info = PD.Info & { accept?: string };

/** Several files. */
export function FilesParam(info?: Info): PD.Value<readonly File[]> {
  return PD.Value<readonly File[]>([], { ...info, isHidden: true });
}

/** One file. */
export function FileParam(info?: Info): PD.Value<File | null> {
  return PD.Value<File | null>(null, { ...info, isHidden: true });
}

/** A command's form for an operation's file params: pickers in their place. */
export function filePickers(params: PD.Params): PD.Params {
  const out: PD.Params = {};
  for (const [k, p] of Object.entries(params)) {
    out[k] =
      p.type === 'value' && Array.isArray(p.defaultValue)
        ? PD.FileList({ label: p.label, description: p.description, multiple: true })
        : p.type === 'value' && p.defaultValue === null
          ? PD.File({ label: p.label, description: p.description })
          : p;
  }
  return out;
}

/** The files a value holds, whether DOM files or Mol*'s file assets. */
export function filesOf(value: unknown): File[] {
  if (!value) return [];
  const list = Array.isArray(value)
    ? value
    : typeof (value as ArrayLike<unknown>).length === 'number'
      ? Array.from(value as ArrayLike<unknown>)
      : [value];
  const out: File[] = [];
  for (const v of list) {
    if (typeof File !== 'undefined' && v instanceof File) out.push(v);
    else if (Asset.isFile(v as Asset) && (v as Asset.File).file) out.push((v as Asset.File).file as File);
  }
  return out;
}

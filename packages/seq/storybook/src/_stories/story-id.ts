/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A story's id as Storybook makes it (the sanitized title, `--`, the sanitized name made of the export's), so a page or
// an index links a story by its title and export, and a test checks the link against the stories that exist.

/** The prefix of every title the repository's Storybook gives this package's stories (`storybook/packages.ts`). */
export const TitlePrefix = 'Seq* (next)';

/** Storybook's `sanitize`: lower case, every run of punctuation and spaces one `-`, none at the ends. */
export function sanitize(s: string): string {
  return s
    .toLowerCase()
    .replace(/[ ’–—―′¿'`~!@#$%^&*()_|+\-=?;:'",.<>{}[\]\\/]/g, '-')
    .replace(/-+/g, '-')
    .replace(/^-+/, '')
    .replace(/-+$/, '');
}

/** An export's story name (`ProteinOnItsGene` is `Protein On Its Gene`, `Chr17` is `Chr 17`: Storybook's `startCase`). */
export function storyNameFromExport(key: string): string {
  return key
    .replace(/([a-z])([A-Z])/g, '$1 $2')
    .replace(/([A-Z]+)([A-Z][a-z])/g, '$1 $2')
    .replace(/([a-zA-Z])(\d)/g, '$1 $2')
    .replace(/(\d)([a-zA-Z])/g, '$1 $2')
    .replace(/[_-]+/g, ' ')
    .trim();
}

/** The id of the export `exportName` of the stories titled `title` (without the package's prefix). */
export function storyId(title: string, exportName: string): string {
  return `${sanitize(`${TitlePrefix}/${title}`)}--${sanitize(storyNameFromExport(exportName))}`;
}

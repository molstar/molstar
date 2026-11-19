# Examples and stories: where they live and how they run

## Intent

The second review found that the quick lane must render every view, but the views would live in the storybook package,
which core's tests cannot import; that nothing said how fixtures are read in Node or in Storybook, or how a story is
structured. This design settles it.

## Decision

- **An examples package**, `packages/seq/examples` (`@molstar/seqx-examples`, private): the datasets' fixtures and one
  module per view, each exporting `build(): Spec` written with the builder. It depends on `@molstar/seqx-core` (and `@molstar/core`) only, so
  the core's and the viewer's tests import it, and the storybook package imports it too; its record script, run by hand
  and never by a test, reads mmCIF through `@molstar/io`, a devDependency ([examples.md](examples.md)).
- **Fixtures** live in `examples/fixtures/<dataset>/`, small files committed (FASTA, `seq` files, a tiny bigWig);
  generated ones are made in code from a seed, never stored ([testing.md](testing.md)).
- **Reading fixtures**: a spec names fixtures by relative location (`fixtures/p53/region.fa`). In Node an assets adapter
  reads them from the package folder; in Storybook a Vite plugin serves the folder (as the previous storybook's did,
  registered in `storybook/packages.ts`). The same spec runs in both.
- **Features and examples are apart** (the author, 2026-11-19): a story that explains how something works (themes and
  their schemes, selection, editing and undo, bookmarks and links, stories, navigators) is `Seq* (next)/Features/<Feature>/<Story>`,
  on whichever dataset shows it best; a dataset's figure is an example. Why: a reader looking for how a feature works
  should not have to guess which dataset holds it, and a dataset's folder shows what its data looks like. A
  dataset's index story links the features shown on its data.
- **Examples load their data from files** (the author, 2026-11-19: "we should have storybook example that load data
  from files too instead of just inlining the values"): a dataset's fixtures are read as sources (FASTA, Stockholm,
  Newick, annotation tables with their transforms, the `seq` format), as a reader of Seq\* would; inline data only
  where the values are few or the inline form is what a story shows.
- **One story per view**: `Seq* (next)/Examples/<Dataset>/<View>`, each composing a viewer from the default spec and
  the view's `build()`, with the canvas and the SVG side by side and the printed spec beside them. A dataset's index
  story lists its views.
- **Tests per view** (the examples package's `_test`): build, validate, resolve and lay out with no unexpected issues,
  render to SVG; the few recorded views compare their draw-list record ([testing.md](testing.md)).
- **Exported files**: each view's JSON spec is written by a script for readers without Seq\* code (the third way a view
  runs, [examples.md](examples.md)); not committed.

## Consequences

- The package joins the inventory, the workspace and the root tsconfig references; its fixtures are excluded from the
  published packages (it is private).
- [packages.md](packages.md) lists it.

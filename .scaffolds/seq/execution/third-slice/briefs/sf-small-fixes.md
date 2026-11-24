# Lane SF: small fixes beside RK and CS

**Goal** (the author, 2026-11-19: "Small fixes"): two leftovers from `followups.md`, each with a test that fails on its
parent:

1. **The size script's `--files`** (HC's leftover): `node packages/seq/_measure/sizes.mjs --files` crashes after its
   first entry (`ReferenceError: entryOutput is not defined`, `sizes.mjs:108`). Make it print each entry's largest input
   files as its usage line says. A check that runs it with `--files` (the script's own, or a quick test if the package
   has a place for one; no network).
2. **`@molstar/ui`'s disabled toolbar item** (C0's followup): a disabled item's color (`--msui-text-subtle`) is barely
   distinct from an enabled one's (`--msui-text-muted`), so Seq\*'s `mountToolbar` dims it with
   `.seqx-toolbar [aria-disabled="true"] { opacity: 0.4; }` (`packages/seq/viewer/src/mount-toolbar.ts:29`). Make the
   disabled state read as disabled in `@molstar/ui` itself (find where its toolbar styles a disabled item; a token or a
   rule, in both themes, contrast measured and reported; WCAG exempts disabled controls, so aim for clearly distinct, not
   4.5:1), then remove Seq\*'s rule. Look at the toolbar in `@molstar/ui`'s Storybook and in a Seq\* story, light and
   dark. A change to a shared token used outside the toolbar: stop and report instead.

**Not in scope**: anything else. RK works in `packages/seq/core` and `packages/seq/viewer`; CS in
`packages/seq/examples` and `packages/seq/storybook`. Your only edit in their packages is the one rule in
`mount-toolbar.ts` (and its test if one asserts the rule).

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/sf-small-fixes
<that commit>`.

**Checks**: `pnpm test:unit packages/ui packages/seq/viewer`; `npx tsc -b packages/ui packages/seq/viewer`; `pnpm -C
packages/seq/viewer typecheck:test`; `pnpm lint`; `pnpm format:check`; the size script with and without `--files`.

**Headers**, **hot files** and **stop when**: as every lane's (see [qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; checks; the contrast before and after in both themes; what you saw; a `.v6/changelogs/` entry's
text for the `@molstar/ui` change.

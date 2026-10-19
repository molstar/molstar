# v6 changelog entries

Pending entries for the [v6 changelog draft](../plans/changelog.md), one file per change. Parallel branches would all
append to the end of `changelog.md` and conflict there. A file of its own per branch never conflicts, and the entries
get folded into the draft now and then.

## Adding an entry

- A change with user-facing effects adds **one** file here instead of editing `../plans/changelog.md`: `<slug>.md`,
  with a slug that is unique and says what the change is (the branch name without its prefix works:
  `perf-marching-cubes.md`, `state-visibility-flags.md`). A follow-up commit on the same branch edits that file; it
  does not add another.
- The content is what the draft's section for this change would say: one `## <Title>` heading, named the way the draft
  names its sections (an area, a feature or a batch: `## Graphics: faster marching cubes`), followed by the bullets.
  If the entries belong to an existing section of the draft, name that section exactly; the reconcile merges them in.
- Bullets follow the draft's style: what changed for a user or an API consumer, terse, `[Breaking]` first when it
  breaks a public API (which is also recorded in `../plans/breaking-v6-changes.md`), and the PR number `(#1234)` at the
  end once it is known. Performance entries give the measured gain and on what (`186k atoms: 12.9 s → 7.6 s`).
- Fixes that also apply to 5.x still go into the `[Unreleased]` section of the root `CHANGELOG.md` as well.

## Reconciling

From time to time, in a commit of its own on master (not inside a feature branch):

1. Fold each file's section into `../plans/changelog.md`: merge it into the section of the same name, or append it as
   a new section at the end in the order the changes landed (`git log --diff-filter=A -- .v6/changelogs/`).
2. Delete the folded files. This `AGENTS.md` stays.

Until then the draft plus the files here are the changelog: whoever writes the release notes or answers "what changed"
reads both.

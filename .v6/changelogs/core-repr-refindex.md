## Core: patches keep their ref index across revisions

- `@molstar/core/repr/patch`: applying a patch derives the new revision's ref index from the old one's, lazily, instead
  of rebuilding it, so a chain of edits on a spec with thousands of refs no longer pays a full index build per edit
  (at 10,000 refs a move's apply about 10 → 6 ms; a one-key set beside 30 story steps about 3 ms).

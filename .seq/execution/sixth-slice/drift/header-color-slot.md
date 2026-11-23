# Drift: `seq:color` is menu-only, not slot 2 (track-header-parts)

**What**: `track-header-parts.md`'s built-ins table gives `seq:color` "action, slot 2" with the `palette` icon. The
code ships it menu only ("Color ▸": one place, never also a slot), so a track's header has one slot (`seq:hide`) by
default, with `header.slots` still 2.

**Where**: `packages/seq/viewer/src/header-entries.ts` (`BuiltInHeaderEntries`); the viewer spec's Behaviors, `headers`,
Built-ins.

**Options**: (a) amend the design's table to "menu only" (the submenu of schemes reads better as a flyout than behind a
slot that opens the menu anyway); (b) give it its slot back; (c) lower `header.slots` to 1 to match.

**Recommendation**: (a), keeping `header.slots` at 2 so a host's entry takes the second column without a setting.

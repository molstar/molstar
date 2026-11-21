# Lane CN: custom navigators

**Task**: build part 2 of [extension-openings](../../../design/extension-openings.md) (cut first): `navigator-picture` registry entries (`<id>:<name>` for a host's), `pictureOf`'s switch becomes the four built-in entries, `content`'s union opened, the window types in the `extension` entry, `NavigatorModel` public as a host-drawn navigator's contract ([navigators](../../../design/navigators.md)); a navigator over several systems stays refused here. Fallback to `axis` when a picture is absent.

**Output**: code, tests (a custom picture drawn and its window driven headless; a picture whose tracks use an extension layer; the built-ins' tests unchanged as the refactor's check), `Features/Navigators` gains a coverage strip and a host-drawn DOM navigator; `_spec` Render (Navigators), viewer Composition.

**Model**: strong (opus).

**Where**: A worktree on a branch `seq/cn-lane` from the base the orchestrator names (after the third slice's close). Hot files (plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a quick-lane test beside the code, failing on the parent for a fix; the look in a visible pane, light and dark, before the report.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface added (names, for the author's acceptance where the brief says so); drift notes; what you saw; the `.v6/changelogs/` entry's text.

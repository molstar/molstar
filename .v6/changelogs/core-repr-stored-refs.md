## Core: stored patches' refs typed by the language

- `@molstar/core/repr/patch`: a stored `set` value (a story step's) is a ref when its key is a ref in the language,
  whatever the spec holds, so it follows renames and moves the same way forward and back, and undoing a patch gives the
  spec back, nested steps included. An edit is no longer refused because it makes one of a step's own ops stop
  applying; a step's ref left dangling or rebound is refused as for any ref.

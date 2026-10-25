## Prism: PyMOL and VMD commands

- New extension `@molstar/command-language-extension` (opt-in; the Prism viewer app lists it): PyMOL and VMD commands
  in the command menu, typed after `pymol `, `vmd ` or `!` (the person's default dialect). A line shows what it will do
  or what does not read, completes commands, arguments, selection keywords, chains and ligands, and runs as one Undo
  step through the viewer's operations. PyMOL's show/hide keep its per-atom meaning on the component tree; VMD's reps
  keep their indices; commands Mol\* does not do fail saying so.
- A Console tab (the rail's foot): runs as a transcript with their status, time and output, a composer that says what
  Enter will do before anything runs, a completion list with the command's signature, history, and script files (`.pml`, `.tcl`, `.vmd`) opened or
  dropped with the files they load; `@file.pml` and `source file.tcl` run included scripts.
- The assistant can run the person's PyMOL or VMD lines (skill `command-language`, tool `run_commands`), one undoable
  step under its approval policy.
- Prism: an extension's tab can have a header shorter than its name (`heading`).
- Prism: extensions can add command-menu prefixes (`ViewerExtension.commandPrefixes`); a prefix row can put text in
  the menu (`run` returns `{ text }`) and post a row's outcome (`CommandPrefixContext.post`); the prefix types moved to
  `commands/types.ts`.
- query-language: PyMOL selections can name objects and named selections through `TranspileOptions.resolveName`.
- Fixed: `not X` in PyMOL and VMD selections selected nothing when `X` selected nothing (`not hydro` in a structure
  without hydrogens).
- Prism: the empty canvas's splash no longer blinks while the window is resized.

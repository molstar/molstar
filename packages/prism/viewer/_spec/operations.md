# Operations

Every operation of `DefaultPrismViewerSpec()` (the core and the built-in features, lazy ones loaded), written from
the code by `src/_test/operations-spec.test.ts` (`UPDATE_SPEC=1 pnpm test:unit packages/prism/viewer/src/_test/operations-spec.test.ts`).
Flags: `read` (no step), `view` (`history: false`), `no undo` (`undoable: false`), the agent's approval (`network`
with `addresses` or always asking), the modules it needs, a `surface` limit. Mol\*'s extensions list theirs in
their own packages (`@molstar/<name>-extension/prism-viewer`).

| Operation | Title | Flags | Params |
| --- | --- | --- | --- |
| `animation.get` | List the animations | read |  |
| `animation.play` | Play an animation | view; needs canvas | name, params |
| `animation.stop` | Stop the animation | view |  |
| `animation.wiggle` | Wiggle |  | mode, amplitude |
| `annotations.apply` | Show what a structure carries | network (addresses) | annotation, structure |
| `annotations.list` | List what the structures carry | read | structure |
| `camera.frame` | Frame | view; needs canvas | target, structure, language, animate |
| `camera.get` | Describe the camera | read |  |
| `camera.orient` | Orient along the principal axes | view; needs canvas | structure, animate |
| `camera.reset` | Reset the camera | view; needs canvas | animate |
| `camera.reset-axes` | Align the view with the x, y and z axes | view; needs canvas | animate |
| `camera.set` | Set the camera | view; needs canvas | position, target, up, animate, distance, radius, sceneRadius, showsAll |
| `camera.zoom` | Zoom | view; needs canvas | factor, animate |
| `canvas.background` | Set the background | needs canvas | color |
| `canvas.camera` | Set the camera | needs canvas | mode, fov |
| `canvas.clipping` | Clip the scene | needs canvas | clip, far |
| `canvas.get` | Describe the canvas | read |  |
| `canvas.lighting` | Change the lighting | needs canvas | occlusion, shadow, outline, fog, illumination |
| `canvas.palette` | Set the data palette |  | palette |
| `canvas.settings` | Change canvas settings | needs canvas; ui-only | section, values |
| `canvas.xr` | Enter or leave XR | view; needs canvas | on |
| `component.add` | Add a component |  | target, label, representation, structure, language |
| `component.list` | List the components | read | structure |
| `component.modify` | Change what a component holds |  | component, target, mode, structure, language |
| `component.options` | Set how representations are made |  | hydrogens, quality, ignoreLight, solidSurface, material, clip, interior, animation |
| `component.remove` | Remove a component | destructive | component, structure |
| `component.rename` | Rename a component |  | component, label, structure |
| `component.visibility` | Show or hide a component |  | component, visible, structure |
| `contents.get` | List the contents | read | structure |
| `contents.visibility` | Show or hide an entity or a chain |  | entity, chain, visible, structure |
| `export.image` | Export an image | view; no undo; needs canvas | size, transparent, crop, box, format, quality, axes, name |
| `export.image-copy` | Copy an image | view; no undo; needs canvas, page | size, transparent, crop, box, format, quality, axes |
| `export.preview` | Preview an image | view; no undo; needs canvas; ui-only | size, transparent, crop, box, format, quality, axes, max |
| `export.viewport` | Export the viewport state | view; no undo; needs canvas | name |
| `focus.again` | Focus again | view | entry |
| `focus.clear` | Clear the focus | view |  |
| `focus.get` | Describe the focus | read |  |
| `focus.history` | List the recent foci | read |  |
| `focus.set` | Focus | view | target, structure, language |
| `highlight.clear` | Clear the highlight | view; ui-only |  |
| `highlight.mark` | Mark residue ranges | view; ui-only | owner, marks, structure, language, dimRest |
| `highlight.set` | Highlight | view; ui-only | target, structure, language, sequence |
| `highlight.unmark` | Clear an owner's marks | view; ui-only | owner |
| `history.list` | List the undoable steps | read |  |
| `history.redo` | Redo | view |  |
| `history.undo` | Undo | view |  |
| `history.undo-turn` | Undo the assistant's turn | view | turn, discardLater |
| `load.action` | Run a loading action | network (asks); ui-only | action, params |
| `load.data` | Load data |  | data, format, label, bytes, models, assembly |
| `load.entry` | Load an entry | network (addresses); needs network | source, id, assembly, representation |
| `load.files` | Open files | needs files | files, visuals |
| `load.mvs` | Load a MolViewSpec | destructive | url, data, format, keepCamera, append, keepCameraOrientation, snapshot, extensions |
| `load.snapshot` | Load a session | destructive; needs network | url, type, id |
| `load.svs` | Load a SeqViewSpec | network (asks) | url, data, file |
| `load.trajectory` | Load a trajectory | network (addresses); needs network | model, modelFormat, modelBinary, coordinates, coordinatesFormat, coordinatesBinary, preset |
| `load.url` | Load from a URL | network (addresses); needs network | url, format, binary, models, assembly |
| `load.viewport` | Load a viewport state | needs canvas | file, json |
| `load.volume` | Load a volume | network (addresses); needs network | url, format, binary, isovalues, entryId, lazy |
| `measure.angle` | Measure an angle |  | a, b, c, text, structure, language |
| `measure.dihedral` | Measure a dihedral angle |  | a, b, c, d, text, structure, language |
| `measure.distance` | Measure a distance |  | a, b, text, structure, language |
| `measure.get` | Describe the measurements | read |  |
| `measure.label` | Add a label |  | a, text, structure, language |
| `measure.options` | Set the measurement options |  | distanceUnitLabel, textColor |
| `measure.orientation` | Show the principal axes |  | target, structure, language |
| `measure.plane` | Show the best-fit plane |  | target, structure, language |
| `measure.remove` | Remove a measurement | destructive | measurement |
| `measure.visibility` | Show or hide a measurement |  | visible, measurement |
| `particles.get` | Describe the particles | read |  |
| `particles.remove` | Remove particles | destructive | particles, representation, node |
| `particles.representation` | Show particles as… |  | type, params, particles |
| `particles.update` | Update a particle representation | ui-only | node, params |
| `particles.visibility` | Show or hide particles |  | visible, particles, representation, node |
| `representation.clipping` | Exclude from clipping |  | target, groups, structure, language |
| `representation.color` | Color |  | target, color, theme, structure, language |
| `representation.emissive` | Set emissive strength |  | target, strength, structure, language |
| `representation.layer-remove` | Remove a style layer |  | layer |
| `representation.layer-update` | Change a style layer |  | layer, value |
| `representation.layers` | List the style layers | read | structure |
| `representation.material` | Set material |  | material, component, structure |
| `representation.remove` | Remove a representation | destructive | component, type, structure, node |
| `representation.reset` | Reset a style layer |  | layer, structure |
| `representation.set` | Show a component as… |  | component, type, params, structure |
| `representation.size` | Set size |  | target, theme, factor, structure, language |
| `representation.transparency` | Set transparency |  | target, alpha, structure, language |
| `representation.update` | Update a representation | ui-only | node, params |
| `representation.visibility` | Show or hide a representation |  | node, visible |
| `select.all` | Select everything | view | structure |
| `select.around` | Select the surroundings | view | radius, asWholeResidues, structure |
| `select.clear` | Clear the selection | view |  |
| `select.expand` | Expand the selection | view | to, structure |
| `select.get` | Describe the selection | read | structure |
| `select.granularity` | Picking level | view | granularity |
| `select.invert` | Invert the selection | view | structure |
| `select.mode` | Selection mode | view | on |
| `select.set` | Select | view | target, mode, structure, language |
| `sequence.colors` | Color the sequence's letters | agent-only | scheme |
| `sequence.follow` | Show the sequence of a structure | view; agent-only | structure |
| `sequence.get` | Describe the sequence view | read |  |
| `sequence.style` | Show the sequence as tracks or text | agent-only | style |
| `session.clear` | Clear the session | destructive |  |
| `session.open` | Open a session | destructive; needs files | file |
| `session.remote.list` | List the remote states | read; network (addresses); needs network | server |
| `session.remote.open` | Open a remote state | destructive; needs network | entry, server |
| `session.remote.remove` | Remove a remote state | no undo; destructive; needs network | entry, server |
| `session.remote.upload` | Upload the state | view; no undo; upload; needs network | name, description, playOnLoad, server |
| `session.save` | Save the session | view; no undo | kind, name, record, options |
| `settings.get` | Read the settings | read | group |
| `settings.set` | Change settings | agent-only | group, values |
| `snapshot.add` | Add a snapshot | no undo | name, description, key, showFor, options |
| `snapshot.apply` | Apply a snapshot |  | entry |
| `snapshot.clear` | Remove all snapshots | no undo; destructive; agent-only |  |
| `snapshot.link` | Follow a link in a description | view; ui-only | href, event |
| `snapshot.list` | List the snapshots | read |  |
| `snapshot.move` | Move a snapshot | no undo; agent-only | entry, position |
| `snapshot.next` | Next snapshot |  |  |
| `snapshot.play` | Play the snapshots | view | delay |
| `snapshot.previous` | Previous snapshot |  |  |
| `snapshot.remove` | Remove a snapshot | no undo; destructive; agent-only | entry |
| `snapshot.replace` | Replace a snapshot with the current view | no undo; agent-only | entry, showFor, options |
| `snapshot.seek` | Show a frame of the snapshot's transition | view; agent-only | time |
| `snapshot.stop` | Stop playing the snapshots | view |  |
| `snapshot.transition` | Play the snapshot's transition | view; needs canvas; agent-only | play |
| `snapshot.update` | Edit a snapshot | no undo; agent-only | entry, name, description, key, showFor |
| `state.action` | Apply an action to a state cell | network (asks); agent-only | node, action, params |
| `state.apply` | Apply a transform | agent-only | node, transformer, params, tags |
| `state.get` | Read the state tree | read | node |
| `state.remove` | Remove state cells | destructive | node |
| `state.select` | Find state cells | read | node |
| `state.update` | Update state cells | agent-only | node, params |
| `state.visibility` | Show or hide state cells |  | node, visible |
| `structure.assembly` | Show an assembly |  | id, structure |
| `structure.get` | Describe the scene | read |  |
| `structure.model` | Show another model | view | index, structure |
| `structure.preset` | Apply a representation preset |  | preset, structure |
| `structure.remove` | Remove a structure | destructive | structure |
| `structure.symmetry` | Show crystal symmetry |  | kind, radius, cells, structure |
| `structure.transform` | Move or rotate a structure |  | translation, axis, angle, center, matrix, compose, structure |
| `structure.transform-reset` | Reset a structure’s transform |  | structure |
| `structure.unit-cell` | Show the unit cell |  | visible, color, opacity, structure |
| `structure.visibility` | Show or hide a structure |  | visible, structure |
| `styles.apply` | Apply a style |  | style |
| `superpose.atoms` | Superpose by atoms |  | language, targets |
| `superpose.chains` | Superpose by chains |  | language, targets, alignSequences, traceOnly |
| `superpose.db` | Superpose by UniProt (SIFTS) |  | structures, traceOnly |
| `superpose.ligands` | Superpose by ligands |  | language, targets |
| `superpose.list` | List the superpositions | read |  |
| `superpose.remove` | Remove a superposition |  | structure |
| `superpose.tm-align` | Superpose by TM-align |  | language, targets, traceOnly |
| `trajectory.frame` | Show a frame | view | structure, index, by |
| `trajectory.get` | Describe the trajectories | read |  |
| `trajectory.play` | Play the frames | view; needs canvas | mode, seconds |
| `trajectory.stop` | Stop playing the frames | view |  |
| `view.close` | Close the open panel and card | view; agent-only |  |
| `view.expand` | Expand the viewport | view | on |
| `view.fullscreen` | Fullscreen | view; needs page | on |
| `view.open` | Open a panel | view; agent-only | tab, card, section |
| `volume.get` | Describe the volumes | read |  |
| `volume.isovalue` | Set the iso value |  | value, kind, volume, representation, node |
| `volume.remove` | Remove a volume | destructive | volume, representation, node |
| `volume.representation` | Show a volume as… |  | type, params, volume |
| `volume.streaming` | Stream density for a structure | network (addresses); needs network | on, structure |
| `volume.streaming-view` | Set how streamed density shows | network (addresses); needs network | view, radius, detail, channels, structure |
| `volume.update` | Update a volume representation | ui-only | node, params |
| `volume.visibility` | Show or hide a volume |  | visible, volume, representation, node |

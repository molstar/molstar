# Plan: What the draw table observes

## Status

- **Paused**: the draw table ([draw-table.md](draw-table.md)) observes every value cell of an object. Step 4 done
  (observing all of 94k members' cells, about 2.8M, hung a merged NPC capsid). Steps 2, 3 and 5 not taken for now:
  measured in the Viewer (36ZA as spacefill unmerged, 2888 objects, WebGPU, heap after commit, collected), observing
  the value cells costs 23 MB of the representation's 550 MB (about 8 KB an object, 4%), no measurable time (4.4 s to
  build and commit either way). Steps 2 and 3 would save under half of it, step 5 most of it (a public `ValueCell`
  change); neither pays for its risk or complexity while scenes of 10k units and more go structure-level and merged
  objects observe their members' marks only. Taken up if many unmerged objects (texture meshes, `renderMerged` off)
  show the observers as a larger share.

## Question

The table observes every value cell of an object (`getValueCells`): 129 for a mesh (52 uniforms, 30 defines, 18
textures, 7 attributes, some 20 others), each with a `Set` of observers. 36ZA unmerged is 2888 objects, about 372k
cells: some 40 ms to observe and some 75 MB of heap (a `Set` a cell, measured in Node; the uniforms' `Set`s exist
already, held by the uniform layout's block observers). A merged object adds its members' marker cells, 2 a member
(`Renderable.memberCells`, step 4). That's the cost of one rule, "any change of an object's values advances its values
version", which no consumer needs in full.

What reads what the table observes:

1. **Rows** (`SlotObserver`): `DrawFactCells`, `boundingSphere`, `uLod` and the state. Unchanged.
2. **Values versions** (`valuesVersions`, `valuesEpoch`): the renderer's plans and recordings (`walkPass`,
   `isChunkInputsCurrent`, `flushDataChanges`, `walkCaps`). A changed version means "draw it again without encoding
   and compare what it binds" (`captureDraws`): what matters is what `render()` writes or binds. On WebGPU
   (`webgpu/render-item.ts`) that is the uniform blocks (material, object, draw), the vertex buffers (attribute
   versions compared in `getVertexBuffer`, written when drawn), the index buffer, the textures and data layers
   (`update`), `drawIndirect`, and the pipeline: a define change replaces the program (`invalidatePrograms`), seen only
   because the define cell advances the values version (`Renderable.programVersion` advances only with
   `setTransparency`). `drawCount` and `instanceCount` are read through `update` too: a row cell, but a changed count
   must flush the draws (not only re-walk). Not read by drawing: the averages and minimum (rows), the `*Type` and
   `meta` values, `matrix`, `transform` and `extraTransform` (inputs of `aTransform` and of culling, which have their
   own versions: `drawVersion`), `hasReflection`, `clipObjects` (`getClipInfo`).
3. **Changed objects** (`takeChanged`, `Scene.updateChanged`): what `update()` reads and changes outside a
   `Scene.update` of the object. The only caller is marking (`canvas3d.mark`): `markerStatus` and `tMarker`, and a
   merged object's members' (`sync` reads `markerStatus` and the textures' versions; every other member change comes
   with an update of the merged object).

## Proposal

1. **Audit** what `render()` and `update()` read on both backends (the list above is WebGPU's) and what changes
   outside a `Scene.update`: visual setters (`setOverpaint`, `setTransparency`, `setEmissive`, `setSubstance`,
   `setClipping`, theme changes) and how each reaches the scene; a member change of a merged object other than a mark.
   The result is the role of each schema entry, recorded in `draw-table.md`.
2. **Defines through the program version**: `Renderable.programVersion` advances with the render item's program
   invalidations (defines, layer pages, replaced textures: `invalidatePrograms`), on both backends, so that a recording
   is re-recorded on a define change without observing the defines.
3. **Values cells by role**: `getValueCells` takes the cells of the kinds drawing reads, by the render item's schema
   (`RenderableSchema`): uniforms, attributes, elements, textures and data layers, and the values `drawCount`,
   `instanceCount`, `drawIndirect`; not the defines (2) nor the other values. About 80 of 129 for a mesh, the uniforms
   among them already observed.
4. **Merged members' marks only** (done): `Renderable.memberCells` are the members' marker cells (`markerStatus`,
   `tMarker`): 2 a member instead of about 30.
5. **A single observer without a `Set`** (`@molstar/core` `ValueCell`): `observers` holds one observer, a `Set` from
   the second on (`observe`, `unobserve`, `notify`). The cells that forward `observers` to a derived cell
   (`affine-rows.ts`, `depth-range.ts`, `oct-normals.ts`) follow. For 372k cells about 8 ms and 15 MB instead of 40 ms
   and 75 MB (Node, a field added after creation; declared in `ValueRef.create`, less). A public API of `@molstar/core`:
   agreed with the maintainers first. The uniforms keep two observers (block and table) and their `Set`s: see the open
   questions.
6. **Measured**: 36ZA unmerged (2888 objects) and merged (Gaussian surface on WebGPU, 2521 members): cells observed,
   heap after adding all objects (snapshot, before and after), time to add and remove them, the hover marking time per
   frame (unchanged or less). Tests:
   - `draw-table.test.ts`: values versions advance for a uniform, an attribute, a texture, `drawCount`, not for a
     define or an average; a merged object is listed by a member's mark;
   - `renderer` recordings: a define change re-records, a uniform, attribute or count change flushes;
   - `value-cell` in core: one, two and none observers, unobserve back to one;
   - images: the marking and the recording cases (`marking/*`, the draw-list cases) identical on both backends.

## Open questions

- Whether the uniforms need the table at all: their block observers (`uniform-layout.ts`) already flag a change; a
  renderable's data version advanced from them (and from `update()` where it writes or replaces data) would let the
  table observe none of the drawing cells, only the rows and the marks. Taken up if 3 and 5 leave the uniforms' `Set`s
  a measurable share.
- Whether the changed list (3) should stay general (every cell `update()` reads) for callers other than marking, at
  the cost of observing the textures and attributes for it too: they are observed for the values versions anyway (3),
  so it costs nothing until those are narrowed further.

# Plan: Two-phase occlusion culling of GPU-culled objects

## Status

- **Current**: built, on by default (`hiZ.twoPhase`), for objects culled on the GPU, the sphere raster included.
- **What it does**: the first phase culls as one phase did, against the previous frame's pyramid, and draws. A pyramid of
  that depth, with the frame's camera, then tests every instance the first left out; those visible are drawn in the same
  frame.

## Why

- **One phase loses what comes into view**: the test projects with the previous frame's camera into its pyramid. Zooming
  or rotating, what the old view hid is missing for a frame.
- **And keeps it lost**: the render loop goes idle once the camera stops, so what the last moving frame left out stays
  missing until something else draws. Zooming into a mesoscale structure in the viewer shows it.
- **Occlusion has to stay**: without it the frames below cost about twice as much.

## Design

- **Exact**: the first phase's depth is part of the frame's final depth, so what the second leaves out is hidden in the
  finished frame. No catch-up draw is needed once the camera stops.
- **No state across frames**: "drawn by the first phase" is rebuilt per dispatch, the top bit of the instance's cell
  (`DrawnBit`, in the word only its own thread reads). So instances do not flip from frame to frame and no third step
  corrects visibility: the earlier design, the previous frame's visible set drawn first with no pyramid, needed both, and
  cost 40-70 % in the opaque pass while moving zoomed in.
- **Second phase** (`MergedCuller.dispatchSecond`): cells and instances against the frame's pyramid, the counters going
  on, so a region holds the first's instances then the second's. The args kernel writes the second's records (`second`:
  count, first instance moved by the first's count) and the sum into the main records, for the frame's later passes and
  still frames.
- **Opaque part 2** (`DrawPass.planOpaqueSecond`): pyramid, second phase, then a pass on the same target, loading. Drawn
  directly, no recording: it exists only while moving. No caps in this part: while a second phase is to come the first part draws caps from whole segments, also where a cap cuts all of them (the object's lists hold the first phase's instances only).
- **Raster** (`SphereRaster.second`): the regions kernel over the second's records appends to the same chunks; words 2 and
  3 of the object's counter mark what is drawn (`uMark`), so args and the raster take the new chunks and rows only. The
  visibility buffer is not cleared, the atomic max keeps the nearest. The resolve runs again; pixels of the first keep
  their depth. Chunk counts are read back after the second phase.
- **Two pyramids a moving frame**: the one within the frame, and the one after it of the complete depth, which is the
  next frame's first phase and the CPU test's (`planFrame` starts no read).
- **Still frames**: no dispatch, no second phase, recordings replay as before. No pyramid in the first phase: it left out
  nothing occluded, no second.
- **Multi-sampling**: both phases in the first sample's pass; later samples draw the sum.

## No occlusion

- **Stereo and XR**: the pyramid is the whole viewport's with the main camera, not an eye's: tested against it, 102,172 of
  921,600 pixels differed on a still camera, for good. So the eyes test nothing and build no levels (their depth would hide what the next mono frame sees): exact, at no occlusion's cost (zoomed
  8.46 ms, 4.13 with the wrong test, 8.61 without occlusion). A pyramid per eye would bring it back.

## Still one phase

- **Unmerged objects drawn by the raster** (`sphereRaster` `on` or `only`, not `auto`): their instances are tested against
  the previous frame's pyramid.
- **Blended to the drawing buffer**: the depth is not in the texture the pyramid reads, so it occludes nothing.

## Results

`occlusion-dolly/*`, `occlusion-lag/*`, 1280×720. Pixels of 921,600 differing from no occlusion by more than 24 levels,
in the frame after the last; the frame's median in ms.

| Case                    | Lost, one phase | Lost, two | One phase | Two  | None |
| ----------------------- | --------------- | --------- | --------- | ---- | ---- |
| PetWorld zoomed, 3 %    | 536             | 0         | 2.21      | 2.83 | 5.11 |
| MG_1184 zoomed, 3 %     | 370             | 0         | 2.46      | 3.50 | 6.08 |
| PetWorld whole, 8 %     | 3,683           | 0         | 1.60      | 2.23 | 4.04 |
| PetWorld zoomed, 2°     | 105             | 0         | 1.26      | 1.67 | 3.10 |
| PetWorld zoomed, 30°    | 53              | 0         | 1.53      | 2.06 | 3.21 |

- **Also exact**: PetWorld zoomed at 1 % (1,601 lost in one phase), MG_1184 whole at 8 % (1,410).
- **Where the cost goes**: moving frames pay 0.4-1.0 ms. The second raster pass 0.4-0.7 ms, the second cull 0.05-0.14,
  the opaque second pass 0.05, submit 0.15-0.35 (the part's indirect records, most of them empty).
- **Still**: `sphere-raster-mg/whole-auto` runs no second phase, its opaque pass replays.
- **Multi-sampling, moving** (`occlusion-view/*`): exact, 9.64 ms against 19.5 without occlusion.
- **By stage**, zoomed, ms of GPU time: the raster kernels 0.26-0.67, growing with the step (the new spheres); the second
  regions pass 0.04-0.20, the second cull 0.05-0.14, the opaque second pass 0.04-0.07, these three by model, not by step.

## Open

- **Skipping the second phase where it finds nothing**: not built. It would save the fixed stages above, 0.13-0.41 ms, and
  the second part's submit, only on frames finding nothing, and nothing tells before the cull ran. Most of the cost is
  drawing what it found.
- **The second part's empty records**: bounds predict them only; a readback cannot help a moving camera.
- **A rare worker crash in the native WebGPU tests, parked**: on Windows, about 1 full run of `webgpu*` in 15, a worker
  exits with `0xC0000005` after its tests passed; cause unknown, and unknown whether it predates this work (16 runs
  without the new test or on the commit before were clean, too few to tell). Not reproduced by Dawn alone (a device,
  in-flight compiles, submitted work, several devices: 245 runs) nor by the files it names alone (35 runs). Next step: a
  crash dump (`Windows Error Reporting\LocalDumps`), not more runs. A shell `timeout` near the suite's 50 s looks alike
  but exits 143.

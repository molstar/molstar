## Graphics: merged units and GPU-culled views

- A merged units visual (`UnitsMergedVisual`): a group a segment of one merged object, its cells written through one
  view over shared values, no visual, values or render object a group: 100k units 5.4 -> 3.6 s, 22.6 -> 14.3 KB a
  unit; a hover at 100k units 0.9 ms. Merged objects built from records and shared values (`MergeRecord`): 24 -> 2.1 KB
  a member, the yeast cell's JS heap 225.8 -> 172 MiB.
- Every merged object culled on the GPU culls a view of its own there (the thickness pass, a ray's pick, a stereo or
  XR eye), not on the CPU: PetWorld zoomed with illumination 27 -> 6.5 ms a frame, MG_1184's CPU submit 2.7 -> 0.9 ms
  (its 46 MB of CPU grids not built).
- In timing mode the graphics timer's passes show on a browser profile's tracks: where they ran on the GPU and when the
  CPU encoded them.

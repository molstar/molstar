## Graphics: faster cartoon geometry

- Cartoon, tube and direction-wedge geometry build faster with identical output: the curve interpolation and the tube
  and sheet builders compute on scalars instead of `Vec3` calls (186k atoms `repr:cartoon`: 169 → 129 ms, allocation
  341 → 212 MB; `createPolymerTraceMesh` 90 → 52 ms)

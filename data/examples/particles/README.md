# Particle fixtures

Small, real particle files, one per reader, for tests and stories. Each is parsed in
`packages/model/src/formats/particles/_test/fixtures.test.ts`. All were downloaded on 2026-10-21; none was edited except
by dropping trailing rows as noted.

| File                      | Format                                 | Size    | Particles              | Pixel size            |
| ------------------------- | -------------------------------------- | ------- | ---------------------- | --------------------- |
| `cryoet-ribosomes.ndjson` | CryoET Data Portal oriented points     | 63.2 KB | 200                    | 6.802 Å (not in file) |
| `relion-subtomo.star`     | RELION STAR (optics + particles)       | 26.8 KB | 300, one tomogram      | 10 Å (optics block)   |
| `dynamo.tbl`              | Dynamo table (35 columns)              | 40.6 KB | 300, tomograms 1 and 2 | not in file           |
| `artiatomi.em`            | Artiatomi/TOM motivelist (EM, float32) | 18.5 KB | 225                    | 13 Å (not in file)    |
| `small.simularium`        | Simularium binary (v3)                 | 4.5 KB  | 5 spheres, 5 frames    | nm units              |

## cryoet-ribosomes.ndjson

- **Source**:
  https://files.cryoetdataportal.cziscience.com/10003/TS_01/Reconstructions/VoxelSpacing6.802/Annotations/100/chloramphenicol_bound_70s_ribosome-1.0_orientedpoint.ndjson
  (annotation metadata: the `chloramphenicol_bound_70s_ribosome-1.0.json` next to it).
- **Entry**: CZ CryoET Data Portal dataset 10003 ("70S ribosome with chloramphenicol", _Mycoplasma pneumoniae_ M129),
  run TS_01, annotation 100 (cytosolic ribosome, GO:0022626; template matching + RELION classification; authors Liang
  Xue, Patrick Cramer, Julia Mahamid). Dataset authors: Dimitry Tegunov, Liang Xue, Christian Dienemann, Patrick Cramer,
  Julia Mahamid; related EMPIAR-10499, EMPIAR-10731; publications doi:10.1038/s41592-020-01054-7,
  doi:10.1038/s41586-022-05255-2.
- **License**: CC0 (the portal's data submission policy: public data is available "subject to a CC0 license"). No
  attribution required; please cite the publications above and the CryoET Data Portal.
- **Trimmed**: the first 200 of 352 lines. Coordinates are in voxels of the 6.802 Å reconstruction.

## relion-subtomo.star

- **Source**: https://zenodo.org/records/12743310/files/01_10.00Apx.star (doi:10.5281/zenodo.12743310)
- **Entry**: "Particle picks of HIV virus-like particles from cryo-electron tomography data", Lorenzo Gaifas (2024),
  modified from Alister Burt, "3D reconstruction of HIV virus-like particles from cryo-electron tomography data"
  (doi:10.5281/zenodo.6504890) to add an optics block with the pixel size. Picks and poses in tomogram 01 (10 Å/px) of
  EMPIAR-10164.
- **License**: CC-BY 4.0. Attribution: "Particle picks of HIV virus-like particles from cryo-electron tomography data,
  Lorenzo Gaifas, Zenodo, doi:10.5281/zenodo.12743310, CC-BY 4.0; based on Alister Burt, doi:10.5281/zenodo.6504890".
  Changes: trimmed.
- **Trimmed**: header, `data_optics` and the `data_particles` loop header kept as they are; the first 300 of 5377
  particle rows. Note the file is RELION-3-style (`rlnCoordinateX/Y/Z` in pixels, `rlnAngleRot/Tilt/Psi`,
  `rlnMicrographName`, no `rlnTomoName`), not a RELION 4/5 tomography project file; it also carries a non-RELION
  `_experiment_id` column.

## dynamo.tbl

- **Source**:
  https://raw.githubusercontent.com/teamtomo/dynamotable/f4198fd1372a1c6f051d3bd269d49d95b57f2924/test/data/t_double.tbl
- **Entry**: the test table of the `dynamotable` Python package (teamtomo), 300 rows over tomograms 1 and 2, 35 columns.
  Its companion `tmap_complex.doc` (same folder) maps tomograms to `tomo_ts01.mrc_15.00Apx.mrc` etc.
- **License**: BSD 3-Clause. Required notice (retain it with the file):

  > Copyright (c) 2020, CCP-EM. All rights reserved.
  >
  > Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
  > following conditions are met: 1. Redistributions of source code must retain the above copyright notice, this list of
  > conditions and the following disclaimer. 2. Redistributions in binary form must reproduce the above copyright
  > notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided
  > with the distribution. 3. Neither the name of the copyright holder nor the names of its contributors may be used to
  > endorse or promote products derived from this software without specific prior written permission.
  >
  > THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
  > INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  > DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  > SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  > SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
  > WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  > OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

- **Trimmed**: no, the whole file.

## artiatomi.em

- **Source**: `mycoplasma_genitalium/demo/motivelist_ribosomes.em` inside
  https://zenodo.org/records/7116187/files/mycoplasma_genitalium.zip (doi:10.5281/zenodo.7116187; the single member read
  by HTTP range requests, the archive is 611 MB).
- **Entry**: "Demo datasets for: ArtiaX: An Electron Tomography Toolbox for the Interactive Handling of Sub-Tomograms in
  UCSF ChimeraX", Utz Heinrich Ermel, Serena M. Arghittu, Achilleas S. Frangakis (2022). Ribosomes in a _Mycoplasma
  genitalium_ tomogram (Seybert et al. 2018, Mol. Microbiol. 108:319–329). The demo's ChimeraX script sets the pixel
  size to 13 Å.
- **License**: CC-BY 4.0. Attribution: "ArtiaX demo dataset, Ermel, Arghittu, Frangakis, Zenodo,
  doi:10.5281/zenodo.7116187, CC-BY 4.0".
- **Trimmed**: no, the whole file (PC byte order, float32, 20 × 225 × 1).

## small.simularium

- **Source**:
  https://raw.githubusercontent.com/simularium/simularium-website/2e665cb45603d04dbdb06f9a0ea163ca783a53f3/test-files/example_plots_binary.simularium
- **Entry**: a test file of the Simularium website (Allen Institute for Cell Science): binary format, trajectory info
  v3, 5 frames of 5 sphere agents (17 agent types declared), plus plot data.
- **License**: Apache License 2.0 (the `simularium-website` repository). Attribution: "From
  simularium/simularium-website, Copyright 2022 Allen Institute for Cell Science, licensed under the Apache License 2.0
  (https://www.apache.org/licenses/LICENSE-2.0)".
- **Trimmed**: no, the whole file. The larger example trajectories in the `aics-simularium-data` S3 bucket carry no
  license, so they were not used.

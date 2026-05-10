/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author ReliaSolve <russ@reliasolve.com>
 */

import { ReaderResult as Result } from '../../../mol-io/reader/result';
import { Task, RuntimeContext } from '../../../mol-task';
import { Kinemage } from './schema';
import KinParser from './kinparser';

async function parseInternal(data: string, ctx: RuntimeContext): Promise<Result<Kinemage[]>> {
  const kinemages: Kinemage[] = [];
  // Split the data into sections based on the '@kinemage' keyword, which indicates one or more kinemages in the file.
  // Handle the case where there is no '@kinemage' keyword by parsing the entire file.
  const kinemageSections = data.split(/@kinemage\s+\d+/); // Split based on '@kinemage' keyword followed by a number

  // If there are one or more @kinemage sections, ignore the portion before the first one.
  // This will either be an empty string (if the first section starts at the beginning of the file)
  // or header data that is not part of a particular kinemage.  This has the effect of removing
  // the header data even in the case where there is a single @kinemage keyword.
  if (kinemageSections.length > 1) {
    kinemageSections.shift();
  }

  for (const section of kinemageSections) {
    if (section.trim()) { // Ignore empty sections
      const NGLParser = new KinParser(section.trim());
      const kinData = NGLParser.kinemage;
      kinemages.push(kinData);
    }
  }

  return Result.success(kinemages);
}

export function parseKin(data: string) {
    return Task.create<Result<Kinemage[]>>('Parse KIN', async ctx => {
        return await parseInternal(data, ctx);
    });
}

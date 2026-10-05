import * as fs from 'fs';

export async function openRead(filename: string) {
  return new Promise<number>((res, rej) => {
    fs.open(filename, 'r', async (err, file) => {
      if (err) {
        rej(err);
        return;
      }

      try {
        res(file);
      } catch (e) {
        fs.closeSync(file);
      }
    });
  });
}

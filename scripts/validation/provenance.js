// SPDX-License-Identifier: GPL-2.0-or-later
import {readFileSync,readdirSync} from 'node:fs';
import {createHash} from 'node:crypto';
export const sha256=path=>createHash('sha256').update(readFileSync(path)).digest('hex');
export function numericalSourceHashes(extra=[]){
  const sources=['potential','inviscid','viscous','euler','numerics','geometry'].flatMap(dir=>
    readdirSync(`src/${dir}`,{recursive:true}).filter(p=>/\.(js|wasm|c|h)$/.test(p)).map(p=>`src/${dir}/${p}`));
  return Object.fromEntries([...new Set([...sources,'scripts/validation/provenance.js',...extra])].sort().map(path=>[path,sha256(path)]));
}
export const changedSources=hashes=>Object.keys(hashes).filter(path=>sha256(path)!==hashes[path]);

// SPDX-License-Identifier: GPL-2.0-or-later
// Rebuild with a WASI-enabled clang/LLVM installation; no Emscripten required.
import {readdirSync,readFileSync,writeFileSync,renameSync,rmSync} from 'node:fs';
import {fileURLToPath} from 'node:url';
import {spawnSync} from 'node:child_process';
import {createHash} from 'node:crypto';
const directory=fileURLToPath(new URL('.',import.meta.url)),source=directory+'source/';
const compiler=process.env.WASM_CC??'clang';
const sources=['KLU','AMD','BTF','COLAMD'].flatMap(part=>readdirSync(source+part+'/Source').filter(p=>p.endsWith('.c')&&p!=='amd_global.c').sort().map(p=>source+part+'/Source/'+p));
sources.push(source+'SuiteSparse_config/SuiteSparse_config.c',directory+'bridge.c');
const options=['--target=wasm32-wasi','-O3','-DNDEBUG','-DNPRINT=','-DNTIMER','-DDINT','-mexec-model=reactor',
  '-Wl,--strip-all','-Wl,-z,stack-size=8388608','-Wl,--initial-memory=16777216','-Wl,--max-memory=2147483648'];
if(process.env.WASI_SYSROOT)options.push('--sysroot='+process.env.WASI_SYSROOT);
const includes=['KLU','AMD','BTF','COLAMD'].map(p=>'-I'+source+p+'/Include').concat('-I'+source+'SuiteSparse_config');
const exports=['malloc','free','mses_klu_factor','mses_klu_factor_with_pivot','mses_klu_solve','mses_klu_free','mses_klu_status','mses_klu_nnz','mses_klu_factor_given','mses_klu_btf','mses_klu_blocks'].map(p=>'-Wl,--export='+p);
const temporary=directory+'klu.next.wasm';
const built=spawnSync(compiler,[...options,...includes,...sources,'-lm','-o',temporary,...exports],{stdio:'inherit'});
if(built.error)throw built.error;if(built.status!==0){rmSync(temporary,{force:true});process.exit(built.status??1);}
const hash=path=>createHash('sha256').update(readFileSync(path)).digest('hex');
const provenance=JSON.parse(readFileSync(directory+'PROVENANCE.json','utf8'));
for(const [path,expected]of Object.entries(provenance.sourceHashes))if(hash(directory+path)!==expected)throw new Error('Upstream source changed: '+path);
const module=await WebAssembly.compile(readFileSync(temporary));
if(WebAssembly.Module.imports(module).length)throw new Error('The sparse math kernel must not require host imports.');
renameSync(temporary,directory+'klu.wasm');
const compilerVersion=spawnSync(compiler,['--version'],{encoding:'utf8'}).stdout.trim();
writeFileSync(directory+'BUILD.json',JSON.stringify({compilerVersion,target:'wasm32-wasi',options:options.filter(s=>!s.startsWith('--sysroot=')&&!s.startsWith('--ld-path=')),
  bridgeSha256:hash(directory+'bridge.c'),buildScriptSha256:hash(directory+'build.js'),wasmSha256:hash(directory+'klu.wasm'),
  librariesNote:'Requires WASI libc/libm and LLVM wasm32 compiler builtins; preserve the supplied copyright notices when redistributing.'},null,2)+'\n');
console.log('Built '+directory+'klu.wasm');

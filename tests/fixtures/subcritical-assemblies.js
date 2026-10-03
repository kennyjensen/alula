// SPDX-License-Identifier: GPL-2.0-or-later
import {naca4,transform} from '../../src/geometry/airfoil.js';

// Deliberately separated thin wakes: these cases test coupled discretization,
// not wake confluence, maximum lift, or an experimental high-lift geometry.
export function slatMainFlap(counts=[120,240,120]){
  return{elements:[
    {name:'Slat',points:transform(naca4('0012',counts[0]),{chord:.2,x:-.24,y:.12,angle:-3})},
    {name:'Main',points:naca4('0012',counts[1])},
    {name:'Flap',points:transform(naca4('0012',counts[2]),{chord:.3,x:1.05,y:-.15,angle:-5})}],
    alpha:0,mach:.2,reynolds:1e6,ncrit:9,trips:[.1,.1],wakeCount:48,wakeLength:2};
}

export function mainFlap(panels=160){
  return{elements:[{name:'Main',points:naca4('0012',panels)},
    {name:'Flap',trips:[.08,.12],points:transform(naca4('0012',panels),{chord:.3,x:1.05,y:-.2})}],
    alpha:2,mach:.2,reynolds:1e6,ncrit:9,trips:[.05,.05],wakeCount:48,wakeLength:2};
}

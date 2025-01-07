function f=upwind(u,v,a)
     f = a/2*(u+v)+abs(a)/2*(u-v);
return

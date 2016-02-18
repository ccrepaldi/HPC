`! unroll fortran do-loops with a desired depth'
`! syntax unroll(index1, nmax, depth, calculation, index2)'
`! index1 and 2 must be different but the output will be with index1 only'
define(`forloop',
       `pushdef(`$1', `$2')_forloop(`$1', `$2', `$3', `$4')popdef(`$1')')
define(`_forloop',
       `$4`'ifelse($1, `$3', ,
                   `define(`$1', incr($1))_forloop(`$1', `$2', `$3', `$4')')')dnl
define(`unroll',`define(`var',$1)dnl
m = mod($2,$3)
do var=1,m
   a(var) = var
end do
dnl
do var = m+1,$2,$3
forloop(`j',0,decr($3),`  define($5, var+`j')$4
')dnl
end do')dnl
define(NMAX, `n')dnl
`! Unroll: depth=1'
`call cpu_time(t1)'
unroll(`i',`n',1,`a(k) = k',`k')
`call cpu_time(t2)'
`print *, "Elapsed time:", t2-t1'

`! Unroll: depth=2'
`call cpu_time(t1)'
unroll(`i',`n',2,`a(k) = k',`k')
`call cpu_time(t2)'
`print *, "Elapsed time:", t2-t1'

`! Unroll: depth=4'
`call cpu_time(t1)'
unroll(`i',`n',4,`a(k) = k',`k')
`call cpu_time(t2)'
`print *, "Elapsed time:", t2-t1'

`! Unroll: depth=8'
`call cpu_time(t1)'
unroll(`i',`n',8,`a(k) = k',`k')
`call cpu_time(t2)'
`print *, "Elapsed time:", t2-t1'

`! Unroll: depth=16'
`call cpu_time(t1)'
unroll(`i',`n',16,`a(k) = k',`k')
`call cpu_time(t2)'
`print *, "Elapsed time:", t2-t1'

`! Unroll: depth=32'
`call cpu_time(t1)'
unroll(`i',`n',32,`a(k) = k',`k')
`call cpu_time(t2)'
`print *, "Elapsed time:", t2-t1'

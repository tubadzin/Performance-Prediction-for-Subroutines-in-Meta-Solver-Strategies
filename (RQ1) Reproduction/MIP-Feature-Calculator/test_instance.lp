\* test_instance.lp *\

Maximize
 obj:  x1 + 2 x2 + 3 x3

Subject To
 c1: x1 + x2 + x3 <= 100
 c2: 2 x1 + x2      >=  50
 c3: x2 + 3 x3      =   80

Bounds
 0 <= x1 <= 40
 2 <= x2 <= 50
 x3 free

End

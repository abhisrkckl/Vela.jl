# pint2vela

Example usage:

```
import pint2vela as p2v
vl = p2v.vl

mv, tv = p2v.read_model_and_toas("pure_rotator.par", "pure_rotator.tim")
lnlike = vl.get_lnlike_parallel_func(mv, tv)

params = np.array((0.1027149062917101, 0.0, -4.19778948625126e-8))
print(lnlike(params))
```

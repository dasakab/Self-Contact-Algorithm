# Verification record and limitations

- Algorithm2 facet-search source compiled under Linux with GNU Fortran, bounds and floating-point checks in the earlier trial.
- Bounded small-step continuation from saved pre-contact increment 35 accepted increment 36 (residual 0.0004248838594) and increment 37 (0.0004742714517), both 54 Newton iterations, both showing zero nonincident surface-triangle crossing pairs in that test; increment 38 was interrupted before acceptance.
- A larger-step Algorithm2 trial converged some increments but yielded 192 triangle crossings at increment 37. Therefore strict force-residual convergence alone cannot be used to certify non-penetration.
- Full t=0 history, full tube collapse, later second contact ring, and native Windows compilation were not demonstrated in that trial. The current packaged fresh-start launcher needs local end-to-end testing.
- Re-run the independent penetration audit on all accepted increments before interpreting stresses or declaring success. The legacy restart state omission and contact pair logic require further work.

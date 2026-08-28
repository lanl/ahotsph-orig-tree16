c © 2026. Triad National Security, LLC. All rights reserved.
c This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

program makeh

        integer
    * 4 i,
    j,
    k integer
        * 4 dim

        parameter(dim = 10000)

            real
        * 8 h(dim),
    rho(dim), r(dim), u(dim), vel(dim), pres(dim), tem(dim) real * 8 dm(dim), tmp(dim),
    ye(dim) real * 8 length, mass,
    convf

    length
    = (6.955d10) !length unit - Rsol - in cm mass
    = (2.d33 * 1.d - 6) !mass unit - 1E-6 Msol - in grams time
    = (1.d2) !time unit - 100s
      - in sec

      !length
    = 6.955d11 !10. Rsol in cm !mass
    = 1.9889d25 !1e-8 msol in g !time
    = 3.600d3 !1 hr in sec

    !length
    = 1.50d12 !0.1 AU in cm !mass
    = 1.9889d25 !1e-8 msol in g !time
    = 3.15576d5 !0.01 yr in sec

    open(10, file = 'model.dat') open(12, file = 'inputmodel.dat')

        i
    = 1 100 continue read(10, '(1p6e12.4)', end = 101) c 1 dm(i),
    tmp(i), c 1 r(i), rho(i), vel(i), tem(i), ye(i) 1 r(i), rho(i), u(i), pres(i), vel(i) 2,
    tem(i)

        write(12, '(1p6e12.4)')(r(i) / length),
    1(rho(i) / mass * (length * length * length)), c 2 1.d30, c 3 1.d30,
    2(u(i) / (length * length) * (time * time)), 3(pres(i) / mass * length * (time * time)),
    4(vel(i) / length * time),
    c 5 1.d30 5 tem(i) i
    = i
      + 1 goto 100 101 continue

      write(*, *) i

      close(10) close(12)

          open(11, file = 'inputh.dat') do j
    = 1,
    i
        - 1 c the coefficient approx.corresponds to the mass of one particle c e.g.9.4 Msol star in
          1M particles ~10 mass code units(1e-6Msol) h(j)
    = 1.0d
      + 0 / (rho(j) / mass * (length * length * length))
            * *(1.d0 / 3.d0) c if (r(i).gt .7.d5 * length.and.r(i).lt .1.5d6 * length) then c h(j)
    = h(j) * 0.5 c endif write(11, '(1p2e21.13)')(r(j) / length),
    (h(j))enddo

    close(11) end

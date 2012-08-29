      doubleprecision function cg (mu, RX, RY, FX, FY, E, B, n)
        doubleprecision mu
        doubleprecision RX
        doubleprecision RY
        doubleprecision FX
        doubleprecision FY
        doubleprecision E
        doubleprecision B
        doubleprecision n
        doubleprecision t13
        doubleprecision t6
        doubleprecision t7
        doubleprecision t20
        doubleprecision t23
        doubleprecision t17
        doubleprecision t19
        doubleprecision t15
        doubleprecision t16
        doubleprecision t5
        doubleprecision t48
        doubleprecision t9
        doubleprecision t10
        doubleprecision t2
        t2 = B ** (-n)
        t5 = dble(2 ** (0.1D1 + n))
        t6 = mu ** 2
        t7 = t6 * E
        t9 = n / 0.2D1 - 0.1D1 / 0.2D1
        t10 = t7 ** t9
        t13 = RX ** 2
        t15 = 0.2D1 * RX * FX
        t16 = FX ** 2
        t17 = RY ** 2
        t19 = 0.2D1 * RY * FY
        t20 = FY ** 2
        t23 = (t13 - t15 + t16 + t17 - t19 + t20 + 0.4D1 * t7) ** t9
        t48 = t13 + t17 + t16 + n * t20 + t20 - 0.2D1 * n * RX * FX - 0.
     #2D1 * n * RY * FY + n * t13 + n * t17 + n * t16 - t15 - t19
        cg = -0.2D1 * (RY - FY) * t2 * (t5 * t10 * t7 - t23 * t13 + 0.2
     #D1 * t23 * RX * FX - t23 * t16 - t23 * t17 + 0.2D1 * t23 * RY * FY
     # - t23 * t20 - 0.4D1 * t23 * t6 * E) / t48
        return
      end